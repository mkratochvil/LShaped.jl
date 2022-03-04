#todo: change this to the arrayid env variable in cluster
#=filestring = string(@__FILE__)
arrayid = 0
n = length(filestring)
if filestring[n-4] == 'c'
    arrayid = parse(Int64,filestring[n-3])
elseif filestring[n-5] == 'c'
    arrayid = parse(Int64,filestring[n-4:n-3])
elseif filestring[n-6] == 'c'
    arrayid = parse(Int64,filestring[n-5:n-3])
elseif filestring[n-7] == 'c'
    arrayid = parse(Int64,filestring[n-6:n-3])
end=#

using LShaped

## make sure these files and their dependencies are in their proper place ##
include("../../FinalProject/parameters.jl")
include("../../FinalProject/get_functions.jl")
loadcsv = CSV.File("../../FinalProject/LOAD.csv");
ptdfdf = DataFrame(CSV.File("../../FinalProject/ptdfsmall.csv"));

include("../../FinalProject/modification_functions_async.jl")
##

## upload these to cluster ##
loaddf = DataFrame(CSV.File("../../FinalProject/loaddata.csv"))
winddf = DataFrame(CSV.File("../../FinalProject/winddata.csv"))
##
## make sure this is uploaded with everything ## 
ptdfdf = DataFrame(CSV.File("../../FinalProject/ptdfsmall.csv"));
##

let 
    lmax = 0.0
    wmax = 0.0
    for i = 1:24
        #println(i, " ", maximum(loaddf[:,2+i]))
        if maximum(loaddf[:,2+i]) > lmax
            lmax = maximum(loaddf[:,2+i])
        end
        #println(maximum(winddf[:,2+i]))
        if maximum(winddf[:,2+i]) > wmax
            wmax = maximum(winddf[:,2+i])
        end
    end
    global lmax = lmax
    global wmax = wmax
end

loaddis = load_distribution_dict(loadcsv);

ptdfdict = Dict()

for i = 1:38
    br = ptdfdf[i,1]
    ptdfdict[br] = Dict()
    for j = 2:25
        bus = parse(Int64,names(ptdfdf)[j])
        ptdfdict[br][bus] = ptdfdf[i,j]
    end
end

lrts = 2850.0
wrts = 713.5

global ercotscens = collect(DataFrame(CSV.File("../../FinalProject/scenarios/shuff133.csv"))[1,:])

#this would be an external variable
infoloc = "./info.csv"

#load in info
info = CSV.File(infoloc) |> Dict

converged = parse(Int64,info["converged"])
tol = parse(Float64,info["tol"])

if converged == 0

    pathloc = info["dir"]

    dataloc = string(pathloc, "data/")

    ssmodelloc = string(pathloc,info["ssscript"])

    #to load in firststage model
    include(ssmodelloc)

    fsfuncname = info["fsmodel"]
    ssfuncname = info["ssmodel"]

    # build first stage model

    model1 = getfield(Main,Symbol(fsfuncname))()

    vardict = getfield(Main,Symbol(info["vardict"]))

    parse(Int64,info["nsubs"])

    nsubs = parse(Int64,info["nsubs"])
    
    firststagevars = Dict()

    for index in 1:length(vardict[1])
        var = vardict[1][index]
        firststagevars[var[1]] = LShaped.FirstStageVariableInfo(var[1], index, var[4], nothing, var[2], var[3])
    end
    
    let
    #model2 = getfield(Main,Symbol(info["ssmodel"]))(1);
        model2 = JuMP.read_from_file("../../FinalProject/storage_expansion_revised/second_stage/noint_PR_exp36p00015_scen_1.mps")
        JuMP.set_optimizer(model2, Gurobi.Optimizer)
        set_optimizer_attribute(model2, "OutputFlag", 0)

        #model2 = getfield(Main,Symbol(info["ssmodel"]))(arrayid);
        for i  = 1:12
            #println("Adjusting to scenario $(ercotscens[i])")
            wind = (1/wmax)*(wrts/100)*collect(winddf[ercotscens[i],3:26])
            load = (1/lmax)*(1.35*lrts/100)*collect(loaddf[ercotscens[i],3:26])

            for bus in buses
                lf = loaddis[bus]
                for ts in timesteps
                    con = get_load_balance(model2, bus, ts)
                    oldval = JuMP.constraint_object(con).set.value
                    lval = load[ts]
                    JuMP.set_normalized_rhs(con, lf*lval)
                    newval = JuMP.constraint_object(con).set.value
                    #println("$(name(con)), $(oldval), $(newval)")
                end
            end

            # change ptdf constraint (remember to run load changes FIRST)
            for ts in timesteps
                for br in branches
                    ptdfcon = get_ptdf_con(model2,br,ts)

                    valold = JuMP.constraint_object(ptdfcon).set.value
                    valnew = 0.0
                    for bus in buses
                        buscon = get_load_balance(model2,bus,ts)

                        loadcon = copy(JuMP.constraint_object(buscon).func)
                        loadval = copy(JuMP.constraint_object(buscon).set.value)

                        valnew -= ptdfdict[br][bus]*loadval

                    end 

                    JuMP.set_normalized_rhs(ptdfcon, valnew)
                    #println("$(JuMP.name(ptdfcon)), $(valold), $(valnew)")

                end
            end

            bus = 122
            for ts in timesteps
                con = get_wind_ub(model2, bus, ts)
                oldval = JuMP.constraint_object(con).set.upper
                wval = wind[ts]
                JuMP.set_normalized_rhs(con, wval)
                newval = JuMP.constraint_object(con).set.upper
                #println("$(name(con)), $(oldval), $(newval)")
            end

            idxtocon = LShaped.IdxToCon(model2)

            #println("computing h for subproblem $(arrayid)...")
            ### rework h as a file to potentially save time
            #h = LShaped.compute_h_new(model2, idxtocon)
            ### 
            h = collect(DataFrame(CSV.File("../../FinalProject/storage_expansion_revised/ercot/h_exp3_scen_$(ercotscens[i]).csv"))[1,:])

            #println("Initializing subproblem $(arrayid)...")
            model2, varstructs, vnametoidx = LShaped.initialize(model2, vardict)

            #println("Creating subprob[$(arrayid)] struct...")
            subprob = LShaped.SubproblemsNew(i#=arrayid=#, model2, #=0.4=#1/12, varstructs, idxtocon, h, 
                    nothing, nothing, vnametoidx, nothing)

            #empty dict temporary until I can figure out how to store things.
            subprobs = Dict([(i#=arrayid=#, subprob)])

            firststage = LShaped.FirstStageInfo(firststagevars, subprobs, dataloc);

            # make dict from var to value

            pathx = string(dataloc, "x.csv")
            x = DataFrame(CSV.File(pathx))
            xvars = String.(names(x))
            curit = size(x,1)
            xvals = collect(x[curit,:])


            xpathdict = Dict()

            for i = 1:length(xvars)
                xpathdict[xvars[i]] = xvals[i]
            end

            #update_second_value equiv
            for var in keys(xpathdict)
                sub = firststage.subproblems[i#=arrayid=#]
                vind = sub.vnametoind[var]
                value = xpathdict[var]
                sub.variableinfo[vind].value = value
            end

            LShaped.solve_sub_and_update!(firststage.subproblems[i#=arrayid=#])
            LShaped.update_second_index!(firststage)

            if size(x,1) == 1
                LShaped.setup_scen_path!(dataloc, i#=arrayid=#)
                LShaped.setup_2nd_paths_regdec!(dataloc, firststage.subproblems[i#=arrayid=#])
            end

            subproblem = firststage.subproblems[i#=arrayid=#]
            #we are set up so that the models have no interval constraints. 
            #LShaped.adjust_h_new!(subproblem) 
            LShaped.compute_Ek_new!(subproblem)
            LShaped.compute_ek_new!(subproblem)
            #println("model1 = ",model1)
            cost = LShaped.get_cost_vector(firststage, model1)

            subproblem.Ek = subproblem.probability*cost + subproblem.Ek
            Ek = subproblem.Ek
            ek = subproblem.ek
            #note I did not check if it was in the same order. It should be?
            w = ek - dot(Ek,xvals)

            #load in theta (only need specific theta)
            patht = string(dataloc, "theta.csv")
            thetadf = DataFrame(CSV.File(patht))
            dfend = size(thetadf,1)
            theta = thetadf[dfend,i#=arrayid=#]

            addcut = 0
            if theta < w - tol
                addcut = 1
            end

            #-compare Ek, and w, and a value whether or not to add the cut
            #-something else....look it up. sorry future (i.e. present) me. Love, past me.
            LShaped.store_Ek_sub!(subproblem, dataloc)
            LShaped.store_ek_sub!(subproblem, dataloc)
            LShaped.store_cut_sub!(addcut, dataloc, i#=arrayid=#)

            objval = subproblem.probability*JuMP.objective_value(subproblem.model)

            LShaped.store_objval!(objval, dataloc, i#=arrayid=#)
        end
    end

else
    println("Yay It worked")
    
end



