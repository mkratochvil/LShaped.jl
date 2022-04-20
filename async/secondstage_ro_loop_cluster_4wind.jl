##todo: change this to the arrayid env variable in cluster
#=filestring = string(@__FILE__)
arrayid = 0
n = length(filestring)
if filestring[n-4] == 'o'
    arrayid = parse(Int64,filestring[n-3])
elseif filestring[n-5] == 'o'
    arrayid = parse(Int64,filestring[n-4:n-3])
elseif filestring[n-6] == 'o'
    arrayid = parse(Int64,filestring[n-5:n-3])
elseif filestring[n-7] == 'o'
    arrayid = parse(Int64,filestring[n-6:n-3])
end=#

using LShaped

K = 24
expid = 5
clid = 24

## make sure these files and their dependencies are in their proper place ##
include("../../FinalProject/parameters.jl")
include("../../FinalProject/get_functions.jl")
loadcsv = CSV.File("../../FinalProject/RTS_data_summary/LOAD.csv");
ptdfdf = DataFrame(CSV.File("../../FinalProject/RTS_data_summary/ptdfsmall.csv"));

include("../../FinalProject/modification_functions.jl")
##

## upload these to cluster ##
loaddf = DataFrame(CSV.File("../../FinalProject/four_wind_turbine_experiments/loaddata.csv"))
#winddf = DataFrame(CSV.File("../../FinalProject/winddata.csv"))
winddfdict = Dict()
for wb in wind_buses
    winddfdict[wb] = DataFrame(CSV.File("../../FinalProject/four_wind_turbine_experiments/winddata$(wb).csv"))
end


#=
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
=#

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

#lrts = 2850.0
#wrts = 713.5

#global ercotscens = collect(DataFrame(CSV.File("../../FinalProject/scenarios/part1458.csv"))[expid,:])
#global ercotscens = collect(DataFrame(CSV.File("../../FinalProject/scenarios/part1458.csv"))[expid,:])
cldict = CSV.File("../../FinalProject/scenarios/four_wind/cldict_$(K)_$(expid).csv") |> Dict

vect = Int64[];
n = length(cldict[clid])
strvals = cldict[clid][2:n-1]
for i = 1:length(strvals)
    if strvals[i] == ','
        if i-1 == 1
            push!(vect, parse(Int64, strvals[i-1]))
        elseif i-2 == 1
            push!(vect, parse(Int64, strvals[i-2:i-1]))
        elseif i-3 == 1
            push!(vect, parse(Int64, strvals[i-3:i-1]))
        elseif i-4 == 1
            push!(vect, parse(Int64, strvals[i-4:i-1]))
        elseif strvals[i-2] == ' ' 
            push!(vect,parse(Int64, strvals[i-1]))
        elseif strvals[i-3] == ' '
            push!(vect,parse(Int64, strvals[i-2:i-1]))
        elseif strvals[i-4] == ' '
            push!(vect,parse(Int64, strvals[i-3:i-1]))
        elseif strvals[i-5] == ' '
            push!(vect,parse(Int64, strvals[i-4:i-1]))
        else
            println("Uh oh!")
        end
    elseif i == n-2
        if i == 1
            push!(vect, parse(Int64, strvals[i]))
        elseif i-1 == 1
            push!(vect, parse(Int64, strvals[i-1:i]))
        elseif i-2 == 1
            push!(vect, parse(Int64, strvals[i-2:i]))
        elseif i-3 == 1
            push!(vect, parse(Int64, strvals[i-3:i]))
        elseif strvals[i-1] == ' ' 
            push!(vect,parse(Int64, strvals[i]))
        elseif strvals[i-2] == ' '
            push!(vect,parse(Int64, strvals[i-1:i]))
        elseif strvals[i-3] == ' '
            push!(vect,parse(Int64, strvals[i-2:i]))
        elseif strvals[i-4] == ' '
            push!(vect,parse(Int64, strvals[i-3:i]))
        else
            println("Uh oh!")
        end

    end
end
global ercotscens = copy(vect)
println()
println("Number of scenarios: ", length(ercotscens))
println


#this would be an external variable
infoloc = "./info.csv"

#load in info
info = CSV.File(infoloc) |> Dict

converged = parse(Int64,info["converged"])

if converged == 0

    print("Second stage: ")
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

    nsubs = length(vect)#parse(Int64,info["nsubs"])

    let
        #model2 = getfield(Main,Symbol(info["ssmodel"]))(1);
        #model2 = JuMP.read_from_file("../../FinalProject/storage_expansion_revised/ercot/PR_exp3_scen_1.mps")
        model2 = JuMP.read_from_file("../../FinalProject/four_wind_turbine_experiments/four_wind/second_stage_model_scen_1.mps")
        JuMP.set_optimizer(model2, Gurobi.Optimizer)
        set_optimizer_attribute(model2, "OutputFlag", 0) 
        
        for i = 1:length(ercotscens)
            print("$(i)..")

            #println("Adjusting to scenario $(ercotscens[i])")
            #wind = (1/wmax)*(wrts/100)*collect(winddf[ercotscens[i],3:26])
            wind = Dict()
            for wb in wind_buses
                wind[wb] = (1/wmaxdict[wb])*(wrtsdict[wb]/100)*collect(winddfdict[wb][ercotscens[i],3:26])
            end
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

            #=bus = 122
            for ts in timesteps
                con = get_wind_ub(model2, bus, ts)
                oldval = JuMP.constraint_object(con).set.upper
                wval = wind[ts]
                JuMP.set_normalized_rhs(con, wval)
                newval = JuMP.constraint_object(con).set.upper
                #println("$(name(con)), $(oldval), $(newval)")
            end=#
            for wb in wind_buses
                for ts in timesteps
                    con = get_wind_ub(model2, wb, ts)
                    oldval = JuMP.constraint_object(con).set.upper
                    wval = wind[wb][ts]
                    JuMP.set_normalized_rhs(con, wval)
                    newval = JuMP.constraint_object(con).set.upper
                    #println("$(name(con)), $(oldval), $(newval)")
                end
            end
            
            idxtocon = LShaped.IdxToCon(model2)

            #println("computing h for subproblem $(arrayid)...") ## unneccessary and waste of time
            #h = LShaped.compute_h_new(model2, idxtocon)
            # do I even need h for robust case?
            h = collect(DataFrame(CSV.File("../../FinalProject/storage_expansion_revised/ercot/h_exp3_scen_$(ercotscens[i]).csv"))[1,:])

            #println("Initializing subproblem $(arrayid)...")
            model2, varstructs, vnametoidx = LShaped.initialize(model2, vardict)

            #println("Creating subprob[$(arrayid)] struct...")
            subprob = LShaped.SubproblemsNew(i#=arrayid=#, model2, 1/12, varstructs, idxtocon, h, 
                    nothing, nothing, vnametoidx, nothing)

            firststagevars = Dict()

            for index in 1:length(vardict[1])
                var = vardict[1][index]
                firststagevars[var[1]] = LShaped.FirstStageVariableInfo(var[1], index, var[4], nothing, var[2], var[3])
            end

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

            subproblem = firststage.subproblems[i#=arrayid=#] 
            objval = subproblem.objective_value
            LShaped.store_objval!(objval, i#=arrayid=#, dataloc)

        end
    end

else
    println("Yay It worked")
    
end



