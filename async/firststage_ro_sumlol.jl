#
## make sure these files and their dependencies are in their proper place ##
include("../../FinalProject/parameters.jl")
include("../../FinalProject/get_functions.jl")
loadcsv = CSV.File("../../FinalProject/RTS_data_summary/LOAD.csv");
ptdfdf = DataFrame(CSV.File("../../FinalProject/RTS_data_summary/ptdfsmall.csv"));

expid = 1

include("../../FinalProject/modification_functions.jl")
##

## upload these to cluster ##
loaddf = DataFrame(CSV.File("../../FinalProject/four_wind_turbine_experiments/loaddata.csv"))
winddfdict = Dict()
for wb in wind_buses
    winddfdict[wb] = DataFrame(CSV.File("../../FinalProject/four_wind_turbine_experiments/winddata$(wb).csv"))
end
##

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

global ercotscens = collect(DataFrame(CSV.File("../../FinalProject/scenarios/part1458.csv"))[expid,:])

#this would be an external variable
infoloc = "./info.csv"

#load in infoa
info = CSV.File(infoloc) |> Dict

converged = parse(Int64,info["converged"])

#this exists so that if the data folder does not exist and converged is 1, we can start over, not manually.
dataloc = string(info["dir"], "data/")

if converged == 0 || ispath(dataloc) == 0

    pathloc = info["dir"]
    tol = parse(Float64,info["tol"])
    nsubs = parse(Int64,info["nsubs"])

    dataloc = string(pathloc, "data/")

    fsmodelloc = string(pathloc,info["fsscript"])

    #to load in firststage model
    include(fsmodelloc)

    fsfuncname = info["fsmodel"]

    # build first stage model
    model = getfield(Main,Symbol(fsfuncname))()

    vardict = getfield(Main,Symbol(info["vardict"]))

    firststagevars = Dict()

    for index in 1:length(vardict[1])
        var = vardict[1][index]
        firststagevars[var[1]] = LShaped.FirstStageVariableInfo(var[1], index, var[4], nothing, var[2], var[3])
    end

    #empty dict temporary until I can figure out how to store things.
    subprob = Dict();

    firststage = LShaped.FirstStageInfo(firststagevars, subprob, dataloc);
    
    curit = 0
    if ispath(string(dataloc))
        
        x = LShaped.get_x_from_file(dataloc)

        #curit, converged, xcur = LShaped.resume_fs!(firststage, model, tol) ## replace with adding theta and models directly?
        # update UB and get new scenario to be added to first stage
        objval = -Inf
        newscen = 0
        for i = 1:nsubs
            subobjdf = DataFrame(CSV.File(string(dataloc,"scen_$(i)/objval.csv")))
            objlength = size(subobjdf,1)
            subobj = subobjdf[objlength,1]
            
            if subobj > objval
                global objval = subobj
                global newscen = i
            end
        end
        println("newscen = $(newscen)")
        
        UBdf = DataFrame(CSV.File(string(dataloc, "ub.csv")))
        UB = UBdf[size(UBdf,1),1]
               
        UB = min(UB, objval)
        
        LShaped.store_UB!(UB, dataloc)
        LShaped.store_scen!(newscen, dataloc)
        
        LBdf = DataFrame(CSV.File(string(dataloc, "lb.csv")))
        LB = LBdf[size(LBdf,1),1]
        
        curit = size(LBdf,1)
        
        println(curit, " ", abs(UB - LB)/(abs(UB)+1.0e-7))
        
        if abs(UB - LB)/(abs(UB)+1.0e-7) < tol
            converged = 1
        else
            # set up for another run
            LShaped.add_theta_to_objective!(model)
            scenset = DataFrame(CSV.File(string(dataloc,"scenarios.csv")))
            
            let
                model2nd = JuMP.read_from_file("../../FinalProject/sum_lol/sum_lol/second_stage_model_lolsum_scen_1.mps")
                
                for i = 1:size(scenset,1)
                    scen = scenset[i,1]
                    #println("Adjusting to scenario $(ercotscens[i])")
                    #wind = (1/wmax)*(wrts/100)*collect(winddf[ercotscens[scen],3:26])
                    wind = Dict()
                    for wb in wind_buses
                        wind[wb] = (1/wmaxdict[wb])*(wrtsdict[wb]/100)*collect(winddfdict[wb][ercotscens[scen],3:26])
                    end
                    load = (1/lmax)*(1.35*lrts/100)*collect(loaddf[ercotscens[scen],3:26])

                    for bus in buses
                        lf = loaddis[bus]
                        for ts in timesteps
                            con = get_load_balance(model2nd, bus, ts)
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
                            ptdfcon = get_ptdf_con(model2nd,br,ts)

                            valold = JuMP.constraint_object(ptdfcon).set.value
                            valnew = 0.0
                            for bus in buses
                                buscon = get_load_balance(model2nd,bus,ts)

                                loadcon = copy(JuMP.constraint_object(buscon).func)
                                loadval = copy(JuMP.constraint_object(buscon).set.value)

                                valnew -= ptdfdict[br][bus]*loadval

                            end 

                            JuMP.set_normalized_rhs(ptdfcon, valnew)
                            #println("$(JuMP.name(ptdfcon)), $(valold), $(valnew)")

                        end
                    end

                    #=
                    bus = 122
                    for ts in timesteps
                        con = get_wind_ub(model2nd, bus, ts)
                        oldval = JuMP.constraint_object(con).set.upper
                        wval = wind[ts]
                        JuMP.set_normalized_rhs(con, wval)
                        newval = JuMP.constraint_object(con).set.upper
                        #println("$(name(con)), $(oldval), $(newval)")
                    end
                    =#
                    for wb in wind_buses
                        for ts in timesteps
                            con = get_wind_ub(model2nd, wb, ts)
                            oldval = JuMP.constraint_object(con).set.upper
                            wval = wind[wb][ts]
                            JuMP.set_normalized_rhs(con, wval)
                            newval = JuMP.constraint_object(con).set.upper
                            #println("$(name(con)), $(oldval), $(newval)")
                        end
                    end

                    LShaped.add_wc_variables_to_fs_model!(model, vardict, scen)
                    v1array = LShaped.get_1st_stage_variable_array(vardict)
                    LShaped.add_wc_constraints!(model, model2nd, scen, v1array)
                    LShaped.add_obj_constraint!(model, model2nd, scen, v1array)
                end
            end
        end
        
    else
        if converged == 1
            converged = 0
            info["converged"] = "0"
            CSV.write(infoloc, info)
        end
        LShaped.setup_ro_paths!(firststage, nsubs) 
    end
    #println("converged = $(converged)")
    
    if converged == 0
        
        curit += 1

        JuMP.optimize!(model)

        LShaped.update_first_value_L!(firststage, model)
        LShaped.store_x!(firststage) 

        if curit > 1 
            LB = JuMP.objective_value(model)
            LShaped.store_LB!(LB, dataloc)
        else
            LShaped.store_LB!(-Inf, dataloc)
            LShaped.store_UB!(Inf, dataloc)
        end
        
    else
        println("Model converged.")
        info["converged"] = "1"
        CSV.write(infoloc, info)
    end

else
    println("Yay it worked")
end


