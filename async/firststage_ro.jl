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
            for i = 1:size(scenset,1)
                scen = scenset[i,1]
                model2nd = getfield(Main,Symbol(info["ssmodel"]))(scen);
                
                LShaped.add_wc_variables_to_fs_model!(model, vardict, scen)
                v1array = LShaped.get_1st_stage_variable_array(vardict)
                LShaped.add_wc_constraints!(model, model2nd, scen, v1array)
                LShaped.add_obj_constraint!(model, model2nd, scen, v1array)
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


