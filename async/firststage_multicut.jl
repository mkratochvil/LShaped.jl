#this would be an external variable
infoloc = "./info.csv"

#load in info
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
        
        #update fs E,e, w, theta here
        ##el = LShaped.get_el_from_file_and_save!(dataloc, nsubs)
        ##El = LShaped.get_El_from_file_and_save!(firststage, model, dataloc, nsubs)
        
        
        x = LShaped.get_x_from_file(dataloc)
        
        ##theta = LShaped.get_theta_from_file(dataloc)
        
        ##w = el - dot(El,x) 
                
        #println("..........w-theta difference: $(w)-$(theta) = $(w-theta).................")
        ##LShaped.store_w_theta!(firststage, w, theta)

        ##curit, converged, xcur = LShaped.resume_fs!(firststage, model, tol)
        model, curit, converged = LShaped.resume_fs_multicut!(firststage, vardict, model, nsubs)
        # + load in each model using iterations that require cuts
        # + track the last iteration, and whether a cut was made, if sum of cuts made is 0, set converged = 1
        # + return current iterations "curit" (files are required to match) and "converged"
    else
        if converged == 1
            converged = 0
            info["converged"] = "0"
            CSV.write(infoloc, info)
        end
        ##LShaped.setup_1st_paths!(firststage)
        LShaped.setup_1st_paths_multicut!(firststage, nsubs)
        # - Create x vars
        # - Create theta vars
    end
    println("converged = $(converged)")
    
    if converged == 0
        
        curit += 1
        println("curit = $(curit)")

        JuMP.optimize!(model)

        LShaped.update_first_value_L!(firststage, model)

        LShaped.store_x!(firststage)
        
        thetavec = [];
        if curit > 1
            for i = 1:nsubs
                push!(thetavec, JuMP.value(JuMP.variable_by_name(model, "theta_$(i)")))
                ##thetavalue = JuMP.value(JuMP.variable_by_name(model, "theta"))
            end
            ##LShaped.store_theta!(firststage, thetavalue)
        else
            ##LShaped.store_theta!(firststage, -Inf)
            for i = 1:nsubs
                push!(thetavec, -Inf)
            end
        end
        LShaped.store_thetak!(firststage, thetavec);

        
    else
        #update converged on info somehow
        println("Model converged.")
        info["converged"] = "1"
        CSV.write(infoloc, info)
    end

else
    println("Yay it worked")
end


