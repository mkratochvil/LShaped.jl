#this would be an external variable
infoloc = "./info.csv"

#load in info
info = CSV.File(infoloc) |> Dict

#add into info
rho = 0.01

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
        
        vdicth = firststage.variables
        nvars = vdicth.count

        header = Array{Union{Nothing, String}}(nothing, nvars)

        for vname in keys(vdicth)

            varinfo = vdicth[vname]
            index = varinfo.index

            header[index] = vname

        end

        ##curit, converged, xcur = LShaped.resume_fs!(firststage, model, tol)
        model, curit, ssobja = LShaped.resume_fs_regdec!(firststage, vardict, model, nsubs, rho, header)
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
        LShaped.setup_1st_paths_regdec!(firststage, nsubs)
        # - Create x vars
        # - Create theta vars
    end
    #println("converged = $(converged)")
    
    if converged == 0
        
        curit += 1
        #println("curit = $(curit)")

        JuMP.optimize!(model)

        LShaped.update_first_value_L!(firststage, model)

        LShaped.store_x!(firststage)
        
        thetavec = [];
        avec = LShaped.get_value_vector(firststage)
        if curit > 1
            
            for i = 1:nsubs
                push!(thetavec, JuMP.value(JuMP.variable_by_name(model, "theta_$(i)")))
                ##thetavalue = JuMP.value(JuMP.variable_by_name(model, "theta"))
            end
            #check for convergence
            fsobj = JuMP.objective_value(model)
            println(curit, " ", fsobj, " ", ssobja, " ", abs(fsobj-ssobja))
            if abs(fsobj - ssobja) < tol

                println("algorithm converged.")
                println("final x = $(avec)")
                
                info["converged"] = "1"
                CSV.write(infoloc, info)
            end
            
        else
            for i = 1:nsubs
                push!(thetavec, -Inf)
            end
            LShaped.store_a!(avec, dataloc)
            
            #ssobja = JuMP.objective_value(model)
            #LShaped.store_ssobja!(ssobja, dataloc)
        end
        
        LShaped.store_thetak!(firststage, thetavec);

        
    else
        #update converged on info somehow
        println("This should never happen.")
        info["converged"] = "1"
        CSV.write(infoloc, info)
    end

else
    println("Yay it worked")
end


