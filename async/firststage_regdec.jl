#this would be an external variable
infoloc = "./info.csv"

#load in info
info = CSV.File(infoloc) |> Dict

#add into info
rho = 1.0
tau = 1e-7
rhomax = 1.0e3#12.8
rhomin = 1.0e-3#.00078125
gamma = 0.01
#avec = [0.01515317295046937, 0.07576586475234684, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.06992444968258622, 0.3496222484129311, 0.024425242996275243, 0.10597209391518082, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.09035704399269796, 0.4517852199634898, 0.0, 0.0, 0.0, 0.0, 0.01246785328231903, 0.06233926641159515, 0.0053185494034406755, 0.026592747017203378, 0.0, 0.0, 0.08532203616288067, 0.4266101808144034, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.2892514173319449, 1.4462570866597244, 0.0, 0.0]
#ssobja = 1507.3174544213448
 

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
        rho = LShaped.get_rho_from_file(dataloc)
        
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
        model, curit, ssobja, rho = LShaped.resume_fs_regdec!(firststage, vardict, model, nsubs, rho, header, rhomin, rhomax, gamma)
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
        
        #=
        ## create regularized problem with initial rho and avec
        vdicth = firststage.variables
        nvars = vdicth.count

        header = Array{Union{Nothing, String}}(nothing, nvars)

        for vname in keys(vdicth)

            varinfo = vdicth[vname]
            index = varinfo.index

            header[index] = vname

        end
        LShaped.add_regularized_decomp_async!(model, avec, header, rho)
        =#
        
    end
    #println("converged = $(converged)")
    
    if converged == 0
        
        curit += 1
        #println("curit = $(curit)")

        JuMP.optimize!(model)
        fsobj = JuMP.objective_value(model)
        ### add this
        LShaped.store_fsobj!(fsobj, dataloc)
        LShaped.store_rho!(rho, dataloc)

        LShaped.update_first_value_L!(firststage, model)

        LShaped.store_x!(firststage)
        
        thetavec = [];
        ## replaced avec here for when an initial avec not decided
        if curit > 1
            ## replaced avec here for when an initial avec not decided
            avec = LShaped.get_value_vector(firststage)
            for i = 1:nsubs
                push!(thetavec, JuMP.value(JuMP.variable_by_name(model, "theta_$(i)")))
                ##thetavalue = JuMP.value(JuMP.variable_by_name(model, "theta"))
            end
            #check for convergence
            println(curit, " ", fsobj, " ", ssobja, " ", abs(fsobj-ssobja)/(tau + abs(ssobja)))
            if abs(fsobj - ssobja) < tol*(tau + abs(ssobja))

                println("algorithm converged.")
                println("final x = $(avec)")
                
                info["converged"] = "1"
                CSV.write(infoloc, info)
            end
            
        else
            for i = 1:nsubs
                push!(thetavec, -Inf)
            end
            #LShaped.store_a!(avec, dataloc)
            #LShaped.store_ssobja!(ssobja, dataloc)
            
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


