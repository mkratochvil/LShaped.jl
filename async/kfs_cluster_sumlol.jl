##using Pkg
#Pkg.activate("..")

#using JuMP
#using Gurobi
#using DataFrames
#using CSV
#using LinearAlgebra

#using LShaped

#this would be an external variable
infoloc = "./info.csv"

#load in info
info = CSV.File(infoloc) |> Dict

#add into info
delta = 1.0
tau = 1e-10
deltamax = 25.#12.8
deltamin = 0.05#.00078125
gamma = 0.01

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
        delta = LShaped.get_delta_from_file(dataloc)
        
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
        model, curit, ssobja, delta, ssobjx, delcrit, numcuts, lb, abscrit = LShaped.resume_fs_trust!(firststage, vardict, model, nsubs, delta, header, deltamin, deltamax, gamma)
        # + load in each model using iterations that require cuts
        # + track the last iteration, and whether a cut was made, if sum of cuts made is 0, set converged = 1
        # + return current iterations "curit" (files are required to match) and "converged"
    else
        if converged == 1
            converged = 0
            info["converged"] = "0"
            CSV.write(infoloc, info)
        end
        LShaped.setup_1st_paths_trust!(firststage, nsubs)
        
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

        JuMP.optimize!(model)
        fsobj = JuMP.objective_value(model)
        LShaped.store_fsobj!(fsobj, dataloc)
        LShaped.store_delta!(delta, dataloc)

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
            println(curit, " ", lb, " ", ssobja, " ", abs(lb-ssobja)/(tau + abs(ssobja)))
            #iteration_summary goes here
            LShaped.store_iteration_summary!(curit, ssobja, lb, ssobjx, abscrit, abs(lb-ssobja)/(tau + abs(ssobja)), delta, delcrit, numcuts, dataloc)
            if abs(lb - ssobja) < tol*(tau + abs(ssobja))

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


