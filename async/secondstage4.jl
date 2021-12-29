#todo: change this to the arrayid env variable in cluster
arrayid = 4

using LShaped

#this would be an external variable
infoloc = "./info.csv"

#load in info
info = CSV.File(infoloc) |> Dict

converged = parse(Int64,info["converged"])

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

    model2 = getfield(Main,Symbol(info["ssmodel"]))(arrayid);

    idxtocon = LShaped.IdxToCon(model2)

    println("computing h for subproblem $(arrayid)...")
    h = LShaped.compute_h_new(model2, idxtocon)

    println("Initializing subproblem $(arrayid)...")
    model2, varstructs, vnametoidx = LShaped.initialize(model2, vardict)

    println("Creating subprob[$(arrayid)] struct...")
    subprob = LShaped.SubproblemsNew(arrayid, model2, 1/12, varstructs, idxtocon, h, 
            nothing, nothing, vnametoidx, nothing)

    firststagevars = Dict()

    for index in 1:length(vardict[1])
        var = vardict[1][index]
        firststagevars[var[1]] = LShaped.FirstStageVariableInfo(var[1], index, var[4], nothing, var[2], var[3])
    end

    #empty dict temporary until I can figure out how to store things.
    subprobs = Dict([(arrayid, subprob)])

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
        sub = firststage.subproblems[arrayid]
        vind = sub.vnametoind[var]
        value = xpathdict[var]
        sub.variableinfo[vind].value = value
    end
    
    ## warm-start attempt here. ##
    #path_id = string(dataloc, "scen_$(arrayid)/")
    #if curit > 1
    #    LShaped.warmstart(firststage.subproblems[arrayid], curit-1, path_id)
    #end

    LShaped.solve_sub_and_update!(firststage.subproblems[arrayid])
    LShaped.update_second_index!(firststage)

    if size(x,1) == 1
        LShaped.setup_scen_path!(dataloc, arrayid)
        LShaped.setup_2nd_paths!(dataloc, firststage.subproblems[arrayid])
    end
    


    subproblem = firststage.subproblems[arrayid]
    LShaped.adjust_h_new!(subproblem) 
    LShaped.compute_Ek_new!(subproblem)
    LShaped.compute_ek_new!(subproblem)

    LShaped.store_Ek_sub!(subproblem, dataloc)
    LShaped.store_ek_sub!(subproblem, dataloc)
    
    #path_id = string(dataloc, "scen_$(arrayid)/")
    #LShaped.save_cur_vars!(path_id, subproblem.model, curit)
    #LShaped.save_cur_duals!(path_id, subproblem, curit)

    #println(JuMP.solution_summary(subproblem.model, verbose=true))
    #println(".........Primal dual difference = 
    #    $(JuMP.objective_value(subproblem.model) - JuMP.dual_objective_value(subproblem.model)).........")
    
else
    println("Yay It worked")
    
end



