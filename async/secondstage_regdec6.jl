#todo: change this to the arrayid env variable in cluster
filestring = string(@__FILE__)
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
end

using LShaped

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

    model2 = getfield(Main,Symbol(info["ssmodel"]))(arrayid);

    idxtocon = LShaped.IdxToCon(model2)

    #println("computing h for subproblem $(arrayid)...")
    ### rework h as a file to potentially save time
    h = LShaped.compute_h_new(model2, idxtocon)
    ### 

    #println("Initializing subproblem $(arrayid)...")
    model2, varstructs, vnametoidx = LShaped.initialize(model2, vardict)

    #println("Creating subprob[$(arrayid)] struct...")
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

    LShaped.solve_sub_and_update!(firststage.subproblems[arrayid])
    LShaped.update_second_index!(firststage)

    if size(x,1) == 1
        LShaped.setup_scen_path!(dataloc, arrayid)
        LShaped.setup_2nd_paths_regdec!(dataloc, firststage.subproblems[arrayid])
    end

    subproblem = firststage.subproblems[arrayid]
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
    theta = thetadf[dfend,arrayid]
    
    addcut = 0
    if theta < w - tol
        addcut = 1
    end
    
    #-compare Ek, and w, and a value whether or not to add the cut
    #-something else....look it up. sorry future (i.e. present) me. Love, past me.
    LShaped.store_Ek_sub!(subproblem, dataloc)
    LShaped.store_ek_sub!(subproblem, dataloc)
    LShaped.store_cut_sub!(addcut, dataloc, arrayid)
    
    objval = subproblem.probability*JuMP.objective_value(subproblem.model)
    
    LShaped.store_objval!(objval, dataloc, arrayid)
    
else
    println("Yay It worked")
    
end



