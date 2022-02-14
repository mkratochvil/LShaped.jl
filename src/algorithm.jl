function update_constraint_values_nac!(model::JuMP.Model, varstructs::Dict{Int64,LocalVariableInfo})
    
    for var in keys(varstructs)
        name = varstructs[var].name
        value = varstructs[var].value
        JuMP.fix(JuMP.variable_by_name(model,name),value)
    end
    
    return
end
    
function update_gradients_nac!(model::JuMP.Model, varstructs::Dict{Int64,LocalVariableInfo})
    
    for var in keys(varstructs)
        name = varstructs[var].name
    
        varstructs[var].gradient = JuMP.dual(JuMP.FixRef(JuMP.variable_by_name(model,name)))
        
    end
    
    return
    
end

function solve_sub_and_update!(subproblem::Union{Subproblems,SubproblemsNew})
        
    varstructs = subproblem.variableinfo
    model = subproblem.model

    update_constraint_values_nac!(model, varstructs)
    
    #println("...Solving subproblem: $(subproblem.id)...")
    optimize!(model)
    update_gradients_nac!(model, varstructs)
    
    subproblem.objective_value = objective_value(model)
    
    return
    
end

function get_cost_vector(firststage::FirstStageInfo, fsmodel::JuMP.Model)
        
    vardict = firststage.variables
    nvars = vardict.count
    
    cost = Vector{Float64}(undef, nvars) 
    
    for vname in keys(vardict)
        
        varinfo = vardict[vname]
        index = varinfo.index
        
        vref = JuMP.variable_by_name(fsmodel, vname)
               
        if vref in keys(JuMP.objective_function(fsmodel).terms)
            cost[index] = JuMP.objective_function(fsmodel).terms[vref]
        else
            cost[index] = 0.0
        end
        
    end
    
    return cost
    
end

function update_first_gradient!(firststage::FirstStageInfo)
        
    for var in keys(firststage.variables)
        grad1 = 0.0;
        # eventually use get functions in parallel.
        for sub = keys(firststage.subproblems)
            prob = firststage.subproblems[sub].probability
            vind = firststage.subproblems[sub].vnametoind[var]
            grad2 = firststage.subproblems[sub].variableinfo[vind].gradient
            grad1 += prob*grad2
        end
        firststage.variables[var].gradient = grad1
    end
    
    return
    
end

#eventually parallelize and put another function inside.
function update_second_value!(firststage::FirstStageInfo) 
    
    for var in keys(firststage.variables)
        for sub = keys(firststage.subproblems)
            subproblem = firststage.subproblems[sub]
            vind = subproblem.vnametoind[var]
            value = firststage.variables[var].value
            subproblem.variableinfo[vind].value = value
        end
    end
    
    return
    
end

function get_grad_vector(firststage::FirstStageInfo)
        
    vardict = firststage.variables
    nvars = vardict.count
    
    grad = Vector{Float64}(undef, nvars) 
    
    for vname in keys(vardict)
        
        varinfo = vardict[vname]
        index = varinfo.index
        gradient = varinfo.gradient
        
        grad[index] = gradient
        
    end
    
    return grad
    
end

#Eventually make this from file.
function get_value_vector(firststage::FirstStageInfo)
        
    vardict = firststage.variables
    nvars = vardict.count
    
    valvec= Vector{Float64}(undef, nvars) 
    
    for vname in keys(vardict)
        
        varinfo = vardict[vname]
        index = varinfo.index
        value = varinfo.value
        
        valvec[index] = value
        
    end
    
    return valvec
    
end

function get_objective_value(firststage::FirstStageInfo)
    
    objval = 0.0;

    subdict = firststage.subproblems

    for sid in keys(subdict)

        subobjval = subdict[sid].objective_value
        vardict = subdict[sid].variableinfo
        prob = subdict[sid].probability

        objval += prob*subobjval

    end
    
    return objval
    
end

#TODO make a dict so we know which constraints to adjust.
function adjust_h(firststage::FirstStageInfo, contoidx::Dict{Any,Int64}, h::Matrix{Float64})
        
    models = firststage.subproblems
    
    ns = models.count
    contoidx.count
        
    for sid in keys(models)
        m = models[sid].model
        
        for (F,S) in list_of_constraint_types(m)
            innercount = 0
            for con in all_constraints(m,F,S)
                innercount += 1
               if occursin("AffExpr",string(F))
                    idx = con.index.value
                    if occursin("EqualTo", string(S))
                    elseif occursin("GreaterThan", string(S))
                    elseif occursin("LessThan", string(S))
                    elseif occursin("Interval", string(S))
                        if dual(con) > 0
                            h[sid, contoidx[(F,S,innercount)]] = constraint_object(con).set.lower
                        else
                            h[sid, contoidx[(F,S,innercount)]] = constraint_object(con).set.upper
                        end
                    else
                        println(con)
                        println("Add ", S, " to hvars.")
                    end
                end
            end
        end
    end
    
    return h
    
end

#TODO since this only applies to interval constraints, classify constraints in subproblem by type, so it
# does not have to parse through every constraint.
function adjust_h_new!(subproblem::Union{Subproblems,SubproblemsNew})

    idxtocon = subproblem.idxtocon
    h = subproblem.h
    
    for idx in keys(idxtocon)
        
        con = idxtocon[idx]
        
        #this is a placeholder
        ctype = string(typeof(con))
        
        #do nothing if not an interval constraint. If something weird is happening say to add constraint type.
        if occursin("EqualTo", ctype)
        elseif occursin("GreaterThan", ctype)
        elseif occursin("LessThan", ctype)
        # my suspition right now is that the below is solver dependent. Works fine in Gurobi, but issues elsewhere.
        elseif occursin("Interval", ctype)
            if dual(con) > 0
                h[idx] = constraint_object(con).set.lower
            else
                h[idx] = constraint_object(con).set.upper
            end
        else
            println(con)
            println("Add ", ctype, " to hvars.")
        end
    end
    
    subproblem.h = h
    
    return
end

function compute_e(firststage::FirstStageInfo, h::Matrix{Float64}, PI::Matrix{Float64})
        
    e_k = 0.0;
    for i in keys(firststage.subproblems)
        prob = firststage.subproblems[i].probability
        e_k += prob*dot(PI[i,:],h[i,:])
    end
    return e_k
end
    

function compute_PI(firststage::FirstStageInfo, contoidx::Dict{Any,Int64})
        
    N = firststage.subproblems.count
    nc = contoidx.count
    PI = Array{Float64}(undef,N, nc)
    for i = 1:N
        subproblem = firststage.subproblems[i]

        m = subproblem.model

        count = 0
        for (F,S) in list_of_constraint_types(m)
            innercount = 0
            for con in all_constraints(m,F,S)
                innercount+=1
               if occursin("AffExpr",string(F))
                    idx = con.index.value
                    dual = JuMP.dual(con)
                    PI[i, contoidx[(F,S,innercount)]] = dual
                count += 1
               end
            end
        end
    end
    return PI
end
;

function add_theta_to_objective!(fs::JuMP.Model)
        
    @variable(fs, theta)
    obj_old = objective_function(fs)
    @objective(fs, Min, obj_old + theta)
    
    return fs
    
end

function add_theta_to_objective!(fs::JuMP.Model, nscen::Int64)
        
    for sid = 1:nscen
        var = @variable(fs, base_name = "theta_$(sid)")
        obj_old = objective_function(fs)
        @objective(fs, Min, obj_old + var)
    end
    
    return fs
    
end
    

function add_constraint_to_objective!(fs::JuMP.Model, E::Array{Float64}, e_k::Float64, v_dict::Dict{Int64,Array{Any}})
        
    nvar = length(E)
    
    @constraint(fs, sum(E[i]*variable_by_name(fs, v_dict[1][i][1]) for i = 1:nvar) 
                        + variable_by_name(fs, "theta") >= e_k)
        
    return fs
    
end

function add_constraint_to_objective!(fs::JuMP.Model, E::Array{Float64}, e_k::Float64, v_dict::Dict{Int64,Array{Any}}, sid::Int64)
        
    nvar = length(E)
    
    @constraint(fs, sum(E[i]*variable_by_name(fs, v_dict[1][i][1]) for i = 1:nvar) 
                        + variable_by_name(fs, "theta_$(sid)") >= e_k)
        
    return fs
    
end


function update_first_value_L!(firststage::FirstStageInfo, fs::JuMP.Model)
        
    for vname in keys(firststage.variables)
        
        varinfo = firststage.variables[vname]
        name = varinfo.name
        
        varinfo.value = JuMP.value(JuMP.variable_by_name(fs, name))
    end
    
    return
    
end

#to do, update so that w and theta are completely separate.
function setup_1st_paths!(firststage::FirstStageInfo)
    
    path = firststage.store
    vardict = firststage.variables
    nvars = vardict.count
    
    header = Array{Union{Nothing, String}}(nothing, nvars)
    
    for vname in keys(vardict)
        
        varinfo = vardict[vname]
        index = varinfo.index
        
        header[index] = vname
        
    end
    
    nvars = length(header)
    
    df = DataFrame()
    dfe = DataFrame(econst = Float64[])
    dft = DataFrame(theta = Float64[])
    dfwt = DataFrame(w = Float64[], theta = Float64[])
    
    for i = 1:nvars
        insertcols!(df, i, Symbol(header[i])=>Float64[])
    end
    
    mkdir(path)
    
    xcsv = string(path, "x.csv")
    Ecsv = string(path, "E.csv")
    ecsv = string(path, "ek.csv")
    #tcsv used to store just the theta value in async. Not used in serial/synchronous coding
    tcsv = string(path, "theta.csv")
    wtcsv = string(path, "w_theta.csv")
    
    CSV.write(xcsv, df)
    CSV.write(Ecsv, df)
    CSV.write(ecsv, dfe)
    CSV.write(tcsv, dft)
    CSV.write(wtcsv, dfwt)
    
    return
end

function setup_ro_paths!(firststage::FirstStageInfo, nscen::Int64)
    
    path = firststage.store
    vardict = firststage.variables
    nvars = vardict.count
    
    header = Array{Union{Nothing, String}}(nothing, nvars)
    
    for vname in keys(vardict)
        
        varinfo = vardict[vname]
        index = varinfo.index
        
        header[index] = vname
        
    end
    
    nvars = length(header)
    
    dfx = DataFrame()
    dflb = DataFrame(lower_bound = Float64[])
    dfub = DataFrame(upper_bound = Float64[])
    dfscen = DataFrame(scenario = Int64[])
    dfobjval = DataFrame(objval = Float64[])
    
    for i = 1:nvars
        insertcols!(dfx, i, Symbol(header[i])=>Float64[])
    end
    
    mkdir(path)
    
    xcsv = string(path, "x.csv")
    lbcsv = string(path, "lb.csv")
    ubcsv = string(path, "ub.csv")
    scencsv = string(path, "scenarios.csv")
    
    CSV.write(xcsv, dfx)
    CSV.write(lbcsv, dflb)
    CSV.write(ubcsv, dfub)
    CSV.write(scencsv, dfscen)
    
    for i = 1:nscen
        scenpath = string(path, "scen_$(i)/")
        mkdir(scenpath)
        
        objvalcsv = string(scenpath, "objval.csv")
        
        CSV.write(objvalcsv, dfobjval)
        
    end
    
    return
end

function setup_1st_paths_multicut!(firststage::FirstStageInfo, nscen::Int64)
        
    path = firststage.store
    vardict = firststage.variables
    nvars = vardict.count
    
    header = Array{Union{Nothing, String}}(nothing, nvars)
    
    for vname in keys(vardict)
        
        varinfo = vardict[vname]
        index = varinfo.index
        
        header[index] = vname
        
    end
    
    nvars = length(header)
    
    dfx = DataFrame()
    dft = DataFrame()
    
    for i = 1:nvars
        insertcols!(dfx, i, Symbol(header[i])=>Float64[])
    end
    
    for i = 1:nscen
        insertcols!(dft, i, Symbol(i)=>Float64[])
    end
    
    mkdir(path)
    
    xcsv = string(path, "x.csv")
    tcsv = string(path, "theta.csv")
    
    CSV.write(xcsv, dfx)
    CSV.write(tcsv, dft)
    
    return
end

function setup_1st_paths_regdec!(firststage::FirstStageInfo, nscen::Int64)
        
    path = firststage.store
    vardict = firststage.variables
    nvars = vardict.count
    
    header = Array{Union{Nothing, String}}(nothing, nvars)
    
    for vname in keys(vardict)
        
        varinfo = vardict[vname]
        index = varinfo.index
        
        header[index] = vname
        
    end
    
    nvars = length(header)
    
    dfx = DataFrame()
    dft = DataFrame()
    dfa = DataFrame()
    dfo = DataFrame(objval = Float64[])
    dffs = DataFrame(fsobj = Float64[])
    dfr = DataFrame(rho = Float64[])
    
    for i = 1:nvars
        insertcols!(dfx, i, Symbol(header[i])=>Float64[])
        insertcols!(dfa, i, Symbol(header[i])=>Float64[])
    end
    
    for i = 1:nscen
        insertcols!(dft, i, Symbol(i)=>Float64[])
    end
    
    mkdir(path)
    
    xcsv = string(path, "x.csv")
    tcsv = string(path, "theta.csv")
    acsv = string(path, "a.csv")
    ocsv = string(path, "ssobja.csv")
    fscsv = string(path, "fsobj.csv")
    rcsv = string(path, "rho.csv")
    
    CSV.write(xcsv, dfx)
    CSV.write(tcsv, dft)
    CSV.write(acsv, dfa)
    CSV.write(ocsv, dfo)
    CSV.write(fscsv, dffs)
    CSV.write(rcsv, dfr)
    
    return
end

function setup_1st_paths_trust!(firststage::FirstStageInfo, nscen::Int64)
        
    path = firststage.store
    vardict = firststage.variables
    nvars = vardict.count
    
    header = Array{Union{Nothing, String}}(nothing, nvars)
    
    for vname in keys(vardict)
        
        varinfo = vardict[vname]
        index = varinfo.index
        
        header[index] = vname
        
    end
    
    nvars = length(header)
    
    dfs = DataFrame(ssobja = Float64[], fsobj = Float64[], ssobjx = Float64[], concrit = Float64[], delta = Float64[], delcrit = Float64[]);
    dfx = DataFrame()
    dft = DataFrame()
    dfa = DataFrame()
    dfo = DataFrame(objval = Float64[])
    dffs = DataFrame(fsobj = Float64[])
    dfd = DataFrame(delta = Float64[])
    
    for i = 1:nvars
        insertcols!(dfx, i, Symbol(header[i])=>Float64[])
        insertcols!(dfa, i, Symbol(header[i])=>Float64[])
    end
    
    for i = 1:nscen
        insertcols!(dft, i, Symbol(i)=>Float64[])
    end
    
    mkdir(path)
    
    scsv = string(path, "iteration_summary.csv")
    xcsv = string(path, "x.csv")
    tcsv = string(path, "theta.csv")
    acsv = string(path, "a.csv")
    ocsv = string(path, "ssobja.csv")
    fscsv = string(path, "fsobj.csv")
    dcsv = string(path, "delta.csv")
    
    CSV.write(scsv, dfs)
    CSV.write(xcsv, dfx)
    CSV.write(tcsv, dft)
    CSV.write(acsv, dfa)
    CSV.write(ocsv, dfo)
    CSV.write(fscsv, dffs)
    CSV.write(dcsv, dfd)
    
    return
end

function store_x!(firststage::FirstStageInfo)
    
    path = firststage.store
    xcsv = string(path, "x.csv")
    
    x = get_value_vector(firststage)
    
    sx = string(x)
    n = length(string(x))

    sx = sx[2:n-1]
    open(xcsv, "a") do io
        write(io, "$(sx) \n")
    end
    
    return
    
end

function store_LB!(LB::Float64, path::String)
   
    lbcsv = string(path, "lb.csv")
    
    open(lbcsv, "a") do io
        write(io, "$(LB) \n")
    end
    
    return
end

function store_UB!(UB::Float64, path::String)
   
    ubcsv = string(path, "ub.csv")
    
    open(ubcsv, "a") do io
        write(io, "$(UB) \n")
    end
    
    return
end

function store_scen!(scen::Int64, path::String)
   
    scsv = string(path, "scenarios.csv") ## doublecheck name
    
    open(scsv, "a") do io
        write(io, "$(scen) \n")
    end
    
    return
end

function store_objval!(objval::Float64, sid::Int64, path::String)
   
    scenpath = string(path, "scen_$(sid)/")
    ovcsv = string(scenpath, "objval.csv")
    
    open(ovcsv, "a") do io
        write(io, "$(objval) \n")
    end
    
    return
end
    

function store_a!(a, path::String)
    
    acsv = string(path, "a.csv")
    
    
    sa = string(a)
    n = length(string(a))

    sa = sa[2:n-1]
    open(acsv, "a") do io
        write(io, "$(sa) \n")
    end
    
    return
    
end

function store_iteration_summary!(ssobja, fsobj, ssobjx, concrit, delta, delcrit, path::String)
    
    scsv = string(path, "iteration_summary.csv")
    
    str = string(ssobja, ", ", fsobj, ", ", ssobjx, ", ",concrit, ", ", delta, ", ", delcrit)

    open(scsv, "a") do io
        write(io, "$(str) \n")
    end
    
    return
    
end


function store_ssobja!(ssobja, path::String)
    
    ocsv = string(path, "ssobja.csv")
    
    open(ocsv, "a") do io
        write(io, "$(ssobja) \n")
    end
    
    return
end

function store_fsobj!(fsobj, path::String)
    
    ocsv = string(path, "fsobj.csv")
    
    open(ocsv, "a") do io
        write(io, "$(fsobj) \n")
    end
    
    return
end

function store_rho!(rho, path::String)
    
    rcsv = string(path, "rho.csv")
    
    open(rcsv, "a") do io
        write(io, "$(rho) \n")
    end
    
    return
end

function store_delta!(delta, path::String)
    
    dcsv = string(path, "delta.csv")
    
    open(dcsv, "a") do io
        write(io, "$(delta) \n")
    end
    
    return
end

function store_thetak!(firststage::FirstStageInfo, thetavec)
    
    path = firststage.store
    tcsv = string(path, "theta.csv")
        
    st = string(thetavec)
    n = length(string(thetavec))

    st = st[5:n-1]
    open(tcsv, "a") do io
    write(io, "$(st) \n")
    end
    
    return
end

function store_E!(firststage::FirstStageInfo, E::Array{Float64})
    
    path = firststage.store
    Ecsv = string(path, "E.csv")
        
    sE = string(E)
    n = length(string(E))

    sE = sE[2:n-1]
    open(Ecsv, "a") do io
        write(io, "$(sE) \n")
    end
    
    return
    
end

function store_e!(firststage::FirstStageInfo, e_k::Float64)
    
    path = firststage.store
    ecsv = string(path, "ek.csv")
        
    open(ecsv, "a") do io
        write(io, "$(e_k) \n")
    end
    
    return
    
end

function store_w_theta!(firststage::FirstStageInfo, w::Float64, theta::Float64)
    
    path = firststage.store
    wtcsv = string(path, "w_theta.csv")
    
    open(wtcsv, "a") do io
        write(io, "$(w), $(theta) \n")
    end
    
    return
    
end

function store_theta!(firststage::FirstStageInfo, theta::Float64)
    
    path = firststage.store
    tcsv = string(path, "theta.csv")
    
    open(tcsv, "a") do io
        write(io, "$(theta) \n")
    end
    
    return
    
end

function store_w!(firststage::FirstStageInfo, w::Float64)
    
    path = firststage.store
    wcsv = string(path, "w.csv")
    
    open(wcsv, "a") do io
        write(io, "$(w) \n")
    end
    
    return
    
end

function setup_2nd_paths!(path::String, subproblem::Union{Subproblems,SubproblemsNew})
    
    vardict = subproblem.variableinfo
    sid = subproblem.id
    nvars = vardict.count
    
    header = Array{Union{Nothing, String}}(nothing, nvars)
    
    for vid in keys(vardict)
        
        vname = vardict[vid].name
        index = vardict[vid].index
        
        header[index] = vname
        
    end
    
    nvars = length(header)
    
    df = DataFrame()
    dfe = DataFrame(econst = Float64[])
    
    for i = 1:nvars
        insertcols!(df, i, Symbol(header[i])=>Float64[])
    end
        
    Ecsv = string(path, "scen_$(sid)/E.csv")
    ekcsv = string(path, "scen_$(sid)/ek.csv")
    
    CSV.write(Ecsv, df)
    CSV.write(ekcsv, dfe)   
    
    return
end

function setup_2nd_paths!(path::String, subproblem::Union{Subproblems,SubproblemsNew}, multicut)
    
    vardict = subproblem.variableinfo
    sid = subproblem.id
    nvars = vardict.count
    
    header = Array{Union{Nothing, String}}(nothing, nvars)
    
    for vid in keys(vardict)
        
        vname = vardict[vid].name
        index = vardict[vid].index
        
        header[index] = vname
        
    end
    
    nvars = length(header)
    
    df = DataFrame()
    dfe = DataFrame(econst = Float64[])
    dfc = DataFrame(addcut = Int64[])
    
    for i = 1:nvars
        insertcols!(df, i, Symbol(header[i])=>Float64[])
    end
        
    Ecsv = string(path, "scen_$(sid)/E.csv")
    ekcsv = string(path, "scen_$(sid)/ek.csv")
    cutcsv = string(path, "scen_$(sid)/addcut.csv")
    
    CSV.write(Ecsv, df)
    CSV.write(ekcsv, dfe)   
    CSV.write(cutcsv, dfc)
    
    return
end

function setup_2nd_paths_regdec!(path::String, subproblem::Union{Subproblems,SubproblemsNew})
    
    vardict = subproblem.variableinfo
    sid = subproblem.id
    nvars = vardict.count
    
    header = Array{Union{Nothing, String}}(nothing, nvars)
    
    for vid in keys(vardict)
        
        vname = vardict[vid].name
        index = vardict[vid].index
        
        header[index] = vname
        
    end
    
    nvars = length(header)
    
    df = DataFrame()
    dfe = DataFrame(econst = Float64[])
    dfc = DataFrame(addcut = Int64[])
    dfo = DataFrame(objval = Float64[])
    
    for i = 1:nvars
        insertcols!(df, i, Symbol(header[i])=>Float64[])
    end
        
    Ecsv = string(path, "scen_$(sid)/E.csv")
    ekcsv = string(path, "scen_$(sid)/ek.csv")
    cutcsv = string(path, "scen_$(sid)/addcut.csv")
    objcsv = string(path, "scen_$(sid)/ssobjx.csv")
    
    CSV.write(Ecsv, df)
    CSV.write(ekcsv, dfe)   
    CSV.write(cutcsv, dfc)
    CSV.write(objcsv, dfo)
    
    return
end

function setup_scen_path!(path::String, sid::Int64)
    
    scenpath = string(path, "scen_$(sid)/")
    
    mkdir(scenpath)
    
    return 
    
end

function store_Es_es!(path::String, subproblem::Union{Subproblems,SubproblemsNew}, pi_k::Array{Float64}, h_k::Array{Float64})
           
    prob = subproblem.probability
    vardict = subproblem.variableinfo
    sid = subproblem.id
    nvars = vardict.count
    
    # initialize E vector
    Evec = Vector{Float64}(undef, nvars) 
    
    for var in keys(vardict)
        vind = vardict[var].index
        grad2 = vardict[var].gradient
        
        # this is using the NAC method
        Evec[vind] = -prob*grad2
        
    end
    
    #TODO have this not based on firststage stuff
    e_k = prob*dot(pi_k,h_k)
    
    Ecsv = string(path, "scen_$(sid)/E.csv")
    ekcsv = string(path, "scen_$(sid)/ek.csv")
    
    sE = string(Evec)
    n = length(string(Evec))

    sE = sE[2:n-1]
    open(Ecsv, "a") do io
        write(io, "$(sE) \n")
    end 
            
    open(ekcsv, "a") do io
        write(io, "$(e_k) \n")
    end
    
    
    return
    
end


function iterate_L(firststage::FirstStageInfo, fs::JuMP.Model, contoidx::Dict, h::Matrix{Float64}, v_dict::Dict{Int64,Array{Any}}, addtheta::Int64, tol::Float64, niter::Int64, verbose::Int64)
        
    x = 0
    
    cost = get_cost_vector(firststage, fs)
    if verbose == 1
        println("cost = $(cost)")
    end
    println("cost = $(cost)")
    
    for i = 1:niter

        # step 1 set v = v+1 and solve first stage problem.
        println("Iteration $(i)")

        optimize!(fs)

        update_first_value_L!(firststage, fs)
        
        if firststage.store != nothing
        
            if i == 1
                setup_1st_paths!(firststage)
            end
            
            store_x!(firststage)
            
        end

        # step 3 * update x-variables in second stage
        update_second_value!(firststage)

        #        * solve second stage problems
        # done separately to eventually parallelize 
        
        for sid in keys(firststage.subproblems)  
            solve_sub_and_update!(firststage.subproblems[sid])
            
            if firststage.store!= nothing
                
                #I will have to store this locally at some point...
                path = firststage.store
                
                if i == 1
                    setup_scen_path!(path, sid)
                    
                    setup_2nd_paths!(path, firststage.subproblems[sid])
                end
                                                   
            end
        end

        #        * get simplex multipliers, update E and e
        # update E
        update_first_gradient!(firststage)
        
        grad = get_grad_vector(firststage)
        if verbose == 1
            println("grad = $(grad)")
        end
                
        h = adjust_h(firststage, contoidx, h)

        E = cost - grad
        if firststage.store != nothing
            store_E!(firststage, E)
        end
        if verbose == 1
            println("E = $(E)")
        end

        #update PI
        PI = compute_PI(firststage, contoidx)
                
        if firststage.store != nothing
            for sid in keys(firststage.subproblems)
                store_Es_es!(firststage.store, firststage.subproblems[sid], PI[sid,:], h[sid,:])
            end
        end

        #update e_k
        e_k = compute_e(firststage,h,PI)
        if verbose == 1
            println("e = $(e_k)")
        end
        
        if firststage.store != nothing
            store_e!(firststage, e_k)
        end


        x = get_value_vector(firststage)
        if verbose == 1
            println("x = $(x)")
        end

        if addtheta == 1
            theta = JuMP.value(JuMP.variable_by_name(fs, "theta"))
            println("theta = $(theta)")
        end

        w = e_k - dot(E,x)
        
        if firststage.store != nothing
            if addtheta == 0
                store_w_theta!(firststage, w, -Inf)
            elseif addtheta == 1
                store_w_theta!(firststage, w, theta)
            end
        end
        
        println("w = $(w)")
        if addtheta == 1
            if theta >= w - tol
                println("algorithm converged.")
                println("final x = $(x)")
                return x, firststage, fs, contoidx, h, i
            end
        end
        
        if addtheta == 0
            fs = add_theta_to_objective!(fs)
            addtheta = 1
        end

        fs = add_constraint_to_objective!(fs, E, e_k, v_dict)

    end
    
    println("L-Shaped Algorithm Failed to converge in $(niter) iterations.")
    
    return x, firststage, fs, contoidx, h, niter
    
end

function compute_Ek_new!(subproblem::Union{Subproblems,SubproblemsNew})
    
    prob = subproblem.probability
    vardict = subproblem.variableinfo
    nvars = vardict.count
    
    # initialize E vector
    Evec = Vector{Float64}(undef, nvars) 
    
    for var in keys(vardict)
        vind = vardict[var].index
        grad2 = vardict[var].gradient
        
        # this is using the NAC method
        Evec[vind] = -prob*grad2
        
    end
    
    subproblem.Ek = Evec
    
    return
    
end

function compute_ek_new!(subproblem::Union{Subproblems,SubproblemsNew})
    
    prob = subproblem.probability
    idxtocon = subproblem.idxtocon
    h = subproblem.h
    
    ek = 0.0;
    
    for idx in keys(idxtocon)
        
        con = idxtocon[idx]
        dual = JuMP.dual(con)
        
        ek += h[idx]*dual
        
    end

    subproblem.ek = prob*ek
    
    return
    
end

function store_Ek_sub!(subproblem::Union{Subproblems,SubproblemsNew}, path::String)
    
    sid = subproblem.id
    Evec = subproblem.Ek
    
    Ecsv = string(path, "scen_$(sid)/E.csv")
    
    sE = string(Evec)
    n = length(string(Evec))

    sE = sE[2:n-1]
    open(Ecsv, "a") do io
        write(io, "$(sE) \n")
    end 
    
    return
    
end

function store_ek_sub!(subproblem::Union{Subproblems,SubproblemsNew}, path::String)
    
    sid = subproblem.id
    ek = subproblem.ek
    
    ecsv = string(path, "scen_$(sid)/ek.csv")
    
    open(ecsv, "a") do io
        write(io, "$(ek) \n")
    end
    
    return
    
end

function store_cut_sub!(addcut::Int64, path::String, sid::Int64)
    
    ccsv = string(path, "scen_$(sid)/addcut.csv")
    
    open(ccsv, "a") do io
        write(io, "$(addcut) \n")
    end
    
    return
    
end

function store_objval!(objval::Float64, path::String, sid::Int64)
    
    ocsv = string(path, "scen_$(sid)/ssobjx.csv")
    
    open(ocsv, "a") do io
        write(io, "$(objval) \n")
    end
    
    return
    
end
    
#as opposed to get_El_from_file to be made
function get_El_from_sub!(firststage::FirstStageInfo)
    
    subproblems = firststage.subproblems
    nvars = firststage.variables.count
    
    El = zeros(nvars)
    
    for sid in keys(subproblems)
        El += subproblems[sid].Ek
    end
    
    return El
    
end

#possibly inefficient. WIll fix if get to be too slow.
# put back in firststage once in algorithms.
function get_El_from_file_and_save!(firststage::FirstStageInfo, model::JuMP.Model, path::String, nsubs::Int64)
    
    vardict = firststage.variables
    nvars = vardict.count
    
    El = LShaped.get_cost_vector(firststage, model)
    
    for sid = 1:nsubs
        Epath = string(path, "scen_$(sid)/E.csv")
        Edf = DataFrame(CSV.File(Epath))
                
        Esub = collect(Edf[size(Edf,1),:])
        
        El += Esub
        
    end
    
    Ecsv = string(path, "E.csv")
        
    sE = string(El)
    n = length(string(El))

    sE = sE[2:n-1]
    open(Ecsv, "a") do io
        write(io, "$(sE) \n")
    end
    
    return El
    
end


function get_el_from_sub!(firststage::FirstStageInfo)
    
    subproblems = firststage.subproblems
    
    el = 0.0
    
    for sid in keys(subproblems)
        el += subproblems[sid].ek
    end
    
    return el
    
end

# path should be path that contains the scen_# directories
function get_el_from_file_and_save!(path::String, nsubs::Int64)
    
    el = 0.0
    
    for i = 1:nsubs
        epath = string(path, "scen_$(i)/ek.csv")
        
        #open as a dataframe, as you are using DataFrames package anyway. More efficient to load in last row, though.
        edf = DataFrame(CSV.File(epath))
        
        nrows = size(edf,1)
        
        el += edf[nrows, 1]
        
    end
    
    ecsv = string(path, "ek.csv")
        
    open(ecsv, "a") do io
        write(io, "$(el) \n")
    end
    
    return el
    
end

function load_current_fs!(model::JuMP.Model, E::DataFrame, ek::DataFrame, xvars::Vector{String}, curit::Int64)
    
    add_theta_to_objective!(model)
    
    nvars = length(xvars)
    
    for it = 1:curit
        
        El = collect(E[it,:])
        el = ek[it,1]
        
        @constraint(model, sum(El[i]*variable_by_name(model, xvars[i]) for i = 1:nvars) 
                        + variable_by_name(model, "theta") >= el)
    
    end
    
    return
    
end

function get_x_from_file(path::String)
    
    xpath = string(path, "x.csv")
    
    xdf = DataFrame(CSV.File(xpath))
    
    x = collect(xdf[size(xdf,1),:])
    
    return x
    
end

function get_rho_from_file(path::String)
    
    rpath = string(path, "rho.csv")
    
    rdf = DataFrame(CSV.File(rpath))
    
    rho = rdf[size(rdf,1),1]
    
    return rho
end

function get_delta_from_file(path::String)
    
    dpath = string(path, "delta.csv")
    
    ddf = DataFrame(CSV.File(dpath))
    
    delta = ddf[size(ddf,1),1]
    
    return delta
end

function get_theta_from_file(path::String)
    
    tpath = string(path, "theta.csv")
    
    tdf = DataFrame(CSV.File(tpath))
    
    theta = tdf[size(tdf,1),1]
    
    return theta
    
end

#TODO: add fs to firststage struct.
function resume_fs!(firststage::FirstStageInfo, model::JuMP.Model, tol::Float64)
    
    converged = 0
    #model = firststage.model
    path = firststage.store
    pathE = string(path, "E.csv")
    pathe = string(path, "ek.csv")
    pathx = string(path, "x.csv")
    pathwt = string(path, "w_theta.csv")
    
    E = DataFrame(CSV.File(pathE))
    ek = DataFrame(CSV.File(pathe))

    w_theta = DataFrame(CSV.File(pathwt))
    #perhaps temporary.
    x = DataFrame(CSV.File(pathx))
    
    #double check command to exit the program
    curit=0
    if size(ek,1) == size(E,1) 
        curit = size(ek,1)
    else 
        println("Error: Look into this dataset. The number of rows in x, Ek, and ek should match.")
        return 0, 0, collect(x[curit,:])
    end
    
    #check convergence
    if w_theta[curit,2] >= w_theta[curit,1] - tol
        converged = 1
        return curit, converged, collect(x[curit,:])
    end
    
    xvars = String.(names(E))
    
    load_current_fs!(model, E, ek, xvars, curit)
    
    return curit, converged, collect(x[curit,:])
end

function resume_fs_multicut!(firststage::FirstStageInfo, v_dict, model::JuMP.Model, nscen::Int64)
    
    converged = 0
    path = firststage.store
    ##pathE = string(path, "E.csv")
    ##pathe = string(path, "ek.csv")
    pathx = string(path, "x.csv")
    ##pathwt = string(path, "w_theta.csv")
    
    ##E = DataFrame(CSV.File(pathE))
    ##ek = DataFrame(CSV.File(pathe))

    ##w_theta = DataFrame(CSV.File(pathwt))
    x = DataFrame(CSV.File(pathx))
    
    #double check command to exit the program
    curit= size(x,1)

    #check convergence
    numcuts = 0
    
    add_theta_to_objective!(model, nscen)
    
    for sid = 1:nscen
        pathE = string(path, "scen_$(sid)/E.csv")
        pathe = string(path, "scen_$(sid)/ek.csv")
        pathc = string(path, "scen_$(sid)/addcut.csv")
        
        E = DataFrame(CSV.File(pathE))
        ek = DataFrame(CSV.File(pathe))
        cuts = DataFrame(CSV.File(pathc))
        
        #this determines if a cut was generated in the most recent iteration
        n = size(cuts,1)
        numcuts += cuts[n,1]
        
        for j = 1:n
            if cuts[j,1] == 1
                add_constraint_to_objective!(model, collect(E[j,:]), ek[j,1], v_dict, sid)
            end
        end        
    end
    
    if numcuts == 0
        converged = 1
        return model, curit, converged
    end
        
    
    return model, curit, converged
end

function resume_fs_regdec!(firststage::FirstStageInfo, v_dict, model::JuMP.Model, nscen::Int64, rho::Float64, header, rhomin, rhomax, gamma)
    
    converged = 0
    path = firststage.store
    pathx = string(path, "x.csv")
    patha = string(path, "a.csv")
    ssobjacsv = string(path, "ssobja.csv")
    fsobjcsv = string(path, "fsobj.csv")
    rhocsv = string(path, "rho.csv")

    xdf = DataFrame(CSV.File(pathx))
    adf = DataFrame(CSV.File(patha))
    fsobjdf = DataFrame(CSV.File(fsobjcsv))
    rhodf = DataFrame(CSV.File(rhocsv))
    
    na = size(adf,1)
    nx = size(xdf,1)
    nfs = size(fsobjdf,1)
    nrho = size(rhodf,1)
    
    curit= size(xdf,1)
    ssobjadf = DataFrame()
    ssobja = 0.0
    if curit > 1
        ssobjadf = DataFrame(CSV.File(ssobjacsv))
        ssobja = ssobjadf[size(ssobjadf,1),1]
        a = collect(adf[na,:])
    end
    
    x = collect(xdf[nx,:])
    fsobj = fsobjdf[nfs,1]
    
    #check convergence
    numcuts = 0
    
    add_theta_to_objective!(model, nscen)
    
    ssobjx = 0.0
    
    
    for sid = 1:nscen
        pathE = string(path, "scen_$(sid)/E.csv")
        pathe = string(path, "scen_$(sid)/ek.csv")
        pathc = string(path, "scen_$(sid)/addcut.csv")
        pathsobj = string(path, "scen_$(sid)/ssobjx.csv")
        
        E = DataFrame(CSV.File(pathE))
        ek = DataFrame(CSV.File(pathe))
        cuts = DataFrame(CSV.File(pathc))
        ssobjxdf = DataFrame(CSV.File(pathsobj))
        
        ssobjx += ssobjxdf[size(ssobjxdf,1),1]
        
        #this determines if a cut was generated in the most recent iteration
        n = size(cuts,1)
        numcuts += cuts[n,1]
        
        for j = 1:n
            if cuts[j,1] == 1
                add_constraint_to_objective!(model, collect(E[j,:]), ek[j,1], v_dict, sid)
            end
        end        
    end
    
    #initial ssobja = ssobjx since an initial a comes from first (linear) solve
    if curit == 1
        ssobja = ssobjx
        LShaped.store_ssobja!(ssobja, path)
    end
    
    print(curit, " ", ssobja, " ", ssobjx, " ", numcuts, " ")
    if numcuts == 0
        a = x #from last iter
        ssobja = ssobjx
        store_a!(a, path)
        store_ssobja!(ssobja, path)
    else
        if ssobjx <= ssobja
            a = x
            ssobja = ssobjx
            store_a!(a, path)
            store_ssobja!(ssobja, path)
        else
            store_a!(a, path)
            store_ssobja!(ssobja, path)
        end
    end
    
    ##### put alg stuff here #####
    if ssobja <= ssobjx - gamma*(ssobjx - fsobj)
        rho = max(rhomin, 0.5*rho)
        #if rho == rhomin
        #    rho = 0.0
        #end
    else
        #if rho > 0.0
            rho = min(rhomax, 2.0*rho)
        #else
        #    rho = rhomin
        #end
    end
    
    println(rho)

    ## load in model, x or a vector, and header to adjust the model
    model = add_regularized_decomp_async!(model, a, header, rho)
    
    return model, curit, ssobja, rho
end

function resume_fs_trust!(firststage::FirstStageInfo, v_dict, model::JuMP.Model, nscen::Int64, delta::Float64, header, deltamin, deltamax, gamma)
    
    converged = 0
    path = firststage.store
    pathx = string(path, "x.csv")
    patha = string(path, "a.csv")
    ssobjacsv = string(path, "ssobja.csv")
    fsobjcsv = string(path, "fsobj.csv")
    deltacsv = string(path, "delta.csv")

    xdf = DataFrame(CSV.File(pathx))
    adf = DataFrame(CSV.File(patha))
    fsobjdf = DataFrame(CSV.File(fsobjcsv))
    deltadf = DataFrame(CSV.File(deltacsv))
    
    na = size(adf,1)
    nx = size(xdf,1)
    nfs = size(fsobjdf,1)
    ndelta = size(deltadf,1)
    
    curit= size(xdf,1)
    ssobjadf = DataFrame()
    ssobja = 0.0
    if curit > 1
        ssobjadf = DataFrame(CSV.File(ssobjacsv))
        ssobja = ssobjadf[size(ssobjadf,1),1]
        a = collect(adf[na,:])
    end
    
    x = collect(xdf[nx,:])
    fsobj = fsobjdf[nfs,1]
    
    #check convergence
    numcuts = 0
    
    add_theta_to_objective!(model, nscen)
    
    ssobjx = 0.0
    
    for sid = 1:nscen
        pathE = string(path, "scen_$(sid)/E.csv")
        pathe = string(path, "scen_$(sid)/ek.csv")
        pathc = string(path, "scen_$(sid)/addcut.csv")
        pathsobj = string(path, "scen_$(sid)/ssobjx.csv")
        
        E = DataFrame(CSV.File(pathE))
        ek = DataFrame(CSV.File(pathe))
        cuts = DataFrame(CSV.File(pathc))
        ssobjxdf = DataFrame(CSV.File(pathsobj))
        
        ssobjx += ssobjxdf[size(ssobjxdf,1),1]
        
        #this determines if a cut was generated in the most recent iteration
        n = size(cuts,1)
        numcuts += cuts[n,1]
        
        for j = 1:n
            if cuts[j,1] == 1
                add_constraint_to_objective!(model, collect(E[j,:]), ek[j,1], v_dict, sid)
            end
        end        
    end
    
    #initial ssobja = ssobjx since an initial a comes from first (linear) solve
    if curit == 1
        ssobja = ssobjx
        a = x
        LShaped.store_ssobja!(ssobja, path)
        LShaped.store_a!(a, path)
    end
       
    ##### put alg stuff here #####
    if (ssobjx - fsobj) != 0
        println(-(ssobja-ssobjx)/(ssobjx-fsobj))
    end
    delcrit = -(ssobja-ssobjx)/(1.0e-10+(ssobjx-fsobj))
    if ssobja < ssobjx - 0.15*(ssobjx - fsobj)
        delta = min(deltamax, 2.0*delta)
        #if rho == rhomin
        #    rho = 0.0
        #end
    elseif ssobja > ssobjx - 0.05*(ssobjx - fsobj)
        #if rho > 0.0
            delta = max(deltamin, 0.25*delta)
        #else
        #    rho = rhomin
        #end
    end
    
    #add in lower bound to problem (ssobja is upper bound)
    JuMP.optimize!(model)
    lb = JuMP.objective_value(model)
    
    print(curit, " ", ssobja, " ", ssobjx, " ", numcuts, " ")
    if numcuts == 0
        a = x #from last iter
        ssobja = ssobjx
        store_a!(a, path)
        store_ssobja!(ssobja, path)
    else
        if ssobjx < ssobja
            a = x
            ssobja = ssobjx
            store_a!(a, path)
            store_ssobja!(ssobja, path)
        else
            store_a!(a, path)
            store_ssobja!(ssobja, path)
        end
    end
    concrit = (fsobj-ssobja)/(1.0e-10+abs(ssobja))
    println(delta)

    ## load in model, x or a vector, and header to adjust the model
    model = add_trust_async!(model, a, header, delta)
    
    store_iteration_summary!(ssobja, fsobj, ssobjx, concrit, delta, delcrit, path)
    
    return model, curit, ssobja, delta, lb
end


function add_regularized_decomp_async!(model, a, varnames, rho)
    
    n = length(a)
    
    linobj = JuMP.objective_function(model)
    
    newfunc = JuMP.GenericQuadExpr{Float64,VariableRef}()
    
    for i = 1:n
        
        var = JuMP.variable_by_name(model, varnames[i])
               
        JuMP.add_to_expression!(newfunc, rho/2*(var - a[i])^2 )
    end
    
    JuMP.@objective(model, Min, linobj + newfunc)
        
    return model
end

function add_trust_async!(model, a, varnames, delta)
    
    n = length(a)
            
    for i = 1:n
        
        var = JuMP.variable_by_name(model, varnames[i])
        
        JuMP.@constraint(model, var-a[i] >= -delta, base_name = "trlb_$(var)")
        JuMP.@constraint(model, var-a[i] <= delta, base_name = "trub_$(var)")
               
    end
            
    return model
end

function save_cur_duals!(path::String, subproblem::SubproblemsNew, curit::Int64)
    
    idxtocon = subproblem.idxtocon
    
    #the indexes should start from 1 and go to end
    count = idxtocon.count
    
    connums = [];
    dualvals = [];
    
    for i = 1:count
        
        push!(connums, i)
        
        con = idxtocon[i]
        push!(dualvals, JuMP.dual(con))
        
    end
        
    ddf = DataFrame(conidx=connums, dualval=dualvals)
        
    conpath = string(path, "con_$(curit).csv")
    CSV.write(conpath, ddf)
    
    return
end

# stopped here 11/19 3:21 PM
# path should be path ..scen_(sid)/ for second stage problems
# path should be ..data/ for for firststage problems
function save_cur_vars!(path::String, model::JuMP.Model, curit::Int64)
        
    vars = JuMP.all_variables(model)
    vals = Vector{Float64}(undef, size(vars,1))
    
    for i = 1:size(vars,1)
        vals[i] = JuMP.value(vars[i])
    end
    
    varpath = string(path, "var_$(curit).csv")
    
    vdf = DataFrame(variables = vars, values = vals)
    
    CSV.write(varpath, vdf)
    
    return
end

function warmstart(subproblem::SubproblemsNew, curit::Int64, path::String)
    
    model = subproblem.model
    idxtocon = subproblem.idxtocon
    
    varcsvname = string(path, "var_$(curit).csv")
    concsvname = string(path, "con_$(curit).csv")

    vdf = DataFrame(CSV.File(varcsvname))
    #cdf = DataFrame(CSV.File(concsvname))
    
    for i = 1:size(vdf,1)
        
        #when saving to file, type is String31 for some reason.
        vname = ""
        for s = 1:length(vdf[i,1])
            vname = string(vname,vdf[i,1][s])
        end
        
        vref = JuMP.variable_by_name(model, vname)
        JuMP.set_start_value(vref, vdf[i,2])
        
    end
    #=
    for i = 1:size(cdf,1)
        
        cref = idxtocon[cdf[i,1]]
        JuMP.set_dual_start_value(cref, cdf[i,2])
        
    end
    =#
    
    return
    
end

#to do: determine the first stage values and subtract off 1st stage part before inputting a lower bound.
function add_lower_bound_to_first_stage!(model::JuMP.Model, lowerbound::Float64)
   
    obj_func = JuMP.objective_function(model)
    
    JuMP.@constraint(model, obj_func >= lowerbound)
    
    return
        
end

function iterate_L_new(firststage::FirstStageInfo, fs::JuMP.Model, v_dict::Dict{Int64,Array{Any}}, addtheta::Int64, tol::Float64, niter::Int64, verbose::Int64, resume::Int64, lowerbound::Union{Float64,Nothing})
        
    x = 0
    
    cost = get_cost_vector(firststage, fs)
    if verbose == 1
        println("cost = $(cost)")
    end
    
    curit = 0
    if resume > 0
        println("!!!!!!!!!!!!!! Resuming from iteration $(resume) using path $(firststage.store) !!!!!!!!!!!!!!")
        addtheta = 1
        #rebuild based on values of E and e
        curit, converged, xcur = resume_fs!(firststage, fs, tol)
        
        if converged > 0
            println("Model has already converged.")
            
            # 0 is placeholder for a better x or just having firststage hold everything
            return xcur, firststage, fs, curit
        end
    end
    if typeof(lowerbound) != Nothing
        
        fs = add_theta_to_objective!(fs)
        addtheta = 1
                
        #todo create this function
        add_lower_bound_to_first_stage!(fs, lowerbound)
    end
    
    for i = 1:niter
        
        itnum = i+curit
        itmax = curit+niter

        # step 1 set v = v+1 and solve first stage problem.
        println("Iteration $(itnum)/$(itmax)")

        optimize!(fs)

        if verbose == 1
            println("Updating first value...")
        end
        update_first_value_L!(firststage, fs)
        
        if firststage.store != nothing
        
            if itnum == 1
                if verbose == 1
                    println("Setting up First Stage paths...")
                end
                setup_1st_paths!(firststage)
            end
            
            if verbose == 1
                println("Storing x...")
            end
            store_x!(firststage)
            
        end

        # step 3 * update x-variables in second stage
        if verbose == 1
            println("Updating second values...")
        end
        update_second_value!(firststage)

        #        * solve second stage problems
        # done separately to eventually parallelize 
        
        for sid in keys(firststage.subproblems)  
            println("Solving subproblems and updating...")
            solve_sub_and_update!(firststage.subproblems[sid])
            
            if firststage.store!= nothing
                
                #I will have to store this locally at some point...
                path = firststage.store
                
                if itnum == 1
                    println("Setting up second stage paths...")
                    setup_scen_path!(path, sid)
                    
                    setup_2nd_paths!(path, firststage.subproblems[sid])
                end
                                                   
            end
        end

        #        * get simplex multipliers, update E and e
        # update E
        println("Updating first stage gradients...")
        update_first_gradient!(firststage)

        grad = get_grad_vector(firststage)
        if verbose == 1
            println("grad = $(grad)")
        end
        
        for sid in keys(firststage.subproblems)  
            println("For subproblem $(sid)..")
            subproblem = firststage.subproblems[sid]
            println("...adjusting h...")
            #temporarily removing since there should not be interval constraints
            #adjust_h_new!(subproblem) 
            #see current Ek_ek folder)
            println("...computing Ek...")
            compute_Ek_new!(subproblem) 
            println("...computing ek...")
            compute_ek_new!(subproblem) 
            
            if firststage.store != nothing
                println("Storing Ek and ek...")
                store_Ek_sub!(subproblem, firststage.store) 
                store_ek_sub!(subproblem, firststage.store) 
                
                println("...saving primal variables...")
                sid = subproblem.id
                path_id = string(firststage.store, "scen_$(sid)/")
                save_cur_vars!(path_id, subproblem.model, itnum)
                
                println("...saving dual variables...")
                save_cur_duals!(path_id, subproblem, itnum)
            end
        end
                    
        #get it from subproblems. In async make get_Ek_from_file function
        println("Updating First stage stuff...")
        #TODO change this to gradient, as two-stage problems have a certain first stage cost ONLY
        El = get_El_from_sub!(firststage) #done
        El = cost + El
        if verbose == 1
            println("El = $(El)")
        end
        el = get_el_from_sub!(firststage) #done
        if verbose == 1
            println("el = $(el)")
        end
        
        x = get_value_vector(firststage)
        if verbose == 1
            println("x = $(x)")
        end
    
        if addtheta == 1
            theta = JuMP.value(JuMP.variable_by_name(fs, "theta"))
            println("theta = $(theta)")
        end

        w = el - dot(El,x)
        
        if firststage.store != nothing
            store_E!(firststage, El)
            store_e!(firststage, el)
            if addtheta == 0
                store_w_theta!(firststage, w, -Inf)
            elseif addtheta == 1
                store_w_theta!(firststage, w, theta)
            end
        end
        
        println("w = $(w)")
        if addtheta == 1
            if theta >= w - tol
                println("algorithm converged.")
                println("final x = $(x)")
                return x, firststage, fs, i
            end
        end
        
        if addtheta == 0
            fs = add_theta_to_objective!(fs)
            addtheta = 1
        end

        fs = add_constraint_to_objective!(fs, El, el, v_dict)

    end
    
    println("L-Shaped Algorithm Failed to converge in $(niter) iterations.")
    
    return x, firststage, fs, niter
    
end

function iterate_L_multicut(firststage::FirstStageInfo, fs::JuMP.Model, v_dict::Dict{Int64,Array{Any}}, addtheta::Int64, tol::Float64, niter::Int64, verbose::Int64, resume::Int64, lowerbound::Union{Float64,Nothing})
        
    x = 0
    
    addthetak = zeros(Int64, firststage.subproblems.count)
    
    cost = get_cost_vector(firststage, fs)
    if verbose == 1
        println("cost = $(cost)")
    end
    
    curit = 0
    if resume > 0
        println("!!!!!!!!!!!!!! Resuming from iteration $(resume) using path $(firststage.store) !!!!!!!!!!!!!!")
        addtheta = 1
        #rebuild based on values of E and e
        ### NEEDS TO CHANGE
        curit, converged, xcur = resume_fs!(firststage, fs, tol)
        ###
        
        if converged > 0
            println("Model has already converged.")
            
            # 0 is placeholder for a better x or just having firststage hold everything
            return xcur, firststage, fs, curit
        end
    end
    if typeof(lowerbound) != Nothing
        ### NEEDS TO CHANGE
        fs = add_theta_to_objective!(fs)
        addtheta = 1
        ###
                
        #todo create this function
        ### (potentially) NEEDS TO CHANGE
        add_lower_bound_to_first_stage!(fs, lowerbound)
        ###
    end
    
    for i = 1:niter
        
        itnum = i+curit
        itmax = curit+niter

        # step 1 set v = v+1 and solve first stage problem.
        println("Iteration $(itnum)/$(itmax)")

        optimize!(fs)

        if verbose == 1
            println("Updating first value...")
        end
        update_first_value_L!(firststage, fs)
        
        if firststage.store != nothing
        
            if itnum == 1
                if verbose == 1
                    println("Setting up First Stage paths...")
                end
                setup_1st_paths!(firststage)
            end
            
            if verbose == 1
                println("Storing x...")
            end
            store_x!(firststage)
            
        end

        # step 3 * update x-variables in second stage
        if verbose == 1
            println("Updating second values...")
        end
        update_second_value!(firststage)

        #        * solve second stage problems
        # done separately to eventually parallelize 
        
        for sid in keys(firststage.subproblems)  
            println("Solving subproblems and updating...")
            solve_sub_and_update!(firststage.subproblems[sid])
            
            if firststage.store!= nothing
                
                #I will have to store this locally at some point...
                path = firststage.store
                
                if itnum == 1
                    println("Setting up second stage paths...")
                    setup_scen_path!(path, sid)
                    
                    setup_2nd_paths!(path, firststage.subproblems[sid])
                end
                                                   
            end
        end

        #        * get simplex multipliers, update E and e
        # update E
        println("Updating first stage gradients...")
        update_first_gradient!(firststage)

        grad = get_grad_vector(firststage)
        if verbose == 1
            println("grad = $(grad)")
        end
        
        for sid in keys(firststage.subproblems)  
            println("For subproblem $(sid)..")
            subproblem = firststage.subproblems[sid]
            println("...adjusting h...")
            #temporarily removing, since there should not be interval constraints
            #adjust_h_new!(subproblem) 
            #see current Ek_ek folder)
            println("...computing Ek...")
            compute_Ek_new!(subproblem) 
            println("...computing ek...")
            compute_ek_new!(subproblem) 
            
            if firststage.store != nothing
                println("Storing Ek and ek...")
                store_Ek_sub!(subproblem, firststage.store) 
                store_ek_sub!(subproblem, firststage.store) 
                
                println("...saving primal variables...")
                sid = subproblem.id
                path_id = string(firststage.store, "scen_$(sid)/")
                save_cur_vars!(path_id, subproblem.model, itnum)
                
                println("...saving dual variables...")
                save_cur_duals!(path_id, subproblem, itnum)
            end
        end
                    
        #get it from subproblems. In async make get_Ek_from_file function
        println("Updating First stage stuff...")
        #TODO change this to gradient, as two-stage problems have a certain first stage cost ONLY
        x = get_value_vector(firststage)
        if verbose == 1
            println("x = $(x)")
        end
        
        cutcount = 0
        if addtheta == 0
            fs = add_theta_to_objective!(fs, firststage.subproblems.count)
            addtheta = 1
        end
        
        
        for sid in keys(firststage.subproblems)
            prob = firststage.subproblems[sid].probability
            Ek = prob*cost + firststage.subproblems[sid].Ek
            ek = firststage.subproblems[sid].ek
            
            println(sid, " ", Ek, " ", ek)
            #make this function
            
            w = ek - dot(Ek,x)
            
            if addthetak[sid] == 1
                thetak = JuMP.value(JuMP.variable_by_name(fs, "theta_$(sid)"))
                if thetak < w - tol
                    cutcount += 1
                    println(sid, " ", thetak, " ", w)
                    println("...subproblem $(sid) adding cut...")
                    fs = add_constraint_to_objective!(fs, Ek, ek, v_dict, sid)
                end
            else
                addthetak[sid] = 1
                cutcount += 1
                println("...subproblem $(sid) adding cut...")
                fs = add_constraint_to_objective!(fs, Ek, ek, v_dict, sid)
            end
                    
        end
                              
        ### NEEDS TO CHANGE
        if firststage.store != nothing
            store_E!(firststage, El)
            store_e!(firststage, el)
            if addtheta == 0
                store_w_theta!(firststage, w, -Inf)
            elseif addtheta == 1
                store_w_theta!(firststage, w, theta)
            end
        end
        ### 
        
        if cutcount == 0
            println("algorithm converged.")
            println("final x = $(x)")
            return x, firststage, fs, i
        end
        
        #may need to move for storing stuff
        #=
        if addtheta == 0
            fs = add_theta_to_objective!(fs, firststage.subproblems.count)
            addtheta = 1
        end
        =#


    end
    
    println("L-Shaped Algorithm Failed to converge in $(niter) iterations.")
    
    return x, firststage, fs, niter
    
end

function add_regularized_decomp!(firststage, model, linobj, rho)
    
    vardict = firststage.variables
    
    newfunc = JuMP.GenericQuadExpr{Float64,VariableRef}()
    
    for vname in keys(vardict)
        varinfo = vardict[vname]
        index = varinfo.index
        value = varinfo.value
        
        var = JuMP.variable_by_name(model, vname)
               
        JuMP.add_to_expression!(newfunc, rho/2*(var - value)^2 )
    end
    
    JuMP.@objective(model, Min, linobj + newfunc)
        
    return model
end
    
    

function iterate_L_regularized_decomp(firststage::FirstStageInfo, fs::JuMP.Model, v_dict::Dict{Int64,Array{Any}}, addtheta::Int64, tol::Float64, niter::Int64, verbose::Int64, resume::Int64, lowerbound::Union{Float64,Nothing}, rho, rhomin, rhomax, gamma)
    
    println(rho)
    println(rhomin)
    println(rhomax)
    println(gamma)
    
    if rho > rhomax
        rho = rhomax
    elseif rho < rhomin
        rho = rhomin
    end
        
    tau = 1e-8
    x = 0
    
    addthetak = zeros(Int64, firststage.subproblems.count)
    
    cost = get_cost_vector(firststage, fs)
    if verbose == 1
        println("cost = $(cost)")
    end
    
    curit = 0
    if resume > 0
        println("!!!!!!!!!!!!!! Resuming from iteration $(resume) using path $(firststage.store) !!!!!!!!!!!!!!")
        addtheta = 1
        #rebuild based on values of E and e
        ### NEEDS TO CHANGE
        curit, converged, xcur = resume_fs!(firststage, fs, tol)
        ###
        
        if converged > 0
            println("Model has already converged.")
            
            # 0 is placeholder for a better x or just having firststage hold everything
            return xcur, firststage, fs, curit
        end
    end
    if typeof(lowerbound) != Nothing
        ### NEEDS TO CHANGE
        fs = add_theta_to_objective!(fs)
        addtheta = 1
        ###
                
        #todo create this function
        ### (potentially) NEEDS TO CHANGE
        add_lower_bound_to_first_stage!(fs, lowerbound)
        ###
    end
    
    linobj = GenericAffExpr{Float64,VariableRef}()
    ssobja = 0.0;
    
    for i = 1:niter
        
        println("rho = ", rho)
        itnum = i+curit
        itmax = curit+niter

        # step 1 set v = v+1 and solve first stage problem.
        #println("Iteration $(itnum)/$(itmax)")
        
        optimize!(fs)
        
        avec = get_value_vector(firststage)
        
        fsobj = JuMP.objective_value(fs)        

        if verbose == 1
            println("Updating first value...")
        end
        update_first_value_L!(firststage, fs)
        
        if firststage.store != nothing
        
            if itnum == 1
                if verbose == 1
                    println("Setting up First Stage paths...")
                end
                setup_1st_paths!(firststage)
            end
            
            if verbose == 1
                println("Storing x...")
            end
            store_x!(firststage)
            
        end
        
        if itnum > 1
            println(itnum, " ", fsobj, " ", ssobja, " ", abs(fsobj-ssobja)/(tau + abs(ssobja)))
            if abs(fsobj - ssobja) < tol*(tau + abs(ssobja))

                println("algorithm converged.")
                println("final x = $(avec)")
                return avec, firststage, fs, i
            end
        end

        # step 3 * update x-variables in second stage
        if verbose == 1
            println("Updating second values...")
        end
        update_second_value!(firststage)

        #        * solve second stage problems
        # done separately to eventually parallelize 
        
        ssobjx = 0.0
        for sid in keys(firststage.subproblems)  
            #println("Solving subproblems and updating...")
            solve_sub_and_update!(firststage.subproblems[sid])
            
            ssobjx += firststage.subproblems[sid].probability*firststage.subproblems[sid].objective_value
            
            if firststage.store!= nothing
                
                #I will have to store this locally at some point...
                path = firststage.store
                
                if itnum == 1
                    #println("Setting up second stage paths...")
                    setup_scen_path!(path, sid)
                    
                    setup_2nd_paths!(path, firststage.subproblems[sid])
                end
                                                   
            end
        end
        if itnum == 1
            println(itnum, " ", fsobj, " ", ssobjx, " ", tol*(tau + abs(ssobjx))-abs(fsobj-ssobjx))
            if abs(fsobj - ssobjx) < tol*(tau + abs(ssobjx))

                println("algorithm converged.")
                println("final x = $(avec)")
                return avec, firststage, fs, i
            end
        end
        
        
        #        * get simplex multipliers, update E and e
        # update E
        #println("Updating first stage gradients...")
        update_first_gradient!(firststage)

        grad = get_grad_vector(firststage)
        if verbose == 1
            println("grad = $(grad)")
        end
        
        for sid in keys(firststage.subproblems)  
          #  println("For subproblem $(sid)..")
            subproblem = firststage.subproblems[sid]
           # println("...adjusting h...")
            #temporarily removing, since there should not be interval constraints
            #adjust_h_new!(subproblem) 
            #see current Ek_ek folder)
            #println("...computing Ek...")
            compute_Ek_new!(subproblem) 
            #println("...computing ek...")
            compute_ek_new!(subproblem) 
            
            if firststage.store != nothing
                println("Storing Ek and ek...")
                store_Ek_sub!(subproblem, firststage.store) 
                store_ek_sub!(subproblem, firststage.store) 
                
                println("...saving primal variables...")
                sid = subproblem.id
                path_id = string(firststage.store, "scen_$(sid)/")
                save_cur_vars!(path_id, subproblem.model, itnum)
                
                println("...saving dual variables...")
                save_cur_duals!(path_id, subproblem, itnum)
            end
        end
                    
        #get it from subproblems. In async make get_Ek_from_file function
        #println("Updating First stage stuff...")
        #TODO change this to gradient, as two-stage problems have a certain first stage cost ONLY
        x = get_value_vector(firststage)
        if verbose == 1
            println("x = $(x)")
        end
        
        cutcount = 0
        if addtheta == 0
            fs = add_theta_to_objective!(fs, firststage.subproblems.count)
            linobj = JuMP.objective_function(fs)
        end
        
        if addtheta == 0
            # x = a here
            fs = add_regularized_decomp!(firststage, fs, linobj, rho)
            ssobja = ssobjx
            addtheta = 1
        end
        
        
        for sid in keys(firststage.subproblems)
            prob = firststage.subproblems[sid].probability
            Ek = prob*cost + firststage.subproblems[sid].Ek
            ek = firststage.subproblems[sid].ek
            
            #println(sid, " ", Ek, " ", ek)
            #make this function
            
            w = ek - dot(Ek,x)
            
            if addthetak[sid] == 1
                thetak = JuMP.value(JuMP.variable_by_name(fs, "theta_$(sid)"))
                if thetak < w - tol
                    cutcount += 1
                    #println(sid, " ", thetak, " ", w)
                    #println("...subproblem $(sid) adding cut...")
                    fs = add_constraint_to_objective!(fs, Ek, ek, v_dict, sid)
                end
            else
                addthetak[sid] = 1
                cutcount += 1
                #println("...subproblem $(sid) adding cut...")
                fs = add_constraint_to_objective!(fs, Ek, ek, v_dict, sid)
            end
                    
        end
                              
        ### NEEDS TO CHANGE
        if firststage.store != nothing
            store_E!(firststage, El)
            store_e!(firststage, el)
            if addtheta == 0
                store_w_theta!(firststage, w, -Inf)
            elseif addtheta == 1
                store_w_theta!(firststage, w, theta)
            end
        end
        ### 
        
        ##### THIS WHOLE BLOCK SHOULD BE ADJUSTED ACCORDING TO WHETHER ADAPTIVE RHO IS BEING USED #####
        println(itnum, " ", ssobja, " ", ssobjx, " ", cutcount)
        changea = 0
        if cutcount == 0
            #fs = add_regularized_decomp!(firststage, fs, linobj, rho)
            changea = 1
            ssobja = ssobjx
        elseif ssobjx <= ssobja
            #fs = add_regularized_decomp!(firststage, fs, linobj, rho)
            changea = 1
            ssobja = ssobjx
        end

        if ssobja <= ssobjx - gamma*(ssobjx - fsobj)
            rho = max(rhomin, 0.5*rho)
        else
            rho = min(rhomax, 2.0*rho)
        end
        
        ### this is WRONG in the case where a is not changed but rho is. Perhaps add function to make it right (async's?)
        if changea == 1
            fs = add_regularized_decomp!(firststage, fs, linobj, rho)
        else
            fs = change_rho_only!(firststage, fs, linobj, avec, rho)
        end
        #####

    end
    
    println("L-Shaped Algorithm Failed to converge in $(niter) iterations.")
    
    return x, firststage, fs, niter
    
end

function change_rho_only!(firststage, model, linobj, avec, rho)
    
    vardict = firststage.variables
    
    newfunc = JuMP.GenericQuadExpr{Float64,VariableRef}()
    
    for vname in keys(vardict)
        varinfo = vardict[vname]
        index = varinfo.index
        value = avec[index]
        
        var = JuMP.variable_by_name(model, vname)
               
        JuMP.add_to_expression!(newfunc, rho/2*(var - value)^2 )
    end
    
    JuMP.@objective(model, Min, linobj + newfunc)
        
    return model
end

function iterate_CaCG(firststage::FirstStageInfo, fs::JuMP.Model, v_dict::Dict{Int64,Array{Any}}, addtheta::Int64, tol::Float64, niter::Int64, verbose::Int64, resume::Int64)
        
    x = 0
    UB = Inf
    LB = -Inf
    itnum = 0
    
    cost = get_cost_vector(firststage, fs)
    if verbose == 1
        println("cost = $(cost)")
    end
    
    curit = 0
    if resume > 0
        println("!!!!!!!!!!!!!! Resuming from iteration $(resume) using path $(firststage.store) !!!!!!!!!!!!!!")
        addtheta = 1
        #rebuild based on values of E and e
        curit, converged, xcur = resume_fs!(firststage, fs, tol)
        
        if converged > 0
            println("Model has already converged.")
            
            # 0 is placeholder for a better x or just having firststage hold everything
            return xcur, firststage, fs, curit
        end
    end
    
    wcscens = [];
    for i = 1:niter       
        
        itnum = i+curit
        itmax = curit+niter
        
        if abs(UB - LB)/UB < tol
            println(abs(UB - LB)/UB)
            x = get_value_vector(firststage)
            println("algorithm converged.")
            println("final x = $(x)")
            return x, firststage, fs, itnum-1, wcscens
        end

        # step 1 set v = v+1 and solve first stage problem.
        println("Iteration $(itnum)/$(itmax)")

        optimize!(fs)
        
        if itnum > 1
            LB = JuMP.objective_value(fs)
        end
        println("New LB is $(LB)")

        if verbose == 1
            println("Updating first value...")
        end
        update_first_value_L!(firststage, fs)
        
        if firststage.store != nothing
        
            if itnum == 1
                if verbose == 1
                    println("Setting up First Stage paths...")
                end
                setup_1st_paths!(firststage)
            end
            
            if verbose == 1
                println("Storing x...")
            end
            store_x!(firststage)
            
        end

        # step 3 * update x-variables in second stage
        if verbose == 1
            println("Updating second values...")
        end
        update_second_value!(firststage)

        #        * solve second stage problems
        # done separately to eventually parallelize 
        wcval = -Inf
        wcscen = 0
        for sid in keys(firststage.subproblems)  
            #println("Solving subproblems and updating...")
            solve_sub_and_update!(firststage.subproblems[sid])
            objval = JuMP.objective_value(firststage.subproblems[sid].model)
            println(sid, " ", objval)
            if  objval > wcval
                wcval = objval 
                wcscen = sid
            end
            if firststage.store!= nothing
                
                #I will have to store this locally at some point...
                path = firststage.store
                
                if itnum == 1
                    println("Setting up second stage paths...")
                    setup_scen_path!(path, sid)
                    
                    setup_2nd_paths!(path, firststage.subproblems[sid])
                end
                                                   
            end
        end
        push!(wcscens, wcscen)
        
        UB = min(UB,wcval)
        println("New UB is $(UB)")
        println("New WC scens are $(wcscens)")
        
        if addtheta == 0
            fs = add_theta_to_objective!(fs)
            addtheta = 1
        end
        
        if abs(UB - LB)/UB >= tol
            println(abs(UB - LB)/UB)
            n = length(wcscens)
            if n > 1
                if wcscens[n] in wcscens[1:n-1]
                else
                    add_wc_variables_to_fs_model!(fs, v_dict, wcscen)
                    v1array = get_1st_stage_variable_array(v_dict)
                    add_wc_constraints!(fs, firststage.subproblems[wcscen].model, wcscen, v1array)
                    add_obj_constraint!(fs, firststage.subproblems[wcscen].model, wcscen, v1array)
                end
            else
                add_wc_variables_to_fs_model!(fs, v_dict, wcscen)
                v1array = get_1st_stage_variable_array(v_dict)
                add_wc_constraints!(fs, firststage.subproblems[wcscen].model, wcscen, v1array)
                add_obj_constraint!(fs, firststage.subproblems[wcscen].model, wcscen, v1array)                
            end
        end
                
    end
    x = get_value_vector(firststage)
    println("L-Shaped Algorithm Failed to converge in $(niter) iterations.")
    
    return x, firststage, fs, itnum, wcscens
    
end

function add_obj_constraint!(exmodel::JuMP.Model, submodel::JuMP.Model, scen::Int64, v1array)
    
    subobj = JuMP.objective_function(submodel)
    terms = subobj.terms
    theta = JuMP.variable_by_name(exmodel, "theta")
    objfunc = AffExpr()
    for (subvar, coeff) in terms
        subvarname = JuMP.name(subvar)
        if subvarname in v1array
        else
            exvar = JuMP.variable_by_name(exmodel, string(subvarname,"_$(scen)"))
            JuMP.add_to_expression!(objfunc, coeff, exvar)
        end
    end
    JuMP.@constraint(exmodel, objfunc <= theta)
    
    return 
    
end
            
            
# this needs to be edited before building
function add_wc_constraints!(exmodel::JuMP.Model, submodel::JuMP.Model, scen::Int64, v1array)
    
    for (F,S) in list_of_constraint_types(submodel)
        for con in all_constraints(submodel,F,S)
            if occursin("AffExpr",string(F))
                newfunc = JuMP.AffExpr()
                terms = JuMP.constraint_object(con).func.terms
                if occursin("EqualTo", string(S))
                    value = constraint_object(con).set.value
                    for (subvar, coeff) in terms
                        subvarname = JuMP.name(subvar)
                        if subvarname in v1array 
                            exvar = JuMP.variable_by_name(exmodel, subvarname)
                        else
                            exvar = JuMP.variable_by_name(exmodel, string(subvarname,"_$(scen)"))
                        end
                        JuMP.add_to_expression!(newfunc, coeff, exvar)
                    end
                    JuMP.@constraint(exmodel, newfunc == value)
                elseif occursin("GreaterThan", string(S))
                    value = constraint_object(con).set.lower
                    count=0
                    for (subvar, coeff) in terms
                        subvarname = JuMP.name(subvar)
                        if subvarname in v1array 
                            exvar = JuMP.variable_by_name(exmodel, subvarname)
                        else
                            exvar = JuMP.variable_by_name(exmodel, string(subvarname,"_$(scen)"))
                        end
                        JuMP.add_to_expression!(newfunc, coeff, exvar)
                    end
                    JuMP.@constraint(exmodel, newfunc >= value)
                elseif occursin("LessThan", string(S))
                    value = constraint_object(con).set.upper
                    count = 0
                    for (subvar, coeff) in terms
                        subvarname = JuMP.name(subvar)
                        if subvarname in v1array 
                            exvar = JuMP.variable_by_name(exmodel, subvarname)
                        else
                            exvar = JuMP.variable_by_name(exmodel, string(subvarname,"_$(scen)"))
                        end
                        JuMP.add_to_expression!(newfunc, coeff, exvar)
                    end
                    JuMP.@constraint(exmodel, newfunc <= value)
                elseif occursin("Interval", string(S))
                    valu = constraint_object(con).set.upper
                    vall = constraint_object(con).set.lower
                    count = 0
                    for (subvar, coeff) in terms
                        subvarname = JuMP.name(subvar)
                        if subvarname in v1array 
                            exvar = JuMP.variable_by_name(exmodel, subvarname)
                        else
                            exvar = JuMP.variable_by_name(exmodel, string(subvarname,"_$(scen)"))
                        end
                        JuMP.add_to_expression!(newfunc, coeff, exvar)
                    end
                    JuMP.@constraint(exmodel, vall <= newfunc <= valu)
                else
                    println(con)
                    println("Add ", S, " to constraints.")
                end
            else
                subvar = JuMP.constraint_object(con).func
                if occursin("EqualTo", string(S))
                    value = constraint_object(con).set.value
                    subvarname = JuMP.name(subvar)
                    if subvarname in v1array 
                        # ignore fixed first stage variable constants
                        #exvar = JuMP.variable_by_name(exmodel, subvarname)
                    else
                        #println(con)
                        exvar = JuMP.variable_by_name(exmodel, string(subvarname,"_$(scen)"))
                        JuMP.@constraint(exmodel, exvar == value)
                    end
                elseif occursin("GreaterThan", string(S))
                    value = constraint_object(con).set.lower
                    subvarname = JuMP.name(subvar)
                    if subvarname in v1array 
                        exvar = JuMP.variable_by_name(exmodel, subvarname)
                    else
                        exvar = JuMP.variable_by_name(exmodel, string(subvarname,"_$(scen)"))
                    end
                    JuMP.@constraint(exmodel, exvar >= value)
                elseif occursin("LessThan", string(S))
                    value = constraint_object(con).set.upper
                    subvarname = JuMP.name(subvar)
                    if subvarname in v1array 
                        exvar = JuMP.variable_by_name(exmodel, subvarname)
                    else
                        exvar = JuMP.variable_by_name(exmodel, string(subvarname,"_$(scen)"))
                    end
                    JuMP.@constraint(exmodel, exvar <= value)
                elseif occursin("Interval", string(S))
                    valu = constraint_object(con).set.upper
                    vall = constraint_object(con).set.lower
                    subvarname = JuMP.name(subvar)
                    if subvarname in v1array 
                        exvar = JuMP.variable_by_name(exmodel, subvarname)
                    else
                        exvar = JuMP.variable_by_name(exmodel, string(subvarname,"_$(scen)"))
                    end
                    JuMP.@constraint(exmodel, vall <= exvar <= valu)
                else
                    println(con)
                    println("Add ", S, " to constraints.")
                end
            end
        end
    end
    
    return
end

function get_1st_stage_variable_array(v_dict)
    
    n = length(v_dict[1])
    v1array = [];
    for i = 1:n
        push!(v1array, v_dict[1][i][1])
    end
    
    return v1array
end

function add_wc_variables_to_fs_model!(model::JuMP.Model, v_dict, scen::Int64)
    
    n = length(v_dict[2])
    
    for i = 1:n
        JuMP.@variable(model, base_name = string(v_dict[2][i],"_$(scen)"))
    end
        
    return
end

