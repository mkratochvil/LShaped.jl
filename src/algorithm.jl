function update_constraint_values!(varstructs, linkedcons)
    
    for con in keys(linkedcons)
        curval = linkedcons[con].initvalue
        for var in linkedcons[con].variables
            curval -= varstructs[var].conval[con]*varstructs[var].value
        end
        settype = string(linkedcons[con].set)
        functype = string(linkedcons[con].func)
        if occursin("EqualTo", settype) || occursin("EqualTo", functype)
            MOI.set(linkedcons[con].ref.model, MOI.ConstraintSet(), linkedcons[con].ref, MOI.EqualTo(curval))
        elseif occursin("LessThan", settype) || occursin("LessThan", functype)
            MOI.set(linkedcons[con].ref.model, MOI.ConstraintSet(), linkedcons[con].ref, MOI.LessThan(curval))
        elseif occursin("GreaterThan", settype) || occursin("GreaterThan", functype)
            MOI.set(linkedcons[con].ref.model, MOI.ConstraintSet(), linkedcons[con].ref, MOI.GreaterThan(curval))
        else
            println(settype)
            println(functype)
            println("Something went wrong.")
        end
        linkedcons[con].curvalue = curval
    end
    
    return
end

function update_constraint_values_nac!(model, varstructs)
    
    for var in keys(varstructs)
        name = varstructs[var].name
        value = varstructs[var].value
        JuMP.fix(JuMP.variable_by_name(model,name),value)
    end
    
    return
end

function update_gradients!(varstructs, linkedcons) 
    
    for var in keys(varstructs)
        grad = varstructs[var].cost
        for con in keys(varstructs[var].conval)
            #dual*coeff*var
            dual = JuMP.dual(linkedcons[con].ref)
            if occursin("LessThan", string(linkedcons[con].func))
                if dual > 0
                    #println("Positive dual $(dual) at constraint $(con)")
                    dual = 0
                end
            elseif occursin("GreaterThan", string(linkedcons[con].func))
                if dual < 0
                    #println("Negative dual $(dual) at constraint $(con)")
                    dual = 0
                end
            end
            coeff = varstructs[var].conval[con]
            grad -= dual*coeff
        end
        varstructs[var].gradient = grad
    end
    
    return
    
end
    
function update_gradients_nac!(model, varstructs)
    
    for var in keys(varstructs)
        name = varstructs[var].name
        
        cost = varstructs[var].cost
        #println("$(name), $(cost)")
    
        varstructs[var].gradient = JuMP.dual(JuMP.FixRef(JuMP.variable_by_name(model,name)))
        
    end
    
    return
    
end

function solve_sub_and_update!(subproblem)
        
    varstructs = subproblem.variableinfo
    #linkedcons = subproblem.linkedconstraintinfo
    model = subproblem.model

    #update_constraint_values!(varstructs,linkedcons)
    update_constraint_values_nac!(model, varstructs)
    
    println("solving subproblem: $(subproblem.id)")
    optimize!(model)
    update_gradients_nac!(model, varstructs)
    
    subproblem.objective_value = objective_value(model)
    
    return
    
end

function get_cost_vector(firststage, fsmodel)
        
    vardict = firststage.variables
    nvars = vardict.count
    
    cost = Vector{Float64}(undef, nvars) 
    
    for vname in keys(vardict)
        
        varinfo = vardict[vname]
        index = varinfo.index
        
        vref = JuMP.variable_by_name(fsmodel, vname)
               
        #this might not work. there is in fact a second stage cost at the moment.
        if vref in keys(JuMP.objective_function(fsmodel).terms)
            cost[index] = JuMP.objective_function(fsmodel).terms[vref]
        else
            cost[index] = 0.0
        end
        
    end
    
    return cost
    
end

function update_first_gradient!(firststage)
        
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
        #variable = firststage.variables[var].name
    end
    
    return
    
end

# maybe adjust alpha eventually.

function update_first_value!(firststage, alpha)
        
    for var in keys(firststage.variables)
        old_value = firststage.variables[var].value
        grad = firststage.variables[var].gradient
        new_value = old_value - alpha*grad
        firststage.variables[var].value = new_value
        variable = firststage.variables[var].name
    end
    
    return
    
end

#eventually parallelize and put another function inside.
function update_second_value!(firststage) 
    
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

### this uses direct definition and not nac constraints ###
function iterate!(firststage) 
        
    # done separately to eventually parallelize 
    for i in keys(firststage.subproblems)  
        solve_sub_and_update!(firststage.subproblems[i])
    end

    update_first_gradient!(firststage)

    update_first_value!(firststage, 0.1)

    update_second_value!(firststage)
    
    return
    
end

function get_grad_vector(firststage)
        
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

function get_value_vector(firststage)
        
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

function get_objective_value(firststage)
    
    #change this. this could lead to issues.
    objval = 0.0;

    subdict = firststage.subproblems

    for sid in keys(subdict)

        subobjval = subdict[sid].objective_value
        
        vardict = subdict[sid].variableinfo
        
        for var in keys(vardict)
            
            varinfo = vardict[var]
            
            value = varinfo.value
            cost = varinfo.cost
            
            subobjval += cost*value
            
        end
            
            
        prob = subdict[sid].probability

        objval += prob*subobjval

    end
    
    return objval
    
end

function adjust_h(firststage, contoidx, h)
        
    models = firststage.subproblems
    
    ns = models.count
    nc = models[1].ncons#contoidx.count
        
    for sid in keys(models)
        m = models[sid].model
        
        for (F,S) in list_of_constraint_types(m)
            innercount = 0
            for con in all_constraints(m,F,S)
                innercount += 1
               if occursin("AffExpr",string(F))
                    idx = con.index.value
                    if occursin("EqualTo", string(S))
                        #val = constraint_object(con).set.value
                        #h[sid, contoidx[idx]] = val
                    elseif occursin("GreaterThan", string(S))
                        #val = constraint_object(con).set.lower
                        #h[sid, contoidx[idx]] = val
                    elseif occursin("LessThan", string(S))
                        #val = constraint_object(con).set.upper
                        #h[sid, contoidx[idx]] = val
                    elseif occursin("Interval", string(S))
                        #val = constraint_object(con).set.upper
                        #for key in keys(constraint_object(con).func.terms)
                        #    #println("$(key), $(JuMP.value(key))")
                        #    varval
                        #end
                        #println("hval = ", h[sid, contoidx[idx]])
                        #println("dual = ", dual(con))
                        #h[sid, contoidx[idx]] = val
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
function adjust_h_new!(subproblem)

    idxtocon = subproblem.idxtocon
    h = subproblem.h
    
    for idx in keys(idxtocon)
        
        con = idxtocon[idx]
        
        #this is a placeholder, as below is obviously poor coding practice.
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

function compute_e(firststage, h, PI)
        
    e_k = 0.0;
    for i in keys(firststage.subproblems)
        prob = firststage.subproblems[i].probability
        e_k += prob*dot(PI[i,:],h[i,:])
    end
    return e_k
end
    

function compute_PI(firststage, contoidx)
        
    N = firststage.subproblems.count
    #nc = firststage.subproblems[1].ncons
    nc = contoidx.count
    println(nc)
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
                    #f = MOI.get(moi_backend, MOI.ConstraintFunction(), con.index)
                    #conidxtoref[index] = (F,S,f,con, con.index.value)
                    idx = con.index.value
                    dual = JuMP.dual(con)
                    PI[i, contoidx[(F,S,innercount)]] = dual
                count += 1
               end
            end
        end
        println(count)
    end
    return PI
end
;

function add_theta_to_objective!(fs)
        
    @variable(fs, theta)
    
    obj_old = objective_function(fs)
    
    @objective(fs, Min, obj_old + theta)
    
    #fsnew = Model();
    
    #MathOptInterface.copy_to(fsnew.moi_backend.model_cache.model, fs.moi_backend.model_cache.model)
    
    return fs#new
    
end
    

function add_constraint_to_objective!(fs, E, e_k, v_dict)
        
    nvar = length(E)
    
    @constraint(fs, sum(E[i]*variable_by_name(fs, v_dict[1][i][1]) for i = 1:nvar) 
                        + variable_by_name(fs, "theta") >= e_k)
        
    return fs
    
end


function update_first_value_L!(firststage, fs)
        
    for vname in keys(firststage.variables)
        #println(vname)
        varinfo = firststage.variables[vname]
        #println(varinfo)
        name = varinfo.name
        
        varinfo.value = JuMP.value(JuMP.variable_by_name(fs, name))
        #println(name, " ", varinfo.value)
    end
    
    return
    
end

function setup_1st_paths!(firststage)
    
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
    dfwt = DataFrame(w = Float64[], theta = Float64[])
    
    for i = 1:nvars
        insertcols!(df, i, Symbol(header[i])=>Float64[])
    end
    
    mkdir(path)
    
    xcsv = string(path, "x.csv")
    Ecsv = string(path, "E.csv")
    ecsv = string(path, "ek.csv")
    wtcsv = string(path, "w_theta.csv")
    
    CSV.write(xcsv, df)
    CSV.write(Ecsv, df)
    CSV.write(ecsv, dfe)
    CSV.write(wtcsv, dfwt)
    
    return
end

function store_x!(firststage)
    
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

function store_E!(firststage, E)
    
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

function store_e!(firststage, e_k)
    
    path = firststage.store
    ecsv = string(path, "ek.csv")
        
    open(ecsv, "a") do io
        write(io, "$(e_k) \n")
    end
    
    return
    
end

function store_w_theta!(firststage, w, theta)
    
    path = firststage.store
    wtcsv = string(path, "w_theta.csv")
    
    open(wtcsv, "a") do io
        write(io, "$(w), $(theta) \n")
    end
    
    return
    
end

function setup_2nd_paths!(path, subproblem)
    
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

function setup_scen_path!(path::String, sid::Int64)
    
    scenpath = string(path, "scen_$(sid)/")
    
    mkdir(scenpath)
    
    return 
    
end

function store_Es_es!(path, subproblem, pi_k, h_k)
           
    prob = subproblem.probability
    vardict = subproblem.variableinfo
    sid = subproblem.id
    nvars = vardict.count
    
    # initialize E vector
    Evec = Vector{Float64}(undef, nvars) 
    
    for var in keys(vardict)
        vind = vardict[var].index
        grad2 = vardict[var].gradient
        cost = vardict[var].cost
        
        # this is using the NAC method
        Evec[vind] = prob*(cost - grad2)
        
    end
    
    # now it's e's turn
    #=m = subproblem.model
    nc = contoidx.count
    PIk = Array{Float64}(undef, nc)

    for (F,S) in list_of_constraint_types(m)
        for con in all_constraints(m,F,S)
           if occursin("AffExpr",string(F))
                #f = MOI.get(moi_backend, MOI.ConstraintFunction(), con.index)
                #conidxtoref[index] = (F,S,f,con, con.index.value)
                idx = con.index.value
                dual = JuMP.dual(con)
                PIk[contoidx[idx]] = dual
           end
        end
    end=#
    
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


function iterate_L(firststage, fs, contoidx, h, v_dict, addtheta = 0, tol = 1e-6, niter = 10)
        
    x = 0
    
    cost = get_cost_vector(firststage, fs)
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
        println("grad = $(grad)")
                
        h = adjust_h(firststage, contoidx, h)
        println("h = $(h)")

        E = cost - grad
        if firststage.store != nothing
            store_E!(firststage, E)
        end
        println("E = $(E)")


        #update PI
        PI = compute_PI(firststage, contoidx)
        println("PI = $(PI)")
        
        #this needs to be redesigned. e_k in the second stage needs to be recorded first in the fashion of E
        #this requires restructuring how contoidx, pi, and h are constructed in the second stage.
        #this will probably take a full day of work but will be needed to do the async stuff.
        
        if firststage.store != nothing
            for sid in keys(firststage.subproblems)
                store_Es_es!(firststage.store, firststage.subproblems[sid], PI[sid,:], h[sid,:])
            end
        end

        #update e_k
        e_k = compute_e(firststage,h,PI)
        println("e = $(e_k)")
        
        if firststage.store != nothing
            store_e!(firststage, e_k)
        end


        x = get_value_vector(firststage)
        println("x = $(x)")
    

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

function compute_Ek_new!(subproblem)
    
    prob = subproblem.probability
    vardict = subproblem.variableinfo
    nvars = vardict.count
    
    # initialize E vector
    Evec = Vector{Float64}(undef, nvars) 
    
    for var in keys(vardict)
        vind = vardict[var].index
        grad2 = vardict[var].gradient
        cost = vardict[var].cost
        println(var, " ", cost)
        
        # this is using the NAC method
        Evec[vind] = prob*(cost - grad2)
        
    end
    
    subproblem.Ek = Evec
    
    return
    
end

function compute_ek_new!(subproblem)
    
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

function store_Ek_sub!(subproblem, path)
    
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

function store_ek_sub!(subproblem, path)
    
    sid = subproblem.id
    ek = subproblem.ek
    
    ecsv = string(path, "scen_$(sid)/ek.csv")
    
    open(ecsv, "a") do io
        write(io, "$(ek) \n")
    end
    
    return
    
end
    
    
#=function store_Es_es!(path, subproblem, pi_k, h_k)
           
    prob = subproblem.probability
    vardict = subproblem.variableinfo
    sid = subproblem.id
    nvars = vardict.count
    
    # initialize E vector
    Evec = Vector{Float64}(undef, nvars) 
    
    for var in keys(vardict)
        vind = vardict[var].index
        grad2 = vardict[var].gradient
        cost = vardict[var].cost
        
        # this is using the NAC method
        Evec[vind] = prob*(cost - grad2)
        
    end
    
    # now it's e's turn
    #=m = subproblem.model
    nc = contoidx.count
    PIk = Array{Float64}(undef, nc)

    for (F,S) in list_of_constraint_types(m)
        for con in all_constraints(m,F,S)
           if occursin("AffExpr",string(F))
                #f = MOI.get(moi_backend, MOI.ConstraintFunction(), con.index)
                #conidxtoref[index] = (F,S,f,con, con.index.value)
                idx = con.index.value
                dual = JuMP.dual(con)
                PIk[contoidx[idx]] = dual
           end
        end
    end=#
    
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
    
end=#

#as opposed to get_El_from_file to be made
function get_El_from_sub!(firststage)
    
    subproblems = firststage.subproblems
    nvars = firststage.variables.count
    
    El = zeros(nvars)
    
    for sid in keys(subproblems)
        El += subproblems[sid].Ek
    end
    
    return El
    
end

function get_el_from_sub!(firststage)
    
    subproblems = firststage.subproblems
    
    el = 0.0
    
    for sid in keys(subproblems)
        el += subproblems[sid].ek
    end
    
    return el
    
end

function iterate_L_new(firststage, fs, v_dict, addtheta = 0, tol = 1e-6, niter = 10)
        
    x = 0
    
    cost = get_cost_vector(firststage, fs)
    println("cost = $(cost)")
    
    for i = 1:niter

        # step 1 set v = v+1 and solve first stage problem.
        println("Iteration $(i)")

        optimize!(fs)

        println("Updating first value...")
        update_first_value_L!(firststage, fs)
        
        if firststage.store != nothing
        
            if i == 1
                println("Setting up First Stage paths...")
                setup_1st_paths!(firststage)
            end
            
            println("Storing x...")
            store_x!(firststage)
            
        end

        # step 3 * update x-variables in second stage
        println("Updating second values...")
        update_second_value!(firststage)

        #        * solve second stage problems
        # done separately to eventually parallelize 
        
        for sid in keys(firststage.subproblems)  
            println("Solving subproblems and updating...")
            solve_sub_and_update!(firststage.subproblems[sid])
            
            if firststage.store!= nothing
                
                #I will have to store this locally at some point...
                path = firststage.store
                
                if i == 1
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
        println("grad = $(grad)")
        
        for sid in keys(firststage.subproblems)  
            println("For subproblem $(sid)..")
            subproblem = firststage.subproblems[sid]
            println("...adjusting h...")
            adjust_h_new!(subproblem) #done
            #see current Ek_ek folder)
            println("...computing Ek...")
            compute_Ek_new!(subproblem) #done
            println("...computing ek...")
            compute_ek_new!(subproblem) #done
            
            if firststage.store != nothing
                println("Storing Ek and ek...")
                store_Ek_sub!(subproblem, firststage.store) #done
                store_ek_sub!(subproblem, firststage.store) #done
            end
        end
                    
        #E = cost - grad
        #get it from subproblems. In async make get_Ek_from_file function
        println("Updating First stage stuff...")
        #TODO change this to gradient, as two-stage problems have a certain first stage cost ONLY
        El = get_El_from_sub!(firststage) #done
        El = cost + El
        println("El = $(El)")
        el = get_el_from_sub!(firststage) #done
        println("el = $(el)")
        
       

        x = get_value_vector(firststage)
        println("x = $(x)")
    
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








