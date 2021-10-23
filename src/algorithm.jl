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
    
        varstructs[var].gradient = JuMP.dual(JuMP.FixRef(JuMP.variable_by_name(model,name)))
        
    end
    
    return
    
end

function update_variable_values!(varstructs, alpha) 
    for var in keys(varstructs)
        value = varstructs[var].value
        grad = varstructs[var].gradient

        value -= alpha*grad

        varstructs[var].value = value

    end
    
    return
    
end

#=
function iterate!(model, varstruct, linkedcons, alpha, niter=1)
    
    for i = 1:niter
        
        update_constraint_values!(varstructs, linkedcons)

        optimize!(model)

        update_gradients!(varstructs, linkedcons)

        update_variable_values!(varstructs, alpha) 

        println("x = ", varstructs[1].value)

    end
    
    return
    
end
=#

function solve_sub_and_update!(subproblem)
    
    varstructs = subproblem.variableinfo
    linkedcons = subproblem.linkedconstraintinfo
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
               
        cost[index] = JuMP.objective_function(fsmodel).terms[vref]
        
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

function update_first_gradient_active_set!(firststage)
    for var in keys(firststage.variables)
        value = firststage.variables[var].value
        grad1 = 0.0;
        # eventually use get functions in parallel.
        for sub = keys(firststage.subproblems)
            prob = firststage.subproblems[sub].probability
            vind = firststage.subproblems[sub].vnametoind[var]
            grad2 = firststage.subproblems[sub].variableinfo[vind].gradient
            grad1 += prob*grad2
        end
        firststage.variables[var].gradient = grad1
        variable = firststage.variables[var].name
        
        lowerbound = firststage.variables[var].lowerbound
        upperbound = firststage.variables[var].upperbound
        
        if value == lowerbound
            firststage.variables[var].gradient = min(0.0, grad1)
        elseif value == upperbound
            firststage.variables[var].gradient = max(0.0, grad1)
        end
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

function update_first_value_active_set!(firststage, alpha)
    for var in keys(firststage.variables)
        lowerbound = firststage.variables[var].lowerbound
        upperbound = firststage.variables[var].upperbound
        
        old_value = firststage.variables[var].value
        grad = firststage.variables[var].gradient
        new_value = old_value - alpha*grad
        if new_value > lowerbound
            if new_value < upperbound
                firststage.variables[var].value = new_value
            else
                firststage.variables[var].value = upperbound
            end
        else
            firststage.variables[var].value = lowerbound
        end
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


### Define separately as iterate and comment out the old one. ###
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

### Define separately as iterate and comment out the old one. ###
function iterate_active_set!(firststage)
    
    update_second_value!(firststage)
    
    # done separately to eventually parallelize 
    for i in keys(firststage.subproblems)  
        solve_sub_and_update!(firststage.subproblems[i])
    end

    update_first_gradient_active_set!(firststage)

    update_first_value_active_set!(firststage, 0.1)
    
    return
    
end

### First stage stuff below ###

# put in algorithm.jl

function make_active_set(firststage)

    lbset = [];
    ubset = [];
    naset = [];

    for var in keys(firststage.variables)

        a = 10^-5
        b = 10^-5 

        varinfo = firststage.variables[var]

        name = varinfo.name
        lowerbound = varinfo.lowerbound
        upperbound = varinfo.upperbound
        value = varinfo.value
        gradient = varinfo.gradient

        #println("$(name), $(value), $(gradient)")

        index = varinfo.index
        
        if value - a*gradient <= lowerbound
            #println("Hopefully Lowerbound")
            push!(lbset, index)
            varinfo.status='l'
        elseif value - b*gradient >= upperbound
            #println("Hopefully Upperbound")
            push!(ubset, index)
            varinfo.status='u'
        else
            #println("Hopefully non-active")
            push!(naset, index)
            varinfo.status='n'
        end

    end
    
    sort!(lbset)
    sort!(ubset)
    sort!(naset)
    
    nasdict = Dict()
    
    for i=1:length(naset)
        nasdict[naset[i]] = i
    end
    
    return lbset, ubset, naset, nasdict
    
end

#nvars is firststage.variables.count

#put in algorithm.jl

function initialize_H(nvars, theta=0.1)
    return theta*Matrix(I, nvars,nvars)
end

# create Z based on nonactive set "naset" knowing the number of first stage variables "nvars"

# put in algorithm.jl

function make_Z_matrix(naset, nvars)

    nna = length(naset)
    Z = Array{Float64}(undef, nvars, 0)
    for i = 1:nna
        vector = zeros(nvars)
        vector[i] = 1

        Z = [Z vector]
    end
    
    if length(Z) == 0
        println("Warning: Z is empty. This may cause problems.")
    end

    return Z
    
end

# put in algorithm.jl

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

# put in algorithm.jl

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

# put in algorithm.jl

function bound_value_differences(firststage, d, g)

    alpha_max = 1.0;

    vardict = firststage.variables

    for vname in keys(vardict)

        varinfo = vardict[vname]

        index = varinfo.index
        status = varinfo.status
        lb = varinfo.lowerbound
        ub = varinfo.upperbound
        value = varinfo.value

        if status == 'l'
            #println("l, $(lb-value), $(d[index])")
            d[index] = lb-value
            #if lb-value == 0
            #    println("$(index), 0")
            #else
            #    println("$(index), $(lb-value)")
            #end
        elseif status == 'u'
            #println("u, $(ub-value)")
            d[index] = ub-value
        else
            upper = ub - value
            lower = lb - value
            dval = d[index]
            
            if abs(dval) < 10^-6
                println("skipping, dval near 0")
                println(dval)
            elseif dval < 0
                
                #if lower/dval < 1e-7
                #    println("$(index), $(lower/dval), $(g[index]), $(lower), $(dval)")
                #end
                if lower/dval < 1e-6
                    d[index] = lower
                else
                    alpha_max = min(alpha_max, lower/dval)
                end
                #alpha_max = min(alpha_max, lower/dval)
            elseif dval > 0
                alpha_max = min(alpha_max, upper/dval)
            else
                println(dval)
                println("something is wrong")
            end
        end
    end

    for vname in keys(vardict)
        
        varinfo = vardict[vname]
        
        index = varinfo.index
        status = varinfo.status
        value = varinfo.value
        
        if status == 'n'
            d[index] = alpha_max*d[index]
        end
    end

    #println("alpha_max = $(alpha_max)")
    return d, alpha_max

end


# updates value based on x vector computed prior. 
# simply plugs into the struct.

# put in algorithm.jl
function update_first_value_bfgs!(firststage, x)
    
    for vname in keys(firststage.variables)
        
        varinfo = firststage.variables[vname]
        index = varinfo.index
        
        varinfo.value = x[index]
        
    end
    
    return
    
end

# put in algorithm.jl

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

# put in algorithm.jl

function update_second_value_bfgs!(firststage, x)
    for sub in keys(firststage.subproblems)
        subproblem = firststage.subproblems[sub]
        
        vardict = subproblem.variableinfo
        for var in keys(vardict)
            varinfo = vardict[var]
            index = varinfo.index
            varinfo.value = x[index]
        end
    end
    return
end

# put in algorithm.jl

function update_first_value_bfgs!(firststage, x)
    
    for var in keys(firststage.variables)
        index = firststage.variables[var].index
        firststage.variables[var].value = x[index]
    end
    
    return
    
end

function adjust_h(firststage, contoidx, h)
# make link
    
    models = firststage.subproblems
    
    ns = models.count
    nc = contoidx.count
        
    for sid in keys(models)
        m = models[sid].model
        
        for (F,S) in list_of_constraint_types(m)
            for con in all_constraints(m,F,S)
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
                            h[sid, contoidx[idx]] = constraint_object(con).set.lower
                        else
                            h[sid, contoidx[idx]] = constraint_object(con).set.upper
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
    nc = contoidx.count
    PI = Array{Float64}(undef,N, nc)
    for i = 1:N
        subproblem = firststage.subproblems[i]

        m = subproblem.model

        for (F,S) in list_of_constraint_types(m)
            for con in all_constraints(m,F,S)
               if occursin("AffExpr",string(F))
                    #f = MOI.get(moi_backend, MOI.ConstraintFunction(), con.index)
                    #conidxtoref[index] = (F,S,f,con, con.index.value)
                    idx = con.index.value
                    dual = JuMP.dual(con)
                    PI[i, contoidx[idx]] = dual
               end
            end
        end
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
    

function iterate_L(firststage, fs, contoidx, h, v_dict, addtheta = 0, tol = 1e-6, niter = 10)
    
    x = 0
    
    cost = get_cost_vector(firststage, fs)
    
    for i = 1:niter

        # step 1 set v = v+1 and solve first stage problem.
        println("Iteration $(i)")

        optimize!(fs)

        update_first_value_L!(firststage, fs)

        # step 3 * update x-variables in second stage
        update_second_value!(firststage)

        #        * solve second stage problems
        # done separately to eventually parallelize 
        for i in keys(firststage.subproblems)  
            solve_sub_and_update!(firststage.subproblems[i])
        end

        #        * get simplex multipliers, update E and e
        # update E
        update_first_gradient!(firststage)

        grad = get_grad_vector(firststage)
                
        h = adjust_h(firststage, contoidx, h)

        E = cost - grad
        #println("E = $(E)")


        #update PI
        PI = compute_PI(firststage, contoidx)

        #update e_k
        e_k = compute_e(firststage,h,PI)
        println("e = $(e_k)")


        x = get_value_vector(firststage)
    

        if addtheta == 1
            theta = JuMP.value(JuMP.variable_by_name(fs, "theta"))
            println("theta = $(theta)")
        end

        w = e_k - dot(E,x)
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








