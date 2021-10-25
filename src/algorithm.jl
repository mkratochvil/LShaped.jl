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
    
    println("function alg 264 iterate!")
    
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








