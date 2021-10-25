#### ALGORITHM #####

function update_variable_values!(varstructs, alpha) 
    println("function alg 84 update_variable_values!")
    
    for var in keys(varstructs)
        value = varstructs[var].value
        grad = varstructs[var].gradient

        value -= alpha*grad

        varstructs[var].value = value

    end
    
    return
    
end

#=function update_first_gradient_active_set!(firststage)
    
    println("function alg 169 update_first_gradient_active_set!")
    
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
    
end=#

#=function update_first_value_active_set!(firststage, alpha) 
    
    println("function alg 218 update_first_value_active_set!")
    
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
    
end=#

### Define separately as iterate and comment out the old one. ###
#=function iterate_active_set!(firststage)
    
    println("function alg 286 iterate_active_set!")
    
    update_second_value!(firststage)
    
    # done separately to eventually parallelize 
    for i in keys(firststage.subproblems)  
        solve_sub_and_update!(firststage.subproblems[i])
    end

    update_first_gradient_active_set!(firststage)

    update_first_value_active_set!(firststage, 0.1)
    
    return
    
end=#

### First stage stuff below ###

# put in algorithm.jl

#=function make_active_set(firststage)
    
    println("function alg 309 make_active_set")

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
    
end=#

#nvars is firststage.variables.count

#put in algorithm.jl

#=function initialize_H(nvars, theta=0.1)
    
    println("function alg 368 initialize_H")
    
    return theta*Matrix(I, nvars,nvars)
end=#

# create Z based on nonactive set "naset" knowing the number of first stage variables "nvars"

# put in algorithm.jl

#=function make_Z_matrix(naset, nvars)
    
    println("function alg 379 make_Z_matrix")

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
    
end=#

# put in algorithm.jl

# put in algorithm.jl

#=function bound_value_differences(firststage, d, g)
    
    println("function 452 bound_value_differences")

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
    
    println("function alg 532 update_first_value_bfgs!")
    
    for vname in keys(firststage.variables)
        
        varinfo = firststage.variables[vname]
        index = varinfo.index
        
        varinfo.value = x[index]
        
    end
    
    return
    
end=#

# put in algorithm.jl

# put in algorithm.jl

#=function update_second_value_bfgs!(firststage, x)
    
    println("function alg 589 update_second_value_bfgs!")
    
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
end=#

#=function update_first_value_bfgs!(firststage, x)
    
    println("function alg 608 update_first_value_bfgs!")
    
    for var in keys(firststage.variables)
        index = firststage.variables[var].index
        firststage.variables[var].value = x[index]
    end
    
    return
    
end=#






