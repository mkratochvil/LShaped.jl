function varRefToIdx(m)
        
    varrefs = all_variables(m)
    varreftoidx = Dict()
    
    for var in varrefs
        varreftoidx[var] = var.index.value
    end
    
    return varreftoidx

end
    
function varIdxToRef(m)
        
    varrefs = all_variables(m)
    varidxtoref = Dict()

    for var in varrefs
        varidxtoref[var.index.value] = var
    end

    return varidxtoref
end

function varIdxToCost(m)
        
    varidxtocost = Dict()
    
    obj= MOI.get(m, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}())
    
    for term in obj.terms
        varidxtocost[term.variable_index.value] = term.coefficient
    end
    
    return varidxtocost
    
end

function conIdxToRef(m)
        
    conidxtoref = Dict()
    moi_backend = backend(m)
    index = 0
    for (F,S) in list_of_constraint_types(m)
        for con in all_constraints(m,F,S)
            index += 1
           if occursin("AffExpr",string(F))
                f = MOI.get(moi_backend, MOI.ConstraintFunction(), con.index)
                conidxtoref[index] = (F,S,f,con, con.index.value)
           end
        end
    end
    
    return conidxtoref
    
end

function make_dicts_and_arrays(m)
        
    vr = varRefToIdx(m)
    vi = varIdxToRef(m)
    vc = varIdxToCost(m)
    cr = conIdxToRef(m)
    #ae, al, ag, ie, il, ig = create_eq_lt_gt_arrays(m,cr,vi)
    
    #return vr, vi, vc, cr, ae, al, ag, ie, il, ig
    return vr, vi, vc, cr #ae, al, ag, ie, il, ig
    
end


function make_VariableInfo(vardict, vartocost, condict)
    
    varstructs = Dict()
    
    for vname in keys(vardict[1])
        vid = vardict[1][vname]
        conval=Dict();
        # note vartocost may be obsolete. directly implement in class?
        cost = vartocost[vid]
        # note this may change. I seek out the constraints in which variables appear iteratively
        # better to do this by constraint or by variable?
        for cid in keys(condict)
            if occursin("MathOptInterface.ScalarAffineFunction{Float64}", string(typeof(condict[cid][3])))
                for term in 1:length(condict[cid][3].terms)
                    if vid == condict[cid][3].terms[term].variable_index.value
                        concoeff = condict[cid][3].terms[term].coefficient
                        #tuple = (cid, concoeff)
                        #push!(conval,tuple)
                        conval[cid] = concoeff
                    end
                end
                #SingleVariable
            elseif occursin("MathOptInterface.ScalarAffineTerm", string(typeof(condict[cid][3])))
                if vid == condict[cid][3].variable.value
                    println("it appears! :O")
                    # this is not added to the variable class, as it is assumed bounds will be in the first stage
                    #println(typeof(condict[cid][3]),cid)
                end
            else
                #println(typeof(condict[cid][3]))
            end
        end
        varstructs[vid] = LocalVariableInfo(vid, vname, nothing, cost, conval, nothing, nothing)
    end

    return varstructs
end

function make_LinkedConstraintInfo(varstructs, condict)
        
    linkedcons = Dict()
    for var in keys(varstructs)
        convals = varstructs[var].conval
        for (cid, cval) in convals
            if cid in keys(linkedcons)
                #needs to be pushed to add for other variables
                #println(cid)
                #unchecked
                push!(linkedcons[cid].variables, var)
            else
                if occursin("MathOptInterface.EqualTo{Float64}", string(constraint_object(condict[cid][4]).set))
                    val = constraint_object(condict[cid][4]).set.value
                    linkedcons[cid] = LinkedConstraintInfo(condict[cid][5], 
                                                            condict[cid][1], 
                                                            condict[cid][2], 
                                                            condict[cid][4], [var], val, val)
                elseif occursin("MathOptInterface.GreaterThan{Float64}",string(constraint_object(condict[cid][4]).set))
                    val = constraint_object(condict[cid][4]).set.lower
                    linkedcons[cid] = LinkedConstraintInfo(condict[cid][5], 
                                                                condict[cid][2], condict[cid][3], condict[cid][4], [var], val, val)
                elseif occursin("MathOptInterface.LessThan{Float64}",string(constraint_object(condict[cid][4]).set))
                    val = constraint_object(condict[cid][4]).set.upper
                    linkedcons[cid] = LinkedConstraintInfo(condict[cid][5], 
                                                                condict[cid][2], condict[cid][3], condict[cid][4], [var], val, val)
                else
                    # change this to a warning eventually
                    #println(condict[cid][4])
                    println("Add ", string(constraint_object(condict[cid][4]).set), " to initialization. con: $(condict[cid][4])")
                end
            end
        end
    end
    
    return linkedcons
end

function stage_name_idx(model::JuMP.Model, vardict)
        
    vdict = Dict()
    vdict[1]=Dict()
    vdict[2]=Dict()
    #vdict[3]=[]
    vardict_array = [];
    for varinfo in vardict[1]
        push!(vardict_array, varinfo[1])
    end
    for var in JuMP.all_variables(model)
        name = JuMP.name(var)
        if  name in vardict_array #occursin("x", name)
            val = JuMP.variable_by_name(model, name).index.value
            vdict[1][name] = val
        else #if occursin("Thsp", name)
            val = JuMP.variable_by_name(model, name).index.value
            vdict[2][name] = val
        #else
        #    push!(vdict[3], name)
        end
    end
    return vdict
end

function initialize(model, vnames)
        
    vardict = stage_name_idx(model, vnames)
    
    #vrm1, vindtoref, varcost, condict, Ae, Al, Ag, Ie, Il, Ig = make_dicts_and_arrays(model);
    vrm1, vindtoref, varcost, condict = make_dicts_and_arrays(model);
    
    varstructs = make_VariableInfo(vardict, varcost, condict);
    
    #arrays = create_eq_lt_gt_arrays(m,condict,vindtoref, varstructs)
    
    linkedcons = make_LinkedConstraintInfo(varstructs, condict)
    
    #for i in 1:length(x_init)
    #    varstructs[i].value = x_init[i]
    #end
    
    #for name in keys(vardict[1])
    #    var = JuMP.variable_by_name(model, name)
    #    delete(model, var)
    #end
    
    #arrays = Arrays(Ae, Al, Ag, Ie, Il, Ig)
    
    return model, varstructs, linkedcons, vardict[1], 0
    
end

function update_second_index!(firststage)
    
    for var in keys(firststage.variables)
        for sub = keys(firststage.subproblems)
            subproblem = firststage.subproblems[sub]
            vind = subproblem.vnametoind[var]
            index = firststage.variables[var].index
            subproblem.variableinfo[vind].index = index
        end
    end
    
    return
    
end


# huge assumption that the subproblems will have the same constraint order.
# this should be true for my problems.
# May cause issues with WSGEP
function ConToIdx(m)
    
    contoidx = Dict()

    count = 0

    for (F,S) in list_of_constraint_types(m)
        for con in all_constraints(m,F,S)
           if occursin("AffExpr",string(F))
                count += 1
                # this may have to change because of how JuMP stores things. we'll see.
                contoidx[con.index.value] = count
           end
        end
    end
    
    return contoidx
    
end

#creates a constraint number and lists the constraint in it.
function IdxToCon(m)
    
    idxtocon = Dict()

    count = 0

    for (F,S) in list_of_constraint_types(m)
        for con in all_constraints(m,F,S)
           if occursin("AffExpr",string(F))
                count += 1
                idxtocon[count] = con
           end
        end
    end
    
    return idxtocon
    
end

function compute_h_new(model,idxtocon)

    nc = idxtocon.count
    
    h = Array{Float64}(undef, nc)
    
    for idx in keys(idxtocon)
        
        con = idxtocon[idx]
        
        #this is a placeholder, as below is obviously poor coding practice.
        ctype = string(typeof(con))
        println(con)
        
        if occursin("EqualTo", ctype)
            val = constraint_object(con).set.value
                h[idx] = val
        elseif occursin("GreaterThan", ctype)
            val = constraint_object(con).set.lower
            h[idx] = val
        elseif occursin("LessThan", ctype)
            val = constraint_object(con).set.upper
            h[idx] = val
        elseif occursin("Interval", ctype)
            val = constraint_object(con).set.upper
            h[idx] = val
        else
            println(con)
            println("Add ", ctype, " to hvars.")
        end
    end
    
    return h
end

function compute_h(models, contoidx)
    
    ns = models.count
    nc = contoidx.count
    
    h = Array{Float64}(undef, ns, nc)
    
    for sid in keys(models)
        m = models[sid]
        
        for (F,S) in list_of_constraint_types(m)
            for con in all_constraints(m,F,S)
               if occursin("AffExpr",string(F))
                    idx = con.index.value
                    if occursin("EqualTo", string(S))
                        val = constraint_object(con).set.value
                        h[sid, contoidx[idx]] = val
                    elseif occursin("GreaterThan", string(S))
                        val = constraint_object(con).set.lower
                        h[sid, contoidx[idx]] = val
                    elseif occursin("LessThan", string(S))
                        val = constraint_object(con).set.upper
                        h[sid, contoidx[idx]] = val
                    elseif occursin("Interval", string(S))
                        #println(con)
                        #if idx == 3455
                        #    println("Success!")
                        #    println(fnto(constraint_object(con).set))
                        #end
                        val = constraint_object(con).set.upper
                        h[sid, contoidx[idx]] = val
                        #println("$(idx) Add ", S, " to initialization.")
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

function make_two_stage_setup_L(subproblem_generator, v_dict, N, probs = 1/N*ones(N), store = nothing)
    
    models = Dict()

    for i = 1:N
        models[i] = subproblem_generator(i);
    end

    contoidx = ConToIdx(models[1])
    
    h = compute_h(models, contoidx)
    
    subprob = Dict()

    for i = 1:N
        model, varstructs, linkedcons, vnametoidx, arrays = initialize(models[i], v_dict)

        subprob[i] = Subproblems(i, model, probs[i], varstructs, linkedcons, vnametoidx, arrays, nothing)

    end

    firststagevars = Dict()

    for index in 1:length(v_dict[1])
        var = v_dict[1][index]
        # very temporary
        #firststagevars[var[1]] = FirstStageVariableInfo(var[1], index, x_init[1], var[4], var[2], var[3], nothing)
        firststagevars[var[1]] = FirstStageVariableInfo(var[1], index, var[4], nothing, var[2], var[3], nothing)

    end

    firststage = FirstStageInfo(firststagevars, subprob, store);
    
    #temporary? ideally this would be when the second stage is made.
    update_second_index!(firststage)
    
    return firststage, contoidx, h
    
end


function make_two_stage_setup_L_new(subproblem_generator, v_dict, N, probs = 1/N*ones(N), store = nothing)
    
    subprob = Dict()

    for i = 1:N
        model = subproblem_generator(i);
        
        idxtocon = IdxToCon(model)
    
        h = compute_h_new(model, idxtocon)
    
        model, varstructs, linkedcons, vnametoidx, arrays = initialize(model, v_dict)

        #add descriptions to these.
        subprob[i] = SubproblemsNew(i, model, probs[i], varstructs, linkedcons, idxtocon, 
                                            h, nothing, nothing, vnametoidx, arrays, nothing)

    end

    firststagevars = Dict()

    for index in 1:length(v_dict[1])
        var = v_dict[1][index]
        # very temporary
        #firststagevars[var[1]] = FirstStageVariableInfo(var[1], index, x_init[1], var[4], var[2], var[3], nothing)
        firststagevars[var[1]] = FirstStageVariableInfo(var[1], index, var[4], nothing, var[2], var[3], nothing)

    end

    firststage = FirstStageInfo(firststagevars, subprob, store);
    
    #temporary? ideally this would be when the second stage is made.
    update_second_index!(firststage)
    
    return firststage
    
end
