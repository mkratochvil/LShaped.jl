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
        varstructs[vid] = LocalVariableInfo(vid, vname, nothing, 0.0, nothing, 0.0, nothing)
    end

    return varstructs
end


function stage_name_idx(model::JuMP.Model, vardict)
        
    vdict = Dict()
    vdict[1]=Dict()
    vdict[2]=Dict()
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
        end
    end
    return vdict
end

function initialize(model, vnames)
        
    println("...making stage_name_idx...")
    vardict = stage_name_idx(model, vnames)
    
    println("...Making dicts and arrays...")
    vrm1, vindtoref, varcost, condict = make_dicts_and_arrays(model);
    
    println("...making variable_info...") #TODO fix the bottleneck here. 
    varstructs = make_VariableInfo(vardict, varcost, condict);
    
    return model, varstructs, vardict[1], 0
    
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
        innercount = 0
        for con in all_constraints(m,F,S)
            #exists to make sure nac constraints do not get added.
           if occursin("AffExpr",string(F))
                innercount +=1
                count += 1
                contoidx[(F,S,innercount)] = count
           end
        end
    end
    println("contoidx og count = $(count).")
    
    return contoidx, count
    
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

function compute_h(models, contoidx, count)
    
    ns = models.count
    nc = count
    
    h = Array{Float64}(undef, ns, nc)
    
    for sid in keys(models)
        m = models[sid]
        
        for (F,S) in list_of_constraint_types(m)
            innercount = 0
            for con in all_constraints(m,F,S)
                innercount += 1
               if occursin("AffExpr",string(F))
                    idx = con.index.value
                    if occursin("EqualTo", string(S))
                        val = constraint_object(con).set.value
                        h[sid, contoidx[(F,S,innercount)]] = val
                    elseif occursin("GreaterThan", string(S))
                        val = constraint_object(con).set.lower
                        h[sid, contoidx[(F,S,innercount)]] = val
                    elseif occursin("LessThan", string(S))
                        val = constraint_object(con).set.upper
                        h[sid, contoidx[(F,S,innercount)]] = val
                    elseif occursin("Interval", string(S))
                        val = constraint_object(con).set.upper
                        h[sid, contoidx[(F,S,innercount)]] = val
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

function make_two_stage_setup_L(subproblem_generator, v_dict, N, probs, store, verbose)
    
    models = Dict()

    for i = 1:N
        models[i] = subproblem_generator(i);
    end
    
    contoidx, count = ConToIdx(models[1])
    
    h = compute_h(models, contoidx, count)
    
    subprob = Dict()

    for i = 1:N
        
        model, varstructs, vnametoidx, arrays = initialize(models[i], v_dict)
                
        update_constraint_values_nac!(model, varstructs)
        
        subprob[i] = Subproblems(i, model, probs[i], varstructs, count, nothing, vnametoidx, arrays, nothing)

    end
    
    #contoidx, count = ConToIdx(subprob[1].model)
    
    #h = compute_h(models, contoidx, count)

    firststagevars = Dict()

    for index in 1:length(v_dict[1])
        var = v_dict[1][index]
        firststagevars[var[1]] = FirstStageVariableInfo(var[1], index, var[4], nothing, var[2], var[3])
    end

    firststage = FirstStageInfo(firststagevars, subprob, store);
    
    #temporary? ideally this would be when the second stage is made.
    update_second_index!(firststage)
    
    return firststage, contoidx, h
    
end


function make_two_stage_setup_L_new(subproblem_generator, v_dict, N, probs, store, verbose)
    
    subprob = Dict()

    for i = 1:N
        model = subproblem_generator(i);
        
        idxtocon = IdxToCon(model)
    
        println("computing h for subproblem $(i)...")
        h = compute_h_new(model, idxtocon)
    
        println("Initializing subproblem $(i)...")
        model, varstructs, vnametoidx, arrays = initialize(model, v_dict)

        #add descriptions to these.
        println("Creating subprob[$(i)] struct...")
        subprob[i] = SubproblemsNew(i, model, probs[i], varstructs, nothing, idxtocon, 
                                            h, nothing, nothing, vnametoidx, arrays, nothing)

    end

    firststagevars = Dict()

    for index in 1:length(v_dict[1])
        var = v_dict[1][index]
        firststagevars[var[1]] = FirstStageVariableInfo(var[1], index, var[4], nothing, var[2], var[3])
    end

    firststage = FirstStageInfo(firststagevars, subprob, store);
    
    #temporary? ideally this would be when the second stage is made.
    update_second_index!(firststage)
    
    return firststage
    
end
