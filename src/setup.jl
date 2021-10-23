using JuMP
using LinearAlgebra 
using MathOptInterface

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
#=
function conIdxToRef(m)
    
    conidxtoref = Dict()
    moi_backend = backend(m)
    
    for (F,S) in list_of_constraint_types(m)
        for con in all_constraints(m,F,S)
            f = MOI.get(moi_backend, MOI.ConstraintFunction(), con.index)
            conidxtoref[con.index.value] = (F,S,f,con)
        end
    end
    
    return conidxtoref
    
end
=#

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

function create_eq_lt_gt_arrays(m,conidxtoref,vr)

    #nvar = num_variables(m) # fails at index 11 since there are now 10 variables. make max value for keys?
    nvar = maximum(keys(vr))

    Aeq = Array{Float64}(undef, 0, nvar)
    Ieq = Dict()
    counteq = 0

    Ageq = Array{Float64}(undef, 0, nvar)
    Igeq = Dict()
    countgeq = 0
    
    Aleq = Array{Float64}(undef, 0, nvar)
    Ileq = Dict()
    countleq = 0


    for i in keys(conidxtoref)
        row = zeros(nvar)
        if(occursin("LessThan",string(conidxtoref[i][2]))) #find something better than this.
            if(occursin("ScalarAffine", string(conidxtoref[i][3])))
                for term in conidxtoref[i][3].terms
                    row[term.variable_index.value] = term.coefficient
                end
            else
                row[conidxtoref[i][3].variable.value] = 1.0;
            end
            countleq += 1
            Ileq[countleq] = i
            Aleq = [Aleq; row']
        end
        if(occursin("GreaterThan",string(conidxtoref[i][2]))) #find something better than this.
            if(occursin("ScalarAffine", string(conidxtoref[i][3])))
                for term in conidxtoref[i][3].terms
                    row[term.variable_index.value] = term.coefficient
                end
            else
                row[conidxtoref[i][3].variable.value] = 1.0;
            end
            countgeq += 1
            Igeq[countgeq] = i
            Ageq = [Ageq; row']
        end
        if(occursin("EqualTo",string(conidxtoref[i][2]))) #find something better than this.
            if(occursin("ScalarAffine", string(conidxtoref[i][3])))
                for term in conidxtoref[i][3].terms
                    row[term.variable_index.value] = term.coefficient
                end
            else
                row[conidxtoref[i][3].variable.value] = 1.0;
            end
            counteq += 1
            Ieq[counteq] = i
            Aeq = [Aeq; row']
        end
    end
    
    return Aeq, Aleq, Ageq, Ieq, Ileq, Igeq
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
# make link
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

function make_arrays(subproblem)

    m = subproblem.model
    lci = subproblem.linkedconstraintinfo
    nvar1 = subproblem.variableinfo.count
    nvar2 = num_variables(m)

    reindexvar1 = Dict()
    reindexvar2 = Dict()
    
    count = 0
    for i in keys(subproblem.variableinfo)
        count += 1
        reindexvar1[i] = count
    end

    count = 0
    allvar = all_variables(m)
    for i in 1:nvar2
        count += 1
        reindexvar2[allvar[i].index.value] = count
    end

    
    Te = Array{Float64}(undef, 0, nvar1)
    We = Array{Float64}(undef, 0, nvar2)
    Ce = Dict()
    counteq = 0

    Tg = Array{Float64}(undef, 0, nvar1)
    Wg = Array{Float64}(undef, 0, nvar2)
    Cg = Dict()
    countgeq = 0

    Tl = Array{Float64}(undef, 0, nvar1)
    Wl = Array{Float64}(undef, 0, nvar2)
    Cl = Dict()
    countleq = 0

    for i in keys(lci)
        row1 = zeros(nvar1)
        row2 = zeros(nvar2)

        if(occursin("EqualTo",string(lci[i].set))) #find something better than this.
            
            for vnum in lci[i].variables
                coeff = subprob.variableinfo[vnum].conval[lci[i].id]
                row1[reindexvar1[vnum]] = coeff
            end
            
            Te = [Te; row1']
            
            entries = MOI.get(backend(m), MOI.ConstraintFunction(), lci[i].ref.index).terms
            for entry in entries
                row2[reindexvar2[entry.variable_index.value]] = entry.coefficient
            end
            We = [We; row2']
            
            counteq += 1
            Ce[counteq] = i
            
        elseif(occursin("LessThan",string(lci[i].set))) #find something better than this.
            
            for vnum in lci[i].variables
                coeff = subprob.variableinfo[vnum].conval[lci[i].id]
                row1[reindexvar1[vnum]] = coeff
            end
            
            Tl = [Tl; row1']
            
            entries = MOI.get(backend(m), MOI.ConstraintFunction(), lci[i].ref.index).terms
            for entry in entries
                row2[reindexvar2[entry.variable_index.value]] = entry.coefficient
            end
            Wl = [Wl; row2']
            
            countleq += 1
            Cl[countleq] = i
        elseif(occursin("GreaterThan",string(lci[i].set))) #find something better than this.
            
            for vnum in lci[i].variables
                coeff = subprob.variableinfo[vnum].conval[lci[i].id]
                row1[reindexvar1[vnum]] = coeff
            end
            
            Tg = [Tg; row1']
            
            entries = MOI.get(backend(m), MOI.ConstraintFunction(), lci[i].ref.index).terms
            for entry in entries
                row2[reindexvar2[entry.variable_index.value]] = entry.coefficient
            end
            Wg = [Wg; row2']
            
            countgeq += 1
            Cg[countgeq] = i
        end
    end

    return Te, We, Ce, Tl, Wl, Cl, Tg, Wl, Cg

end

# make_dict_and_arrays needs to be updated as varcost and condict is all that is needed.
### x_init not needed apparently
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

# make into setup functions

function make_two_stage_setup(subproblem_generator, v_dict, N)

    models = Dict()

    for i = 1:N
        models[i] = subproblem_generator(i);
    end

    subprob = Dict()

    for i = 1:N
        model, varstructs, linkedcons, vnametoidx, arrays = initialize(models[i], v_dict)

        subprob[i] = Subproblems(i, model, 1/N, varstructs, linkedcons, vnametoidx, arrays, nothing)

    end

    firststagevars = Dict()

    for index in 1:length(v_dict[1])
        var = v_dict[1][index]
        # very temporary
        #firststagevars[var[1]] = FirstStageVariableInfo(var[1], index, x_init[1], var[4], var[2], var[3], nothing)
        firststagevars[var[1]] = FirstStageVariableInfo(var[1], index, var[4], nothing, var[2], var[3], nothing)

    end

    firststage = FirstStageInfo(firststagevars, subprob);
    
    #temporary? ideally this would be when the second stage is made.
    update_second_index!(firststage)
    
    return firststage
    
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

function compute_h(models, contoidx)
# make link
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
                        #println("$(idx) Add ", S, " to fsasdinitialization.")
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

function make_two_stage_setup_L(subproblem_generator, v_dict, N, probs = 1/N*ones(N))

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

    firststage = FirstStageInfo(firststagevars, subprob);
    
    #temporary? ideally this would be when the second stage is made.
    update_second_index!(firststage)
    
    return firststage, contoidx, h
    
end
