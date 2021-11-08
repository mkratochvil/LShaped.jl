function make_VariableInfo(vardict::Dict{Int64}{Dict{String,Int64}})
    
    varstructs = Dict{Int64,LocalVariableInfo}()
    
    for vname in keys(vardict[1])
        vid = vardict[1][vname]
        varstructs[vid] = LocalVariableInfo(vid, vname, nothing, 0.0, nothing)
    end

    return varstructs
end

function stage_name_idx(model::JuMP.Model, vardict::Dict{Int64,Array{Any}})
        
    vdict = Dict{Int64}{Dict{String,Int64}}()
    vdict[1]=Dict{String,Int64}()
    vdict[2]=Dict{String,Int64}()
    
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

function initialize(model::JuMP.Model, vnames::Dict{Int64,Array{Any}})
        
    println("...making stage_name_idx...")
    vardict = stage_name_idx(model, vnames)
    
    println("...making variable_info...") #TODO fix the bottleneck here. 
    varstructs = make_VariableInfo(vardict)
    
    return model, varstructs, vardict[1]
    
end

function update_second_index!(firststage::FirstStageInfo)
    
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
# this should be true for my problems, but still will not use when in parallel.
function ConToIdx(m::JuMP.Model)
    
    contoidx = Dict{Any,Int64}()

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
    
    return contoidx
    
end

function IdxToCon(m::JuMP.Model)
    
    idxtocon = Dict{Int64,JuMP.ConstraintRef}()

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

function compute_h_new(model::JuMP.Model,idxtocon::Dict{Int64,JuMP.ConstraintRef})

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

function compute_h(models::Dict{Int64,JuMP.Model}, contoidx::Dict{Any,Int64})
    
    ns = models.count
    nc = contoidx.count
    
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

function make_two_stage_setup_L(subproblem_generator::Function, v_dict::Dict{Int64,Array{Any}}, N::Int64, probs::Array{Float64}, store::Union{String,Nothing}, verbose::Int64)
    
    models = Dict{Int64,JuMP.Model}()

    for i = 1:N
        models[i] = subproblem_generator(i);
    end
    
    contoidx = ConToIdx(models[1])
    
    h = compute_h(models, contoidx)
    println(typeof(h))
    
    subprob = Dict()

    for i = 1:N
        
        model, varstructs, vnametoidx = initialize(models[i], v_dict)
                
        update_constraint_values_nac!(model, varstructs)
        
        subprob[i] = Subproblems(i, model, probs[i], varstructs, vnametoidx, nothing)

    end
        
    firststagevars = Dict()

    for index in 1:length(v_dict[1])
        var = v_dict[1][index]
        firststagevars[var[1]] = FirstStageVariableInfo(var[1], index, var[4], nothing, var[2], var[3])
    end

    firststage = FirstStageInfo(firststagevars, subprob, store);
    
    update_second_index!(firststage)
    
    return firststage, contoidx, h
    
end


function make_two_stage_setup_L_new(subproblem_generator::Function, v_dict::Dict{Int64,Array{Any}}, N::Int64, probs::Array{Float64}, store::Union{String,Nothing}, verbose::Int64)
    
    subprob = Dict()

    for i = 1:N
        model = subproblem_generator(i);
        
        idxtocon = IdxToCon(model)
    
        println("computing h for subproblem $(i)...")
        h = compute_h_new(model, idxtocon)
    
        println("Initializing subproblem $(i)...")
        model, varstructs, vnametoidx = initialize(model, v_dict)

        println("Creating subprob[$(i)] struct...")
        subprob[i] = SubproblemsNew(i, model, probs[i], varstructs, idxtocon, h, nothing, nothing, vnametoidx, nothing)

    end

    firststagevars = Dict()

    for index in 1:length(v_dict[1])
        var = v_dict[1][index]
        firststagevars[var[1]] = FirstStageVariableInfo(var[1], index, var[4], nothing, var[2], var[3])
    end

    firststage = FirstStageInfo(firststagevars, subprob, store);
    
    update_second_index!(firststage)
    
    return firststage
    
end
