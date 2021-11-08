using JuMP
using Gurobi
using CSV
using DataFrames
using OrderedCollections

using LShaped

function create_model_var_obj(dfvar)
    
    model = Model(with_optimizer(Gurobi.Optimizer; OutputFlag = 0))

    cexpr = JuMP.GenericAffExpr(0.0, OrderedDict{VariableRef,Float64}())
    for i = 1:size(dfvar)[1]
        
        if i % 1000 == 0
            println("Creating variable $(i)...")
        end

        vname31 = dfvar[i,2]
        vlb = dfvar[i,3]
        vub = dfvar[i,4]
        cost = dfvar[i,5]
        
	#attempt to change vname::String31 to String
	vname = ""
	for s = 1:length(vname31)
	    vname = string(vname,vname31[s])
	end

	#println("$(vname), $(typeof(vname))")

        if vlb > -Inf
            if vub < Inf
                vref = @variable(model, base_name = vname, lower_bound = vlb, upper_bound = vub)
            else
                vref = @variable(model, base_name = vname, lower_bound = vlb)
            end
        else
            if vub < Inf
                vref = @variable(model, base_name = vname, upper_bound = vub)
            else
                vref = @variable(model, base_name = string(vname))
            end
        end

        add_to_expression!(cexpr, cost, vref)

    end

    @objective(model, Min, cexpr)
    
    return model
    
end

function add_model_cons!(model, dfvar, dfcon, dfmat)
    
    #println("Creating Sparse Matrix...")
    #con_mat = sparse(dfmat[:,1], dfmat[:,2], dfmat[:,3])
    
    count = 1
    
    for connum = 1:size(dfcon)[1]
        
        if connum % 1000 == 0
            println("Creating Constraint $(connum)...")
        end
        
        conexpr = JuMP.GenericAffExpr(0.0, OrderedDict{VariableRef,Float64}())
        
        varids = Int64[];
        vals = Float64[];
                
        t1 = time_ns()
        for i = count:size(dfmat,1)

            # temporary. probably should do something like "find nonzeros for this row"
            #(varids, vals) = findnz(con_mat[rownum,:])
            
            if connum == dfmat[i,1]
                push!(varids, dfmat[i,2])
                push!(vals, dfmat[i,3])
            else
                count = i
                break
            end
            
        end
        #=
        if connum < 20 || connum == size(dfcon)[1]
            println("varids = $(varids)")
            println("vals = $(vals)")
            println("count = $(count)")
        end
        =#
        t2 = time_ns()
            
            
        for  i = 1:size(varids,1)

            coeff = vals[i] #dfcon[rownum,colnum]

	    vname_old = dfvar[varids[i],2]
	    vname = ""
	    for s = 1:length(vname_old)
		vname = string(vname,vname_old[s])
	    end

            vref = variable_by_name(model, vname)

            add_to_expression!(conexpr, coeff, vref)

        end
        t3 = time_ns()

        contype = dfcon[connum, 3]

        if occursin("EqualTo", contype)

            value = dfcon[connum,4]

            @constraint(model, conexpr == value)

        elseif occursin("GreaterThan", contype)

            value = dfcon[connum,4]

            @constraint(model, conexpr >= value)

         elseif occursin("LessThan", contype)

            value = dfcon[connum,4]

            @constraint(model, conexpr <= value)

         elseif occursin("Interval", contype)

            lb = dfcon[connum,4]
            ub = dfcon[connum,5]

            @constraint(model, lb <= conexpr <= ub)

         else 

            println("constraint missing from model. This is not good.")

        end
            t4 = time_ns()

            if connum % 1000 == 0
                nztime = (t2-t1)/1e9
                exptime = (t3-t2)/1e9
                contime = (t4-t3)/1e9
                tottime = (t4-t1)/1e9
                println("For this iteration: nzt =$(nztime), expt = $(exptime), cont = $(contime), tot = $(tottime)")
            end
        end
    
    return
    
end


# creates the model whose info is stored in "dir/scenid"
function wsgep2(scenid)
    
    dir = "../ScenarioPrimal/ScenarioPrimal/jumpmodels/wsgep/twostage/ts3/"
    
    varfile = string(dir, scenid, "/vars_eff.csv")
    confile = string(dir, scenid, "/cons_eff.csv")
    matfile = string(dir, scenid, "/mat_eff.csv")

    println("Loading in variable dataframe...")
    dfvar = DataFrame(CSV.File(varfile));
    println("Loading in constraint dataframe...")
    dfcon = DataFrame(CSV.File(confile));
    println("Loading in constraint matrix...")
    dfmat = DataFrame(CSV.File(matfile));
    
    println("Creating model...")
    model = create_model_var_obj(dfvar);
    
    add_model_cons!(model, dfvar, dfcon, dfmat)
    
    return model
    
end


dir = "../ScenarioPrimal/ScenarioPrimal/jumpmodels/wsgep/twostage/ts3/"

scenid = 1

varfile = string(dir, scenid, "/vars_eff.csv")
#confile = string(dir, scenid, "/cons.csv")

dfvar = DataFrame(CSV.File(varfile));
#dfcon = DataFrame(CSV.File(confile));

dfvar[size(dfvar,1),1]

x_init = ones(150);

function variable_dict_primal(varlist, x)
    vdict = Dict{Int64,Array{Any}}()
    vdict[1]=[]
    vdict[2]=[]
    i = 0
    for name in varlist
        vname = ""
        for s = 1:length(name)
            vname = string(vname,name[s])
        end
        if occursin("m_EE", name) || occursin("r_WE", name)
            
            i+=1
            
            # for some reason, the lower bound on these variables is not showing...I'm not sure why.
            lb = 0
            ub = 1000.0
            #if has_lower_bound(var) == 1
            #    lb = lower_bound(var)
            #end
            #if has_upper_bound(var) == 1
            #    ub = upper_bound(var)
            #end
            

            push!(vdict[1], (vname, lb, ub, x[i]))
        else
            push!(vdict[2], vname)
        end
    end
    return vdict
end

wsgepv = variable_dict_primal(dfvar[:,2], x_init);

function wsgep1()
    
    names = ["r_WE[101_WIND_2]"
"r_WE[102_WIND_2]"
"r_WE[103_WIND_2]"
"r_WE[104_WIND_2]"
"r_WE[105_WIND_2]"
"r_WE[106_WIND_2]"
"r_WE[107_WIND_2]"
"r_WE[108_WIND_2]"
"r_WE[109_WIND_2]"
"r_WE[110_WIND_2]"
"r_WE[111_WIND_2]"
"r_WE[112_WIND_2]"
"r_WE[113_WIND_2]"
"r_WE[114_WIND_2]"
"r_WE[115_WIND_2]"
"r_WE[116_WIND_2]"
"r_WE[117_WIND_2]"
"r_WE[118_WIND_2]"
"r_WE[119_WIND_2]"
"r_WE[120_WIND_2]"
"r_WE[121_WIND_2]"
"r_WE[122_WIND_1]"
"r_WE[122_WIND_2]"
"r_WE[123_WIND_2]"
"r_WE[124_WIND_2]"
"r_WE[201_WIND_2]"
"r_WE[202_WIND_2]"
"r_WE[203_WIND_2]"
"r_WE[204_WIND_2]"
"r_WE[205_WIND_2]"
"r_WE[206_WIND_2]"
"r_WE[207_WIND_2]"
"r_WE[208_WIND_2]"
"r_WE[209_WIND_2]"
"r_WE[210_WIND_2]"
"r_WE[211_WIND_2]"
"r_WE[212_WIND_2]"
"r_WE[213_WIND_2]"
"r_WE[214_WIND_2]"
"r_WE[215_WIND_2]"
"r_WE[216_WIND_2]"
"r_WE[217_WIND_2]"
"r_WE[218_WIND_2]"
"r_WE[219_WIND_2]"
"r_WE[220_WIND_2]"
"r_WE[221_WIND_2]"
"r_WE[222_WIND_2]"
"r_WE[223_WIND_2]"
"r_WE[224_WIND_2]"
"r_WE[301_WIND_2]"
"r_WE[302_WIND_2]"
"r_WE[303_WIND_1]"
"r_WE[303_WIND_2]"
"r_WE[304_WIND_2]"
"r_WE[305_WIND_2]"
"r_WE[306_WIND_2]"
"r_WE[307_WIND_2]"
"r_WE[308_WIND_2]"
"r_WE[309_WIND_1]"
"r_WE[309_WIND_2]"
"r_WE[310_WIND_2]"
"r_WE[311_WIND_2]"
"r_WE[312_WIND_2]"
"r_WE[313_WIND_2]"
"r_WE[314_WIND_2]"
"r_WE[315_WIND_2]"
"r_WE[316_WIND_2]"
"r_WE[317_WIND_1]"
"r_WE[317_WIND_2]"
"r_WE[318_WIND_2]"
"r_WE[319_WIND_2]"
"r_WE[320_WIND_2]"
"r_WE[321_WIND_2]"
"r_WE[322_WIND_2]"
"r_WE[323_WIND_2]"
"r_WE[324_WIND_2]"
"r_WE[325_WIND_2]"
"m_EE[101_STORAGE_1]"
"m_EE[102_STORAGE_1]"
"m_EE[103_STORAGE_1]"
"m_EE[104_STORAGE_1]"
"m_EE[105_STORAGE_1]"
"m_EE[106_STORAGE_1]"
"m_EE[107_STORAGE_1]"
"m_EE[108_STORAGE_1]"
"m_EE[109_STORAGE_1]"
"m_EE[110_STORAGE_1]"
"m_EE[111_STORAGE_1]"
"m_EE[112_STORAGE_1]"
"m_EE[113_STORAGE_1]"
"m_EE[114_STORAGE_1]"
"m_EE[115_STORAGE_1]"
"m_EE[116_STORAGE_1]"
"m_EE[117_STORAGE_1]"
"m_EE[118_STORAGE_1]"
"m_EE[119_STORAGE_1]"
"m_EE[120_STORAGE_1]"
"m_EE[121_STORAGE_1]"
"m_EE[122_STORAGE_1]"
"m_EE[123_STORAGE_1]"
"m_EE[124_STORAGE_1]"
"m_EE[201_STORAGE_1]"
"m_EE[202_STORAGE_1]"
"m_EE[203_STORAGE_1]"
"m_EE[204_STORAGE_1]"
"m_EE[205_STORAGE_1]"
"m_EE[206_STORAGE_1]"
"m_EE[207_STORAGE_1]"
"m_EE[208_STORAGE_1]"
"m_EE[209_STORAGE_1]"
"m_EE[210_STORAGE_1]"
"m_EE[211_STORAGE_1]"
"m_EE[212_STORAGE_1]"
"m_EE[213_STORAGE_1]"
"m_EE[214_STORAGE_1]"
"m_EE[215_STORAGE_1]"
"m_EE[216_STORAGE_1]"
"m_EE[217_STORAGE_1]"
"m_EE[218_STORAGE_1]"
"m_EE[219_STORAGE_1]"
"m_EE[220_STORAGE_1]"
"m_EE[221_STORAGE_1]"
"m_EE[222_STORAGE_1]"
"m_EE[223_STORAGE_1]"
"m_EE[224_STORAGE_1]"
"m_EE[301_STORAGE_1]"
"m_EE[302_STORAGE_1]"
"m_EE[303_STORAGE_1]"
"m_EE[304_STORAGE_1]"
"m_EE[305_STORAGE_1]"
"m_EE[306_STORAGE_1]"
"m_EE[307_STORAGE_1]"
"m_EE[308_STORAGE_1]"
"m_EE[309_STORAGE_1]"
"m_EE[310_STORAGE_1]"
"m_EE[311_STORAGE_1]"
"m_EE[312_STORAGE_1]"
"m_EE[313_STORAGE_1]"
"m_EE[314_STORAGE_1]"
"m_EE[315_STORAGE_1]"
"m_EE[316_STORAGE_1]"
"m_EE[317_STORAGE_1]"
"m_EE[318_STORAGE_1]"
"m_EE[319_STORAGE_1]"
"m_EE[320_STORAGE_1]"
"m_EE[321_STORAGE_1]"
"m_EE[322_STORAGE_1]"
"m_EE[323_STORAGE_1]"
"m_EE[324_STORAGE_1]"
"m_EE[325_STORAGE_1]"
];

    fs = Model(with_optimizer(Gurobi.Optimizer, OutputFlag = 1));

    x = Vector{VariableRef}()
    for i = 1:150
        push!(x, @variable(fs, base_name = names[i], lower_bound = 0.0, upper_bound = 1000.0))
    end
    
    @objective(fs, Min, sum(0.9040815000000001*x[i] for i = 1:77)+sum(1.1059811266666666*x[i] for i = 78:150))
    #@objective(fs, Min, 0.0)
        
    return fs
    
end


