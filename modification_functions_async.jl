# load in helper functions
#include("./parameters.jl")
#include("./get_functions.jl")

# load distribution factors for RTS-GMLC Data Set csv file
#loadcsv = CSV.File("./LOAD.csv");

# ptdf matrix (calculated elsewhere) for Zone A of RTS-GMLC Data Set.
#ptdfdf = DataFrame(CSV.File("./ptdfsmall.csv"));


function load_distribution_dict(loadcsv; size="small")
    
    loaddis = Dict{Int64,Float64}()

    if size == "small"
        for i = 1:24
            loaddis[loadcsv[i][1]] = loadcsv[i][2]
        end
    elseif size  == "large"
        for i = 1:73
            loaddis[loadcsv[i][1]] = loadcsv[i][2]
        end
    else
        println("Only set up for 24 bus (size = small) or 73 bus (size = large) cases.")
        return nothing
    end
    
    return loaddis
    
end

#note this is currently only set up for the small case.
function ptdf_dict(ptdfdf)
    
    ptdfdict = Dict()

    for i = 1:38
        br = ptdfdf[i,1]
        ptdfdict[br] = Dict()
        for j = 2:25
            bus = parse(Int64,names(ptdfdf)[j])
            ptdfdict[br][bus] = ptdfdf[i,j]
        end
    end
    
    return ptdfdict
    
end

function load_in_models()
    
    return JuMP.Model(), JuMP.read_from_file("./new_store_exp_Z1_newptdf_newobj.mps")
    
end

function add_firststage_variables!(exmodel::JuMP.Model, submodel::JuMP.Model)
    
    #add first stage variables to extensive form

    for bus in buses        
        subnamePR = JuMP.name(get_PR_variable(submodel, bus))
        subnameER = JuMP.name(get_ER_variable(submodel, bus))

        JuMP.@variable(exmodel, base_name = string(subnamePR,"_0"))
        JuMP.@variable(exmodel, base_name = string(subnameER,"_0"))
    end
    
    return
    
end

function add_ER_lb!(exmodel::JuMP.Model, submodel::JuMP.Model)
    
    outercount = 0
    for bus in buses
        outercount += 1
        subcon = get_ER_lb(submodel, bus)
        subconname = JuMP.name(subcon)
        exconname = string(subconname, "_0")

        terms = JuMP.constraint_object(subcon).func.terms
        value = JuMP.constraint_object(subcon).set.lower

        innercount = 0
        for (subvar, coeff) in terms
            exvar = JuMP.variable_by_name(exmodel, string(JuMP.name(subvar),"_0"))
            if innercount == 0
                JuMP.@constraint(exmodel, coeff*exvar >= value, 
                                            base_name = exconname)
                innercount += 1
            else
                excon = JuMP.constraint_by_name(exmodel, exconname)
                JuMP.set_normalized_coefficient(excon, exvar, coeff)
            end
        end
    end

    return
    
end

function add_PR_lb!(exmodel::JuMP.Model, submodel::JuMP.Model)
        
    #add PR_lb constraint
    outercount = 0
    for bus in buses
        outercount += 1
        subcon = get_PR_lb(submodel, bus)
        subconname = JuMP.name(subcon)
        exconname = string(subconname, "_0")

        terms = JuMP.constraint_object(subcon).func.terms
        value = JuMP.constraint_object(subcon).set.lower

        innercount = 0
        for (subvar, coeff) in terms
            exvar = JuMP.variable_by_name(exmodel, string(JuMP.name(subvar),"_0"))
            if innercount == 0
                JuMP.@constraint(exmodel, coeff*exvar >= value, 
                                            base_name = exconname)
                innercount += 1
            else
                excon = JuMP.constraint_by_name(exmodel, exconname)
                JuMP.set_normalized_coefficient(excon, exvar, coeff)
            end
        end
    end
    
    return
    
end

#add ER_ub constraint
function add_ER_ub!(exmodel::JuMP.Model, submodel::JuMP.Model)
    outercount = 0
    for bus in buses
        outercount += 1
        subcon = get_ER_ub(submodel, bus)
        subconname = JuMP.name(subcon)
        exconname = string(subconname, "_0")

        terms = JuMP.constraint_object(subcon).func.terms
        value = JuMP.constraint_object(subcon).set.upper

        innercount = 0
        for (subvar, coeff) in terms
            exvar = JuMP.variable_by_name(exmodel, string(JuMP.name(subvar),"_0"))
            if innercount == 0
                JuMP.@constraint(exmodel, coeff*exvar <= value, 
                                            base_name = exconname)
                innercount += 1
            else
                excon = JuMP.constraint_by_name(exmodel, exconname)
                JuMP.set_normalized_coefficient(excon, exvar, coeff)
            end
        end
    end
    
    return
end

#add expansion budget
function add_expansion_budget!(exmodel::JuMP.Model, submodel::JuMP.Model)
    
    subcon = JuMP.constraint_by_name(submodel, "expansion_budget")
    subconname = JuMP.name(subcon)
    exconname = string(subconname, "_0")

    terms = JuMP.constraint_object(subcon).func.terms
    value = JuMP.constraint_object(subcon).set.upper

    innercount = 0
    for (subvar, coeff) in terms
        exvar = JuMP.variable_by_name(exmodel, string(JuMP.name(subvar),"_0"))
        if innercount == 0
            JuMP.@constraint(exmodel, coeff*exvar <= value, 
                                        base_name = exconname)
            innercount += 1
        else
            excon = JuMP.constraint_by_name(exmodel, exconname)
            JuMP.set_normalized_coefficient(excon, exvar, coeff)
        end
    end
    
    return
end


### LOOPING FUNCTIONS ###

#adjust load balance values (remember to change load)
function adjust_load_balance_constraint!(submodel::JuMP.Model, loadvec, loaddis)

    for bus in buses
        lf = loaddis[bus]
        for ts in timesteps
            con = get_load_balance(submodel, bus, ts)
            oldval = JuMP.constraint_object(con).set.value
            lval = loadvec[ts]
            JuMP.set_normalized_rhs(con, lf*lval)
            newval = JuMP.constraint_object(con).set.value
            #println("$(name(con)), $(oldval), $(newval)")
        end
    end
    
    return
end


# change ptdf constraint (remember to run load changes FIRST)
function adjust_ptdf_constraint!(exmodel::JuMP.Model, submodel::JuMP.Model, ptdfdict, loadvec, loaddis)
    for ts in timesteps
        for br in branches
            ptdfcon = get_ptdf_con(submodel,br,ts)

            valold = JuMP.constraint_object(ptdfcon).set.value
            valnew = 0.0
            
            for bus in buses
                lf = loaddis[bus]
                lval = loadvec[ts]
                valnew -= ptdfdict[br][bus]*lf*lval
            end
            #=
            for bus in buses
                buscon = get_load_balance(submodel,bus,ts)

                #loadcon = copy(JuMP.constraint_object(buscon).func)
                loadval = copy(JuMP.constraint_object(buscon).set.value)

                valnew -= ptdfdict[br][bus]*loadval

            end 
            =#

            JuMP.set_normalized_rhs(ptdfcon, valnew)
            #println("$(JuMP.name(ptdfcon)), $(valold), $(valnew)")

        end
    end
    
    return
    
end

#adjust wind ub values (remember to change wind)
function adjust_wind_ub!(submodel::JuMP.Model, windvec)
    bus = 122
    for ts in timesteps
        con = get_wind_ub(submodel, bus, ts)
        oldval = JuMP.constraint_object(con).set.upper
        wval = windvec[ts]
        JuMP.set_normalized_rhs(con, wval)
        newval = JuMP.constraint_object(con).set.upper
        #println("$(name(con)), $(oldval), $(newval)")
    end
    return
end

# add thermal vars to extensive form
function add_thermal_variables!(exmodel::JuMP.Model, submodel::JuMP.Model, scen::Int64)
    for ts in timesteps
        for gen in gens
            name = JuMP.name(get_thermal_variable(submodel, gen, ts))
            JuMP.@variable(exmodel, base_name = string(name, "_$(scen)"))
        end
    end
    return
end

#add branch vars to extensive form
function add_branch_variables!(exmodel::JuMP.Model, submodel::JuMP.Model, scen::Int64)
    for ts in timesteps

        for br in branches
            name = JuMP.name(get_line_variable(submodel, br, ts))

            JuMP.@variable(exmodel, base_name = string(name, "_$(scen)"))
        end
    end  
    return
end

#add wind vars to extensive form
function add_wind_variables!(exmodel::JuMP.Model, submodel::JuMP.Model, scen::Int64)
    for ts in timesteps

        bus = 122
        name = JuMP.name(get_wind_variable(submodel, bus, ts))

        JuMP.@variable(exmodel, base_name = string(name, "_$(scen)"))
    end  
    return
end

# add charging vars to extensive form
function add_charging_variables!(exmodel::JuMP.Model, submodel::JuMP.Model, scen::Int64)
    for ts in timesteps

        for bus in buses

            name = JuMP.name(get_charging_variable(submodel, bus, ts))

            JuMP.@variable(exmodel, base_name = string(name, "_$(scen)"))
        end
    end
    return
end


# add discharging vars to extensive form
function add_discharging_variables!(exmodel::JuMP.Model, submodel::JuMP.Model, scen::Int64)
    for ts in timesteps

        for bus in buses

            name = JuMP.name(get_discharging_variable(submodel, bus, ts))

            JuMP.@variable(exmodel, base_name = string(name, "_$(scen)"))
        end
    end
    return
end


# add storage vars to extensive form
function add_storage_variables!(exmodel::JuMP.Model, submodel::JuMP.Model, scen::Int64)
    for ts in timesteps

        for bus in buses

            name = JuMP.name(get_stored_variable(submodel, bus, ts))

            JuMP.@variable(exmodel, base_name = string(name, "_$(scen)"))
        end
    end
    return
end



# add lossofload vars to extensive form
function add_lossofload_variables!(exmodel::JuMP.Model, submodel::JuMP.Model, scen::Int64)
    for ts in timesteps

        for bus in buses

            name = JuMP.name(get_lossofload_variable(submodel, bus, ts))

            JuMP.@variable(exmodel, base_name = string(name, "_$(scen)"))
        end
    end
    return
end



# add overload vars to extensive form
function add_overload_variables!(exmodel::JuMP.Model, submodel::JuMP.Model, scen::Int64)
    for ts in timesteps

        for bus in buses

            name = JuMP.name(get_overload_variable(submodel, bus, ts))

            JuMP.@variable(exmodel, base_name = string(name, "_$(scen)"))
        end
    end
    return
end


# add storage balance constraint

function add_storage_balance_constraint!(exmodel::JuMP.Model, submodel::JuMP.Model, scen::Int64)
    for ts in timesteps
        for bus in buses
            ### change constraint here
            subcon = get_storage_balance(submodel, bus, ts)
            ###
            subconname = JuMP.name(subcon)
            exconname = string(subconname, "_$(scen)")

            terms = JuMP.constraint_object(subcon).func.terms
            ### change sense here
            value = JuMP.constraint_object(subcon).set.value
            ###

            count = 0
            for (subvar, coeff) in terms
                if JuMP.name(subvar)[1:2] == "PR" || JuMP.name(subvar)[1:2] == "ER"
                    exvar = JuMP.variable_by_name(exmodel, string(JuMP.name(subvar),"_0"))
                else
                    exvar = JuMP.variable_by_name(exmodel, string(JuMP.name(subvar),"_$(scen)"))
                end
                if count == 0
                    ### change objective type here
                    JuMP.@constraint(exmodel, coeff*exvar == value, 
                                                base_name = exconname)
                    ###
                    count += 1
                else
                    excon = JuMP.constraint_by_name(exmodel, exconname)
                    JuMP.set_normalized_coefficient(excon, exvar, coeff)
                end
            end
            #println(JuMP.constraint_by_name(exmodel, exconname))
        end
    end
    return
end

# add load balance constraint
function add_load_balance_constraint!(exmodel::JuMP.Model, submodel::JuMP.Model, scen::Int64)
    for ts in timesteps
        for bus in buses
            ### change constraint here
            subcon = get_load_balance(submodel, bus, ts)
            ###
            subconname = JuMP.name(subcon)
            exconname = string(subconname, "_$(scen)")

            terms = JuMP.constraint_object(subcon).func.terms
            ### change sense here
            value = JuMP.constraint_object(subcon).set.value
            ###

            count = 0
            for (subvar, coeff) in terms
                if JuMP.name(subvar)[1:2] == "PR" || JuMP.name(subvar)[1:2] == "ER"
                    exvar = JuMP.variable_by_name(exmodel, string(JuMP.name(subvar),"_0"))
                else
                    exvar = JuMP.variable_by_name(exmodel, string(JuMP.name(subvar),"_$(scen)"))
                end
                if count == 0
                    ### change objective type here
                    JuMP.@constraint(exmodel, coeff*exvar == value, 
                                                base_name = exconname)
                    ###
                    count += 1
                else
                    excon = JuMP.constraint_by_name(exmodel, exconname)
                    JuMP.set_normalized_coefficient(excon, exvar, coeff)
                end
            end
            #println(JuMP.constraint_by_name(exmodel, exconname))
        end
    end
    return
end


# add ptdf constraint
function add_ptdf_constraint!(exmodel::JuMP.Model, submodel::JuMP.Model, scen::Int64)
    for ts in timesteps
        for br in branches
            ### change constraint here
            subcon = get_ptdf_con(submodel, br, ts)
            ###
            subconname = JuMP.name(subcon)
            exconname = string(subconname, "_$(scen)")

            terms = JuMP.constraint_object(subcon).func.terms
            ### change sense here
            value = JuMP.constraint_object(subcon).set.value
            ###

            count = 0
            for (subvar, coeff) in terms
                if JuMP.name(subvar)[1:2] == "PR" || JuMP.name(subvar)[1:2] == "ER"
                    exvar = JuMP.variable_by_name(exmodel, string(JuMP.name(subvar),"_0"))
                else
                    exvar = JuMP.variable_by_name(exmodel, string(JuMP.name(subvar),"_$(scen)"))
                end
                if count == 0
                    ### change objective type here
                    JuMP.@constraint(exmodel, coeff*exvar == value, 
                                                base_name = exconname)
                    ###
                    count += 1
                else
                    excon = JuMP.constraint_by_name(exmodel, exconname)
                    JuMP.set_normalized_coefficient(excon, exvar, coeff)
                end
            end
            #println(JuMP.constraint_by_name(exmodel, exconname))
        end
    end
    return
end


# add wind lb constraint
function add_wind_lb_constraint!(exmodel::JuMP.Model, submodel::JuMP.Model, scen::Int64)
    for ts in timesteps
        bus = 122
        ### change constraint here
        subcon = get_wind_lb(submodel, bus, ts)
        ###
        subconname = JuMP.name(subcon)
        exconname = string(subconname, "_$(scen)")

        terms = JuMP.constraint_object(subcon).func.terms
        ### change sense here
        value = JuMP.constraint_object(subcon).set.lower
        ###

        count = 0
        for (subvar, coeff) in terms
            if JuMP.name(subvar)[1:2] == "PR" || JuMP.name(subvar)[1:2] == "ER"
                exvar = JuMP.variable_by_name(exmodel, string(JuMP.name(subvar),"_0"))
            else
                exvar = JuMP.variable_by_name(exmodel, string(JuMP.name(subvar),"_$(scen)"))
            end
            if count == 0
                ### change objective type here
                JuMP.@constraint(exmodel, coeff*exvar >= value, 
                                            base_name = exconname)
                ###
                count += 1
            else
                excon = JuMP.constraint_by_name(exmodel, exconname)
                JuMP.set_normalized_coefficient(excon, exvar, coeff)
            end
        end
        #println(JuMP.constraint_by_name(exmodel, exconname))
    end
    return
end


# add ramping lb constraint
function add_ramping_lb_constraint!(exmodel::JuMP.Model, submodel::JuMP.Model, scen::Int64)
    for ts in timesteps
        if ts > 1
            for gen in gens
                ### change constraint here
                subcon = get_ramp_lb(submodel, gen, ts)
                ###
                subconname = JuMP.name(subcon)
                exconname = string(subconname, "_$(scen)")

                terms = JuMP.constraint_object(subcon).func.terms
                ### change sense here
                value = JuMP.constraint_object(subcon).set.lower
                ###

                count = 0
                for (subvar, coeff) in terms
                    if JuMP.name(subvar)[1:2] == "PR" || JuMP.name(subvar)[1:2] == "ER"
                        exvar = JuMP.variable_by_name(exmodel, string(JuMP.name(subvar),"_0"))
                    else
                        exvar = JuMP.variable_by_name(exmodel, string(JuMP.name(subvar),"_$(scen)"))
                    end
                    if count == 0
                        ### change objective type here
                        JuMP.@constraint(exmodel, coeff*exvar >= value, 
                                                    base_name = exconname)
                        ###
                        count += 1
                    else
                        excon = JuMP.constraint_by_name(exmodel, exconname)
                        JuMP.set_normalized_coefficient(excon, exvar, coeff)
                    end
                end
                #println(JuMP.constraint_by_name(exmodel, exconname))
            end
        end
    end
    return
end


# add charging lb constraint
function add_charging_lb_constraint!(exmodel::JuMP.Model, submodel::JuMP.Model, scen::Int64)
    for ts in timesteps
        for bus in buses
            ### change constraint here
            subcon = get_charge_lb(submodel, bus, ts)
            ###
            subconname = JuMP.name(subcon)
            exconname = string(subconname, "_$(scen)")

            terms = JuMP.constraint_object(subcon).func.terms
            ### change sense here
            value = JuMP.constraint_object(subcon).set.lower
            ###

            count = 0
            for (subvar, coeff) in terms
                if JuMP.name(subvar)[1:2] == "PR" || JuMP.name(subvar)[1:2] == "ER"
                    exvar = JuMP.variable_by_name(exmodel, string(JuMP.name(subvar),"_0"))
                else
                    exvar = JuMP.variable_by_name(exmodel, string(JuMP.name(subvar),"_$(scen)"))
                end
                if count == 0
                    ### change objective type here
                    JuMP.@constraint(exmodel, coeff*exvar >= value, 
                                                base_name = exconname)
                    ###
                    count += 1
                else
                    excon = JuMP.constraint_by_name(exmodel, exconname)
                    JuMP.set_normalized_coefficient(excon, exvar, coeff)
                end
            end
            #println(JuMP.constraint_by_name(exmodel, exconname))
        end
    end
    return
end


# add discharging lb constraint
function add_discharging_lb_constraint!(exmodel::JuMP.Model, submodel::JuMP.Model, scen::Int64)
    for ts in timesteps
        for bus in buses
            ### change constraint here
            subcon = get_discharge_lb(submodel, bus, ts)
            ###
            subconname = JuMP.name(subcon)
            exconname = string(subconname, "_$(scen)")

            terms = JuMP.constraint_object(subcon).func.terms
            ### change sense here
            value = JuMP.constraint_object(subcon).set.lower
            ###

            count = 0
            for (subvar, coeff) in terms
                if JuMP.name(subvar)[1:2] == "PR" || JuMP.name(subvar)[1:2] == "ER"
                    exvar = JuMP.variable_by_name(exmodel, string(JuMP.name(subvar),"_0"))
                else
                    exvar = JuMP.variable_by_name(exmodel, string(JuMP.name(subvar),"_$(scen)"))
                end
                if count == 0
                    ### change objective type here
                    JuMP.@constraint(exmodel, coeff*exvar >= value, 
                                                base_name = exconname)
                    ###
                    count += 1
                else
                    excon = JuMP.constraint_by_name(exmodel, exconname)
                    JuMP.set_normalized_coefficient(excon, exvar, coeff)
                end
            end
            #println(JuMP.constraint_by_name(exmodel, exconname))
        end
    end
    return
end


# add storage lb constraint
function add_storage_lb_constraint!(exmodel::JuMP.Model, submodel::JuMP.Model, scen::Int64)
    for ts in timesteps
        for bus in buses
            ### change constraint here
            subcon = get_storage_lb(submodel, bus, ts)
            ###
            subconname = JuMP.name(subcon)
            exconname = string(subconname, "_$(scen)")

            terms = JuMP.constraint_object(subcon).func.terms
            ### change sense here
            value = JuMP.constraint_object(subcon).set.lower
            ###

            count = 0
            for (subvar, coeff) in terms
                if JuMP.name(subvar)[1:2] == "PR" || JuMP.name(subvar)[1:2] == "ER"
                    exvar = JuMP.variable_by_name(exmodel, string(JuMP.name(subvar),"_0"))
                else
                    exvar = JuMP.variable_by_name(exmodel, string(JuMP.name(subvar),"_$(scen)"))
                end
                if count == 0
                    ### change objective type here
                    JuMP.@constraint(exmodel, coeff*exvar >= value, 
                                                base_name = exconname)
                    ###
                    count += 1
                else
                    excon = JuMP.constraint_by_name(exmodel, exconname)
                    JuMP.set_normalized_coefficient(excon, exvar, coeff)
                end
            end
            #println(JuMP.constraint_by_name(exmodel, exconname))
        end
    end
    return
end

# add loss of load lb constraint
function add_lossofload_lb_constraint!(exmodel::JuMP.Model, submodel::JuMP.Model, scen::Int64)
    for ts in timesteps
        for bus in buses
            ### change constraint here
            subcon = get_loss_of_load_lb(submodel, bus, ts)
            ###
            subconname = JuMP.name(subcon)
            exconname = string(subconname, "_$(scen)")

            terms = JuMP.constraint_object(subcon).func.terms
            ### change sense here
            value = JuMP.constraint_object(subcon).set.lower
            ###

            count = 0
            for (subvar, coeff) in terms
                if JuMP.name(subvar)[1:2] == "PR" || JuMP.name(subvar)[1:2] == "ER"
                    exvar = JuMP.variable_by_name(exmodel, string(JuMP.name(subvar),"_0"))
                else
                    exvar = JuMP.variable_by_name(exmodel, string(JuMP.name(subvar),"_$(scen)"))
                end
                if count == 0
                    ### change objective type here
                    JuMP.@constraint(exmodel, coeff*exvar >= value, 
                                                base_name = exconname)
                    ###
                    count += 1
                else
                    excon = JuMP.constraint_by_name(exmodel, exconname)
                    JuMP.set_normalized_coefficient(excon, exvar, coeff)
                end
            end
            #println(JuMP.constraint_by_name(exmodel, exconname))
        end
    end
    return
end

# add overload lb constraint
function add_overload_lb_constraint!(exmodel::JuMP.Model, submodel::JuMP.Model, scen::Int64)
    for ts in timesteps
        for bus in buses
            ### change constraint here
            subcon = get_overload_lb(submodel, bus, ts)
            ###
            subconname = JuMP.name(subcon)
            exconname = string(subconname, "_$(scen)")

            terms = JuMP.constraint_object(subcon).func.terms
            ### change sense here
            value = JuMP.constraint_object(subcon).set.lower
            ###

            count = 0
            for (subvar, coeff) in terms
                if JuMP.name(subvar)[1:2] == "PR" || JuMP.name(subvar)[1:2] == "ER"
                    exvar = JuMP.variable_by_name(exmodel, string(JuMP.name(subvar),"_0"))
                else
                    exvar = JuMP.variable_by_name(exmodel, string(JuMP.name(subvar),"_$(scen)"))
                end
                if count == 0
                    ### change objective type here
                    JuMP.@constraint(exmodel, coeff*exvar >= value, 
                                                base_name = exconname)
                    ###
                    count += 1
                else
                    excon = JuMP.constraint_by_name(exmodel, exconname)
                    JuMP.set_normalized_coefficient(excon, exvar, coeff)
                end
            end
            #println(JuMP.constraint_by_name(exmodel, exconname))
        end
    end
    return
end

# add wind ub constraint
function add_wind_ub_constraint!(exmodel::JuMP.Model, submodel::JuMP.Model, scen::Int64)
    for ts in timesteps
        bus = 122
        ### change constraint here
        subcon = get_wind_ub(submodel, bus, ts)
        ###
        subconname = JuMP.name(subcon)
        exconname = string(subconname, "_$(scen)")

        terms = JuMP.constraint_object(subcon).func.terms
        ### change sense here
        value = JuMP.constraint_object(subcon).set.upper
        ###

        count = 0
        for (subvar, coeff) in terms
            if JuMP.name(subvar)[1:2] == "PR" || JuMP.name(subvar)[1:2] == "ER"
                exvar = JuMP.variable_by_name(exmodel, string(JuMP.name(subvar),"_0"))
            else
                exvar = JuMP.variable_by_name(exmodel, string(JuMP.name(subvar),"_$(scen)"))
            end
            if count == 0
                ### change objective type here
                JuMP.@constraint(exmodel, coeff*exvar <= value, 
                                            base_name = exconname)
                ###
                count += 1
            else
                excon = JuMP.constraint_by_name(exmodel, exconname)
                JuMP.set_normalized_coefficient(excon, exvar, coeff)
            end
        end
        #println(JuMP.constraint_by_name(exmodel, exconname))
    end
    return
end

# add ramping ub constraint
function add_ramping_ub_constraint!(exmodel::JuMP.Model, submodel::JuMP.Model, scen::Int64)
    for ts in timesteps
        if ts > 1
            for gen in gens
                ### change constraint here
                subcon = get_ramp_ub(submodel, gen, ts)
                ###
                subconname = JuMP.name(subcon)
                exconname = string(subconname, "_$(scen)")

                terms = JuMP.constraint_object(subcon).func.terms
                ### change sense here
                value = JuMP.constraint_object(subcon).set.upper
                ###

                count = 0
                for (subvar, coeff) in terms
                    if JuMP.name(subvar)[1:2] == "PR" || JuMP.name(subvar)[1:2] == "ER"
                        exvar = JuMP.variable_by_name(exmodel, string(JuMP.name(subvar),"_0"))
                    else
                        exvar = JuMP.variable_by_name(exmodel, string(JuMP.name(subvar),"_$(scen)"))
                    end
                    if count == 0
                        ### change objective type here
                        JuMP.@constraint(exmodel, coeff*exvar <= value, 
                                                    base_name = exconname)
                        ###
                        count += 1
                    else
                        excon = JuMP.constraint_by_name(exmodel, exconname)
                        JuMP.set_normalized_coefficient(excon, exvar, coeff)
                    end
                end
                #println(JuMP.constraint_by_name(exmodel, exconname))
            end
        end
    end
    return
end

# add charging ub constraint
function add_charging_ub_constraint!(exmodel::JuMP.Model, submodel::JuMP.Model, scen::Int64)
    for ts in timesteps
        for bus in buses
            ### change constraint here
            subcon = get_charge_ub(submodel, bus, ts)
            ###
            subconname = JuMP.name(subcon)
            exconname = string(subconname, "_$(scen)")

            terms = JuMP.constraint_object(subcon).func.terms
            ### change sense here
            value = JuMP.constraint_object(subcon).set.upper
            ###

            count = 0
            for (subvar, coeff) in terms
                if JuMP.name(subvar)[1:2] == "PR" || JuMP.name(subvar)[1:2] == "ER"
                    exvar = JuMP.variable_by_name(exmodel, string(JuMP.name(subvar),"_0"))
                else
                    exvar = JuMP.variable_by_name(exmodel, string(JuMP.name(subvar),"_$(scen)"))
                end
                if count == 0
                    ### change objective type here
                    JuMP.@constraint(exmodel, coeff*exvar <= value, 
                                                base_name = exconname)
                    ###
                    count += 1
                else
                    excon = JuMP.constraint_by_name(exmodel, exconname)
                    JuMP.set_normalized_coefficient(excon, exvar, coeff)
                end
            end
            #println(JuMP.constraint_by_name(exmodel, exconname))
        end
    end
    return
end

# add discharging ub constraint
function add_discharging_ub_constraint!(exmodel::JuMP.Model, submodel::JuMP.Model, scen::Int64)
    for ts in timesteps
        for bus in buses
            ### change constraint here
            subcon = get_discharge_ub(submodel, bus, ts)
            ###
            subconname = JuMP.name(subcon)
            exconname = string(subconname, "_$(scen)")

            terms = JuMP.constraint_object(subcon).func.terms
            ### change sense here
            value = JuMP.constraint_object(subcon).set.upper
            ###

            count = 0
            for (subvar, coeff) in terms
                if JuMP.name(subvar)[1:2] == "PR" || JuMP.name(subvar)[1:2] == "ER"
                    exvar = JuMP.variable_by_name(exmodel, string(JuMP.name(subvar),"_0"))
                else
                    exvar = JuMP.variable_by_name(exmodel, string(JuMP.name(subvar),"_$(scen)"))
                end
                if count == 0
                    ### change objective type here
                    JuMP.@constraint(exmodel, coeff*exvar <= value, 
                                                base_name = exconname)
                    ###
                    count += 1
                else
                    excon = JuMP.constraint_by_name(exmodel, exconname)
                    JuMP.set_normalized_coefficient(excon, exvar, coeff)
                end
            end
            #println(JuMP.constraint_by_name(exmodel, exconname))
        end
    end
    return
end

# add storage ub constraint
function add_storage_ub_constraint!(exmodel::JuMP.Model, submodel::JuMP.Model, scen::Int64)
    for ts in timesteps
        for bus in buses
            ### change constraint here
            subcon = get_storage_ub(submodel, bus, ts)
            ###
            subconname = JuMP.name(subcon)
            exconname = string(subconname, "_$(scen)")

            terms = JuMP.constraint_object(subcon).func.terms
            ### change sense here
            value = JuMP.constraint_object(subcon).set.upper
            ###

            count = 0
            for (subvar, coeff) in terms
                if JuMP.name(subvar)[1:2] == "PR" || JuMP.name(subvar)[1:2] == "ER"
                    exvar = JuMP.variable_by_name(exmodel, string(JuMP.name(subvar),"_0"))
                else
                    exvar = JuMP.variable_by_name(exmodel, string(JuMP.name(subvar),"_$(scen)"))
                end
                if count == 0
                    ### change objective type here
                    JuMP.@constraint(exmodel, coeff*exvar <= value, 
                                                base_name = exconname)
                    ###
                    count += 1
                else
                    excon = JuMP.constraint_by_name(exmodel, exconname)
                    JuMP.set_normalized_coefficient(excon, exvar, coeff)
                end
            end
            #println(JuMP.constraint_by_name(exmodel, exconname))
        end
    end
    return
end


    # add thermal interval constraint
function add_thermal_interval_constraint!(exmodel::JuMP.Model, submodel::JuMP.Model, scen::Int64)
    for ts in timesteps
        for gen in gens
            ### change constraint here
            subcon = get_thermal_interval(submodel, gen, ts)
            ###
            subconname = JuMP.name(subcon)
            exconname = string(subconname, "_$(scen)")

            terms = JuMP.constraint_object(subcon).func.terms
            ### change sense here
            upper = JuMP.constraint_object(subcon).set.upper
            lower = JuMP.constraint_object(subcon).set.lower
            ###

            count = 0
            for (subvar, coeff) in terms
                if JuMP.name(subvar)[1:2] == "PR" || JuMP.name(subvar)[1:2] == "ER"
                    exvar = JuMP.variable_by_name(exmodel, string(JuMP.name(subvar),"_0"))
                else
                    exvar = JuMP.variable_by_name(exmodel, string(JuMP.name(subvar),"_$(scen)"))
                end
                if count == 0
                    ### change objective type here
                    JuMP.@constraint(exmodel, lower <= coeff*exvar <= upper, 
                                                base_name = exconname)
                    ###
                    count += 1
                else
                    excon = JuMP.constraint_by_name(exmodel, exconname)
                    JuMP.set_normalized_coefficient(excon, exvar, coeff)
                end
            end
            #println(JuMP.constraint_by_name(exmodel, exconname))
        end
    end
    return
end

# add branch interval constraint
function add_branch_interval_constraint!(exmodel::JuMP.Model, submodel::JuMP.Model, scen::Int64)
    for ts in timesteps
        for br in branches
            ### change constraint here
            subcon = get_branch_interval(submodel, br, ts)
            ###
            subconname = JuMP.name(subcon)
            exconname = string(subconname, "_$(scen)")

            terms = JuMP.constraint_object(subcon).func.terms
            ### change sense here
            upper = JuMP.constraint_object(subcon).set.upper
            lower = JuMP.constraint_object(subcon).set.lower
            ###

            count = 0
            for (subvar, coeff) in terms
                if JuMP.name(subvar)[1:2] == "PR" || JuMP.name(subvar)[1:2] == "ER"
                    exvar = JuMP.variable_by_name(exmodel, string(JuMP.name(subvar),"_0"))
                else
                    exvar = JuMP.variable_by_name(exmodel, string(JuMP.name(subvar),"_$(scen)"))
                end
                if count == 0
                    ### change objective type here
                    JuMP.@constraint(exmodel, lower <= coeff*exvar <= upper, 
                                                base_name = exconname)
                    ###
                    count += 1
                else
                    excon = JuMP.constraint_by_name(exmodel, exconname)
                    JuMP.set_normalized_coefficient(excon, exvar, coeff)
                end
            end
            #println(JuMP.constraint_by_name(exmodel, exconname))
        end
    end
    return
end

# change objective function for extensive form

function add_to_objective_function!(exmodel::JuMP.Model,submodel::JuMP.Model, scen::Int64)

    exobjfunc = JuMP.objective_function(exmodel)

    for (subvar, coeff) in JuMP.objective_function(submodel).terms

        if JuMP.name(subvar)[1:2] == "PR" || JuMP.name(subvar)[1:2] == "ER"
            exvar = JuMP.variable_by_name(exmodel, string(JuMP.name(subvar),"_0"))
        else
            exvar = JuMP.variable_by_name(exmodel, string(JuMP.name(subvar),"_$(scen)"))
        end

        JuMP.add_to_expression!(exobjfunc, coeff, exvar)

    end
    
    JuMP.@objective(exmodel, Min, exobjfunc);
    
    return
end

### END LOOPING FUNCTIONS ###

# set extensive form objective to 1/nscen times itself
function set_ef_objective_function!(exmodel::JuMP.Model, nscen::Int64)
    exobjfunc = JuMP.objective_function(exmodel)

    JuMP.@objective(exmodel, Min, 1/nscen*exobjfunc)

    #println(JuMP.objective_function(exmodel))
    return
end


### main extensive form function
function construct_ef_model(exmodel, submodel, load, wind, nscen)
    
    
    println("Adding first stage info.")
    add_firststage_variables!(exmodel, submodel)
    add_ER_lb!(exmodel, submodel)
    add_PR_lb!(exmodel, submodel)
    add_ER_ub!(exmodel, submodel)
    add_expansion_budget!(exmodel, submodel)
    
    # initialize objective function
    JuMP.@objective(exmodel, Min, 0)
    
    println("Adding second stage info...")
    for scen = 1:nscen
        
        println("...Scenario $(scen)...")
        
        loadvec = load[scen]
        windvec = wind[scen]
        
        println("......Adjusting uncertainty...")
        println(".........load balance...")
        adjust_load_balance_constraint!(submodel, loadvec, loaddis)  
        println(".........ptdf...")
        adjust_ptdf_constraint!(exmodel, submodel, ptdfdict, loadvec, loaddis) 
        println(".........wind upper bound...")
        adjust_wind_ub!(submodel, windvec)
        
        println("......Adding variables...")
        add_thermal_variables!(exmodel, submodel, scen)
        add_branch_variables!(exmodel, submodel, scen)
        add_wind_variables!(exmodel, submodel, scen)
        add_charging_variables!(exmodel, submodel, scen)
        add_discharging_variables!(exmodel, submodel, scen)
        add_storage_variables!(exmodel, submodel, scen)
        add_lossofload_variables!(exmodel, submodel, scen)
        add_overload_variables!(exmodel, submodel, scen)
        
        println("......Adding constraints...")
        add_storage_balance_constraint!(exmodel, submodel, scen)
        add_load_balance_constraint!(exmodel, submodel, scen)
        add_ptdf_constraint!(exmodel, submodel, scen)
        add_wind_lb_constraint!(exmodel, submodel, scen)
        add_ramping_lb_constraint!(exmodel, submodel, scen)
        add_charging_lb_constraint!(exmodel, submodel, scen)
        add_discharging_lb_constraint!(exmodel, submodel, scen)
        add_storage_lb_constraint!(exmodel, submodel, scen)
        add_lossofload_lb_constraint!(exmodel, submodel, scen)
        add_overload_lb_constraint!(exmodel, submodel, scen)
        add_wind_ub_constraint!(exmodel, submodel, scen)
        add_ramping_ub_constraint!(exmodel, submodel, scen)
        add_charging_ub_constraint!(exmodel, submodel, scen)
        add_discharging_ub_constraint!(exmodel, submodel, scen)
        add_storage_ub_constraint!(exmodel, submodel, scen)
        add_thermal_interval_constraint!(exmodel, submodel, scen)
        add_branch_interval_constraint!(exmodel, submodel, scen)
        
        println("......Adding to objective function...")
        add_to_objective_function!(exmodel,submodel,scen)
        
    end
    
    println("Setting objective function.")
    
    set_ef_objective_function!(exmodel, nscen)
    
    println("Done.")
    
    return exmodel, submodel
end
        
    

