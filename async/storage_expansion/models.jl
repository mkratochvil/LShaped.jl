include("../../../FinalProject/parameters.jl")
include("../../../FinalProject/get_functions.jl")
#include("../../../FinalProject/modification_functions.jl")

function first_func()
    
    model = JuMP.read_from_file("../../FinalProject/storage_expansion_revised/first_stage/first_stage_model_PR_exp1B.mps")
    JuMP.set_optimizer(model, Gurobi.Optimizer)
    set_optimizer_attribute(model, "OutputFlag", 0) 
    JuMP.set_optimizer_attribute(model, "Method", 1)  
    
    return model
end

function second_func(i)
    
    scen = ercotscens[i]
    #model = JuMP.read_from_file("../../FinalProject/storage_expansion_revised/second_stage/noint_PR_exp3_scen_$(scen).mps")
    model = JuMP.read_from_file("../../FinalProject/storage_expansion_revised/ercot/PR_exp3_scen_$(scen).mps")
    JuMP.set_optimizer(model, Gurobi.Optimizer)
    set_optimizer_attribute(model, "OutputFlag", 0) 
    JuMP.set_optimizer_attribute(model, "Method", 1) 

    return model
end

 sedict = Dict{Int64,Array{Any}}()

#modelvars = JuMP.read_from_file("../../FinalProject/storage_expansion_revised/second_stage/noint_PR_exp5p0015_scen_1.mps")
modelvars = JuMP.read_from_file("../../FinalProject/storage_expansion_revised/ercot/PR_exp3_scen_1.mps")

firstvars = []
for bus in buses
    push!(firstvars,(string(get_PR_variable(modelvars,bus)), 0.0, Inf, 0.0))
    #push!(firstvars,(string(get_ER_variable(modelvars,bus)), 0.0, Inf, 0.0))
end

 sedict[1] = firstvars

secondvars = [];

for ts in timesteps
    
    for gen in gens
        push!(secondvars,string(get_thermal_variable(modelvars,gen,ts)))
    end
    for br in branches
        push!(secondvars,string(get_line_variable(modelvars,br,ts)))
    end
    
    push!(secondvars, string(get_wind_variable(modelvars,122,ts)))
    
    for bus in buses
        push!(secondvars, string(get_charging_variable(modelvars,bus,ts)))
        push!(secondvars, string(get_discharging_variable(modelvars,bus,ts)))
        push!(secondvars, string(get_stored_variable(modelvars,bus,ts)))
        push!(secondvars, string(get_lossofload_variable(modelvars,bus,ts)))
        push!(secondvars, string(get_overload_variable(modelvars,bus,ts)))
    end
end

 sedict[2] = secondvars