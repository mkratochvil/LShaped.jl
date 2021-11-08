mutable struct LocalVariableInfo
    id::Int64
    name::String
    index::Union{Int64,Nothing}
    value::Union{Float64,Nothing}
    gradient::Union{Float64,Nothing}
end

mutable struct Subproblems
    id::Int64
    model::JuMP.Model
    probability::Float64
    variableinfo::Dict{Int64,LocalVariableInfo}
    vnametoind::Dict{String,Int64}
    objective_value::Union{Float64,Nothing}
end

mutable struct SubproblemsNew
    id::Int64
    model::JuMP.Model
    probability::Float64
    variableinfo::Dict{Int64,LocalVariableInfo}
    idxtocon::Dict{Int64,JuMP.ConstraintRef}
    h::Union{Array{Float64},Nothing}
    Ek::Union{Array{Float64},Nothing}
    ek::Union{Float64,Nothing}
    vnametoind::Dict{String,Int64}
    objective_value::Union{Float64,Nothing}
end

mutable struct FirstStageInfo
    variables::Dict 
    subproblems::Dict
    store::Union{String,Nothing}
end

mutable struct FirstStageVariableInfo
    name::String
    index::Int64
    value::Union{Float64,Nothing}
    gradient::Union{Float64,Nothing}
    lowerbound::Union{Float64,Nothing}
    upperbound::Union{Float64,Nothing}
end