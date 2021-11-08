mutable struct LocalVariableInfo
    id::Int64
    name::String
    index::Any
    value::Union{Float64,Nothing}
    gradient::Union{Float64,Nothing}
end

mutable struct Subproblems
    id::Int64
    model::Any
    probability::Float64
    variableinfo::Any
    vnametoind::Dict
    objective_value::Any
end

mutable struct SubproblemsNew
    id::Int64
    model::Any
    probability::Float64
    variableinfo::Any
    idxtocon::Dict
    h::Any
    Ek::Any
    ek::Any
    vnametoind::Dict
    objective_value::Any
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