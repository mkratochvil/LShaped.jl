mutable struct LocalVariableInfo
    id::Int64
    name::String
    index::Any
    cost::Float64
    conval::Any
    value::Union{Float64,Nothing}
    gradient::Union{Float64,Nothing}
end

mutable struct LinkedConstraintInfo
    id::Int64
    func::Any
    set::Any
    ref::Any
    variables::Array
    initvalue::Float64
    curvalue::Float64
end

mutable struct Subproblems
    id::Int64
    model::Any
    probability::Float64
    variableinfo::Any
    linkedconstraintinfo::Any
    vnametoind::Dict
    arrays::Any
    objective_value::Any
end

mutable struct Arrays
    Ae::Any
    Al::Any
    Ag::Any
    Ie::Any
    Il::Any
    Ig::Any
end

# More stuff stored here like convergence criteria, etc stored here
mutable struct FirstStageInfo
    variables::Dict # dict from variable string name to FirstStageVariableInfo 
    subproblems::Dict
    store::Union{String,Nothing}
    #probability::Float64
    #variableinfo::Any
    #linkedconstraintinfo::Any
end

mutable struct FirstStageVariableInfo
    #id::Int64
    name::String
    index::Int64
    #cost::Float64
    #conval::Any
    value::Union{Float64,Nothing}
    gradient::Union{Float64,Nothing}
    lowerbound::Union{Float64,Nothing}
    upperbound::Union{Float64,Nothing}
    status::Any
end