using Pkg
Pkg.activate("..")

using JuMP
using Gurobi
using DataFrames
using CSV
using LinearAlgebra

using LShaped


for i = 1:5
    
    include("firststage_ro_cluster.jl")
    
    #replace with asynchronous for loop with arrayid
    include("secondstage_ro_loop_cluster.jl")
    #=include("secondstage_ro1.jl")
    include("secondstage_ro2.jl")
    include("secondstage_ro3.jl")
    include("secondstage_ro4.jl")
    include("secondstage_ro5.jl")
    include("secondstage_ro6.jl")
    include("secondstage_ro7.jl")
    include("secondstage_ro8.jl")
    include("secondstage_ro9.jl")
    include("secondstage_ro10.jl")
    include("secondstage_ro11.jl")
    include("secondstage_ro12.jl") =#   
end
