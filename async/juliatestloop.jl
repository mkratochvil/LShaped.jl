using Pkg
Pkg.activate("..")

using JuMP
using Gurobi
using DataFrames
using CSV
using LinearAlgebra

using LShaped


for i = 1:200
    
    include("firststage.jl")
    
    #replace with asynchronous for loop with arrayid
    include("secondstage1.jl")
    include("secondstage2.jl")
    include("secondstage3.jl")
    include("secondstage4.jl")
    include("secondstage5.jl")
    include("secondstage6.jl")
    include("secondstage7.jl")
    include("secondstage8.jl")
    include("secondstage9.jl")
    include("secondstage10.jl")
    include("secondstage11.jl")
    include("secondstage12.jl")    
end
