using Pkg
Pkg.activate("..")

using JuMP
using Gurobi
using DataFrames
using CSV
using LinearAlgebra

using LShaped


for i = 1:10
    
    include("firststage.jl")
    
    #replace with asynchronous for loop with arrayid
    include("secondstage1.jl")
    include("secondstage2.jl")
    
end