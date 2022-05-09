using Pkg
Pkg.activate("..")

using JuMP
using Gurobi
using DataFrames
using CSV
using LinearAlgebra

using LShaped


for i = 1:25
    
    include("firststage_ro_sumlol.jl")
    
    include("secondstage_ro_loop_sumlol.jl") 
end
