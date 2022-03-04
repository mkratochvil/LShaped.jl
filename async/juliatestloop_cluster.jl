using Pkg
Pkg.activate("..")

using JuMP
using Gurobi
using DataFrames
using CSV
using LinearAlgebra

using LShaped

for i = 1:400
    
    include("kfs_cluster.jl")
    include("kss_cluster.jl")

end