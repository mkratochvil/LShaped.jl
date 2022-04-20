using Pkg
Pkg.activate("..")

using JuMP
using Gurobi
using DataFrames
using CSV
using LinearAlgebra

using LShaped

for i = 1:400
    
    include("kfs_cluster_4wind.jl")
    include("kss_cluster_4wind.jl")

end