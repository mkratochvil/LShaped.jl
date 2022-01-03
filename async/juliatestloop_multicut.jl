using Pkg
Pkg.activate("..")

using JuMP
using Gurobi
using DataFrames
using CSV
using LinearAlgebra

using LShaped


for i = 1:250
    
    include("firststage_multicut.jl")
    
    #replace with asynchronous for loop with arrayid
    include("secondstage_multicut1.jl")
    include("secondstage_multicut2.jl")
    include("secondstage_multicut3.jl")
    include("secondstage_multicut4.jl")
    include("secondstage_multicut5.jl")
    include("secondstage_multicut6.jl")
    include("secondstage_multicut7.jl")
    include("secondstage_multicut8.jl")
    include("secondstage_multicut9.jl")
    include("secondstage_multicut10.jl")
    include("secondstage_multicut11.jl")
    include("secondstage_multicut12.jl")
    #=
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
    =#
end
