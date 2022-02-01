using Pkg
Pkg.activate("..")

using JuMP
using Gurobi
using DataFrames
using CSV
using LinearAlgebra

using LShaped

for i = 1:150
    
    include("firststage_trust.jl")
    
    #replace with asynchronous for loop with arrayid
    include("secondstage_regdec_loop1.jl")
    include("secondstage_regdec_loop2.jl")
    include("secondstage_regdec_loop3.jl")
    include("secondstage_regdec_loop4.jl")
    #=include("secondstage_regdec2.jl")
    include("secondstage_regdec3.jl")
    include("secondstage_regdec4.jl")
    include("secondstage_regdec5.jl")
    include("secondstage_regdec6.jl")
    include("secondstage_regdec7.jl")
    include("secondstage_regdec8.jl")
    include("secondstage_regdec9.jl")
    include("secondstage_regdec10.jl")
    include("secondstage_regdec11.jl")
    include("secondstage_regdec12.jl")=#

end