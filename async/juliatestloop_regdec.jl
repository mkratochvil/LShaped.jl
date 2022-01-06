using Pkg
Pkg.activate("..")

using JuMP
using Gurobi
using DataFrames
using CSV
using LinearAlgebra

using LShaped

for i = 1:100
    
    include("firststage_regdec.jl")
    
    #replace with asynchronous for loop with arrayid
    include("secondstage_regdec1.jl")
    include("secondstage_regdec2.jl")
    include("secondstage_regdec3.jl")
    include("secondstage_regdec4.jl")
    include("secondstage_regdec5.jl")
    include("secondstage_regdec6.jl")
    include("secondstage_regdec7.jl")
    include("secondstage_regdec8.jl")
    include("secondstage_regdec9.jl")
    include("secondstage_regdec10.jl")
    include("secondstage_regdec11.jl")
    include("secondstage_regdec12.jl")

end
