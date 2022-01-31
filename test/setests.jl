using Pkg
#make sure to run from directly inside the LShaped.jl package until you figure out how tests actually work.
Pkg.activate(".")

using Suppressor
using LinearAlgebra
using DataFrames
using CSV
using LShaped

global ercotscens = collect(DataFrame(CSV.File("../FinalProject/scenarios/strat133.csv"))[1,:])

function get_sol_info(xvec, model)

    PRsum = 0.; 
    ERsum = 0.; 
    for i = 0:47 
        println(1+i/2, " ", xm[i+1]) 
        if i % 2 == 0 
            PRsum += xm[i+1] 
        else 
            ERsum += xm[i+1] 
        end 
    end 
    println("PRcost = $(PRsum)") 
    println("ERcost = $(ERsum)")
    println("objective_value = $(JuMP.objective_value(model))")

    return
end

function get_sol_info_PR(xvec, model)

    PRsum = 0.; 
    for i = 1:24
        println(100+i, " ", xm[i]) 
        PRsum += xm[i]  
    end 
    println("PRcost = $(PRsum)") 
    println("objective_value = $(JuMP.objective_value(model))")

    return
end

include("./storage_expansion_models.jl")
include("./bl_ch5_ex1.jl")

testsol =  zeros(48);
testsol[27] = 0.2970924180099013;
testsol[28] = 1.4854620900495066;
testsol[31] = 1.6240825598620396;
testsol[32] = 8.120412799310198;
testsol[33] = 1.469783051206798;
testsol[34] = 7.348915256033989;
testsol[43] = 0.08005897890460051;
testsol[44] = 0.40029489452300254;

tol = 1e-8

blxtest = [46.66666666666667, 36.25]


#blx, blfstruct, blfsmodel = LShaped.L_Shaped_Algorithm_new(bl2, blv, 2, bl1, 1e-6, 10, [0.4, 0.6]);

#blxm, blfstructm, blfsmodelm = LShaped.L_Shaped_Algorithm_new(bl2, blv, 2, bl1, 1e-6, 10, [0.4, 0.6]; multicut = 1);

#blxm, blfstructm, blfsmodelm = LShaped.L_Shaped_Algorithm_new(bl2, blv, 2, bl1, 1e-6, 10, [0.4, 0.6]; multicut = 3);

#xm, probstructm, fsmodelm, ittimem, niterm = LShaped.L_Shaped_Algorithm_new(second_func, sedict, 12, first_func, 1e-8, 1000; multicut=2, rho = 1.0, rhomin = .0001, rhomax = 1000., gamma = 0.01);

xm, probstructm, fsmodelm, ittimem, niterm = LShaped.L_Shaped_Algorithm_new(second_func, sedict, 12, first_func, 1e-6, 5; multicut=3);
#for i = 0:47 println(1+i/2, " ", xm[i+1]) end
#println(JuMP.objective_value(fsmodelm))
get_sol_info_PR(xm, fsmodelm)
println("$(niterm)")

println(" ")

#x, probstruct, fsmodel, ittime, niter = LShaped.L_Shaped_Algorithm_new(second_func, sedict, 12, first_func, 1e-6, 200);
