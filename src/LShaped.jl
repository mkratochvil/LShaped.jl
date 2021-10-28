module LShaped

using JuMP
using LinearAlgebra 
using MathOptInterface
using DataFrames
using CSV
const MOI = MathOptInterface

include("structs.jl")
include("setup.jl")
include("algorithm.jl")



"""
    L_Shaped_Algorithm(subproblem_generator, 
	v_dict, 
	N, 
	master_generator, 
	tol = 1e-6, 
	maxiter=10, 
	probs = 1/N*ones(N))

**Arguments**

* `subproblem_generator` : function that generates your subproblem. Should have a single argument: scenario id number
* `v_dict` : variable dictionary. keys are stage id (e.g. 1, 2) and values should be 4ples: (vname::string, lb, ub, init)
* `N` : Total number of scenarios.
* `tol` : tolerance for convergence.
* `maxiter` : maximum number of iterations of L-Shaped Method that should be implemented.
* `probs` : vector of probabilities for each scenario (indexed by id number). Defaults to 1/N.
"""
function L_Shaped_Algorithm(subproblem_generator, 
                            v_dict, 
                            N, 
                            master_generator, 
                            tol = 1e-6, 
                            maxiter=10, 
                            probs = 1/N*ones(N);
                            store = nothing)
    
    println(store)
    
    firststage, contoidx, h = make_two_stage_setup_L(subproblem_generator, v_dict, N, probs, store);
    
    fs = master_generator()
    
    ittime = @elapsed x, firststage, fs, contoidx, h, niter = iterate_L(firststage, fs, contoidx, h, v_dict, 0, tol, maxiter)
    
    return x, firststage, fs, contoidx, h, ittime, niter
    
end

end # module
