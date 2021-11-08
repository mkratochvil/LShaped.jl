using JuMP
using Gurobi
using CSV
using DataFrames

using LShaped

# implementation of Birge, Louveax Ch. 5 Ex 1 Second Stage

function bl2(sid)
    
    q1 = [-24, -28]
    q2 = [-28, -32]
    d1 = [500, 300]
    d2 = [100, 300]
    
    model = Model(with_optimizer(Gurobi.Optimizer; OutputFlag=0))
    
    @variable(model, x1 )# >= 40.0)
    @variable(model, x2 )# >= 20.0)
    @variable(model, y1[sid] >= 0)
    @variable(model, y2[sid] >= 0)
    
    @objective(model, Min, 100*x1 + 150*x2 + q1[sid]*y1[sid] + q2[sid]*y2[sid])
    
   # @constraint(model, x1 + x2 <= 120)
    
    @constraint(model, 6*y1[sid] + 10*y2[sid] <= 60*x1)
    @constraint(model, 8*y1[sid] + 5*y2[sid] <= 80*x2)
    
    @constraint(model, y1[sid] <= d1[sid])
    @constraint(model, y2[sid] <= d2[sid])
    
    return model
    
end


blv = Dict{Int64,Array{Any}}()

blv[1] = [("x1", 40, Inf, 40), ("x2", 20, Inf, 20)]
blv[2] = ["y1", "y2"]

function bl1()
  
    fs = Model(with_optimizer(Gurobi.Optimizer, OutputFlag=0));
    
    @variable(fs, x1 >= 40)
    @variable(fs, x2 >= 20)
    
    #@objective(fs, Min, 100*x1 + 150*x2)
    @objective(fs, Min, 0.0)
    
    @constraint(fs, x1 + x2 <= 120)
    
    return fs
    
end



