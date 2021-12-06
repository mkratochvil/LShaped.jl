blv = Dict{Int64,Array{Any}}()

blv[1] = [("x1", 40, Inf, 40), ("x2", 20, Inf, 20)]
blv[2] = ["y1", "y2"]

function bl1()
  
    fs = Model(with_optimizer(Gurobi.Optimizer, OutputFlag=0));
    
    @variable(fs, x1 >= 40)
    @variable(fs, x2 >= 20)
    
    @objective(fs, Min, 100*x1 + 150*x2)
    #@objective(fs, Min, 0.0)
    
    @constraint(fs, x1 + x2 <= 120)
    
    return fs
    
end
