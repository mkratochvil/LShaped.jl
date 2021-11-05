using JuMP
using Gurobi
using CSV
using DataFrames

using LShaped

# implementation of Birge, Louveax Ch. 5 Ex 1 Second Stage

function subproblem_constructor(sid)
    
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


v_dict = Dict()

v_dict[1] = [("x1", 40, Inf, 40), ("x2", 20, Inf, 20)]
v_dict[2] = ["y1", "y2"]

function create_first_stage()
  
    fs = Model(with_optimizer(Gurobi.Optimizer, OutputFlag=0));
    
    @variable(fs, x1 >= 40)
    @variable(fs, x2 >= 20)
    
    @objective(fs, Min, 100*x1 + 150*x2)
    
    @constraint(fs, x1 + x2 <= 120)
    
    return fs
    
end


xn, firststage, fs = LShaped.L_Shaped_Algorithm(subproblem_constructor, 
                                        v_dict, 2, create_first_stage, 1e-6, 10, [0.4, 0.6])#; store="./bl_data/");

xnn, firststagen, fsn = LShaped.L_Shaped_Algorithm_new(subproblem_constructor, 
                                        v_dict, 2, create_first_stage, 1e-6, 10, [0.4, 0.6])#; store="./bl_data_new/");

#rm("bl_data",recursive=true)
