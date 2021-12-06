using Pkg
#make sure to run from directly inside the LShaped.jl package until you figure out how tests actually work.
Pkg.activate(".")

using Suppressor
using LinearAlgebra

tol = 1e-8

println(pwd())
include("./bl_ch5_ex1.jl")
include("./wsgep_2_288.jl")

blxtest = [46.66666666666667, 36.25]

wsgepxtest = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1.757008472,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,35.3516819,0,0,0,
0,0,0,0,0,0,0,0.254070409,0,0,0,0,0,0,3.565270862,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];

blx, blfstruct, blfsmodel = LShaped.L_Shaped_Algorithm(bl2, blv, 2, bl1, 1e-6, 10, [0.4, 0.6]; store="./bl_data/");

blt = norm(blxtest-blx)
if blt < tol
    println("BL Old Test passed. $(blt)")
else
    println("BL Old Test FAILED!!! $(blt)")
end

blxn, blfstructn, blfsmodeln = LShaped.L_Shaped_Algorithm_new(bl2, blv, 2, bl1, 1e-6, 3, [0.4, 0.6]; store="./bl_data_new/", verbose=1, resume=0);

blxn, blfstructn, blfsmodeln = LShaped.L_Shaped_Algorithm_new(bl2, blv, 2, bl1, 1e-6, 10, [0.4, 0.6]; store="./bl_data_new/", verbose=1, resume=1);

bltn = norm(blxtest-blxn)
if bltn < tol
	println("BL New Test passed. $(bltn)")
else
	println("BL New Test FAILED!!! $(bltn)")
end

wsgepx, wsgepfsstruct, wsgepfsmodel = @suppress LShaped.L_Shaped_Algorithm(wsgep2, wsgepv, 2, wsgep1, 1e-6, 100; store="./ws_data/");

wsgept = norm(wsgepxtest-wsgepx)
if wsgept < tol
    println("WSGEP Old Test passed. $(wsgept)")
else
    println("WSGEP Old Test FAILED!!! $(wsgept)")
end

wsgepxn, wsgepfsstructn, wsgepfsmodeln = @suppress LShaped.L_Shaped_Algorithm_new(wsgep2, wsgepv, 2, wsgep1, 1e-6, 100; store="./ws_data_new/");

wsgeptn = norm(wsgepxtest-wsgepxn)
if wsgeptn < tol
    println("WSGEP New Test passed. $(wsgeptn)")
else
    println("WSGEP New Test FAILED!!! $(wsgeptn)")
end


xpause = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.040444786119284, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 20.9746459833302, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.264501879667551, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.79094469682893, 0.0, 0.0, 0.0, 0.0, 3.5420215765336325, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

xpauset, wsgepfsstructn, wsgepfsmodeln = LShaped.L_Shaped_Algorithm_new(wsgep2, wsgepv, 2, wsgep1, 1e-6, 25; store="./ws_data_new_pause/");

xpt = norm(xpause-xpauset)
if xpt < tol
    println("WSGEP Pause Test passed. $(xpt)")
else
    println("WSGEP Pause Test FAILED!!! $(xpt)")
end

wsgepxn, wsgepfsstructn, wsgepfsmodeln = LShaped.L_Shaped_Algorithm_new(wsgep2, wsgepv, 2, wsgep1, 1e-6, 100; store="./ws_data_new_pause/", resume=1);

wsgeptn = norm(wsgepxtest-wsgepxn)
if wsgeptn < tol
    println("WSGEP New Test passed. $(wsgeptn)")
else
    println("WSGEP New Test FAILED!!! $(wsgeptn)")
end
