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

blx, blfstruct, blfsmodel =  @suppress LShaped.L_Shaped_Algorithm(bl2, blv, 2, bl1, 1e-6, 10, [0.4, 0.6])#; store="./bl_data/");

blxn, blfstructn, blfsmodeln = @suppress LShaped.L_Shaped_Algorithm_new(bl2, blv, 2, bl1, 1e-6, 10, [0.4, 0.6])#; store="./bl_data_new/");

wsgepx, wsgepfsstruct, wsgepfsmodel = @suppress LShaped.L_Shaped_Algorithm(wsgep2, wsgepv, 2, wsgep1, 1e-6, 100)#; store="./bl_data/");

wsgepxn, wsgepfsstructn, wsgepfsmodeln = @suppress LShaped.L_Shaped_Algorithm_new(wsgep2, wsgepv, 2, wsgep1, 1e-6, 100)#; store="./bl_data_new/");

blt = norm(blxtest-blx)
if blt < tol
	println("BL Old Test passed. $(blt)")
else
	println("BL Old Test FAILED!!! $(blt)")
end

bltn = norm(blxtest-blxn)
if bltn < tol
	println("BL New Test passed. $(bltn)")
else
	println("BL New Test FAILED!!! $(bltn)")
end

wsgept = norm(wsgepxtest-wsgepx)
if wsgept < tol
	println("WSGEP Old Test passed. $(wsgept)")
else
	println("WSGEP Old Test FAILED!!! $(wsgept)")
end

wsgeptn = norm(wsgepxtest-wsgepxn)
if wsgeptn < tol
	println("WSGEP New Test passed. $(wsgeptn)")
else
	println("WSGEP New Test FAILED!!! $(wsgeptn)")
end

