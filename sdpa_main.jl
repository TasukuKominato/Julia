include("sdpa_file.jl")
include("sdpa02_tsk.jl")
# include("sdpa_newton_test.jl")

# m, nBlock, blockStruct, b, F0, F = sdpa_read("example1.dat-s")
# m, nBlock, blockStruct, b, F0, F = sdpa_read("control1.dat-s")
# m, nBlock, blockStruct, b, F0, F = sdpa_read("arch0.dat-s")
# m, nBlock, blockStruct, b, F0, F = sdpa_read("arch2.dat-s")
# m, nBlock, blockStruct, b, F0, F = sdpa_read("arch8.dat-s")
m, nBlock, blockStruct, b, F0, F = sdpa_read("example1.dat-s")
sdpa_newton(m, nBlock, blockStruct, -b, -F0, -F)
