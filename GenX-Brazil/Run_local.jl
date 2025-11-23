case = dirname(@__FILE__);

genx_benders_path = dirname(dirname(case))*"/GenX-Benders/"

import Pkg
Pkg.activate(genx_benders_path)
include(genx_benders_path*"src/GenX-Benders.jl")

run_genx_case!(case)