case = dirname(@__FILE__);

genx_benders_path = dirname(dirname(case))*"/GenX-Benders/"

import Pkg
Pkg.activate(genx_benders_path)
include(genx_benders_path*"src/GenX-Benders.jl")

using Distributed,ClusterManagers
ntasks = parse(Int, ENV["SLURM_NTASKS"]);
cpus_per_task = parse(Int, ENV["SLURM_CPUS_PER_TASK"]);
addprocs(SlurmManager(ntasks);exeflags=["-t $cpus_per_task"])

@everywhere genx_benders_path = dirname(dirname(case))*"/GenX-Benders/"

@everywhere begin
    import Pkg
    Pkg.activate(genx_benders_path)
end

println("Number of procs: ", nprocs())
println("Number of workers: ", nworkers())

@everywhere include(genx_benders_path*"src/GenX-Benders.jl")

run_genx_case!(case)