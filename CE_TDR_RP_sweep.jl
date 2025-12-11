# CE_TDR_RP_sweep.jl
#
# Sweep over different representative periods per year for the Brazil Benders case

using YAML
using CSV, DataFrames

const CASE = joinpath(@__DIR__, "GenX-Brazil")
const TDR_SETTINGS_PATH = joinpath(CASE, "Settings", "time_domain_reduction_settings.yml")
const RESULTS_ROOT = joinpath(CASE, "Results_Benders_Int_LevelSet")
const CE_RESULTS_ARCHIVE_ROOT = joinpath(CASE, "CE_TDR_RP_results")
const OUT_SUMMARY_CSV = joinpath(CASE, "CE_TDR_RP_summary.csv")

rep_periods = [1, 5, 10, 15, 20]

isdir(CE_RESULTS_ARCHIVE_ROOT) || mkpath(CE_RESULTS_ARCHIVE_ROOT)

println("=== CE + TDR Representative Period Sweep ===")
println("Case path: $CASE")

function read_cost_for_run(results_root::AbstractString)
    candidate_paths = [
        joinpath(results_root, "master_sol_costs.csv"),
        joinpath(results_root, "Master_sol_costs.csv")
    ]
    summary_path = first(filter(isfile, candidate_paths); init = "")
    if summary_path == ""
        @warn "master_sol_costs.csv not found"
        return missing
    end
    df = CSV.read(summary_path, DataFrame)
    needed = [:FixCost, :NetworkExpCost, :UnmetCapReqCost]
    missing_cols = setdiff(needed, names(df))
    if !isempty(missing_cols)
        @warn "missing columns in cost file"
        return missing
    end
    row = df[1, :]
    return float(row[:FixCost]) + float(row[:NetworkExpCost]) + float(row[:UnmetCapReqCost])
end

function archive_ce_results(rp::Integer)
    if !isdir(RESULTS_ROOT)
        @warn "Results root missing, nothing archived"
        return
    end
    rp_dir = joinpath(CE_RESULTS_ARCHIVE_ROOT, "RP_$(rp)")
    isdir(rp_dir) && rm(rp_dir; recursive = true)
    mkpath(rp_dir)
    dst = joinpath(rp_dir, "Results_Benders_Int_LevelSet")
    cp(RESULTS_ROOT, dst; force = true)
    println("  Archived CE results to $dst")
end

tdr_settings = YAML.load_file(TDR_SETTINGS_PATH)
include("src/GenX-Benders.jl")

results_df = DataFrame(
    Rep_Periods_per_Year = Int[],
    TDR_Runtime_s        = Float64[],
    CE_Runtime_s         = Float64[],
    Total_Runtime_s      = Float64[],
    Total_System_Cost    = Union{Missing,Float64}[],
)

for (idx, rp) in pairs(rep_periods)
    println("\n=== Sweep $(idx)/$(length(rep_periods))  RP_per_year = $rp ===")

    println("  Updating TDR settings...")
    tdr_settings["MinPeriods"] = rp
    tdr_settings["MaxPeriods"] = rp
    YAML.write_file(TDR_SETTINGS_PATH, tdr_settings)

    println("  Clearing old TDR outputs...")
    tdr_top = joinpath(CASE, "TDR_Results")
    isdir(tdr_top) && rm(tdr_top; recursive = true)
    scen_dir = joinpath(CASE, "Scenarios")
    if isdir(scen_dir)
        for scen in filter(name -> startswith(name, "Scenario_"), readdir(scen_dir))
            tdr_scen = joinpath(scen_dir, scen, "TDR_Results")
            isdir(tdr_scen) && rm(tdr_scen; recursive = true)
        end
    end

    println("  Running run_timedomainreduction!() ...")
    t_tdr = @elapsed run_timedomainreduction!(CASE)
    println("   TDR time = $(round(t_tdr, digits=1)) s")

    println("  Clearing previous CE results...")
    isdir(RESULTS_ROOT) && rm(RESULTS_ROOT; recursive = true)

    println("  Running CE solve via Run_cluster.jl...")
    run_script = joinpath(CASE, "Run_cluster.jl")
    t_ce = @elapsed include(run_script)
    println("   CE time = $(round(t_ce, digits=1)) s")

    total_runtime = t_tdr + t_ce
    println("   Total runtime = $(round(total_runtime, digits=1)) s")

    println("  Extracting cost...")
    cost = read_cost_for_run(RESULTS_ROOT)
    println("   Total cost = ", cost)

    println("  Archiving CE results...")
    archive_ce_results(rp)

    push!(results_df, (
        Rep_Periods_per_Year = rp,
        TDR_Runtime_s = t_tdr,
        CE_Runtime_s = t_ce,
        Total_Runtime_s = total_runtime,
        Total_System_Cost = cost,
    ))
end

println("\nFinal sweep summary:")
println(results_df)

CSV.write(OUT_SUMMARY_CSV, results_df)
println("\nSummary written to: $OUT_SUMMARY_CSV")
println("Sweep complete.")
