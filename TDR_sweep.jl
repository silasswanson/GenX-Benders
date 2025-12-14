# TDR_sweep.jl
#
# Sweep time-domain reduction over several choices of
# representative periods per year, running:
#
#    RP_per_year ∈ [20, 15, 10, 5, 1]
#
# For each RP:
#   1. Edit time_domain_reduction_settings.yml
#   2. Clear old TDR output in case folder
#   3. Run run_timedomainreduction!(CASE)
#   4. Save / archive outputs into:
#
#        TDR_sweep_results/RP_<RP>/TDR_Results
#        TDR_sweep_results/RP_<RP>/Scenario_<i>/TDR_Results
#
using YAML
include("src/GenX-Benders.jl")   # defines run_timedomainreduction!

# -------- user settings -------- #

const CASE = joinpath(@__DIR__, "GenX-Brazil")

# Representative periods per year to test 
rep_periods = [52,50,40,30,20,15,10,5,1]

# Path to TDR settings file
const TDR_SETTINGS_PATH = joinpath(CASE, "Settings", "time_domain_reduction_settings.yml")

# Root folder for sweep results #### change this based on what you are running####
const SWEEP_OUTDIR = joinpath(CASE, "TDR_sweep_results_global")
isdir(SWEEP_OUTDIR) || mkpath(SWEEP_OUTDIR)

println("Case directory:        $CASE")
println("TDR settings file:     $TDR_SETTINGS_PATH")
println("Sweep output root:     $SWEEP_OUTDIR")
println("Rep periods per year:  ", rep_periods)

# Load settings once
tdr_settings = YAML.load_file(TDR_SETTINGS_PATH)

for (k, rp) in pairs(rep_periods)

    println("\n========================================")
    println(" Sweep $(k) / $(length(rep_periods)): RP_per_year = $rp")
    println("========================================")

    # 1. UPDATE YAML SETTINGS FOR TDR
    tdr_settings["MinPeriods"] = rp
    tdr_settings["MaxPeriods"] = rp
    YAML.write_file(TDR_SETTINGS_PATH, tdr_settings)
    println("  Updated MinPeriods/MaxPeriods in settings to $rp")

    # 2. REMOVE PREVIOUS TDR OUTPUT FROM CASE DIRECTORY

    # Delete the aggregate 8-year TDR_Results folder (so we don't mix runs)
    tdr_top = joinpath(CASE, "TDR_Results")
    if isdir(tdr_top)
        println("  Removing previous 8-year TDR_Results at $tdr_top")
        rm(tdr_top; recursive = true)
    end

    # Delete per-scenario TDR results (fresh sweep)
    scen_dir = joinpath(CASE, "Scenarios")
    if isdir(scen_dir)
        for scen_name in filter(n -> startswith(n, "Scenario_"), readdir(scen_dir))
            tdr_scen = joinpath(scen_dir, scen_name, "TDR_Results")
            if isdir(tdr_scen)
                println("  Removing per-scenario TDR_Results at $tdr_scen")
                rm(tdr_scen; recursive = true)
            end
        end
    end

    # 3. RUN TDR FOR THIS RP (writes into CASE/*)
    println("  Running run_timedomainreduction!(CASE) for RP = $rp")
    t = @elapsed run_timedomainreduction!(CASE)
    println("  Finished TDR call in $(round(t; digits = 1)) seconds")

    # 4. ARCHIVE RESULTS FOR THIS RP

    # Create output folder for this RP (RP_<rp>)
    out_rp_dir = joinpath(SWEEP_OUTDIR, "RP_$(rp)")
    if isdir(out_rp_dir)
        println("  Removing old sweep folder at $out_rp_dir")
        rm(out_rp_dir; recursive = true)
    end
    mkpath(out_rp_dir)

    # 4a. SAVE 8-YEAR AGGREGATED OUTPUT FOR THIS RP
    # Copy the case-level TDR_Results here (Period_map, reduced Load_data, etc.)
    if isdir(tdr_top)
        dst_top = joinpath(out_rp_dir, "TDR_Results")
        println("  SAVING aggregate 8-year TDR results -> $dst_top")
        cp(tdr_top, dst_top; force = true)
    else
        println("  WARNING: expected TDR_Results missing at $tdr_top")
    end

    # 4b. SAVE YEAR-SPECIFIC RESULTS FOR EACH SCENARIO
    if isdir(scen_dir)
        for scen_name in filter(n -> startswith(n, "Scenario_"), readdir(scen_dir))
            src = joinpath(scen_dir, scen_name, "TDR_Results")
            if isdir(src)
                dst = joinpath(out_rp_dir, scen_name, "TDR_Results")
                println("  SAVING per-scenario TDR results -> $dst")
                mkpath(dirname(dst))
                cp(src, dst; force = true)
            end
        end
    end

    # Optional: save runtime info per RP
    rt_file = joinpath(out_rp_dir, "runtime_seconds.txt")
    open(rt_file, "w") do io
        write(io, "Runtime_seconds = $(round(t; digits = 3))\n")
    end
    println("  Saved runtime info at $rt_file")

    println("  Sweep for RP = $rp ARCHIVED under $out_rp_dir")
end

println("\nAll TDR sweeps complete.")
println("Sweep folders written under: $SWEEP_OUTDIR")
