# TDR_plot_RP20.jl
#
# Plot, for each scenario, the mapping from original weeks to
# representative weeks for RP = 20.

using CSV
using DataFrames
using Plots

# ----------------- USER SETTINGS ----------------- #

const CASE = joinpath(@__DIR__, "GenX-Brazil")
const RP   = 20          # number of representative periods used
const NSCEN = 8          # 8 weather years / scenarios

# Base dir where the sweep results live
const SWEEP_DIR = joinpath(CASE, "TDR_sweep_results", "RP_$RP")

# Make a folder for plots
const PLOTS_DIR = joinpath(SWEEP_DIR, "plots")
isdir(PLOTS_DIR) || mkpath(PLOTS_DIR)

println("Reading TDR results from: $SWEEP_DIR")
println("Saving plots in:         $PLOTS_DIR")

plots = Any[]

for scen in 1:NSCEN
    scen_name = "Scenario_$scen"
    pm_path = joinpath(SWEEP_DIR, scen_name, "TDR_Results", "Period_map.csv")

    if !isfile(pm_path)
        @warn "No Period_map.csv for $scen_name at $pm_path – skipping"
        continue
    end

    println("  Reading $pm_path")
    pm = CSV.read(pm_path, DataFrame)

    # --- Figure out which columns are "original week" and "representative week" ---

    # Try to guess column names, with fallbacks
    names_pm = names(pm)

    # heuristic for original period column
    week_col = if :Week in names_pm
        :Week
    elseif :OriginalPeriod in names_pm
        :OriginalPeriod
    else
        names_pm[1]  # fallback: first column
    end

    # heuristic for representative period column
    rep_col = if :RepresentativePeriod in names_pm
        :RepresentativePeriod
    elseif :Rep_Period in names_pm
        :Rep_Period
    else
        names_pm[2]  # fallback: second column
    end

    weeks = pm[!, week_col]
    reps  = pm[!, rep_col]

    # --- Plot mapping: original week index vs representative week index ---

    p = scatter(
        weeks, reps;
        xlabel = "Original week index",
        ylabel = "Representative week",
        title  = "Scenario $scen (RP = $RP)",
        legend = false,
        markersize = 5,
    )

    push!(plots, p)

    # also save individual PNG per scenario
    png_path = joinpath(PLOTS_DIR, "Scenario_$(scen)_RP_$RP")
    println("    Saving $png_path.png")
    savefig(p, png_path)
end

# Optional: combined multi-panel figure (4x2 grid) if we have enough plots
if !isempty(plots)
    grid = plot(plots...; layout = (4, 2), size = (1200, 800))
    grid_path = joinpath(PLOTS_DIR, "All_scenarios_RP_$RP")
    println("Saving combined plot to $grid_path.png")
    savefig(grid, grid_path)
end

println("Done.")
