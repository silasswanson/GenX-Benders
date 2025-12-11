# TDR_sweep_reports.jl
#
# For each RP_* folder under GenX-Brazil/TDR_sweep_results, this script:
#   1. Reconstructs the 8-year total system load from RP_k/TDR_Results
#   2. Computes global + per-year RMSE (saved to CSV in each RP folder)
#   3. Plots 8-year LDC + per-year LDCs (saved in each RP folder)
#   4. Plots the difference LDC (recon - original)
#
using CSV, DataFrames, Statistics, Plots

# ------------------ GLOBAL SETTINGS ------------------ #

const CASE   = joinpath(@__DIR__, "GenX-Brazil")
const NYEARS = 8
const PERIODS_PER_YEAR = 52

const ORIG_LOAD_PATH = joinpath(CASE, "Load_data.csv")
const SWEEP_ROOT     = joinpath(CASE, "TDR_sweep_results")

println("Case directory:   $CASE")
println("Original load:    $ORIG_LOAD_PATH")
println("Sweep root:       $SWEEP_ROOT")

isdir(SWEEP_ROOT) || error("No TDR_sweep_results folder at $SWEEP_ROOT")

# ------------------ 1. Read original total system load ------------------ #

ld = CSV.read(ORIG_LOAD_PATH, DataFrame)

zone_cols = filter(colname -> occursin("Load_MW_z", String(colname)), names(ld))
@assert !isempty(zone_cols) "No Load_MW_z* columns found in original Load_data.csv"

println("Detected zone load columns (original): ", zone_cols)

orig = reduce(+, (ld[!, col] for col in zone_cols))
T_total = length(orig)
periods_total = NYEARS * PERIODS_PER_YEAR

@assert T_total % periods_total == 0 "Total hours ($T_total) not divisible by periods_total ($periods_total)"
H_per_period = T_total ÷ periods_total

println("Total hours       = $T_total")
println("Periods total     = $periods_total")
println("Hours per period  = $H_per_period")

hours_per_year = T_total ÷ NYEARS
@assert hours_per_year * NYEARS == T_total

# ------------------ Helper: list RP_* folders ------------------ #

rp_dirs = sort(filter(name -> startswith(name, "RP_"), readdir(SWEEP_ROOT)))

if isempty(rp_dirs)
    error("No RP_* folders found under $SWEEP_ROOT")
end

println("\nFound RP folders: ", rp_dirs)

# Summary of global RMSE across RPs
summary_rows = DataFrame(RP_per_year = Int[], Global_RMSE = Float64[])

# ------------------ 2. Loop over each RP_* folder ------------------ #

for rp_dir_name in rp_dirs
    println("\n====================================================")
    println(" Processing sweep folder: $rp_dir_name")
    println("====================================================")

    rp_dir = joinpath(SWEEP_ROOT, rp_dir_name)
    results_root = joinpath(rp_dir, "TDR_Results")
    isdir(results_root) || error("Missing TDR_Results in $rp_dir")

    # Paths inside this RP_* set
    tdr_map_path  = joinpath(results_root, "Period_map.csv")
    tdr_load_path = joinpath(results_root, "Load_data.csv")

    # --- read Period_map and determine rp_per_year ---
    pm = CSV.read(tdr_map_path, DataFrame)
    names_pm = String.(names(pm))

    period_idx_name = first(filter(nm -> endswith(nm, "Period_Index"), names_pm))
    rep_idx_name    = first(filter(nm -> endswith(nm, "Rep_Period_Index"), names_pm))

    period_idx = Int.(pm[!, Symbol(period_idx_name)])
    rep_idx    = Int.(pm[!, Symbol(rep_idx_name)])

    @assert length(period_idx) == periods_total "Period_map row count mismatch in $rp_dir_name"
    @assert maximum(period_idx) == periods_total "Unexpected max Period_Index in $rp_dir_name"

    n_rep_total = maximum(rep_idx)
    @assert n_rep_total % NYEARS == 0 "Max Rep_Period_Index ($n_rep_total) not divisible by NYEARS in $rp_dir_name"
    rp_per_year = n_rep_total ÷ NYEARS

    # sanity-check folder name vs inferred RP
    parts = split(rp_dir_name, "_")
    if length(parts) == 2
        folder_rp = try
            parse(Int, parts[2])
        catch
            nothing
        end
        if folder_rp !== nothing && folder_rp != rp_per_year
            println("  WARNING: folder name RP_$folder_rp but Period_map implies RP_$rp_per_year")
        end
    end

    println("  Detected representative periods per year = $rp_per_year (total rep = $n_rep_total)")

    # --- read reduced Load_data for this RP ---
    tdr_ld = CSV.read(tdr_load_path, DataFrame)

    zone_cols_red = filter(colname -> occursin("Load_MW_z", String(colname)), names(tdr_ld))
    @assert !isempty(zone_cols_red) "No Load_MW_z* columns found in reduced Load_data.csv for $rp_dir_name"

    @assert nrow(tdr_ld) % H_per_period == 0 "Reduced Load_data row count $(nrow(tdr_ld)) not divisible by H_per_period in $rp_dir_name"
    @assert nrow(tdr_ld) == n_rep_total * H_per_period "Reduced Load_data row count mismatch in $rp_dir_name"

    reduced_total = reduce(+, (tdr_ld[!, col] for col in zone_cols_red))

    # representative period profiles, k = 1..n_rep_total
    function get_profile(k::Int)
        start_row = (k - 1) * H_per_period + 1
        end_row   = k * H_per_period
        return @view reduced_total[start_row:end_row]
    end

    # --- reconstruct full 8-year series ---
    recon = zeros(Float64, T_total)

    for r in eachindex(period_idx)
        p  = period_idx[r]
        rp = rep_idx[r]

        @assert 1 <= rp <= n_rep_total "Rep_Period_Index out of range: $rp in $rp_dir_name"

        start_hour = (p - 1) * H_per_period + 1
        end_hour   = p * H_per_period
        recon[start_hour:end_hour] .= get_profile(rp)
    end

    # ====================================================
    # 2a. RMSE: global + per year, save to CSV
    # ====================================================

    err = recon .- orig
    rmse_global = sqrt(sum(err .^ 2) / length(err))

    println("  Global RMSE over 8 years = $rmse_global")

    scopes = String[]
    rmses  = Float64[]

    push!(scopes, "Global")
    push!(rmses, rmse_global)

    for y in 1:NYEARS
        hstart = (y - 1) * hours_per_year + 1
        hend   = y * hours_per_year
        err_y = err[hstart:hend]
        rmse_y = sqrt(sum(err_y .^ 2) / length(err_y))
        println("  Year $y RMSE              = $rmse_y")
        push!(scopes, "Year $y")
        push!(rmses, rmse_y)
    end

    rmse_df = DataFrame(
        Rep_Periods_per_Year = fill(rp_per_year, length(scopes)),
        Scope = scopes,
        RMSE  = rmses,
    )

    rmse_csv_path = joinpath(rp_dir, "RMSE_SystemLoad_RP_$(rp_per_year).csv")
    CSV.write(rmse_csv_path, rmse_df)
    println("  Saved RMSE results: $rmse_csv_path")

    push!(summary_rows, (RP_per_year = rp_per_year, Global_RMSE = rmse_global))

    # ====================================================
    # 2b. Global + per-year LDC plots
    # ====================================================

    ldcout_dir_year = joinpath(rp_dir, "LDC_per_year")
    isdir(ldcout_dir_year) || mkpath(ldcout_dir_year)

    ldc_orig  = sort(orig,  rev = true)
    ldc_recon = sort(recon, rev = true)

    p_global = plot(
        ldc_orig,
        label = "Original",
        xlabel = "Hour (sorted)",
        ylabel = "Load (MW)",
        title  = "8-Year Load Duration Curve (RP = $rp_per_year)",
        linewidth = 2,
        size = (900, 500),
    )
    plot!(p_global, ldc_recon, label = "Reconstructed", linewidth = 2, linestyle = :dash)

    png_global = joinpath(rp_dir, "LDC_8yr_System_RP_$rp_per_year.png")
    savefig(p_global, png_global)
    println("  Saved global LDC: $png_global")

    # per-year LDCs
    for y in 1:NYEARS
        hstart = (y - 1) * hours_per_year + 1
        hend   = y * hours_per_year

        orig_y  = orig[hstart:hend]
        recon_y = recon[hstart:hend]

        ldc_orig_y  = sort(orig_y,  rev = true)
        ldc_recon_y = sort(recon_y, rev = true)

        p_y = plot(
            ldc_orig_y,
            label = "Original",
            xlabel = "Hour (sorted)",
            ylabel = "Load (MW)",
            title  = "Year $y LDC (RP = $rp_per_year)",
            linewidth = 2,
            size = (800, 450),
        )
        plot!(p_y, ldc_recon_y, label = "Reconstructed", linewidth = 2, linestyle = :dash)

        png_y = joinpath(ldcout_dir_year, "LDC_Year_$(y)_System_RP_$rp_per_year.png")
        savefig(p_y, png_y)
        println("  Saved per-year LDC: $png_y")
    end

    # ====================================================
    # 2c. LDC difference plot (recon - original)
    # ====================================================

    delta = ldc_recon .- ldc_orig

    pΔ = plot(
        delta,
        xlabel = "Hour (sorted)",
        ylabel = "Δ Load (MW) = recon - original",
        title  = "Difference in 8-Year LDC (RP = $rp_per_year)",
        linewidth = 1.8,
        legend = false,
        size = (900, 500),
    )
    plot!(pΔ, fill(0.0, length(delta)), linewidth = 1, linestyle = :dash)

    png_diff = joinpath(rp_dir, "LDC_8yr_System_Diff_RP_$rp_per_year.png")
    savefig(pΔ, png_diff)
    println("  Saved LDC difference plot: $png_diff")

    println("  Min Δ = $(minimum(delta)), Max Δ = $(maximum(delta))")
end

# ------------------ 3. Save sweep summary ------------------ #

summary_csv = joinpath(SWEEP_ROOT, "RMSE_SystemLoad_summary.csv")
CSV.write(summary_csv, summary_rows)
println("\nSweep summary written to: $summary_csv")

println("\nAll sweep reports complete.")
