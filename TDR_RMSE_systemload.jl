# TDR_RMSE_systemload.jl
#
# Computes RMSE between:
#   - original 8-year total system load (from GenX-Brazil/Load_data.csv)
#   - reconstructed load using TDR_Results (Period_map + reduced Load_data)
#
# Saves results to:
#   GenX-Brazil/TDR_Results/RMSE_SystemLoad_RP_<rp_per_year>.csv
#   GenX-Brazil/TDR_sweep_results/RP_<rp_per_year>/RMSE_SystemLoad_RP_<rp_per_year>.csv
#
# Notes:
# - Uses ONLY Period_map.csv for week -> rep-period mapping.
# - Assumes reduced Load_data.csv has rep periods stacked in order:
#       rep 1: rows 1:168
#       rep 2: rows 169:336
#       ...
# - Sums load across all Load_MW_z* columns.

using CSV, DataFrames, Statistics

# ------------------ SETTINGS ------------------ #

const CASE   = joinpath(@__DIR__, "GenX-Brazil")
const NYEARS = 8
const PERIODS_PER_YEAR = 52   # weeks per year

const ORIG_LOAD_PATH = joinpath(CASE, "Load_data.csv")
const TDR_MAP_PATH   = joinpath(CASE, "TDR_Results", "Period_map.csv")
const TDR_LOAD_PATH  = joinpath(CASE, "TDR_Results", "Load_data.csv")
const OUTDIR         = joinpath(CASE, "TDR_Results")

println("CASE          = $CASE")
println("Original load = $ORIG_LOAD_PATH")
println("TDR map       = $TDR_MAP_PATH")
println("TDR load      = $TDR_LOAD_PATH")

isdir(OUTDIR) || mkpath(OUTDIR)

# ------------------ 1. Original total system load ------------------ #

ld = CSV.read(ORIG_LOAD_PATH, DataFrame)

zone_cols = filter(col -> occursin("Load_MW_z", String(col)), names(ld))
@assert !isempty(zone_cols) "No Load_MW_z* columns found in original Load_data.csv"

println("Detected zone load columns (original): ", zone_cols)

orig = reduce(+, (ld[!, col] for col in zone_cols))  # length T_total

T_total = length(orig)
periods_total = NYEARS * PERIODS_PER_YEAR

@assert T_total % periods_total == 0 "Total hours ($T_total) not divisible by periods_total ($periods_total)"
H_per_period = T_total ÷ periods_total

println("Total hours       = $T_total")
println("Periods total     = $periods_total")
println("Hours per period  = $H_per_period")

# ------------------ 2. Read Period_map (week -> rep-period mapping) ------------------ #

pm = CSV.read(TDR_MAP_PATH, DataFrame)
names_pm = String.(names(pm))
println("Period_map columns: ", names_pm)

# robustly find Period_Index and Rep_Period_Index columns by suffix
period_idx_name = first(filter(n -> endswith(n, "Period_Index"), names_pm))
rep_idx_name    = first(filter(n -> endswith(n, "Rep_Period_Index"), names_pm))

period_idx_col = Symbol(period_idx_name)
rep_idx_col    = Symbol(rep_idx_name)

period_raw = pm[!, period_idx_col]
rep_raw    = pm[!, rep_idx_col]

# Strict: do not proceed if there are missing mapping entries
miss_period = findall(ismissing, period_raw)
miss_rep    = findall(ismissing, rep_raw)

if !isempty(miss_period) || !isempty(miss_rep)
    println("\nERROR: Missing values detected in Period_map.csv:")
    if !isempty(miss_period)
        println("  $(period_idx_name) missing at rows: ", miss_period)
    end
    if !isempty(miss_rep)
        println("  $(rep_idx_name) missing at rows: ", miss_rep)
    end
    error("Fix Period_map.csv before computing RMSE.")
end

period_idx = Int.(period_raw)
rep_idx    = Int.(rep_raw)

@assert length(period_idx) == periods_total "Period_map row count mismatch (got $(length(period_idx)), expected $periods_total)"
@assert maximum(period_idx) == periods_total "Unexpected max Period_Index (got $(maximum(period_idx)), expected $periods_total)"

n_rep_total = maximum(rep_idx)            # e.g., 160
@assert n_rep_total % NYEARS == 0 "Max Rep_Period_Index ($n_rep_total) not divisible by NYEARS ($NYEARS)"
rp_per_year = n_rep_total ÷ NYEARS        # e.g., 20

println("Detected representative periods per year = $rp_per_year (total rep periods = $n_rep_total)")

# ------------------ 3. Read reduced TDR load ------------------ #

tdr_ld = CSV.read(TDR_LOAD_PATH, DataFrame)

zone_cols_red = filter(col -> occursin("Load_MW_z", String(col)), names(tdr_ld))
@assert !isempty(zone_cols_red) "No Load_MW_z* columns found in reduced Load_data.csv"

println("Detected zone load columns (reduced): ", zone_cols_red)

@assert nrow(tdr_ld) % H_per_period == 0 "Reduced Load_data row count ($(nrow(tdr_ld))) not divisible by H_per_period ($H_per_period)"
@assert nrow(tdr_ld) == n_rep_total * H_per_period "Reduced Load_data row count mismatch (got $(nrow(tdr_ld)), expected $(n_rep_total * H_per_period))"

# Total system load in reduced space (length = n_rep_total * H_per_period)
reduced_total = reduce(+, (tdr_ld[!, col] for col in zone_cols_red))

# Helper: profile for representative period k (k = 1..n_rep_total)
get_profile(k) = @view reduced_total[(k-1)*H_per_period + 1 : k*H_per_period]

# ------------------ 4. Reconstruct full 8-year total system load ------------------ #

recon = zeros(Float64, T_total)

for r in eachindex(period_idx)
    p  = period_idx[r]    # original period index 1..416
    rp = rep_idx[r]       # representative period index 1..n_rep_total

    @assert 1 <= rp <= n_rep_total "Rep_Period_Index out of range: $rp"

    start_hour = (p - 1) * H_per_period + 1
    end_hour   = p * H_per_period
    recon[start_hour:end_hour] .= get_profile(rp)
end

# ------------------ 5. RMSE: global + per year ------------------ #

err = recon .- orig

se_total = sum(err .^ 2)
rmse_global = sqrt(se_total / length(err))

println("\n========== SYSTEM LOAD RMSE ==========")
println("Rep periods / year       = $rp_per_year")
println("Global RMSE over 8 years = $rmse_global\n")

hours_per_year = T_total ÷ NYEARS
@assert hours_per_year * NYEARS == T_total

# collect results for CSV
years = String[]
rmses = Float64[]

# global row
push!(years, "Global")
push!(rmses, rmse_global)

for y in 1:NYEARS
    hstart = (y - 1) * hours_per_year + 1
    hend   = y * hours_per_year
    err_y = err[hstart:hend]
    rmse_y = sqrt(sum(err_y .^ 2) / length(err_y))
    println("Year $y RMSE              = $rmse_y")

    push!(years, "Year $y")
    push!(rmses, rmse_y)
end

println("\nDone computing RMSEs.")

# ------------------ 6. Save RMSE results to CSV (two locations) ------------------ #

rmse_df = DataFrame(
    Rep_Periods_per_Year = fill(rp_per_year, length(years)),
    Scope   = years,
    RMSE    = rmses,
)

# Primary GenX output location
out_csv_main = joinpath(OUTDIR, "RMSE_SystemLoad_RP_$(rp_per_year).csv")
CSV.write(out_csv_main, rmse_df)
println("Saved RMSE results to main location: $out_csv_main")

# Sweep tracking folder location
sweep_folder = joinpath(CASE, "TDR_sweep_results", "RP_$(rp_per_year)")
isdir(sweep_folder) || mkpath(sweep_folder)

out_csv_sweep = joinpath(sweep_folder, "RMSE_SystemLoad_RP_$(rp_per_year).csv")
CSV.write(out_csv_sweep, rmse_df)
println("Saved RMSE results to sweep folder: $out_csv_sweep")
