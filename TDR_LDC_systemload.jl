# TDR_LDC_systemload.jl
#
# Plots load duration curves (LDCs) comparing:
#   - original 8-year total system load
#   - reconstructed load from TDR_Results
#
# Saves to BOTH:
#   GenX-Brazil/TDR_Results/
#   GenX-Brazil/TDR_sweep_results/RP_<rp>/
#
using CSV, DataFrames, Statistics, Plots

# ------------------ SETTINGS ------------------ #

const CASE   = joinpath(@__DIR__, "GenX-Brazil")
const NYEARS = 8
const PERIODS_PER_YEAR = 52

const ORIG_LOAD_PATH = joinpath(CASE, "Load_data.csv")
const TDR_MAP_PATH   = joinpath(CASE, "TDR_Results", "Period_map.csv")
const TDR_LOAD_PATH  = joinpath(CASE, "TDR_Results", "Load_data.csv")

println("CASE          = $CASE")
println("Original load = $ORIG_LOAD_PATH")
println("TDR map       = $TDR_MAP_PATH")
println("TDR load      = $TDR_LOAD_PATH")

# Canonical output directories
const OUTDIR_GLOBAL = joinpath(CASE, "TDR_Results")
const OUTDIR_YEARS  = joinpath(CASE, "TDR_Results", "LDC_per_year")
isdir(OUTDIR_GLOBAL) || mkpath(OUTDIR_GLOBAL)
isdir(OUTDIR_YEARS)  || mkpath(OUTDIR_YEARS)

# ------------------ 1. Original total system load ------------------ #

ld = CSV.read(ORIG_LOAD_PATH, DataFrame)

zone_cols = filter(col -> occursin("Load_MW_z", String(col)), names(ld))
@assert !isempty(zone_cols) "No Load_MW_z* columns found in original Load_data.csv"

println("Detected zone load columns (original): ", zone_cols)

orig = reduce(+, (ld[!, col] for col in zone_cols))

T_total = length(orig)
periods_total = NYEARS * PERIODS_PER_YEAR

@assert T_total % periods_total == 0 "Total hours ($T_total) not divisible by periods_total ($periods_total)"
H_per_period = T_total ÷ periods_total

# TDR_LDC_systemload.jl
#
# Plots load duration curves (LDCs) comparing:
#   - original 8-year total system load
#   - reconstructed load from TDR_Results
#
# Saves to BOTH:
#   GenX-Brazil/TDR_Results/
#   GenX-Brazil/TDR_sweep_results/RP_<rp>/
#
using CSV, DataFrames, Statistics, Plots

# ------------------ SETTINGS ------------------ #

const CASE   = joinpath(@__DIR__, "GenX-Brazil")
const NYEARS = 8
const PERIODS_PER_YEAR = 52

const ORIG_LOAD_PATH = joinpath(CASE, "Load_data.csv")
const TDR_MAP_PATH   = joinpath(CASE, "TDR_Results", "Period_map.csv")
const TDR_LOAD_PATH  = joinpath(CASE, "TDR_Results", "Load_data.csv")

println("CASE          = $CASE")
println("Original load = $ORIG_LOAD_PATH")
println("TDR map       = $TDR_MAP_PATH")
println("TDR load      = $TDR_LOAD_PATH")

# Canonical output directories
const OUTDIR_GLOBAL = joinpath(CASE, "TDR_Results")
const OUTDIR_YEARS  = joinpath(CASE, "TDR_Results", "LDC_per_year")
isdir(OUTDIR_GLOBAL) || mkpath(OUTDIR_GLOBAL)
isdir(OUTDIR_YEARS)  || mkpath(OUTDIR_YEARS)

# ------------------ 1. Original total system load ------------------ #

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

# ------------------ 2. Read Period_map ------------------ #

pm = CSV.read(TDR_MAP_PATH, DataFrame)
names_pm = String.(names(pm))
println("Period_map columns: ", names_pm)

period_idx_name = first(filter(nm -> endswith(nm, "Period_Index"), names_pm))
rep_idx_name    = first(filter(nm -> endswith(nm, "Rep_Period_Index"), names_pm))

period_idx_col = Symbol(period_idx_name)
rep_idx_col    = Symbol(rep_idx_name)

period_raw = pm[!, period_idx_col]
rep_raw    = pm[!, rep_idx_col]

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
    error("Fix Period_map.csv before plotting LDCs.")
end

period_idx = Int.(period_raw)
rep_idx    = Int.(rep_raw)

@assert length(period_idx) == periods_total "Period_map row count mismatch"
@assert maximum(period_idx) == periods_total "Unexpected max Period_Index"

n_rep_total = maximum(rep_idx)
@assert n_rep_total % NYEARS == 0 "Max Rep_Period_Index ($n_rep_total) not divisible by NYEARS ($NYEARS)"
rp_per_year = n_rep_total ÷ NYEARS

println("Detected representative periods per year = $rp_per_year (total rep periods = $n_rep_total)")

# Sweep folder for saving copies
sweep_folder      = joinpath(CASE, "TDR_sweep_results", "RP_$(rp_per_year)")
sweep_folder_year = joinpath(sweep_folder, "LDC_per_year")
isdir(sweep_folder)      || mkpath(sweep_folder)
isdir(sweep_folder_year) || mkpath(sweep_folder_year)

# ------------------ 3. Read reduced TDR load ------------------ #

tdr_ld = CSV.read(TDR_LOAD_PATH, DataFrame)

zone_cols_red = filter(colname -> occursin("Load_MW_z", String(colname)), names(tdr_ld))
@assert !isempty(zone_cols_red) "No Load_MW_z* columns found in reduced Load_data.csv"

println("Detected zone load columns (reduced): ", zone_cols_red)

@assert nrow(tdr_ld) % H_per_period == 0 "Reduced Load_data row count $(nrow(tdr_ld)) not divisible by H_per_period ($H_per_period)"
@assert nrow(tdr_ld) == n_rep_total * H_per_period "Reduced Load_data row count mismatch (got $(nrow(tdr_ld)), expected $(n_rep_total * H_per_period))"

reduced_total = reduce(+, (tdr_ld[!, col] for col in zone_cols_red))

# representative period profiles, k = 1..n_rep_total
function get_profile(k::Int)
    start_row = (k - 1) * H_per_period + 1
    end_row   = k * H_per_period
    return @view reduced_total[start_row:end_row]
end

# ------------------ 4. Reconstruct full 8-year total system load ------------------ #

recon = zeros(Float64, T_total)

for r in eachindex(period_idx)
    p  = period_idx[r]
    rp = rep_idx[r]

    @assert 1 <= rp <= n_rep_total "Rep_Period_Index out of range: $rp"

    start_hour = (p - 1) * H_per_period + 1
    end_hour   = p * H_per_period
    recon[start_hour:end_hour] .= get_profile(rp)
end

# ------------------ 5. Global 8-year LDC ------------------ #

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

# Save main canonical copy
png_global = joinpath(OUTDIR_GLOBAL, "LDC_8yr_System_RP_$rp_per_year.png")
savefig(p_global, png_global)
println("\nSaved global LDC: $png_global")

# Save sweep copy
png_global_sweep = joinpath(sweep_folder, "LDC_8yr_System_RP_$rp_per_year.png")
savefig(p_global, png_global_sweep)
println("Saved sweep global LDC: $png_global_sweep")

# ------------------ 6. Per-year LDCs ------------------ #

hours_per_year = T_total ÷ NYEARS
@assert hours_per_year * NYEARS == T_total

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

    # Canonical save
    png_y = joinpath(OUTDIR_YEARS, "LDC_Year_$(y)_System_RP_$rp_per_year.png")
    savefig(p_y, png_y)
    println("Saved per-year LDC: $png_y")

    # Sweep folder save
    png_y_sweep = joinpath(sweep_folder_year, "LDC_Year_$(y)_System_RP_$rp_per_year.png")
    savefig(p_y, png_y_sweep)
    println("Saved sweep per-year LDC: $png_y_sweep")
end

println("\nDone.")
