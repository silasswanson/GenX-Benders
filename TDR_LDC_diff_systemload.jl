# TDR_LDC_diff_systemload.jl
#
# Plots the difference between reconstructed and original
# 8-year system-load LDC:
#   Δ(rank) = LDC_recon(rank) - LDC_orig(rank)

using CSV, DataFrames, Statistics, Plots

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

# Output dirs
const OUTDIR_MAIN  = joinpath(CASE, "TDR_Results")
isdir(OUTDIR_MAIN) || mkpath(OUTDIR_MAIN)

# ------------------ 1. Original total system load ------------------ #

ld = CSV.read(ORIG_LOAD_PATH, DataFrame)
zone_cols = filter(colname -> occursin("Load_MW_z", String(colname)), names(ld))
@assert !isempty(zone_cols) "No Load_MW_z* columns found in original Load_data.csv"

orig = reduce(+, (ld[!, col] for col in zone_cols))

T_total = length(orig)
periods_total = NYEARS * PERIODS_PER_YEAR
@assert T_total % periods_total == 0
H_per_period = T_total ÷ periods_total

# ------------------ 2. Period_map and reduced load ------------------ #

pm = CSV.read(TDR_MAP_PATH, DataFrame)
names_pm = String.(names(pm))

period_idx_name = first(filter(nm -> endswith(nm, "Period_Index"), names_pm))
rep_idx_name    = first(filter(nm -> endswith(nm, "Rep_Period_Index"), names_pm))

period_idx_col = Symbol(period_idx_name)
rep_idx_col    = Symbol(rep_idx_name)

period_idx = Int.(pm[!, period_idx_col])
rep_idx    = Int.(pm[!, rep_idx_col])

n_rep_total = maximum(rep_idx)
@assert n_rep_total % NYEARS == 0
rp_per_year = n_rep_total ÷ NYEARS

tdr_ld = CSV.read(TDR_LOAD_PATH, DataFrame)
zone_cols_red = filter(colname -> occursin("Load_MW_z", String(colname)), names(tdr_ld))
@assert !isempty(zone_cols_red)

# TDR_LDC_diff_systemload.jl
#
# Plots difference between reconstructed and original
# 8-year system-load LDC:
#   Δ(rank) = LDC_recon(rank) - LDC_orig(rank)
#
# Saves output under BOTH:
#   - TDR_Results/
#   - TDR_sweep_results/RP_x/

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

# Canonical output folder
const OUTDIR_MAIN = joinpath(CASE, "TDR_Results")
isdir(OUTDIR_MAIN) || mkpath(OUTDIR_MAIN)

# ------------------ 1. Original total system load ------------------ #

ld = CSV.read(ORIG_LOAD_PATH, DataFrame)

zone_cols = filter(colname -> occursin("Load_MW_z", String(colname)), names(ld))
@assert !isempty(zone_cols) "No Load_MW_z* columns in original Load_data.csv"

orig = reduce(+, (ld[!, col] for col in zone_cols))

T_total = length(orig)
periods_total = NYEARS * PERIODS_PER_YEAR
@assert T_total % periods_total == 0

H_per_period = T_total ÷ periods_total

# ------------------ 2. Read Period_map and reduced load ------------------ #

pm = CSV.read(TDR_MAP_PATH, DataFrame)
names_pm = String.(names(pm))

period_idx_name = first(filter(nm -> endswith(nm, "Period_Index"), names_pm))
rep_idx_name    = first(filter(nm -> endswith(nm, "Rep_Period_Index"), names_pm))

period_idx = Int.(pm[!, Symbol(period_idx_name)])
rep_idx    = Int.(pm[!, Symbol(rep_idx_name)])

n_rep_total = maximum(rep_idx)

@assert n_rep_total % NYEARS == 0
rp_per_year = n_rep_total ÷ NYEARS

println("Detected rep periods per year = $rp_per_year (total rep = $n_rep_total)")

# Sweep folder
sweep_folder = joinpath(CASE, "TDR_sweep_results", "RP_$(rp_per_year)")
isdir(sweep_folder) || mkpath(sweep_folder)

tdr_ld = CSV.read(TDR_LOAD_PATH, DataFrame)

zone_cols_red = filter(colname -> occursin("Load_MW_z", String(colname)), names(tdr_ld))
@assert !isempty(zone_cols_red)

@assert nrow(tdr_ld) == (n_rep_total * H_per_period)

reduced_total = reduce(+, (tdr_ld[!, col] for col in zone_cols_red))

function get_profile(k::Int)
    start_row = (k - 1) * H_per_period + 1
    end_row   = k * H_per_period
    return @view reduced_total[start_row:end_row]
end

# ------------------ 3. Reconstruct ------------------ #

recon = zeros(Float64, T_total)

for r in eachindex(period_idx)
    p  = period_idx[r]
    rp = rep_idx[r]
    start_hour = (p - 1) * H_per_period + 1
    end_hour   = p * H_per_period
    recon[start_hour:end_hour] .= get_profile(rp)
end

# ------------------ 4. Difference LDC ------------------ #

ldc_orig  = sort(orig,  rev = true)
ldc_recon = sort(recon, rev = true)

@assert length(ldc_orig) == length(ldc_recon)

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

# Zero line
plot!(pΔ, fill(0.0, length(delta)), linewidth = 1, linestyle = :dash)

# ---- Save canonical copy ----
png_main = joinpath(OUTDIR_MAIN, "LDC_8yr_System_Diff_RP_$rp_per_year.png")
savefig(pΔ, png_main)
println("Saved canonical difference plot: $png_main")

# ---- Save sweep copy ----
png_sweep = joinpath(sweep_folder, "LDC_8yr_System_Diff_RP_$rp_per_year.png")
savefig(pΔ, png_sweep)
println("Saved sweep difference plot: $png_sweep")

println("\nMin Δ = $(minimum(delta))")
println("Max Δ = $(maximum(delta))")

println("\nDone.")
