# TDR_compare_reconstruction.jl
#
# Apples-to-apples reconstruction RMSE:
#   - Reconstruct full 8-year hourly time series for ALL numeric columns
#   - Compare reconstructed vs original
#   - Works for both local and global sweeps
#
# Run:
#   julia --project=. TDR_compare_reconstruction.jl

using CSV, DataFrames, Statistics, Printf

# ------------------ USER SETTINGS ------------------ #
const CASE = joinpath(@__DIR__, "GenX-Brazil")

const SWEEP_LOCAL  = joinpath(CASE, "TDR_sweep_results_local")
const SWEEP_GLOBAL = joinpath(CASE, "TDR_sweep_results_global")

# If you only want RP_20, set: const RP_LIST = ["RP_20"]
const RP_LIST = nothing  # nothing = auto-detect all RP_* in both sweeps

# Expected structure
const NYEARS = 8
const PERIODS_PER_YEAR = 52
const PERIODS_TOTAL = NYEARS * PERIODS_PER_YEAR
const H_PER_PERIOD = 168  # weekly

# Original input paths
const ORIG_LOAD_PATH = joinpath(CASE, "Load_data.csv")
const ORIG_GV_PATH   = joinpath(CASE, "Generators_variability.csv")
const ORIG_FUEL_PATH = joinpath(CASE, "Fuels_data.csv")

# --------------------------------------------------- #

# ---------- Helpers: column name normalization ----------

"""
Return names(df) as Vector{String}.
Works regardless of whether DataFrames is storing names as Symbol.
"""
names_str(df::DataFrame) = String.(names(df))

"""
Return true if a column exists, case-sensitively, by string.
"""
hascol(df::DataFrame, col::AbstractString) = col in names_str(df)

"""
Fetch column by string name, regardless of Symbol/String internal.
"""
function col(df::DataFrame, colname::AbstractString)
    return df[!, Symbol(colname)]
end

# ---------- Helpers: robust CSV reads ----------

"""
Read a GenX-style CSV.
- Keeps missings (we want to detect them, not silently drop them).
- Does not attempt Symbol/String coercion games; we just access by Symbol(colname).
"""
read_csv(path::AbstractString) = CSV.read(path, DataFrame; silencewarnings=true)

"""
Fuels_data.csv in some cases has an extra first row with Time_Index==0 (or blank),
and then 69888 rows of hourly data. This function:
- ensures we return exactly T rows of hourly data with Time_Index 1..T
- and returns the "header row" separately if present (not used in RMSE)
"""
function read_fuels_hourly(path::AbstractString, T::Int)
    df = read_csv(path)

    # Try to find Time_Index column (some files use Time_Index, some time_index)
    tcol = nothing
    for c in names_str(df)
        if c == "Time_Index"
            tcol = c
            break
        end
    end
    tcol === nothing && error("Fuels file $path has no Time_Index column. Found: $(names_str(df))")

    ti = col(df, tcol)

    # Detect an extra first row with Time_Index == 0 OR missing
    has_header = false
    if length(ti) >= 1
        if (ti[1] isa Missing) || (try ti[1] == 0 catch; false end)
            has_header = true
        end
    end

    if has_header
        header_row = df[1:1, :]
        hourly = df[2:end, :]
    else
        header_row = df[0: -1, :]  # empty
        hourly = df
    end

    n = nrow(hourly)
    if n != T
        error("Hourly fuels rows mismatch: expected $T, got $n at $path (header_row=$(has_header))")
    end

    # Optional sanity: Time_Index should be 1..T (but don’t hard fail if it’s absent/float)
    # We just return hourly as-is, RMSE uses row alignment.
    return header_row, hourly
end

# ---------- Reconstruction from Period_map ----------

"""
Reconstruct a full hourly series (length T_total) from:
- reduced_df: rows = K_total * H_PER_PERIOD
- period_map: rows = PERIODS_TOTAL
using Rep_Period_Index mapping.

Returns recon_df with same columns as reduced_df (numeric columns reconstructed).
"""
function reconstruct_from_period_map(reduced_df::DataFrame, period_map::DataFrame, T_total::Int)
    # Validate Period_map columns exist (by string)
    for required in ["Period_Index", "Rep_Period_Index"]
        hascol(period_map, required) || error("Period_map.csv missing $required. Found: $(names_str(period_map))")
    end

    period_idx = Int.(col(period_map, "Period_Index"))
    rep_idx    = Int.(col(period_map, "Rep_Period_Index"))

    maximum(period_idx) == PERIODS_TOTAL || error("Period_map Period_Index max != $PERIODS_TOTAL")
    length(period_idx) == PERIODS_TOTAL  || error("Period_map rows != $PERIODS_TOTAL")

    # K_total inferred from reduced_df rows
    if nrow(reduced_df) % H_PER_PERIOD != 0
        error("Reduced rows $(nrow(reduced_df)) not divisible by $H_PER_PERIOD")
    end
    K_total = nrow(reduced_df) ÷ H_PER_PERIOD

    maximum(rep_idx) <= K_total || error("Rep_Period_Index max $(maximum(rep_idx)) > K_total $K_total")

    # Identify time-series numeric columns to reconstruct.
    # Exclude metadata fields that are not meaningful per-hour and often contain missings after TDR.
    META_COLS = Set([
        "Voll",
        "Demand_Segment",
        "Cost_of_Demand_Curtailment_per_MW",
        "Max_Demand_Curtailment",
        "\$/MWh",
        "Rep_Periods",
        "Timesteps_per_Rep_Period",
        "Sub_Weights",
        "Time_Index",
    ])

    num_cols = String[]
    for c in names_str(reduced_df)
        if c in META_COLS
            continue
        end
        v = reduced_df[!, Symbol(c)]
        T = eltype(v)
        if T <: Number || (T <: Union{Missing, Number})
            # For time-series columns, we require no missings anywhere.
            # If it has missings, it is either metadata-like or genuinely broken; we skip metadata-like by name.
            if any(ismissing, v)
                # keep this as a hard error for *non-metadata* numeric columns
                error("Missing detected in reduced time-series column '$c' (not in META_COLS). This indicates a real issue.")
            end
            push!(num_cols, c)
        end
    end

    # Allocate output
    recon = DataFrame()
    for c in num_cols
        recon[!, Symbol(c)] = Vector{Float64}(undef, T_total)
    end

    # Helper to view k-th rep period block
    @inline function block_range(k::Int)
        a = (k - 1) * H_PER_PERIOD + 1
        b = k * H_PER_PERIOD
        return a:b
    end

    # Fill
    for r in 1:PERIODS_TOTAL
        p  = period_idx[r]
        k  = rep_idx[r]

        hour_a = (p - 1) * H_PER_PERIOD + 1
        hour_b = p * H_PER_PERIOD

        br = block_range(k)

        for c in num_cols
            src = reduced_df[br, Symbol(c)]
            # Detect missings in numeric time-series columns (don’t “fix” silently)
            recon[hour_a:hour_b, Symbol(c)] .= Float64.(src)
        end
    end

    return recon
end

# ---------- RMSE ----------

rmse(x::AbstractVector{<:Real}, y::AbstractVector{<:Real}) = sqrt(mean((x .- y) .^ 2))

"""
Compute RMSE per column for columns common to orig_df and recon_df.
Returns DataFrame with columns: Variable, RMSE
"""
function rmse_table(orig_df::DataFrame, recon_df::DataFrame)
    common = intersect(names_str(orig_df), names_str(recon_df))

    # Only numeric columns
    vars = String[]
    vals = Float64[]
    for c in common
        xo = orig_df[!, Symbol(c)]
        xr = recon_df[!, Symbol(c)]
        if (eltype(xo) <: Number || eltype(xo) <: Union{Missing, Number}) &&
           (eltype(xr) <: Number || eltype(xr) <: Union{Missing, Number})
            if any(ismissing, xo)
                error("Missing detected in ORIGINAL column '$c' — investigate your original inputs.")
            end
            if any(ismissing, xr)
                error("Missing detected in RECON column '$c' — reconstruction produced missing (should not happen).")
            end
            push!(vars, c)
            push!(vals, rmse(Float64.(xo), Float64.(xr)))
        end
    end
    return DataFrame(Variable = vars, RMSE = vals)
end

# ---------- Main ----------

println("\n================ ORIGINAL INPUT SANITY ================")
@assert isfile(ORIG_LOAD_PATH) "Missing $ORIG_LOAD_PATH"
@assert isfile(ORIG_GV_PATH)   "Missing $ORIG_GV_PATH"
@assert isfile(ORIG_FUEL_PATH) "Missing $ORIG_FUEL_PATH"

orig_load = read_csv(ORIG_LOAD_PATH)
orig_gv   = read_csv(ORIG_GV_PATH)

T_expected = NYEARS * 8736
println("Expected rows per original = $T_expected")
println("nrow(orig_load) = ", nrow(orig_load))
println("nrow(orig_gv)   = ", nrow(orig_gv))
@assert nrow(orig_load) == T_expected
@assert nrow(orig_gv)   == T_expected

_, orig_fuel_hourly = read_fuels_hourly(ORIG_FUEL_PATH, T_expected)
println("nrow(orig_fuel) = ", nrow(orig_fuel_hourly))
println("========================================================\n")

function list_rp_dirs(root::String)
    isdir(root) || return String[]
    return sort(filter(d -> startswith(d, "RP_") && isdir(joinpath(root, d)), readdir(root)))
end

rp_local  = list_rp_dirs(SWEEP_LOCAL)
rp_global = list_rp_dirs(SWEEP_GLOBAL)

rp_all = if RP_LIST === nothing
    sort(union(rp_local, rp_global))
else
    RP_LIST
end

isempty(rp_all) && error("No RP_* folders found. Check SWEEP_LOCAL and SWEEP_GLOBAL paths.")

println("RP folders to compare: ", rp_all)

# Output summary
out_summary = DataFrame(Method=String[], RP=String[], Dataset=String[], Variable=String[], RMSE=Float64[])

function process_one(method_name::String, sweep_root::String, rp::String)
    rpdir = joinpath(sweep_root, rp, "TDR_Results")
    isdir(rpdir) || return nothing

    pm_path = joinpath(rpdir, "Period_map.csv")
    ld_path = joinpath(rpdir, "Load_data.csv")
    gv_path = joinpath(rpdir, "Generators_variability.csv")
    fu_path = joinpath(rpdir, "Fuels_data.csv")

    isfile(pm_path) || error("Missing $pm_path")
    isfile(ld_path) || error("Missing $ld_path")
    isfile(gv_path) || error("Missing $gv_path")
    isfile(fu_path) || error("Missing $fu_path")

    period_map = read_csv(pm_path)
    red_load   = read_csv(ld_path)
    red_gv     = read_csv(gv_path)
    _, red_fuel_hourly = read_fuels_hourly(fu_path, nrow(red_load))  # reduced fuels should match reduced load rows

    # Reconstruct
    recon_load = reconstruct_from_period_map(red_load, period_map, T_expected)
    recon_gv   = reconstruct_from_period_map(red_gv,   period_map, T_expected)
    recon_fuel = reconstruct_from_period_map(red_fuel_hourly, period_map, T_expected)

    # RMSE tables
    tab_load = rmse_table(orig_load, recon_load)
    tab_gv   = rmse_table(orig_gv,   recon_gv)
    tab_fuel = rmse_table(orig_fuel_hourly, recon_fuel)

    # Write per-RP files
    outdir = joinpath(sweep_root, rp)
    CSV.write(joinpath(outdir, "RMSE_reconstruction_LOAD.csv"), tab_load)
    CSV.write(joinpath(outdir, "RMSE_reconstruction_GV.csv"),   tab_gv)
    CSV.write(joinpath(outdir, "RMSE_reconstruction_FUELS.csv"), tab_fuel)

    # Add to summary (long format)
    for row in eachrow(tab_load)
        push!(out_summary, (method_name, rp, "LOAD", row.Variable, row.RMSE))
    end
    for row in eachrow(tab_gv)
        push!(out_summary, (method_name, rp, "GV", row.Variable, row.RMSE))
    end
    for row in eachrow(tab_fuel)
        push!(out_summary, (method_name, rp, "FUELS", row.Variable, row.RMSE))
    end

    return nothing
end

for rp in rp_all
    println("====================================================")
    println("RP: $rp")
    println("====================================================")
    if rp in rp_local
        println("  -> processing LOCAL")
        process_one("LOCAL", SWEEP_LOCAL, rp)
    else
        println("  -> LOCAL missing $rp (skipping)")
    end
    if rp in rp_global
        println("  -> processing GLOBAL")
        process_one("GLOBAL", SWEEP_GLOBAL, rp)
    else
        println("  -> GLOBAL missing $rp (skipping)")
    end
end

# Save combined summary
summary_path = joinpath(CASE, "TDR_compare_reconstruction_summary_long.csv")
CSV.write(summary_path, out_summary)
println("\nWrote: $summary_path")

println("\nDone.")
