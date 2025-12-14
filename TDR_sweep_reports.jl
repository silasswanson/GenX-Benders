# TDR_sweep_reports.jl
#
# Sweep reporter for GenX-Brazil TDR_sweep_results_*.
# - Reconstructs system load (sum of Load_MW_z*) from RP_k/TDR_Results and computes RMSE + LDC plots
# - Aggregates RMSE_timeseries.csv written by run_timedomainreduction:
#     * local sweep: RMS across scenarios per variable
#     * global sweep: uses the single RMSE_timeseries.csv in TDR_Results
# - Writes:
#     RMSE_SystemLoad_summary.csv
#     RMSE_timeseries_summary.csv
#     RMSE_timeseries_aggregated_RP_k.csv per RP folder

using CSV
using DataFrames
using Statistics
using Plots
using Printf
using Measures

# ------------------ SETTINGS ------------------ #

const CASE   = joinpath(@__DIR__, "GenX-Brazil")
const NYEARS = 8
const PERIODS_PER_YEAR = 52

const ORIG_LOAD_PATH = joinpath(CASE, "Load_data.csv")

# IMPORTANT: point this at local or global sweep root
const SWEEP_ROOT = joinpath(CASE, "TDR_sweep_results_local")
 #const SWEEP_ROOT = joinpath(CASE, "TDR_sweep_results_global")

println("Case directory:   $CASE")
println("Original load:    $ORIG_LOAD_PATH")
println("Sweep root:       $SWEEP_ROOT")

isdir(SWEEP_ROOT) || error("No sweep folder at $SWEEP_ROOT")

# ------------------ HELPERS ------------------ #

"""
Return the actual column key that exists in df for a desired column name.
Works whether the df uses String names or Symbol names.
"""
function colkey(df::DataFrame, name::AbstractString)
    if name in names(df)
        return name
    elseif Symbol(name) in names(df)
        return Symbol(name)
    else
        error("Missing column '$name'. Columns found: $(names(df))")
    end
end

"""
RMS aggregator: sqrt(mean(x.^2)); NaN for empty.
"""
rms(x::AbstractVector{<:Real}) = isempty(x) ? NaN : sqrt(mean((Float64.(x)).^2))

"""
Compute category summaries from per-variable RMSE dict.
"""
function summarize_categories(agg_rmse::Dict{String,Float64})
    is_load(v)  = startswith(v, "Load_MW_")
    is_solar(v) = occursin("_solar_pv_", v)
    is_wind(v)  = occursin("_on_wind_", v) || occursin("_off_wind_", v)
    is_hydro(v) = occursin("_hydro_", v)

    total_vals = collect(values(agg_rmse))
    rm_total = rms(total_vals)

    rm_load  = rms([r for (v,r) in agg_rmse if is_load(v)])
    rm_wind  = rms([r for (v,r) in agg_rmse if is_wind(v)])
    rm_solar = rms([r for (v,r) in agg_rmse if is_solar(v)])
    rm_hydro = rms([r for (v,r) in agg_rmse if is_hydro(v)])

    return rm_total, rm_load, rm_wind, rm_solar, rm_hydro
end

# ------------------ READ ORIGINAL LOAD ------------------ #

ld = CSV.read(ORIG_LOAD_PATH, DataFrame)

zone_cols = filter(c -> occursin("Load_MW_z", String(c)), names(ld))
@assert !isempty(zone_cols) "No Load_MW_z* columns found in original Load_data.csv"

println("Detected zone load columns (original): ", String.(zone_cols))

orig = reduce(+, (ld[!, c] for c in zone_cols))

T_total = length(orig)
periods_total = NYEARS * PERIODS_PER_YEAR

@assert T_total % periods_total == 0 "Total hours ($T_total) not divisible by periods_total ($periods_total)"
H_per_period = T_total ÷ periods_total

println("Total hours       = $T_total")
println("Periods total     = $periods_total")
println("Hours per period  = $H_per_period")

hours_per_year = T_total ÷ NYEARS
@assert hours_per_year * NYEARS == T_total

# ------------------ LIST RP FOLDERS ------------------ #

rp_dirs = sort(filter(name -> startswith(name, "RP_"), readdir(SWEEP_ROOT)))
isempty(rp_dirs) && error("No RP_* folders found under $SWEEP_ROOT")
println("\nFound RP folders: ", rp_dirs)

summary_rows = DataFrame(RP_per_year = Int[], Global_RMSE = Float64[])

timeseries_summary = DataFrame(
    RP_per_year = Int[],
    Method      = String[],
    RMSE_Total  = Float64[],
    RMSE_Load   = Float64[],
    RMSE_Wind   = Float64[],
    RMSE_Solar  = Float64[],
    RMSE_Hydro  = Float64[],
)

# store LDC reconstructions for combined plots
ldc_orig_ref = nothing
ldc_recon_by_rp = Dict{Int, Vector{Float64}}()



# ------------------ LOOP RPs ------------------ #

for rp_dir_name in rp_dirs
    println("\n====================================================")
    println(" Processing sweep folder: $rp_dir_name")
    println("====================================================")

    rp_dir = joinpath(SWEEP_ROOT, rp_dir_name)
    isdir(rp_dir) || error("Missing RP folder $rp_dir")

    # local sweep detection
    scenario_dirs = filter(n -> startswith(n, "Scenario_") && isdir(joinpath(rp_dir, n)), readdir(rp_dir))
    is_local = !isempty(scenario_dirs)
    method = is_local ? "local" : "global"

    results_root = joinpath(rp_dir, "TDR_Results")
    isdir(results_root) || error("Missing TDR_Results in $rp_dir")

    # ---- Period_map ----
    tdr_map_path = joinpath(results_root, "Period_map.csv")
    isfile(tdr_map_path) || error("Missing Period_map.csv at $tdr_map_path")

    pm = CSV.read(tdr_map_path, DataFrame)

    k_period = colkey(pm, "Period_Index")
    k_repidx = colkey(pm, "Rep_Period_Index")

    period_idx = Int.(pm[!, k_period])
    rep_idx    = Int.(pm[!, k_repidx])

    @assert length(period_idx) == periods_total "Period_map row count mismatch in $rp_dir_name"
    @assert maximum(period_idx) == periods_total "Unexpected max Period_Index in $rp_dir_name"

    n_rep_total = maximum(rep_idx)
    if n_rep_total % NYEARS != 0
        error("Rep_Period_Index max=$n_rep_total not divisible by NYEARS=$NYEARS in $rp_dir_name")
    end
    rp_per_year = n_rep_total ÷ NYEARS

    println("  Detected representative periods per year = $rp_per_year (total rep = $n_rep_total)")

    # ---- Load_data ----
    tdr_load_path = joinpath(results_root, "Load_data.csv")
    isfile(tdr_load_path) || error("Missing reduced Load_data.csv at $tdr_load_path")

    tdr_ld = CSV.read(tdr_load_path, DataFrame)

    zone_cols_red = filter(c -> occursin("Load_MW_z", String(c)), names(tdr_ld))
    @assert !isempty(zone_cols_red) "No Load_MW_z* columns found in reduced Load_data.csv for $rp_dir_name"
    @assert nrow(tdr_ld) == n_rep_total * H_per_period "Reduced Load_data row count mismatch in $rp_dir_name"

    reduced_total = reduce(+, (tdr_ld[!, c] for c in zone_cols_red))
    get_profile(k::Int) = @view reduced_total[(k-1)*H_per_period+1 : k*H_per_period]

    # reconstruct full series
    recon = zeros(Float64, T_total)
    for r in eachindex(period_idx)
        p  = period_idx[r]
        rp = rep_idx[r]
        start_hour = (p - 1) * H_per_period + 1
        end_hour   = p * H_per_period
        recon[start_hour:end_hour] .= get_profile(rp)
    end


    println("Peak orig MW = ", maximum(orig), "   Peak recon MW = ", maximum(recon))
    println("Min  orig MW = ", minimum(orig), "   Min  recon MW = ", minimum(recon))

    # system-load RMSE
    err = recon .- orig
    rmse_global = sqrt(mean(err .^ 2))
    println("  Global RMSE over 8 years (system load recon) = $rmse_global")
    push!(summary_rows, (RP_per_year = rp_per_year, Global_RMSE = rmse_global))

    # LDC plots
    ldc_orig  = sort(orig,  rev = true)
    ldc_recon = sort(recon, rev = true)

    if ldc_orig_ref === nothing
        global ldc_orig_ref = ldc_orig
    end
    ldc_recon_by_rp[rp_per_year] = ldc_recon


    p_global = plot(ldc_orig,
        label="Original", xlabel="Hour (sorted)", ylabel="Load (MW)",
        title="8-Year Load Duration Curve (RP = $rp_per_year)",
        linewidth=2, size=(900,500)
    )
    plot!(p_global, ldc_recon, label="Reconstructed", linewidth=2, linestyle=:dash)
    savefig(p_global, joinpath(rp_dir, "LDC_8yr_System_RP_$rp_per_year.png"))

    delta = ldc_recon .- ldc_orig
    pΔ = plot(delta,
        xlabel="Hour (sorted)", ylabel="Δ Load (MW) = recon - original",
        title="Difference in 8-Year LDC (RP = $rp_per_year)",
        legend=false, linewidth=1.8, size=(900,500)
    )
    plot!(pΔ, fill(0.0, length(delta)), linewidth=1, linestyle=:dash)
    savefig(pΔ, joinpath(rp_dir, "LDC_8yr_System_Diff_RP_$rp_per_year.png"))

    # ---- RMSE_timeseries aggregation ----
    rmse_per_var = Dict{String, Vector{Float64}}()

    if is_local
        sort!(scenario_dirs, by = s -> parse(Int, split(s, "_")[2]))

        for scen in scenario_dirs
            rmse_path = joinpath(rp_dir, scen, "TDR_Results", "RMSE_timeseries.csv")
            isfile(rmse_path) || error("Missing RMSE_timeseries.csv at $rmse_path")

            df_rm = CSV.read(rmse_path, DataFrame)

            k_var = colkey(df_rm, "Variable")
            k_rm  = colkey(df_rm, "RMSE")

            for row in eachrow(df_rm)
                v = String(row[k_var])
                r = Float64(row[k_rm])
                push!(get!(rmse_per_var, v, Float64[]), r)
            end
        end

        for (v, vals) in rmse_per_var
            if length(vals) != length(scenario_dirs)
                error("Inconsistent RMSE_timeseries across scenarios for '$v' in $rp_dir_name: got $(length(vals)) expected $(length(scenario_dirs))")
            end
        end
    else
        rmse_path = joinpath(results_root, "RMSE_timeseries.csv")
        isfile(rmse_path) || error("Missing RMSE_timeseries.csv at $rmse_path")

        df_rm = CSV.read(rmse_path, DataFrame)

        k_var = colkey(df_rm, "Variable")
        k_rm  = colkey(df_rm, "RMSE")

        for row in eachrow(df_rm)
            v = String(row[k_var])
            r = Float64(row[k_rm])
            push!(get!(rmse_per_var, v, Float64[]), r)
        end
    end

    agg_rmse = Dict{String,Float64}()
    for (v, vals) in rmse_per_var
        agg_rmse[v] = rms(vals)
    end

    rm_total, rm_load, rm_wind, rm_solar, rm_hydro = summarize_categories(agg_rmse)

    push!(timeseries_summary, (
        RP_per_year = rp_per_year,
        Method      = method,
        RMSE_Total  = rm_total,
        RMSE_Load   = rm_load,
        RMSE_Wind   = rm_wind,
        RMSE_Solar  = rm_solar,
        RMSE_Hydro  = rm_hydro,
    ))

    rmse_var_df = DataFrame(Variable = collect(keys(agg_rmse)),
                            RMSE     = collect(values(agg_rmse)))
    sort!(rmse_var_df, :Variable)
    CSV.write(joinpath(rp_dir, "RMSE_timeseries_aggregated_RP_$(rp_per_year).csv"), rmse_var_df)
end



# ------------------ COMBINED LDC PLOTS ------------------ #
using Colors

rp_colors = Dict(
    1  => colorant"#d62728",  # red
    5  => colorant"#ff7f0e",  # orange
    10 => colorant"#2ca02c",  # green
    15 => colorant"#1f77b4",  # blue
    20 => colorant"#9467bd",  # purple
    30 => colorant"#8c564b",  # brown
    40 => colorant"#17becf",  # teal
    50 => colorant"#7f7f7f",  # gray
)



using Printf
tight = 6mm


N = length(ldc_orig_ref)
x = collect(1:N)

# nice tick locations every 10k hours
xticks_vals = 0:10000:N
xticks_labels = string.(xticks_vals)


# Ensure RPs plot in increasing order
RP_PLOT = Set([1, 5, 10, 15, 20, 30, 40, 50])  # no 52
rps = sort([rp for rp in keys(ldc_recon_by_rp) if rp in RP_PLOT])

xt = 0:10000:70000
xformatter = x -> @sprintf("%d", Int(x))

sweep_label = occursin("global", SWEEP_ROOT) ? "TDR sweep global" : "TDR sweep local"

p_ldc = plot(
    x, ldc_orig_ref;
    label="Original",
    xlabel="Hour (sorted)",
    ylabel="Load (MW)",
    title="8-Year Load Duration Curve (all RPs) — $(sweep_label)",
    linewidth=3,
    size=(1000,600),
    legend=:topright,
    xticks=(xticks_vals, xticks_labels),
    left_margin=tight,
    right_margin=4mm,
    top_margin=4mm,
    bottom_margin=6mm,
)

for rp in rps
    plot!(
        p_ldc,
        x,
        ldc_recon_by_rp[rp];
        label="Reconstructed RP=$rp",
        linewidth=2,
        color = rp_colors[rp],
    )
end

savefig(p_ldc, joinpath(SWEEP_ROOT, "LDC_8yr_System_AllRPs.png"))


# Combined ΔLDC: recon - original for each RP
p_dldc = plot(
    x,
    zeros(N);
    label=false,
    xlabel="Hour (sorted)",
    ylabel="Δ Load (MW) = recon − original",
    title="Δ 8-Year LDC (all RPs) — $(sweep_label)",
    size=(1000,600),
    legend=:topright,
    xticks=(xticks_vals, xticks_labels),
    left_margin=tight,
    right_margin=4mm,
    top_margin=4mm,
    bottom_margin=6mm,
)

for rp in rps
    delta = ldc_recon_by_rp[rp] .- ldc_orig_ref
    plot!(
        p_dldc,
        x,
        delta;
        label="RP=$rp",
        linewidth=2,
        color = rp_colors[rp],
    )
end

savefig(p_dldc, joinpath(SWEEP_ROOT, "LDC_8yr_System_Diff_AllRPs.png"))





# ------------------ WRITE SUMMARIES ------------------ #

CSV.write(joinpath(SWEEP_ROOT, "RMSE_SystemLoad_summary.csv"), summary_rows)
println("\nSweep summary written to: ", joinpath(SWEEP_ROOT, "RMSE_SystemLoad_summary.csv"))

CSV.write(joinpath(SWEEP_ROOT, "RMSE_timeseries_summary.csv"), timeseries_summary)
println("Timeseries RMSE summary written to: ", joinpath(SWEEP_ROOT, "RMSE_timeseries_summary.csv"))

println("\nAll sweep reports complete.")


all_delta = vcat([ldc_recon_by_rp[rp] .- ldc_orig_ref for rp in rps]...)
lo = -1000 #quantile(all_delta, 0.01)
hi = 1000 #quantile(all_delta, 0.99)

p_dldc_zoom = plot(p_dldc; ylims=(lo, hi), title="Δ 8-Year LDC (zoomed) — $(sweep_label)",
    left_margin=tight, right_margin=4mm, top_margin=4mm, bottom_margin=6mm,
)

savefig(p_dldc_zoom, joinpath(SWEEP_ROOT, "LDC_8yr_System_Diff_AllRPs_Zoom.png"))

