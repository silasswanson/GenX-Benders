# TDR_sweep_edc_plots.jl
#
# Error Duration Curve (EDC) for system load:
#   err(t) = |recon(t) - orig(t)|, then sort(err, rev=true)
#
# Compares LOCAL vs GLOBAL for each RP in RP_PLOT.
#
# IMPORTANT: This script matches your Period_map.csv structure:
#   Period_Index, Rep_Period, Rep_Period_Index
# and uses Rep_Period_Index (1..K) to select the representative-week blocks.

using CSV, DataFrames
using Statistics
using Plots
using Measures

# ------------------ SETTINGS ------------------ #
const CASE = joinpath(pwd(), "GenX-Brazil")
const ORIG_LOAD = joinpath(CASE, "Load_data.csv")

const SWEEP_GLOBAL = joinpath(CASE, "TDR_sweep_results_global")
const SWEEP_LOCAL  = joinpath(CASE, "TDR_sweep_results_local")

const RP_PLOT = [1, 5, 10, 15, 20, 30, 40, 50]   # exclude 52
const HOURS_PER_PERIOD = 168
const NUM_SCENARIOS = 8

const OUT_ALL  = joinpath(CASE, "EDC_SystemLoad_Local_vs_Global_AllRPs.png")
const OUT_RP20 = joinpath(CASE, "EDC_SystemLoad_Local_vs_Global_RP20.png")

const tight = 8mm

# ------------------ HELPERS ------------------ #

# Detect columns like Load_MW_z1, Load_MW_z2, ...
function detect_zone_load_cols(df::DataFrame)::Vector{String}
    cols = String[]
    for nm in names(df)
        s = String(nm)
        if occursin("Load_MW_z", s)
            push!(cols, s)
        end
    end
    isempty(cols) && error("No Load_MW_z* columns found in Load_data.csv")
    return cols
end

# Read system load (sum of zones) from a Load_data.csv path
function read_system_load(load_path::AbstractString)::Vector{Float64}
    df = CSV.read(load_path, DataFrame)
    loadcols = detect_zone_load_cols(df)
    sys = vec(sum(Matrix(df[:, loadcols]), dims=2))
    return sys
end

# Reconstruct full-hourly system load using:
#   - Period_map.csv (Period_Index, Rep_Period, Rep_Period_Index)
#   - representative Load_data.csv in TDR_Results (K rep periods * 168 hours)
#
# Uses Rep_Period_Index (1..K) to index blocks in representative series.
function reconstruct_system_load_from_rp_root(rp_root::AbstractString)::Vector{Float64}
    tdr = joinpath(rp_root, "TDR_Results")
    map_path = joinpath(tdr, "Period_map.csv")
    rep_load_path = joinpath(tdr, "Load_data.csv")

    isfile(map_path) || error("Missing Period_map.csv at: $map_path")
    isfile(rep_load_path) || error("Missing Load_data.csv at: $rep_load_path")

    pm = CSV.read(map_path, DataFrame)
    rep_sys = read_system_load(rep_load_path)

    ("Rep_Period_Index" in names(pm)) || error("Period_map.csv missing Rep_Period_Index column. Found: $(names(pm))")
    rep_idx = Vector{Int}(pm[!, "Rep_Period_Index"])  # 1..K block index

    # sanity: rep_sys must have exactly K blocks of 168
    K_map = maximum(rep_idx)
    K_rep = length(rep_sys) ÷ HOURS_PER_PERIOD
    (length(rep_sys) % HOURS_PER_PERIOD == 0) || error("Rep Load_data length not divisible by 168: $(length(rep_sys))")
    (K_map == K_rep) || error("Mismatch: Period_map implies K=$K_map reps, but rep Load_data has K=$K_rep blocks")

    T = length(rep_idx) * HOURS_PER_PERIOD
    recon = Vector{Float64}(undef, T)

    for (p, rid) in enumerate(rep_idx)
        start_rep = (rid - 1) * HOURS_PER_PERIOD + 1
        stop_rep  = start_rep + HOURS_PER_PERIOD - 1

        start_out = (p - 1) * HOURS_PER_PERIOD + 1
        stop_out  = start_out + HOURS_PER_PERIOD - 1

        recon[start_out:stop_out] .= rep_sys[start_rep:stop_rep]
    end

    return recon
end

# Local sweep: concatenate Scenario_1..Scenario_8 reconstructions
function reconstruct_system_load_local(rp_root::AbstractString)::Vector{Float64}
    chunks = Vector{Vector{Float64}}()
    for s in 1:NUM_SCENARIOS
        scen_root = joinpath(rp_root, "Scenario_$s")
        isdir(joinpath(scen_root, "TDR_Results")) || error("Missing TDR_Results for local scenario at: $(joinpath(scen_root, "TDR_Results"))")
        push!(chunks, reconstruct_system_load_from_rp_root(scen_root))
    end
    return vcat(chunks...)
end

# Error Duration Curve: sort |error| descending
function edc(orig::Vector{Float64}, recon::Vector{Float64})::Vector{Float64}
    length(orig) == length(recon) || error("Length mismatch: orig=$(length(orig)) recon=$(length(recon))")
    return sort(abs.(recon .- orig), rev=true)
end

# Nice integer x-axis ticks (avoid scientific notation)
function edc_xticks(N::Int)
    vals = 0:10000:N
    labels = string.(vals)
    return (vals, labels)
end

# ------------------ MAIN ------------------ #

println("Case directory: ", CASE)
println("Original load:  ", ORIG_LOAD)
orig_sys = read_system_load(ORIG_LOAD)
N = length(orig_sys)
println("Original hours: ", N)

xt = edc_xticks(N)
x = collect(1:N)

plots = []

for rp in RP_PLOT
    rp_global = joinpath(SWEEP_GLOBAL, "RP_$(rp)")
    rp_local  = joinpath(SWEEP_LOCAL,  "RP_$(rp)")

    println("\nRP = $rp")
    println("  global root: ", rp_global)
    println("  local  root: ", rp_local)

    recon_g = reconstruct_system_load_from_rp_root(rp_global)
    recon_l = reconstruct_system_load_local(rp_local)

    # should match original length
    length(recon_g) == N || error("Global recon length $(length(recon_g)) != original length $N at RP=$rp")
    length(recon_l) == N || error("Local recon length $(length(recon_l)) != original length $N at RP=$rp")

    edc_g = edc(orig_sys, recon_g)
    edc_l = edc(orig_sys, recon_l)

    p = plot(
        x, edc_l;
        label="Local",
        xlabel="Hour rank (by |error|)",
        ylabel="|Error| (MW)",
        title="EDC System Load (RP=$rp)",
        linewidth=2,
        size=(900,600),
        legend=:topright,
        xticks=xt,
        left_margin=tight, right_margin=4mm, top_margin=4mm, bottom_margin=8mm,
    )
    plot!(p, x, edc_g; label="Global", linewidth=2)

    push!(plots, p)
end

# 2 rows x 4 columns
p_all = plot(
    plots...;
    layout=(2, 4),
    size=(1800, 850),
    left_margin=tight, right_margin=4mm, top_margin=6mm, bottom_margin=8mm,
)
savefig(p_all, OUT_ALL)
println("\nSaved: ", OUT_ALL)

# single RP=20 plot
if 20 in RP_PLOT
    rp = 20
    recon_g = reconstruct_system_load_from_rp_root(joinpath(SWEEP_GLOBAL, "RP_$(rp)"))
    recon_l = reconstruct_system_load_local(joinpath(SWEEP_LOCAL,  "RP_$(rp)"))
    edc_g = edc(orig_sys, recon_g)
    edc_l = edc(orig_sys, recon_l)

    p20 = plot(
        x, edc_l;
        label="Local",
        xlabel="Hour rank (by |error|)",
        ylabel="|Error| (MW)",
        title="EDC System Load (RP=20)",
        linewidth=3,
        size=(1100,700),
        legend=:topright,
        xticks=xt,
        left_margin=tight, right_margin=4mm, top_margin=6mm, bottom_margin=10mm,
    )
    plot!(p20, x, edc_g; label="Global", linewidth=3)

    savefig(p20, OUT_RP20)
    println("Saved: ", OUT_RP20)
end

println("\nDone.")
