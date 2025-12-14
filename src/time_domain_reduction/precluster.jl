@doc raw"""Run the GenX time domain reduction on the given case folder

case - folder for the case
stage_id - possibly something to do with MultiStage
verbose - print extra outputs

This function overwrites the time-domain-reduced inputs if they already exist.

"""
function timeseries2scenarios(T::Int64, numscenarios::Int64)
	# Breaks up the timeseries into numscenarios scenarios
	# If scenarios are uneven, make the first scenario longest
	scenario_times = Array{Int64}(undef, numscenarios, 2)
	min_scenario_length = Int64(floor(T/numscenarios))
	remainder = T % numscenarios
	for s in 1:numscenarios
		scenario_times[s,1] = (s-1)*min_scenario_length + 1 + remainder
		scenario_times[s,2] = s*min_scenario_length + remainder
	end
	scenario_times[1,1] = 1
	return scenario_times
end

function list_csv_files(directory::String)
    csv_files = String[]
    for file in readdir(directory)
        if occursin(".csv", file)
            push!(csv_files,  file)
        end
    end
    return csv_files
end

function unpack_scenarios(case::AbstractString,setup::Dict)

	scenario_path = joinpath(case,"Scenarios");
	if isdir(scenario_path)
		# do nothing
		println("")
		println("Scenarios folder already exists. Delete it to re-run the unpacking.")
		return
	else
		scenario_path = mkdir(scenario_path)
	end

	df = load_dataframe(joinpath(case,"Load_data.csv"));
	df_fuels = load_dataframe(joinpath(case,"Fuels_data.csv"));
	df_genvar = load_dataframe(joinpath(case,"Generators_variability.csv"));
    if df.Time_Index[end] % setup["NumScenarios"]==0
    	ScenarioLength = Int(df.Time_Index[end]/setup["NumScenarios"]);
	else
        println("")
        error("Error! The total number of time steps in Load_data.csv has to be divisible by the number of scenarios")
    end
	for s = 1:setup["NumScenarios"]	
		scenario_path_s = joinpath(scenario_path,"Scenario_$s");
		mkdir(scenario_path_s)

        csv_files_in_case = list_csv_files(case);
        for f in csv_files_in_case
            cp(joinpath(case,f),joinpath(scenario_path_s,f))
        end

		df_s = df[(s-1)*ScenarioLength+1:s*ScenarioLength,:];
		df_s.Cost_of_Demand_Curtailment_per_MW[1] = df.Cost_of_Demand_Curtailment_per_MW[1];
		df_s.Max_Demand_Curtailment[1] = df.Max_Demand_Curtailment[1];
		df_s.Time_Index = collect(1:ScenarioLength);
		df_s.Demand_Segment[1] = df.Demand_Segment[1];
		df_s.Rep_Periods[1] = 1;
		df_s.Timesteps_per_Rep_Period[1] = ScenarioLength;
		df_s.Sub_Weights[1] = 1.0;
		df_s.Voll[1] = df.Voll[1];
		CSV.write(joinpath(scenario_path_s, "Load_data.csv"), df_s)

		df_genvar_s = df_genvar[(s-1)*ScenarioLength+1:s*ScenarioLength,:];
		df_genvar_s.Time_Index = collect(1:ScenarioLength);
		CSV.write(joinpath(scenario_path_s, "Generators_variability.csv"), df_genvar_s);

		df_fuels_s = df_fuels[1+(s-1)*ScenarioLength+1:1+s*ScenarioLength,:];
		df_fuels_s.Time_Index = collect(1:ScenarioLength);
		insert!(df_fuels_s,1,df_fuels[1,:])
		CSV.write(joinpath(scenario_path_s, "Fuels_data.csv"), df_fuels_s);


	end


end

function pack_scenarios(case,setup)

    TDRpath = joinpath(case, setup["TimeDomainReductionFolder"]);
    if isdir(TDRpath)
        #do nothing
    else
        mkdir(TDRpath)
    end
    scenario_path = joinpath(joinpath(case,"Scenarios"),"Scenario_1");
    TDRpath_1 = joinpath(scenario_path, setup["TimeDomainReductionFolder"]);
    df = load_dataframe(joinpath(TDRpath_1,"Load_data.csv"));
    df_fuels = load_dataframe(joinpath(TDRpath_1,"Fuels_data.csv"));
    df_genvar = load_dataframe(joinpath(TDRpath_1,"Generators_variability.csv"));
    df_periodmap = load_dataframe(joinpath(TDRpath_1,"Period_map.csv"));
    np = length(df_periodmap.Period_Index);
    for s in 2:setup["NumScenarios"]
        scenario_path = joinpath(joinpath(case,"Scenarios"),"Scenario_$s");
        TDRpath_s = joinpath(scenario_path, setup["TimeDomainReductionFolder"]);
        df_s = load_dataframe(joinpath(TDRpath_s,"Load_data.csv"));
        
        sub_weights_scenario = df_s.Sub_Weights[1:df_s.Rep_Periods[1]];
        num_subperiods_scenarios = df_s.Rep_Periods[1];
        df_s.Cost_of_Demand_Curtailment_per_MW[1]=missing;
		df_s.Max_Demand_Curtailment[1] = missing;
		df_s.Demand_Segment[1] = missing;
		df_s.Rep_Periods[1] = missing;
		df_s.Sub_Weights .= missing;
        df_s.Timesteps_per_Rep_Period[1] = missing;
		df_s.Voll[1] = missing;
        
        append!(df,df_s)
        df.Sub_Weights[df.Rep_Periods[1]+1:df.Rep_Periods[1]+num_subperiods_scenarios] = sub_weights_scenario
        df.Rep_Periods[1] += num_subperiods_scenarios;

	    df_fuels_s = load_dataframe(joinpath(TDRpath_s,"Fuels_data.csv"));
        append!(df_fuels,df_fuels_s[2:end,:])
	    df_genvar_s = load_dataframe(joinpath(TDRpath_s,"Generators_variability.csv"));
        append!(df_genvar,df_genvar_s)

        df_periodmap_s = load_dataframe(joinpath(TDRpath_s,"Period_map.csv"));
        df_periodmap_s.Period_Index = collect((s-1)*np+1:s*np);
        df_periodmap_s.Rep_Period = [(s-1)*np+k for k in df_periodmap_s.Rep_Period];
        append!(df_periodmap,df_periodmap_s)
    end
    df.Time_Index = collect(1:size(df,1));
    df_fuels.Time_Index[2:end] = collect(1:size(df,1));
    df_genvar.Time_Index = collect(1:size(df,1));
    rep_period_all = sort(unique(df_periodmap.Rep_Period));
    df_periodmap.Rep_Period_Index = [findfirst(rep_period_all.== k) for k in df_periodmap.Rep_Period];

    CSV.write(joinpath(TDRpath, "Load_data.csv"), df)
    CSV.write(joinpath(TDRpath, "Fuels_data.csv"),df_fuels)
    CSV.write(joinpath(TDRpath, "Generators_variability.csv"),df_genvar)
    CSV.write(joinpath(TDRpath, "Period_map.csv"),df_periodmap)
end

function run_timedomainreduction_scenarios(case::AbstractString,settings_path::AbstractString,setup::Dict)
 
    unpack_scenarios(case,setup)
    
    for s in 1:setup["NumScenarios"]
        scenario_path = joinpath(joinpath(case,"Scenarios"),"Scenario_$s");
        cluster_inputs(scenario_path, settings_path, setup);
    end

    pack_scenarios(case,setup)

    return
end

# NEW: global TDR across all scenarios
function run_timedomainreduction_scenarios_global(case::AbstractString,
                                                  settings_path::AbstractString,
                                                  setup::Dict)

    # Make sure scenario folders with raw 1-year inputs exist
    unpack_scenarios(case, setup)

    # Run a single global TDR on the full multi-year series in `case`
    # In ScenarioTDRMode = 2, MinPeriods/MaxPeriods/WeightTotal are scaled
    cluster_inputs(case, settings_path, setup)

    # (Optional later: distribute_global_tdr_to_scenarios(case, settings_path, setup))
    return
end


function run_timedomainreduction!(case::AbstractString)
    settings_path = get_settings_path(case)                  # Settings YAML path
    genx_settings = get_settings_path(case, "genx_settings.yml")
    mysetup = configure_settings(genx_settings)

    if mysetup["MultiStage"] == 0
        # NEW: read scenario TDR mode from time_domain_reduction_settings.yml
        tdr_settings = get_settings_path(case, "time_domain_reduction_settings.yml")
        TDRSettingsDict = YAML.load(open(tdr_settings))
        ScenarioTDRMode = get(TDRSettingsDict, "ScenarioTDRMode", 1)

        if mysetup["NumScenarios"] <= 1
            # Single-scenario case: always just cluster once
            cluster_inputs(case, settings_path, mysetup)
        else
            # Multi-scenario cases: choose per-scenario vs global behavior
            if ScenarioTDRMode == 1
                # Legacy behavior: 1-year TDR per scenario + pack
                run_timedomainreduction_scenarios(case, settings_path, mysetup)
            elseif ScenarioTDRMode == 2
                # New behavior: one global TDR across all scenarios
                run_timedomainreduction_scenarios_global(case, settings_path, mysetup)
            else
                error("Unexpected value for 'ScenarioTDRMode' in time_domain_reduction_settings.yml. Expected 1 or 2.")
            end
        end

    elseif mysetup["MultiStage"] >= 1
        run_timedomainreduction_multistage!(case)
    else
        error("Unexpected value for key 'MultiStage' in genx_settings.yml. Expected either 0 or 1.")
    end

    return
end


function run_timedomainreduction_multistage!(case::AbstractString)
    # special multistage version
    settings_path = get_settings_path(case)
    genx_settings = get_settings_path(case, "genx_settings.yml")
    mysetup = configure_settings(genx_settings)
    multistage_settings = get_settings_path(case, "multi_stage_settings.yml")

    mysetup["MultiStageSettingsDict"] = YAML.load(open(multistage_settings))

    tdr_settings = get_settings_path(case, "time_domain_reduction_settings.yml")
    TDRSettingsDict = YAML.load(open(tdr_settings))
    if TDRSettingsDict["MultiStageConcatenate"] == 0
        println("Clustering Time Series Data (Individually)...")
        for stage_id in 1:mysetup["MultiStageSettingsDict"]["NumStages"]
            cluster_inputs(case, settings_path, mysetup, stage_id)
        end
    else
        println("Clustering Time Series Data (Grouped)...")
        cluster_inputs(case, settings_path, mysetup)
    end

    return
end

