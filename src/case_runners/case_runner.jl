function get_settings_path(case::AbstractString)
    return joinpath(case, "Settings")
end

function get_settings_path(case::AbstractString, filename::AbstractString)
    return joinpath(get_settings_path(case), filename)
end

function get_default_output_folder(case::AbstractString)
    return joinpath(case, "Results")
end

@doc raw"""Run the GenX in the given folder
case - folder for the case
"""
function run_genx_case!(case::AbstractString)
    genx_settings = get_settings_path(case, "genx_settings.yml") #Settings YAML file path
    setup = configure_settings(genx_settings) # setup dictionary stores settings and GenX-specific parameters

    if setup["MultiStage"] == 0
        if setup["Benders"] == 0
            run_genx_case_simple!(case, setup)
        else
            setup_benders = configure_settings_benders(get_settings_path(case, "benders_settings.yml"))

            setup = merge(setup,setup_benders);
        
            run_genx_case_benders!(case,setup)
        end
    else
        if setup["Benders"] == 0
            run_genx_case_multistage!(case, setup)
        else
            setup_benders = configure_settings_benders(get_settings_path(case, "benders_settings.yml"))
            setup = merge(setup,setup_benders);
            run_genx_case_multistage_benders!(case, setup)
        end
    end
end

function time_domain_reduced_files_exist(tdrpath)
    tdr_load = isfile(joinpath(tdrpath,"Load_data.csv"))
    tdr_genvar = isfile(joinpath(tdrpath,"Generators_variability.csv"))
    tdr_fuels = isfile(joinpath(tdrpath,"Fuels_data.csv"))
    return (tdr_load && tdr_genvar && tdr_fuels)
end

function run_genx_case_simple!(case::AbstractString, setup::Dict)
    settings_path = get_settings_path(case)

    ### Cluster time series inputs if necessary and if specified by the user
    TDRpath = joinpath(case, setup["TimeDomainReductionFolder"])

    if setup["TimeDomainReduction"] == 1
        if setup["NumScenarios"]<=1
          prevent_doubled_timedomainreduction(case)
        end
        if !time_domain_reduced_files_exist(TDRpath)
             println("Clustering Time Series Data (Grouped)...")
             if setup["NumScenarios"]<=1
                cluster_inputs(case, settings_path, setup)
            else
                unpack_scenarios(case,setup) 
                for s in 1:setup["NumScenarios"]
                    scenario_path = joinpath(joinpath(case,"Scenarios"),"Scenario_$s");
                    cluster_inputs(scenario_path, settings_path, setup);
                end
                pack_scenarios(case,setup)
            end
        else
            println("Time Series Data Already Clustered.")
        end
    end

    ### Configure solver
    println("Configuring Solver")
    OPTIMIZER = configure_solver(setup["Solver"], settings_path)

    #### Running a case

    ### Load inputs
    println("Loading Inputs")
    inputs = load_inputs(setup, case)

    println("Generating the Optimization Model")
    time_elapsed = @elapsed EP = generate_model(setup, inputs, OPTIMIZER)
    println("Time elapsed for model building is")
    println(time_elapsed)

    println("Solving Model")
    EP, solve_time = solve_model(EP, setup)
    inputs["solve_time"] = solve_time # Store the model solve time in inputs

    # Run MGA if the MGA flag is set to 1 else only save the least cost solution
    println("Writing Output")
    outputs_path = get_default_output_folder(case)
    elapsed_time = @elapsed write_outputs(EP, outputs_path, setup, inputs)
    println("Time elapsed for writing is")
    println(elapsed_time)
    if setup["ModelingToGenerateAlternatives"] == 1
        println("Starting Model to Generate Alternatives (MGA) Iterations")
        mga(EP, case, setup, inputs, outputs_path)
    end

    if setup["MethodofMorris"] == 1
        println("Starting Global sensitivity analysis with Method of Morris")
        morris(EP, case, setup, inputs, outputs_path, OPTIMIZER)
    end
end

function run_genx_case_benders!(case::AbstractString, setup::Dict)

    settings_path = get_settings_path(case)

    ### Cluster time series inputs if necessary and if specified by the user
    TDRpath = joinpath(case, setup["TimeDomainReductionFolder"])

    if setup["TimeDomainReduction"] == 1
        if setup["NumScenarios"]<=1
            prevent_doubled_timedomainreduction(case)
        end
        println("Clustering Time Series Data (Grouped)...")
        if !time_domain_reduced_files_exist(TDRpath)
            if setup["NumScenarios"]<=1
                cluster_inputs(case, settings_path, setup)
            else
                unpack_scenarios(case,setup) 
                for s in 1:setup["NumScenarios"]
                    scenario_path = joinpath(joinpath(case,"Scenarios"),"Scenario_$s");
                    cluster_inputs(scenario_path, settings_path, setup);
                end
                pack_scenarios(case,setup)
            end
        else
            println("Time Series Data Already Clustered.")
        end
    end

    
    #### Running a case

    ### Load inputs
    println("Loading Inputs")

    inputs = load_inputs(setup, case);

    inputs_decomp = separate_inputs_periods(inputs);

    setup["settings_path"] = settings_path;

    benders_dict = configure_benders_models(setup,inputs,inputs_decomp);

    master_sol_final, EP_master_final, LB_hist, UB_hist, cpu_time,feasibility_hist = benders_decomposition(setup,benders_dict);

    println("Benders decomposition took $(cpu_time[end]) seconds to run")

    println("Writing Output")
    #outputs_path = get_default_output_folder(case)
    if setup["BD_Stab_Method"]=="int_level_set_dynamic" 
        outputs_path = joinpath(case, "Results_Benders_Int_LevelSetDyn")
    elseif setup["BD_Stab_Method"]=="int_level_set" 
        outputs_path = joinpath(case, "Results_Benders_Int_LevelSet")
    elseif setup["BD_Stab_Method"]=="trust_region" 
        outputs_path = joinpath(case, "Results_Benders_TrustRegion")
    elseif setup["BD_Stab_Method"]=="in_out"
        outputs_path = joinpath(case, "Results_Benders_InOut")
    else
        outputs_path = joinpath(case, "Results_Benders")
    end

    # elapsed_time = @elapsed write_outputs(EP, outputs_path, setup, inputs)
    # println("Time elapsed for writing is")
    # println(elapsed_time)
    if setup["OverwriteResults"] == 1
		# Overwrite existing results if dir exists
		# This is the default behaviour when there is no flag, to avoid breaking existing code
		if !(isdir(outputs_path))
		mkdir(outputs_path)
		end
	else
		# Find closest unused ouput directory name and create it
		outputs_path = choose_output_dir(outputs_path)
		mkdir(outputs_path)
	end
    
    elapsed_time = @elapsed write_benders_output(LB_hist,UB_hist,cpu_time,feasibility_hist,outputs_path,setup,inputs,EP_master_final);
    

end

function run_genx_case_multistage_benders!(case::AbstractString, setup::Dict)
    settings_path = get_settings_path(case)
    multistage_settings = get_settings_path(case, "multi_stage_settings.yml") # Multi stage settings YAML file path
    setup["MultiStageSettingsDict"] = YAML.load(open(multistage_settings))

    ### Cluster time series inputs if necessary and if specified by the user
    tdr_settings = get_settings_path(case, "time_domain_reduction_settings.yml") # Multi stage settings YAML file path
    TDRSettingsDict = YAML.load(open(tdr_settings))

    first_stage_path = joinpath(case, "Inputs", "Inputs_p1")
    TDRpath = joinpath(first_stage_path, setup["TimeDomainReductionFolder"])
    if setup["TimeDomainReduction"] == 1
        prevent_doubled_timedomainreduction(first_stage_path)
        if !time_domain_reduced_files_exist(TDRpath)
            if (setup["MultiStage"] >= 1) && (TDRSettingsDict["MultiStageConcatenate"] == 0)
                println("Clustering Time Series Data (Individually)...")
                for stage_id in 1:setup["MultiStageSettingsDict"]["NumStages"]
                    if setup["NumScenarios"]<=1
                        cluster_inputs(case, settings_path, setup, stage_id)
                    else
                        unpack_scenarios(case,setup, stage_id) 
                        for s in 1:setup["NumScenarios"]
                            scenario_path = joinpath(joinpath(case,"Scenarios"),"Scenario_$s");
                            cluster_inputs(scenario_path, settings_path, setup,stage_id);
                        end
                    end
                end
            else
                println("Clustering Time Series Data (Grouped)...")
                cluster_inputs(case, settings_path, setup)
            end
        else
            println("Time Series Data Already Clustered.")
        end
    end

    setup["settings_path"] = settings_path;

    inputs_d,inputs_decomp =  multistage_benders_load_inputs(case,setup);

    benders_dict = configure_benders_models(setup,inputs_d,inputs_decomp);

    master_sol_final, EP_master_final, LB_hist, UB_hist, cpu_time,feasibility_hist = benders_decomposition(setup,benders_dict);

    println("Benders decomposition took $(cpu_time[end]) seconds to run")

    println("Writing Output")

    if setup["BD_Stab_Method"]=="int_level_set" 
        outputs_path = joinpath(case, "Results_Benders_Int_LevelSet")
    elseif setup["BD_Stab_Method"]=="l2_level_set" 
        outputs_path = joinpath(case, "Results_Benders_L2_LevelSet")
    elseif setup["BD_Stab_Method"]=="trust_region" 
        outputs_path = joinpath(case, "Results_Benders_TrustRegion")
    else
        outputs_path = joinpath(case, "Results_Benders")
    end

	if setup["OverwriteResults"] == 1
		# Overwrite existing results if dir exists
		# This is the default behaviour when there is no flag, to avoid breaking existing code
		if !(isdir(outputs_path))
		mkdir(outputs_path)
		end
	else
		# Find closest unused ouput directory name and create it
		outputs_path = choose_output_dir(outputs_path)
		mkdir(outputs_path)
	end
        
    write_benders_output(LB_hist,UB_hist,cpu_time,feasibility_hist,outputs_path,setup,inputs_d,EP_master_final);
end

function run_genx_case_multistage!(case::AbstractString, setup::Dict)
    settings_path = get_settings_path(case)
    multistage_settings = get_settings_path(case, "multi_stage_settings.yml") # Multi stage settings YAML file path
    setup["MultiStageSettingsDict"] = YAML.load(open(multistage_settings))

    ### Cluster time series inputs if necessary and if specified by the user
    tdr_settings = get_settings_path(case, "time_domain_reduction_settings.yml") # Multi stage settings YAML file path
    TDRSettingsDict = YAML.load(open(tdr_settings))

    first_stage_path = joinpath(case, "Inputs", "Inputs_p1")
    TDRpath = joinpath(first_stage_path, setup["TimeDomainReductionFolder"])
    if setup["TimeDomainReduction"] == 1
        prevent_doubled_timedomainreduction(first_stage_path)
        if !time_domain_reduced_files_exist(TDRpath)
            if (setup["MultiStage"] >= 1) && (TDRSettingsDict["MultiStageConcatenate"] == 0)
                println("Clustering Time Series Data (Individually)...")
                for stage_id in 1:setup["MultiStageSettingsDict"]["NumStages"]
                    cluster_inputs(case, settings_path, setup, stage_id)
                end
            else
                println("Clustering Time Series Data (Grouped)...")
                cluster_inputs(case, settings_path, setup)
            end
        else
            println("Time Series Data Already Clustered.")
        end
    end

    ### Configure solver
    println("Configuring Solver")
    OPTIMIZER = configure_solver(setup["Solver"], settings_path)

    model_dict=Dict()
    inputs_dict=Dict()

    for t in 1:setup["MultiStageSettingsDict"]["NumStages"]

        # Step 0) Set Model Year
        setup["MultiStageSettingsDict"]["CurStage"] = t

        # Step 1) Load Inputs
        inpath_sub = joinpath(case, "Inputs", string("Inputs_p",t))

        inputs_dict[t] = load_inputs(setup, inpath_sub)

        inputs_dict[t] = configure_multi_stage_inputs(inputs_dict,setup["MultiStageSettingsDict"],setup["NetworkExpansion"])

        compute_min_cumulative_retirements!(t,inputs_dict)

        # Step 2) Generate model
        model_dict[t] = generate_model(setup, inputs_dict[t], OPTIMIZER)
    end


    ### Solve model
    println("Solving Model")

    # Step 3) Run DDP Algorithm
    ## Solve Model
    model_dict, mystats_d, inputs_dict = run_ddp(model_dict, setup, inputs_dict)

    # Step 4) Write final outputs from each stage

    outpath = get_default_output_folder(case)

    if setup["OverwriteResults"] == 1
        # Overwrite existing results if dir exists
        # This is the default behaviour when there is no flag, to avoid breaking existing code
        if !(isdir(outpath))
            mkdir(outpath)
        end
    else
        # Find closest unused ouput directory name and create it
        outpath = choose_output_dir(outpath)
        mkdir(outpath)
    end

    for p in 1:setup["MultiStageSettingsDict"]["NumStages"]
        outpath_cur = joinpath(outpath, "Results_p$p")
        write_outputs(model_dict[p], outpath_cur, setup, inputs_dict[p])
    end

    # Step 5) Write DDP summary outputs

    write_multi_stage_outputs(mystats_d, outpath, setup, inputs_dict)
end



