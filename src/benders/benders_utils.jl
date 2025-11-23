function configure_benders_models(setup::Dict,inputs::Dict,inputs_decomp::Dict)

	println("Configuring Benders problems")

	MASTER_OPTIMIZER =  configure_benders_master_solver(setup["Solver"], setup["settings_path"]);

    EP_master, master_vars,master_cons = generate_master_problem(setup, inputs, MASTER_OPTIMIZER);

    if setup["BD_Mode"]=="full"
		SUBPROB_OPTIMIZER =  configure_benders_subprob_solver(setup["Solver"], setup["settings_path"]);
        EP_subprob,master_vars_sub = generate_full_operational_subproblem(setup, inputs, SUBPROB_OPTIMIZER,master_vars,master_cons);
    elseif setup["BD_Mode"]=="serial"
		SUBPROB_OPTIMIZER =  configure_benders_subprob_solver(setup["Solver"],setup["settings_path"]);
        EP_subprob, master_vars_sub = generate_decomp_operational_subproblems(setup,inputs_decomp,SUBPROB_OPTIMIZER,master_vars,master_cons);
    elseif setup["BD_Mode"]=="distributed"
        EP_subprob, master_vars_sub = initialize_dist_helpers(setup,inputs_decomp,master_vars,master_cons);
    end

	benders_dict = Dict();
	benders_dict["EP_master"] = EP_master;
	benders_dict["master_vars"] = master_vars;
	benders_dict["master_cons"] = master_cons;
	###### benders_dict["MASTER_OPTIMIZER"] = MASTER_OPTIMIZER;
	benders_dict["EP_subprob"] = EP_subprob;
	benders_dict["master_vars_sub"] = master_vars_sub;
	###### benders_dict["SUBPROB_OPTIMIZER"] = SUBPROB_OPTIMIZER;

	return benders_dict
end

function get_local_master_vars(helper_local::Vector{Dict{Any,Any}})

    local_vars=Dict();

    for m in helper_local
		w = m["SubPeriod"];
        local_vars[w] = m["master_vars_sub"]
    end

    return local_vars


end


function fix_master_vars!(EP::Model,master_sol::NamedTuple,master_vars_sub::Vector{String})
	for y in master_vars_sub
		vy = variable_by_name(EP,y);
		fix(vy,master_sol.values[y];force=true)
		if is_integer(vy)
			unset_integer(vy)
		elseif is_binary(vy)
			unset_binary(vy)
		end
	end
end



function subperiods2scenarios(W::Int64, numscenarios::Int64)
	# Breaks up the subperiods into numscenarios scenarios
	scenario_times = Array{Int64}(undef, numscenarios, 2)
	min_scenario_length = Int64(floor(W/numscenarios))
	remainder = W % numscenarios
	for s in 1:numscenarios
		scenario_times[s,1] = (s-1)*min_scenario_length + 1 + remainder
		scenario_times[s,2] = s*min_scenario_length + remainder
	end
	scenario_times[1,1] = 1
	return scenario_times
end


function map_from_stage_to_subperiod(inputs,num_stages)

	nW=0;
	nco2cap = Dict();
	nesr = Dict();
	stage_subperiods = Dict()
    for t in 1:num_stages
		stage_subperiods[t] = nW+1:nW+inputs[t]["REP_PERIOD"];
        for p in 1:inputs[t]["REP_PERIOD"]
            nW = nW +1;
			if in("NCO2Cap",keys(inputs[t]))
            	nco2cap[nW] = inputs[t]["NCO2Cap"]
			else
				nco2cap[nW] = 0;
			end
			if in("nESR",keys(inputs[t]))
				nesr[nW] = inputs[t]["nESR"]
			else
				nesr[nW] = 0;
			end
        end

    end

	return nW,stage_subperiods,nco2cap, nesr

end


function separate_inputs_periods(inputs::Dict)

    inputs_all=Dict();
    number_periods = inputs["REP_PERIOD"];
    hours_per_subperiod = inputs["hours_per_subperiod"];
    
    ####### entries_to_be_changed = ["omega","REP_PERIOD","C_Fuel_per_MWh","INTERIOR_SUBPERIODS","START_SUBPERIODS","pP_Max","T","fuel_costs","Weights","pD","C_Start"];

    for w in 1:number_periods
        inputs_all[w] = deepcopy(inputs);
        Tw = (w-1)*hours_per_subperiod+1:w*hours_per_subperiod;
        inputs_all[w]["omega"] = inputs["omega"][Tw];
        inputs_all[w]["REP_PERIOD"]=1;
        inputs_all[w]["C_Fuel_per_MWh"] = inputs["C_Fuel_per_MWh"][:,Tw];
        STARTS = 1:hours_per_subperiod:hours_per_subperiod;
        INTERIORS = setdiff(1:hours_per_subperiod,STARTS);   
        inputs_all[w]["INTERIOR_SUBPERIODS"] = INTERIORS;
        inputs_all[w]["START_SUBPERIODS"] = STARTS;
        inputs_all[w]["pP_Max"] = inputs["pP_Max"][:,Tw];
        inputs_all[w]["T"] = hours_per_subperiod;
        for ks in keys(inputs["fuel_costs"])
            inputs_all[w]["fuel_costs"][ks] = inputs["fuel_costs"][ks][Tw];
        end
        inputs_all[w]["Weights"] = [inputs["Weights"][w]];
        inputs_all[w]["pD"] = inputs["pD"][Tw,:];
        inputs_all[w]["C_Start"] = inputs["C_Start"][:,Tw]; 
        inputs_all[w]["SubPeriod"] = w;
		inputs_all[w]["NumScenarios"]=1;
		if haskey(inputs,"Period_Map")
			inputs_all[w]["SubPeriod_Index"] = inputs["Period_Map"].Rep_Period[findfirst(inputs["Period_Map"].Rep_Period_Index.==w)];
		end
		if haskey(inputs,"scenario_probability")
			scenario_periods = subperiods2scenarios(number_periods, length(inputs["scenario_probability"]));
			sw = findfirst(in(w,scenario_periods[s,1]:scenario_periods[s,2]) for s in 1:length(inputs["scenario_probability"]));
			inputs_all[w]["scenario_probability"] = inputs["scenario_probability"][sw];
		else
			inputs_all[w]["scenario_probability"]=1.0;
		end
    end

    return inputs_all

end

function multistage_benders_load_inputs(case::AbstractString, mysetup::Dict)
    inputs_decomp=Dict()
    inputs_d=Dict()
    w=0;
    for t in 1:mysetup["MultiStageSettingsDict"]["NumStages"]
        mysetup["MultiStageSettingsDict"]["CurStage"] = t

        inpath_sub = joinpath(case, "Inputs", string("Inputs_p",t))

        inputs_d[t] = load_inputs(mysetup, inpath_sub)

        inputs_d[t] = configure_multi_stage_inputs(inputs_d,mysetup["MultiStageSettingsDict"],mysetup["NetworkExpansion"])

		inputs_d[t]["CurStage"] = t;
		inputs_d[t]["SubPeriod"] = t;
        compute_min_cumulative_retirements!(t,inputs_d)

        inputs_decomp_t = separate_inputs_periods(inputs_d[t]);
        for p in 1:inputs_d[t]["REP_PERIOD"]
            w = w +1;
            inputs_decomp[w] = inputs_decomp_t[p];
            inputs_decomp[w]["SubPeriod"] = w;
            inputs_decomp[w]["CurStage"] = t;
			if haskey(inputs_d[t],"Period_Map")
				inputs_decomp[w]["SubPeriod_Index"] = sum(size(inputs_d[s]["Period_Map"],1) for s in 1:t-1; init=0) +  inputs_decomp[w]["SubPeriod_Index"]
			end
        end
    end
    return inputs_d,inputs_decomp
end

function splitfun(x)
	return String(split(x,"[")[1])
end

function configure_benders_master_solver(solver::String, solver_settings_path::String)

	gurobi_settings_path = joinpath(solver_settings_path, "gurobi_benders_master_settings.yml")

	mysettings = convert(Dict{String, Any}, YAML.load(open(gurobi_settings_path)))

	settings = Dict("Crossover"=>0,"Method"=>2);

	attributes = merge(settings, mysettings)
	println("Master Gurobi attributes:")
	display(attributes)

    OPTIMIZER = optimizer_with_attributes(()->Gurobi.Optimizer(GRB_ENV),attributes...)
	#OPTIMIZER = optimizer_with_attributes(Gurobi.Optimizer,attributes...)

	# OPTIMIZER = optimizer_with_attributes(()->Gurobi.Optimizer(GRB_ENV),"Crossover"=>0,"Method"=>2)
	
	return OPTIMIZER
end


function configure_benders_subprob_solver(solver::String, solver_settings_path::String)

	gurobi_settings_path = joinpath(solver_settings_path, "gurobi_benders_subprob_settings.yml")

	mysettings = convert(Dict{String, Any}, YAML.load(open(gurobi_settings_path)))

	settings = Dict("Crossover"=>1,"Method"=>2);

	attributes = merge(settings, mysettings)

	println("Subproblem Gurobi attributes:")
	display(attributes)

    OPTIMIZER = optimizer_with_attributes(()->Gurobi.Optimizer(GRB_ENV),attributes...)
	#OPTIMIZER = optimizer_with_attributes(Gurobi.Optimizer,attributes...)

	# OPTIMIZER = optimizer_with_attributes(()->Gurobi.Optimizer(GRB_ENV),"Crossover"=>1,"Method"=>2)

	return OPTIMIZER
end

function configure_settings_benders(settings_path::String)
	
	println("Configuring Benders Settings")
	model_settings = YAML.load(open(settings_path))

	settings = Dict{Any,Any}("BD_Mode"=>"full",
							"BD_ConvTol"=>1e-3,
							"BD_MaxIter"=>300,
							"BD_MaxCpuTime"=>24*3600,
							"BD_StabParam"=>0.0,
							"BD_Stab_Method"=>"off");

	merge!(settings, model_settings)

return settings
end



function check_negative_capacities(EP::Model)

	neg_cap_bool = false;
	tol = -1e-8;
	if any(value.(EP[:eTotalCap]).< tol) 
			neg_cap_bool = true;
	elseif haskey(EP,:eTotalCapEnergy)
		if any(value.(EP[:eTotalCapEnergy]).< tol)
			neg_cap_bool = true;
		end
	elseif haskey(EP,:eTotalCapCharge)
		if any(value.(EP[:eTotalCapCharge]).< tol)
			neg_cap_bool = true;
		end
	elseif haskey(EP,:eAvail_Trans_Cap)
		if any(value.(EP[:eAvail_Trans_Cap]).< tol)
			neg_cap_bool = true;
		end
	end
	return neg_cap_bool
	
end

function write_benders_output(LB_hist::Vector{Float64},UB_hist::Vector{Float64},cpu_time::Vector{Float64},feasibility_hist::Vector{Float64},outpath::AbstractString, setup::Dict,inputs::Dict,EP_master::Model)
	
	dfConv = DataFrame(Iter = 1:length(LB_hist),CPU_Time = cpu_time, LB = LB_hist, UB  = UB_hist, Gap = (UB_hist.-LB_hist)./LB_hist,Feasibility=feasibility_hist)

	if !has_values(EP_master)
		optimize!(EP_master)
	end

	if setup["MultiStage"]==1
		write_master_solution_multistage(outpath, inputs, setup, EP_master)
	else
		write_master_solution_singlestage(outpath, inputs, setup, EP_master)
	end

	CSV.write(joinpath(outpath, "benders_convergence.csv"),dfConv)

	YAML.write_file(joinpath(outpath, "run_settings.yml"),setup)

end

function write_master_solution_multistage(path::AbstractString, inputs_d::Dict, setup::Dict, EP_master::Model)

    number_of_periods = length(inputs_d);

    cap = zeros(number_of_periods,length(inputs_d[1]["RESOURCES"]))
    retcap = zeros(number_of_periods,length(inputs_d[1]["RESOURCES"]))
    existingcap = value.(EP_master[:eExistingCap])
	totcap = value.(EP_master[:eTotalCap]);

    capenergy = zeros(number_of_periods,length(inputs_d[1]["RESOURCES"]))
    retcapenergy = zeros(number_of_periods,length(inputs_d[1]["RESOURCES"]))
    existingcapenergy = zeros(length(inputs_d[1]["RESOURCES"]))
	totcapenergy = zeros(number_of_periods,length(inputs_d[1]["RESOURCES"]))

    capcharge = zeros(number_of_periods,length(inputs_d[1]["RESOURCES"]))
    retcapcharge = zeros(number_of_periods,length(inputs_d[1]["RESOURCES"]))
    existingcapcharge = zeros(length(inputs_d[1]["RESOURCES"]))
	totcapcharge = zeros(number_of_periods,length(inputs_d[1]["RESOURCES"]))


	capsize = inputs_d[1]["dfGen"][!,:Cap_Size];

    for i in inputs_d[1]["NEW_CAP"]
        for p in 1:number_of_periods
			cap[p,i] = value(EP_master[:vCAP][p,i]);
        end
    end
    
    for i in inputs_d[1]["RET_CAP"]
        for p in 1:number_of_periods
			retcap[p,i] = value(EP_master[:vRETCAP][p,i]);
        end
    end
    
    for i in inputs_d[1]["STOR_ASYMMETRIC"]
        existingcapcharge[i] = value(EP_master[:eExistingCapCharge][i]);
        for p in 1:number_of_periods
			totcapcharge[p,i] = value(EP_master[:eTotalCapCharge][p,i])
            if i in inputs_d[1]["NEW_CAP_CHARGE"]
                capcharge[p,i] = value(EP_master[:vCAPCHARGE][p,i])
            end
            if i in inputs_d[1]["RET_CAP_CHARGE"]
                retcapcharge[p,i] = value(EP_master[:vRETCAPCHARGE][p,i])
            end
        end
    end
    
    for i in inputs_d[1]["STOR_ALL"]
        existingcapenergy[i] = value(EP_master[:eExistingCapEnergy][i]);
        for p in 1:number_of_periods
			totcapenergy[p,i] = value(EP_master[:eTotalCapEnergy][p,i])
            if i in inputs_d[1]["NEW_CAP_ENERGY"]
				capenergy[p,i] = value(EP_master[:vCAPENERGY][p,i])
            end
            if i in inputs_d[1]["RET_CAP_ENERGY"]
				retcapenergy[p,i] = value(EP_master[:vRETCAPENERGY][p,i])
            end
        end
    end

    scale_factor = setup["ParameterScale"] == 1 ? ModelScalingFactor : 1

	capsize = capsize.*scale_factor;
	existingcap = existingcap.*scale_factor;
	existingcapcharge = existingcapcharge.*scale_factor;
	existingcapenergy = existingcapenergy.*scale_factor;
	totcap = totcap.*scale_factor;
	totcapcharge = totcapcharge.*scale_factor;
	totcapenergy = totcapenergy.*scale_factor;

    auxnames = ["Resource";"Zone";"ExistingCap_MW";["EndCap_MW_p$t" for t in 1:number_of_periods];"Cap_Size_MW";["RetCap_p$t" for t in 1:number_of_periods];["NewCap_p$t" for t in 1:number_of_periods]]
    dfCap = DataFrame([inputs_d[1]["RESOURCES"] inputs_d[1]["dfGen"][!,:Zone] existingcap reduce(hcat,[totcap[p,:] for p in 1:number_of_periods]) capsize reduce(hcat,[retcap[p,:] for p in 1:number_of_periods]) reduce(hcat,[cap[p,:] for p in 1:number_of_periods])],auxnames)
    CSV.write(joinpath(path, "master_sol_capacity.csv"), dfCap)

	if !isempty(inputs_d[1]["STOR_ASYMMETRIC"])
		auxnames_charge = ["Resource";"Zone";"ExistingCapCharge_MW";["EndCapCharge_MW_p$t" for t in 1:number_of_periods];"Cap_Size_MW";["RetCapCharge_p$t" for t in 1:number_of_periods];["NewCapCharge_p$t" for t in 1:number_of_periods]]
    	dfCapCharge = DataFrame([inputs_d[1]["RESOURCES"] inputs_d[1]["dfGen"][!,:Zone] existingcapcharge reduce(hcat,[totcapcharge[p,:] for p in 1:number_of_periods]) capsize reduce(hcat,[retcapcharge[p,:] for p in 1:number_of_periods]) reduce(hcat,[capcharge[p,:] for p in 1:number_of_periods])],auxnames_charge)
    	CSV.write(joinpath(path, "master_sol_capacity_charge.csv"), dfCapCharge)
	end

	if !isempty(inputs_d[1]["STOR_ALL"])
		auxnames_energy = ["Resource";"Zone";"ExistingCapEnergy_MWh";["EndCapEnergy_MWh_p$t" for t in 1:number_of_periods];"Cap_Size_MWh";["RetCapEnergy_p$t" for t in 1:number_of_periods];["NewCapEnergy_p$t" for t in 1:number_of_periods]]
    	dfCapEnergy = DataFrame([inputs_d[1]["RESOURCES"] inputs_d[1]["dfGen"][!,:Zone] existingcapenergy reduce(hcat,[totcapenergy[p,:] for p in 1:number_of_periods]) capsize reduce(hcat,[retcapenergy[p,:] for p in 1:number_of_periods]) reduce(hcat,[capenergy[p,:] for p in 1:number_of_periods])],auxnames_energy)
		CSV.write(joinpath(path, "master_sol_capacity_energy.csv"), dfCapEnergy)
	end

	newtranscap = zeros(number_of_periods,inputs_d[1]["L"])
    existingtranscap = inputs_d[1]["pTrans_Max"]*scale_factor;
	Line_Size = inputs_d[1]["pLine_Size"]*scale_factor;
	if setup["NetworkExpansion"]==1
		for l in inputs_d[1]["EXPANSION_LINES"]
			for p in 1:number_of_periods
				newtranscap[p,l] = value(EP_master[:vNEW_TRANS_CAP][p,l]);
			end
		end
	end

	tottranscap = value.(EP_master[:eAvail_Trans_Cap])*scale_factor;

	auxnames = ["Line";"ExistingTransCap_MW";["EndTransCap_MW_p$t" for t in 1:number_of_periods];"Line_Size_MW";["NewTransCap_p$t" for t in 1:number_of_periods]]
	dfTransCap = DataFrame([ ["L$i" for i in 1: inputs_d[1]["L"]] existingtranscap reduce(hcat,[tottranscap[p,:] for p in 1:number_of_periods]) Line_Size reduce(hcat,[newtranscap[p,:] for p in 1:number_of_periods])],auxnames)

    CSV.write(joinpath(path, "master_sol_trans_capacity.csv"), dfTransCap)

	### Writing costs
	cFix =zeros(number_of_periods)
	for t in 1:number_of_periods
		cFix[t] = value(EP_master[:eTotalCFix][t]) + (!isempty(inputs_d[1]["STOR_ALL"]) ? value(EP_master[:eTotalCFixEnergy][t]) : 0.0) + (!isempty(inputs_d[1]["STOR_ASYMMETRIC"]) ? value(EP_master[:eTotalCFixCharge][t]) : 0.0)
	end
	cNetworkExp = value.(EP_master[:eTotalCNetworkExp])

	cUnmetCapReq = zeros(number_of_periods);
	for t in 1:number_of_periods
		cUnmetCapReq[t] = (setup["MinCapReq"]==1 ? value(EP_master[:eTotalCMinCapSlack][t]) : 0.0) +  (setup["MaxCapReq"]==1 ? value(EP_master[:eTotalCMaxCapSlack][t]) : 0.0) 
	end

 	# Conversion from Million$ to $
	cFix *= scale_factor^2;
	cNetworkExp *= scale_factor^2;
	cUnmetCapReq *= scale_factor^2;

	dfCost = DataFrame([1:number_of_periods cFix cNetworkExp cUnmetCapReq],["Period";"FixCost";"NetworkExpCost";"UnmetCapReqCost"])
	
	CSV.write(joinpath(path, "master_sol_costs.csv"), dfCost)


end


function write_master_solution_singlestage(path::AbstractString, inputs::Dict, setup::Dict, EP_master::Model)


    cap = zeros(length(inputs["RESOURCES"]))
    retcap = zeros(length(inputs["RESOURCES"]))
    existingcap = value.(EP_master[:eExistingCap])
	totcap = value.(EP_master[:eTotalCap]);

    capenergy = zeros(length(inputs["RESOURCES"]))
    retcapenergy = zeros(length(inputs["RESOURCES"]))
    existingcapenergy = zeros(length(inputs["RESOURCES"]))
	totcapenergy = zeros(length(inputs["RESOURCES"]))

    capcharge = zeros(length(inputs["RESOURCES"]))
    retcapcharge = zeros(length(inputs["RESOURCES"]))
    existingcapcharge = zeros(length(inputs["RESOURCES"]))
	totcapcharge = zeros(length(inputs["RESOURCES"]))

	capsize = inputs["dfGen"][!,:Cap_Size];
	
    for i in inputs["NEW_CAP"]
		cap[i] = value(EP_master[:vCAP][i]);
    end
    
    for i in inputs["RET_CAP"]
		retcap[i] = value(EP_master[:vRETCAP][i]);
    end
    
    for i in inputs["STOR_ASYMMETRIC"]
        existingcapcharge[i] = value(EP_master[:eExistingCapCharge][i]);
		totcapcharge[i] = value(EP_master[:eTotalCapCharge][i])
        if i in inputs["NEW_CAP_CHARGE"]
            capcharge[i] = value(EP_master[:vCAPCHARGE][i])
        end
        if i in inputs["RET_CAP_CHARGE"]
            retcapcharge[i] = value(EP_master[:vRETCAPCHARGE][i])
        end
    end
    
    for i in inputs["STOR_ALL"]
        existingcapenergy[i] = value(EP_master[:eExistingCapEnergy][i]);
		totcapenergy[i] = value(EP_master[:eTotalCapEnergy][i])
        if i in inputs["NEW_CAP_ENERGY"]
            capenergy[i] = value(EP_master[:vCAPENERGY][i])
        end
        if i in inputs["RET_CAP_ENERGY"]
            retcapenergy[i] = value(EP_master[:vRETCAPENERGY][i])
        end
    end

	scale_factor = setup["ParameterScale"] == 1 ? ModelScalingFactor : 1

	capsize = capsize.*scale_factor;
	existingcap = existingcap.*scale_factor;
	existingcapcharge = existingcapcharge.*scale_factor;
	existingcapenergy = existingcapenergy.*scale_factor;
	totcap = totcap.*scale_factor;
	totcapcharge = totcapcharge.*scale_factor;
	totcapenergy = totcapenergy.*scale_factor;

    auxnames = ["Resource";"Zone";"ExistingCap_MW";"EndCap_MW";"Cap_Size_MW";"RetCap";"NewCap"]
    dfCap = DataFrame([inputs["RESOURCES"] inputs["dfGen"][!,:Zone] existingcap totcap capsize retcap cap],auxnames)
    CSV.write(joinpath(path, "master_sol_capacity.csv"), dfCap)

	if !isempty(inputs["STOR_ASYMMETRIC"])
		auxnames_charge = ["Resource";"Zone";"ExistingCapCharge_MW";"EndCapCharge_MW";"Cap_Size_MW";"RetCapCharge";"NewCapCharge"];
    	dfCapCharge = DataFrame([inputs["RESOURCES"] inputs["dfGen"][!,:Zone] existingcapcharge totcapcharge capsize retcapcharge capcharge],auxnames_charge)
    	CSV.write(joinpath(path, "master_sol_capacity_charge.csv"), dfCapCharge)
	end

	if !isempty(inputs["STOR_ALL"])
		auxnames_energy = ["Resource";"Zone";"ExistingCapEnergy_MWh";"EndCapEnergy_MWh";"Cap_Size_MWh";"RetCapEnergy";"NewCapEnergy"]
    	dfCapEnergy = DataFrame([inputs["RESOURCES"] inputs["dfGen"][!,:Zone] existingcapenergy totcapenergy capsize retcapenergy capenergy],auxnames_energy)
		CSV.write(joinpath(path, "master_sol_capacity_energy.csv"), dfCapEnergy)
	end

	newtranscap = zeros(inputs["L"])
    existingtranscap = inputs["pTrans_Max"]*scale_factor;
	Line_Size = inputs["pLine_Size"]*scale_factor;
	if setup["NetworkExpansion"]==1
		for l in inputs["EXPANSION_LINES"]
			newtranscap[l] = value(EP_master[:vNEW_TRANS_CAP][l]);
		end
	end
	tottranscap = value.(EP_master[:eAvail_Trans_Cap])*scale_factor;

	auxnames = ["Line";"ExistingTransCap_MW";"EndTransCap_MW";"Line_Size_MW";"NewTransCap"]
	dfTransCap = DataFrame([ ["L$i" for i in 1: inputs["L"]] existingtranscap tottranscap Line_Size newtranscap],auxnames)

    CSV.write(joinpath(path, "master_sol_trans_capacity.csv"), dfTransCap)

	### Writing costs
	cFix = value(EP_master[:eTotalCFix]) + (!isempty(inputs["STOR_ALL"]) ? value(EP_master[:eTotalCFixEnergy]) : 0.0) + (!isempty(inputs["STOR_ASYMMETRIC"]) ? value(EP_master[:eTotalCFixCharge]) : 0.0)

	cNetworkExp = value.(EP_master[:eTotalCNetworkExp])

	cUnmetCapReq = (setup["MinCapReq"]==1 ? value(EP_master[:eTotalCMinCapSlack]) : 0.0) +  (setup["MaxCapReq"]==1 ? value(EP_master[:eTotalCMaxCapSlack]) : 0.0) 
	
 	# Conversion from Million$ to $
	cFix *= scale_factor^2;
	cNetworkExp *= scale_factor^2;
	cUnmetCapReq *= scale_factor^2;

	dfCost = DataFrame([cFix cNetworkExp cUnmetCapReq],["FixCost";"NetworkExpCost";"UnmetCapReqCost"])
	
	CSV.write(joinpath(path, "master_sol_costs.csv"), dfCost)

end

