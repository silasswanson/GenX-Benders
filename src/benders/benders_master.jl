function generate_master_problem(setup::Dict,inputs::Dict,OPTIMIZER::MOI.OptimizerWithAttributes)

	if setup["MultiStage"]==1
		EP,varnames,consnames = generate_multistage_master(setup,inputs,OPTIMIZER)
	else
		EP,varnames,consnames = generate_singlestage_master(setup,inputs,OPTIMIZER)
	end

    return EP,varnames,consnames
end


function generate_singlestage_master(setup::Dict,inputs::Dict,OPTIMIZER::MOI.OptimizerWithAttributes)

	## Start pre-solve timer
	master_generation_start_time = time()

	# Generate Master Model
	EP = Model(OPTIMIZER)

	# Introduce dummy variable fixed to zero to ensure that expressions like eTotalCap,
	# eTotalCapCharge, eTotalCapEnergy and eAvail_Trans_Cap all have a JuMP variable
	@variable(EP, vZERO == 1);

	# Initialize Objective Function Expression
	@expression(EP, eObj, 0)

	if (setup["MinCapReq"] == 1)
		@expression(EP, eMinCapRes[mincap = 1:inputs["NumberOfMinCapReqs"]], 0)
	end

	if setup["MaxCapReq"] == 1
		@expression(EP, eMaxCapRes[maxcap = 1:inputs["NumberOfMaxCapReqs"]], 0)
	end

	investment_discharge!(EP, inputs, setup)

    investment_transmission!(EP,inputs,setup)

	if !isempty(inputs["STOR_ALL"])
		investment_energy!(EP, inputs, setup)
	end

    if !isempty(inputs["STOR_ASYMMETRIC"])
		investment_charge!(EP, inputs, setup)
	end

    # Model constraints, variables, expression related to retrofit technologies
	if !isempty(inputs["RETRO"])
		EP = retrofit(EP, inputs)
	end

    if (setup["MinCapReq"] == 1)
		minimum_capacity_requirement!(EP, inputs, setup)
	end

	
	if setup["MaxCapReq"] == 1
		maximum_capacity_requirement!(EP, inputs, setup)
	end


	if setup["IntegerInvestments"]==1
		integer_investments!(EP,inputs,setup)
	end

	if setup["BD_Mode"]=="full"
        #### Single benders cut model
		@variable(EP,vTHETA[1:1]>=0)
	else
		
		if inputs["REP_PERIOD"] > 1 && !isempty(inputs["STOR_LONG_DURATION"])
			long_duration_storage_master!(EP, inputs,setup)
		end

		if inputs["REP_PERIOD"] > 1 && !isempty(inputs["STOR_HYDRO_LONG_DURATION"])
			hydro_inter_period_linkage_master!(EP, inputs,setup)
		end

		@variable(EP,vTHETA[1:inputs["REP_PERIOD"]]>=0)
		
		if setup["CO2Cap"]>=1
			co2_cap_master!(EP,inputs,setup)
		end
		
		if setup["EnergyShareRequirement"] >= 1
			energy_share_requirement_master!(EP, inputs, setup)
		end
	end

	@objective(EP,Min, EP[:eObj] + sum(vTHETA))
	 
    ## Record pre-solver time
	master_generation_time = time() - master_generation_start_time

    println("Master problem generation took $master_generation_time seconds")

	varnames = name.(setdiff(all_variables(EP),[EP[:vZERO];EP[:vTHETA]]));
	consnames = name.(all_constraints(EP,include_variable_in_set_constraints=false));
	set_silent(EP);

    return EP,varnames,consnames
end


function generate_multistage_master(setup::Dict,inputs::Dict,OPTIMIZER::MOI.OptimizerWithAttributes)
	# OPEX multiplier scales fixed costs to account for multiple years between two model stages
	# We have already accounted for the opex multipliers on fixed costs in configure_multi_stage_inputs.jl

	DF = Dict([t => 1 / (1 + setup["MultiStageSettingsDict"]["WACC"])^(setup["MultiStageSettingsDict"]["StageLengths"][t] * (t - 1))  for t in 1:setup["MultiStageSettingsDict"]["NumStages"]]);

	## Start pre-solve timer
	master_generation_start_time = time()

	# Generate Master Model
	EP = Model(OPTIMIZER)
 
	# Introduce dummy variable fixed to zero to ensure that expressions like eTotalCap,
	# eTotalCapCharge, eTotalCapEnergy and eAvail_Trans_Cap all have a JuMP variable
	@variable(EP, vZERO == 1);

   	# Initialize Objective Function Expression
	@expression(EP, eObj, 0)

	if (setup["MinCapReq"] == 1)
		@expression(EP, eMinCapRes[t = 1:setup["MultiStageSettingsDict"]["NumStages"], mincap = 1:inputs[t]["NumberOfMinCapReqs"]], 0)
	end

	if setup["MaxCapReq"] == 1
		@expression(EP, eMaxCapRes[t = 1:setup["MultiStageSettingsDict"]["NumStages"], maxcap = 1:inputs[t]["NumberOfMaxCapReqs"]], 0)
	end

	multistage_investment_discharge!(EP, inputs, setup)
	EP[:eObj] += sum(DF[t]*EP[:eTotalCFix][t] for t in 1:setup["MultiStageSettingsDict"]["NumStages"]);

	multistage_endogenous_retirement_discharge!(EP, inputs, setup)

	multistage_investment_transmission!(EP,inputs,setup)
	EP[:eObj] += sum(DF[t]*EP[:eTotalCNetworkExp][t] for t in 1:setup["MultiStageSettingsDict"]["NumStages"]);

	varnames = [name.(EP[:vCAP].data[:]);name.(EP[:vRETCAP].data[:]);name.(EP[:vNEW_TRANS_CAP].data[:])];

	if !isempty(inputs[1]["STOR_ALL"])
		multistage_investment_energy!(EP, inputs, setup)
		append!(varnames,name.(EP[:vCAPENERGY].data[:]))
		append!(varnames,name.(EP[:vRETCAPENERGY].data[:]))
		EP[:eObj] += sum(DF[t]*EP[:eTotalCFixEnergy][t] for t in 1:setup["MultiStageSettingsDict"]["NumStages"])
		multistage_endogenous_retirement_energy!(EP, inputs,setup)
	end

	if !isempty(inputs[1]["STOR_ASYMMETRIC"])
		multistage_investment_charge!(EP, inputs, setup)
		append!(varnames,name.(EP[:vCAPCHARGE].data[:]))
		append!(varnames,name.(EP[:vRETCAPCHARGE].data[:]))
		EP[:eObj] += sum(DF[t]*EP[:eTotalCFixCharge][t] for t in 1:setup["MultiStageSettingsDict"]["NumStages"]);
		multistage_endogenous_retirement_charge!(EP, inputs,setup)
	end

	if !isempty(inputs[1]["RETRO"])
		multistage_retrofit!(EP, inputs, setup)
	end

	if (setup["MinCapReq"] == 1)
		multistage_minimum_capacity_requirement!(EP, inputs, setup)
		if haskey(EP,:eTotalCMinCapSlack)
			EP[:eObj] += sum(EP[:eTotalCMinCapSlack][t] for t in 1:setup["MultiStageSettingsDict"]["NumStages"])
		end
	end

	if setup["MaxCapReq"] == 1
		multistage_maximum_capacity_requirement!(EP, inputs, setup)
		if haskey(EP,:eTotalCMaxCapSlack)
			EP[:eObj] += sum(EP[:eTotalCMaxCapSlack][t] for t in 1:setup["MultiStageSettingsDict"]["NumStages"])
		end
	end

	if setup["IntegerInvestments"]==1
		integer_investments!(EP,inputs,setup)
	end
	
	if setup["BD_Mode"]=="full"
       	#### Single benders cut model for each stage
		@variable(EP,vTHETA[1:setup["MultiStageSettingsDict"]["NumStages"]]>=0)
	else
		nW = sum(inputs[t]["REP_PERIOD"] for t in 1:setup["MultiStageSettingsDict"]["NumStages"])
		@variable(EP,vTHETA[1:nW]>=0)

		if inputs[1]["REP_PERIOD"] > 1 && !isempty(inputs[1]["STOR_LONG_DURATION"])
			multistage_long_duration_storage_master!(EP, inputs,setup)
			append!(varnames,name.(collect(values(EP[:vSOCw].data))))
			append!(varnames,name.(collect(values(EP[:vdSOC].data))))
		end

		if inputs[1]["REP_PERIOD"] > 1 && !isempty(inputs[1]["STOR_HYDRO_LONG_DURATION"])
			multistage_hydro_inter_period_linkage_master!(EP, inputs,setup)
			append!(varnames,name.(collect(values(EP[:vSOC_HYDROw].data))))
			append!(varnames,name.(collect(values(EP[:vdSOC_HYDRO].data))))
		end

		if setup["CO2Cap"]>=1
			multistage_co2_cap_master!(EP,inputs,setup)
			append!(varnames,name.(collect(values(EP[:vCO2budget].data))))
		end
		if setup["EnergyShareRequirement"] >= 1
			multistage_energy_share_requirement_master!(EP, inputs, setup)
			append!(varnames,name.(collect(values(EP[:vESRbudget].data))))
		end
	end

	@objective(EP,Min, EP[:eObj] + sum(vTHETA))
	 
   	## Record pre-solver time
	master_generation_time = time() - master_generation_start_time

   	println("Master problem generation took $master_generation_time seconds")
	
	consnames = name.(all_constraints(EP,include_variable_in_set_constraints=false));
	set_silent(EP);

	return EP,varnames,consnames
end



function load_built_charge_cap!(EP::Model,setup::Dict,inputs::Dict)

	cur_stage = inputs["CurStage"];

	STOR_ASYMMETRIC = inputs["STOR_ASYMMETRIC"] # Set of storage resources with asymmetric (separte) charge/discharge capacity components

	NEW_CAP_CHARGE = inputs["NEW_CAP_CHARGE"] # Set of asymmetric charge/discharge storage resources eligible for new charge capacity
	RET_CAP_CHARGE = inputs["RET_CAP_CHARGE"] # Set of asymmetric charge/discharge storage resources eligible for charge capacity retirements

	### Variables ###

	## Storage capacity built and retired for storage resources with independent charge and discharge power capacities (STOR=2)

	# New installed charge capacity of resource "y"
	@variable(EP, vCAPCHARGE[t in 1:cur_stage, y in NEW_CAP_CHARGE] >= 0)

	# Retired charge capacity of resource "y" from existing capacity
	@variable(EP, vRETCAPCHARGE[t in 1:cur_stage, y in RET_CAP_CHARGE] >= 0)

	### Expressions ###
	@expression(EP, eExistingCapCharge[y in STOR_ASYMMETRIC], inputs["dfGen"][y,:Existing_Charge_Cap_MW]*EP[:vZERO])

	@expression(EP, eTotalCapCharge[y in STOR_ASYMMETRIC],
	if (y in intersect(NEW_CAP_CHARGE, RET_CAP_CHARGE))
		eExistingCapCharge[y] + inputs["dfGen"][y,:Cap_Size]*sum(EP[:vCAPCHARGE][s,y] - EP[:vRETCAPCHARGE][s,y] for s in 1:cur_stage)
	elseif (y in setdiff(NEW_CAP_CHARGE, RET_CAP_CHARGE))
		eExistingCapCharge[y] + inputs["dfGen"][y,:Cap_Size]*sum(EP[:vCAPCHARGE][s,y] for s in 1:cur_stage)
	elseif (y in setdiff(RET_CAP_CHARGE, NEW_CAP_CHARGE))
		eExistingCapCharge[y] - inputs["dfGen"][y,:Cap_Size]*sum(EP[:vRETCAPCHARGE][s,y] for s in 1:cur_stage)
	else
		eExistingCapCharge[y]
	end
	)

end

function load_built_energy_cap!(EP::Model,setup::Dict,inputs::Dict)
	cur_stage = inputs["CurStage"];
	STOR_ALL = inputs["STOR_ALL"] # Set of all storage resources
	NEW_CAP_ENERGY = inputs["NEW_CAP_ENERGY"] # Set of all storage resources eligible for new energy capacity
	RET_CAP_ENERGY = inputs["RET_CAP_ENERGY"] # Set of all storage resources eligible for energy capacity 

	### Variables ###

	## Energy storage reservoir capacity (MWh capacity) built/retired for storage with variable power to energy ratio (STOR=1 or STOR=2)

	# New installed energy capacity of resource "y"
	@variable(EP, vCAPENERGY[t in 1:cur_stage, y in NEW_CAP_ENERGY] >= 0)

	# Retired energy capacity of resource "y" from existing capacity
	@variable(EP, vRETCAPENERGY[t in 1:cur_stage, y in RET_CAP_ENERGY] >= 0)

	### Expressions ###
	@expression(EP, eExistingCapEnergy[y in STOR_ALL], inputs["dfGen"][y,:Existing_Cap_MWh]*EP[:vZERO])

	@expression(EP, eTotalCapEnergy[y in STOR_ALL],
	if (y in intersect(NEW_CAP_ENERGY, RET_CAP_ENERGY))
		eExistingCapEnergy[y] + inputs["dfGen"][y,:Cap_Size]*sum(EP[:vCAPENERGY][s,y] - EP[:vRETCAPENERGY][s,y] for s in 1:cur_stage)
	elseif (y in setdiff(NEW_CAP_ENERGY, RET_CAP_ENERGY))
		eExistingCapEnergy[y] + inputs["dfGen"][y,:Cap_Size]*sum(EP[:vCAPENERGY][s,y] for s in 1:cur_stage)
	elseif (y in setdiff(RET_CAP_ENERGY, NEW_CAP_ENERGY))
		eExistingCapEnergy[y] - inputs["dfGen"][y,:Cap_Size]*sum(EP[:vRETCAPENERGY][s,y] for s in 1:cur_stage)
	else
		eExistingCapEnergy[y]
	end
	)

end

function load_built_transmission_cap!(EP::Model,setup::Dict,inputs::Dict)

	L = inputs["L"]     # Number of transmission lines
	NetworkExpansion = setup["NetworkExpansion"]
	if NetworkExpansion == 1
		# Network lines and zones that are expandable have non-negative maximum reinforcement inputs
		EXPANSION_LINES = inputs["EXPANSION_LINES"]
	end
	#num_stages = setup["MultiStageSettingsDict"]["NumStages"]
	cur_stage = inputs["CurStage"];
	if NetworkExpansion == 1
		# Transmission network capacity reinforcements per line
		@variable(EP, vNEW_TRANS_CAP[t in 1:cur_stage, l in EXPANSION_LINES] >= 0)
	end

	@expression(EP, eExistingTransMax[l=1:L], inputs["pTrans_Max"][l]*EP[:vZERO])
	
	## Transmission power flow and loss related expressions:
	# Total availabile maximum transmission capacity is the sum of existing maximum transmission capacity plus new transmission capacity
	if NetworkExpansion == 1
		@expression(EP, eAvail_Trans_Cap[l=1:L],
			if l in EXPANSION_LINES
				eExistingTransMax[l] + inputs["pLine_Size"][l]*sum(vNEW_TRANS_CAP[s,l] for s in 1:cur_stage)
			else
				eExistingTransMax[l]
			end
		)
	else
		@expression(EP, eAvail_Trans_Cap[l=1:L], eExistingTransMax[l])
	end


end

function load_built_discharge_cap!(EP::Model,setup::Dict,inputs::Dict)

	#num_stages = setup["MultiStageSettingsDict"]["NumStages"];
	
	cur_stage = inputs["CurStage"];
	
	G = inputs["G"] # Number of resources (generators, storage, DR, and DERs)

	NEW_CAP = inputs["NEW_CAP"] # Set of all resources eligible for new capacity
	RET_CAP = inputs["RET_CAP"] # Set of all resources eligible for capacity retirements
	
	# Retired capacity of resource "y" from existing capacity
	@variable(EP, vRETCAP[t in 1:cur_stage, y in RET_CAP] >= 0);

	# New installed capacity of resource "y"
	@variable(EP, vCAP[t in 1:cur_stage, y in NEW_CAP] >= 0);

	@expression(EP, eExistingCap[y in 1:G], inputs["dfGen"][y,:Existing_Cap_MW]*EP[:vZERO])

	@expression(EP, eTotalCap[y in 1:G],
	if y in intersect(NEW_CAP, RET_CAP) # Resources eligible for new capacity and retirements
		eExistingCap[y] + inputs["dfGen"][y,:Cap_Size]*sum(EP[:vCAP][s,y] - EP[:vRETCAP][s,y] for s in 1:cur_stage)
	elseif y in setdiff(NEW_CAP, RET_CAP) # Resources eligible for only new capacity
		eExistingCap[y] + inputs["dfGen"][y,:Cap_Size]*sum(EP[:vCAP][s,y] for s in 1:cur_stage)
	elseif y in setdiff(RET_CAP, NEW_CAP) # Resources eligible for only capacity retirements
		eExistingCap[y] - inputs["dfGen"][y,:Cap_Size]*sum(EP[:vRETCAP][s,y] for s in 1:cur_stage)
	else # Resources not eligible for new capacity or retirements
		eExistingCap[y]
	end
	)

end


function co2_cap_master!(EP::Model, inputs::Dict, setup::Dict)

	@variable(EP,vCO2budget[w=1:inputs["REP_PERIOD"],cap=1:inputs["NCO2Cap"]])

	if setup["CO2Cap"] == 1
		if inputs["NumScenarios"]<=1
			@constraint(EP, cCO2Emissions_systemwide_master[cap=1:inputs["NCO2Cap"]], sum(vCO2budget[w,cap] for w in 1:inputs["REP_PERIOD"]) == sum(inputs["dfMaxCO2"][z,cap] for z=findall(x->x==1, inputs["dfCO2CapZones"][:,cap])))
		else
			scenario_periods = subperiods2scenarios(inputs["REP_PERIOD"], inputs["NumScenarios"])
			@constraint(EP,cCO2Emissions_systemwide_master[cap=1:inputs["NCO2Cap"],s=1:inputs["NumScenarios"]], sum(vCO2budget[w,cap] for w in scenario_periods[s,1]:scenario_periods[s,2]) == sum(inputs["dfMaxCO2"][z,cap] for z=findall(x->x==1, inputs["dfCO2CapZones"][:,cap])))
		end
	elseif setup["CO2Cap"] == 2
		@constraint(EP, cCO2Emissions_systemwide_master[cap=1:inputs["NCO2Cap"]],
			sum(vCO2budget[w,cap] for w in 1:inputs["REP_PERIOD"]) == sum(inputs["dfMaxCO2Rate"][z,cap] * sum(inputs["omega"][t]*inputs["pD"][t,z] for t=1:inputs["T"]) for z = findall(x->x==1, inputs["dfCO2CapZones"][:,cap])))
	elseif setup["CO2Cap"]==3
		@constraint(EP, cCO2Emissions_systemwide_master[cap=1:inputs["NCO2Cap"]],sum(vCO2budget[w,cap] for w in 1:inputs["REP_PERIOD"])==0)
	end
end

function energy_share_requirement_master!(EP::Model, inputs::Dict, setup::Dict)
	T = inputs["T"]     # Number of time steps (hours)

	W = inputs["REP_PERIOD"]

	weighted_demand = [sum(inputs["omega"][t]*inputs["pD"][t,z] for t in 1:T) for z in 1:inputs["Z"]];

	@variable(EP,vESRbudget[w=1:W, ESR=1:inputs["nESR"]])

	@constraint(EP,cESRShare_master[ESR=1:inputs["nESR"]],sum(vESRbudget[w,ESR] for w in 1:inputs["REP_PERIOD"]) ==0)
	# @constraint(EP, cESRShare_master[ESR=1:inputs["nESR"]], sum(vESRbudget[ESR,w] for w in 1:inputs["REP_PERIOD"]) == sum(inputs["dfESR"][z,ESR]*weighted_demand[z] for z=findall(x->x>0,inputs["dfESR"][:,ESR])))

end




function long_duration_storage_master!(EP::Model,inputs::Dict,setup::Dict)

	REP_PERIOD = inputs["REP_PERIOD"]     # Number of representative periods

	STOR_LONG_DURATION = inputs["STOR_LONG_DURATION"]

	dfPeriodMap = inputs["Period_Map"] # Dataframe that maps modeled periods to representative periods
	NPeriods = size(inputs["Period_Map"])[1] # Number of modeled periods

	MODELED_PERIODS_INDEX = 1:NPeriods
	### Variables ###

	# Variables to define inter-period energy transferred between modeled periods

	# State of charge of storage at beginning of each modeled period n
	@variable(EP, vSOCw[y in STOR_LONG_DURATION, n in MODELED_PERIODS_INDEX] >= 0)

	# Build up in storage inventory over each representative period w
	# Build up inventory can be positive or negative
	@variable(EP, vdSOC[y in STOR_LONG_DURATION, w=1:REP_PERIOD])

	### Constraints ###

	# Storage at beginning of period w = storage at beginning of period w-1 + storage built up in period w (after n representative periods)
	## Multiply storage build up term from prior period with corresponding weight
	# @constraint(EP, cSoCBalLongDurationStorage[y in STOR_LONG_DURATION, r in MODELED_PERIODS_INDEX],
	# 				vSOCw[y, mod1(r+1, NPeriods)] == vSOCw[y,r] + vdSOC[y,dfPeriodMap[r,:Rep_Period_Index]])

	if inputs["NumScenarios"]<=1
		@constraint(EP, cSoCBalLongDurationStorage[y in STOR_LONG_DURATION, r in MODELED_PERIODS_INDEX],
		vSOCw[y, mod1(r+1, NPeriods)] == vSOCw[y,r] + vdSOC[y,dfPeriodMap[r,:Rep_Period_Index]])
	else
		### Note, function circular_index is defined in utility.jl
		@constraint(EP, cSoCBalLongDurationStorage[y in STOR_LONG_DURATION, r in MODELED_PERIODS_INDEX],
		vSOCw[y, circular_index(r,MODELED_PERIODS_INDEX,inputs["NumScenarios"],-1)] == vSOCw[y,r] + vdSOC[y,dfPeriodMap[r,:Rep_Period_Index]])
	end

	# Storage at beginning of each modeled period cannot exceed installed energy capacity
	@constraint(EP, cSoCBalLongDurationStorageUpper[y in STOR_LONG_DURATION, r in MODELED_PERIODS_INDEX],
	vSOCw[y,r] <= EP[:eTotalCapEnergy][y])

end


function hydro_inter_period_linkage_master!(EP::Model,inputs::Dict,setup::Dict)

	dfGen = inputs["dfGen"]

	REP_PERIOD = inputs["REP_PERIOD"]     # Number of representative periods

	STOR_HYDRO_LONG_DURATION = inputs["STOR_HYDRO_LONG_DURATION"]


	dfPeriodMap = inputs["Period_Map"] # Dataframe that maps modeled periods to representative periods
	NPeriods = size(inputs["Period_Map"])[1] # Number of modeled periods

	MODELED_PERIODS_INDEX = 1:NPeriods

	### Variables ###

	# Variables to define inter-period energy transferred between modeled periods

	# State of charge of storage at beginning of each modeled period n
	@variable(EP, vSOC_HYDROw[y in STOR_HYDRO_LONG_DURATION, n in MODELED_PERIODS_INDEX] >= 0)

	# Build up in storage inventory over each representative period w
	# Build up inventory can be positive or negative
	@variable(EP, vdSOC_HYDRO[y in STOR_HYDRO_LONG_DURATION, w=1:REP_PERIOD])

	### Constraints ###

	# Storage at beginning of period w = storage at beginning of period w-1 + storage built up in period w (after n representative periods)
	## Multiply storage build up term from prior period with corresponding weight
	# @constraint(EP, cHydroReservoirLongDurationStorage[y in STOR_HYDRO_LONG_DURATION, r in MODELED_PERIODS_INDEX],
	# 				vSOC_HYDROw[y, mod1(r+1, NPeriods)] == vSOC_HYDROw[y,r] + vdSOC_HYDRO[y,dfPeriodMap[r,:Rep_Period_Index]])
	if inputs["NumScenarios"]<=1
		@constraint(EP, cHydroReservoirLongDurationStorage[y in STOR_HYDRO_LONG_DURATION, r in MODELED_PERIODS_INDEX],
					vSOC_HYDROw[y, mod1(r+1, NPeriods)] == vSOC_HYDROw[y,r] + vdSOC_HYDRO[y,dfPeriodMap[r,:Rep_Period_Index]])
	else
		### Note, function circular_index is defined in utility.jl
		@constraint(EP, cHydroReservoirLongDurationStorage[y in STOR_HYDRO_LONG_DURATION, r in MODELED_PERIODS_INDEX],
					vSOC_HYDROw[y, circular_index(r,MODELED_PERIODS_INDEX,inputs["NumScenarios"],-1)] == vSOC_HYDROw[y,r] + vdSOC_HYDRO[y,dfPeriodMap[r,:Rep_Period_Index]])
	end

	# Storage at beginning of each modeled period cannot exceed installed energy capacity
	@constraint(EP, cHydroReservoirLongDurationStorageUpper[y in STOR_HYDRO_LONG_DURATION, r in MODELED_PERIODS_INDEX],
					vSOC_HYDROw[y,r] <= dfGen[y,:Hydro_Energy_to_Power_Ratio]*EP[:eTotalCap][y])

end

function investment_transmission!(EP::Model, inputs::Dict, setup::Dict)

	L = inputs["L"]     # Number of transmission lines

	NetworkExpansion = setup["NetworkExpansion"]

	if NetworkExpansion == 1
		# Network lines and zones that are expandable have non-negative maximum reinforcement inputs
		EXPANSION_LINES = inputs["EXPANSION_LINES"]
	end

	### Variables ###

	if NetworkExpansion == 1
		# Transmission network capacity reinforcements per line
		@variable(EP, vNEW_TRANS_CAP[l in EXPANSION_LINES] >= 0)
	end

	### Expressions ###

	@expression(EP, eTransMax[l=1:L], inputs["pTrans_Max"][l]*EP[:vZERO])
	
	## Transmission power flow and loss related expressions:
	# Total availabile maximum transmission capacity is the sum of existing maximum transmission capacity plus new transmission capacity
	if NetworkExpansion == 1
		@expression(EP, eAvail_Trans_Cap[l=1:L],
			if l in EXPANSION_LINES
				eTransMax[l] + vNEW_TRANS_CAP[l]*inputs["pLine_Size"][l] 
			else
				eTransMax[l]
			end
		)
	else
		@expression(EP, eAvail_Trans_Cap[l=1:L], eTransMax[l])
	end

	## Objective Function Expressions ##

	if NetworkExpansion == 1
		@expression(EP, eTotalCNetworkExp, sum(vNEW_TRANS_CAP[l]*inputs["pLine_Size"][l]*inputs["pC_Line_Reinforcement"][l] for l in EXPANSION_LINES))
		EP[:eObj] += eTotalCNetworkExp
	end

	## End Objective Function Expressions ##

	### Constraints ###

  	## Power flow and transmission (between zone) loss related constraints

	# If network expansion is used:
	if NetworkExpansion == 1
		# Constrain maximum single-stage line capacity reinforcement for lines eligible for expansion
		@constraint(EP, cMaxLineReinforcement[l in EXPANSION_LINES], vNEW_TRANS_CAP[l]*inputs["pLine_Size"][l] <= inputs["pMax_Line_Reinforcement"][l])
		if haskey(inputs, "dfMaxPathTrans")
			# Constrain maximum single-stage line capacity reinforcement for transmission paths
			dfMaxPathTrans = inputs["dfMaxPathTrans"];
			constrained_paths = dfMaxPathTrans[!,:transmission_path_name];
			@constraint(EP,cMaxPathReinforcement[i in 1:length(constrained_paths)], sum(vNEW_TRANS_CAP[l]*inputs["pLine_Size"][l] for l in intersect(EXPANSION_LINES,findall(inputs["TransmissionPath"].==constrained_paths[i]))) <= dfMaxPathTrans[i,:Path_Max_Reinforcement])
		end
	end
	#END network expansion contraints


end



function multistage_endogenous_retirement_charge!(EP::Model, inputs::Dict, setup::Dict)

	multi_stage_settings = setup["MultiStageSettingsDict"]
	num_stages = num_stages = multi_stage_settings["NumStages"];
	stage_lens = multi_stage_settings["StageLengths"]
	STOR_ASYMMETRIC = inputs[1]["STOR_ASYMMETRIC"] # Set of storage resources with asymmetric (separte) charge/discharge capacity components

	NEW_CAP_CHARGE = inputs[1]["NEW_CAP_CHARGE"] # Set of asymmetric charge/discharge storage resources eligible for new charge capacity
	RET_CAP_CHARGE = inputs[1]["RET_CAP_CHARGE"] # Set of asymmetric charge/discharge storage resources eligible for charge capacity retirements


	println("Endogenous Retirement (Charge) Module")

	### Expressions ###

	@expression(EP, eNewCapCharge[t in 1:num_stages, y in RET_CAP_CHARGE],
		if y in NEW_CAP_CHARGE
			EP[:vCAPCHARGE][t,y]
		else
			0*EP[:vZERO]
		end
	)

	# Construct and add the endogenous retirement constraint expressions
	@expression(EP, eRetCapTrackCharge[t in 1:num_stages, y in RET_CAP_CHARGE], sum(EP[:vRETCAPCHARGE][s,y] for s=1:t))

	@expression(EP, eNewCapTrackCharge[t in 1:num_stages, y in RET_CAP_CHARGE], sum(eNewCapCharge[s,y] for s=1:get_retirement_stage(t, inputs[t]["dfGen"][y,:Lifetime], stage_lens)))

	@expression(EP, eMinRetCapTrackCharge[t in 1:num_stages, y in RET_CAP_CHARGE], sum((inputs[s]["dfGen"][y,:Min_Retired_Charge_Cap_MW]/inputs[s]["dfGen"][y,:Cap_Size]) for s=1:t))

	### Constratints ###

	@constraint(EP, cLifetimeRetCharge[t in 1:num_stages, y in RET_CAP_CHARGE], eNewCapTrackCharge[t,y] + eMinRetCapTrackCharge[t,y]  <= eRetCapTrackCharge[t,y])

end


function multistage_endogenous_retirement_energy!(EP::Model, inputs::Dict, setup::Dict)

	multi_stage_settings = setup["MultiStageSettingsDict"]
	num_stages = num_stages = multi_stage_settings["NumStages"];
	stage_lens = multi_stage_settings["StageLengths"]

	STOR_ALL = inputs[1]["STOR_ALL"] # Set of all storage resources
	NEW_CAP_ENERGY = inputs[1]["NEW_CAP_ENERGY"] # Set of all storage resources eligible for new energy capacity
	RET_CAP_ENERGY = inputs[1]["RET_CAP_ENERGY"] # Set of all storage resources eligible for energy capacity 

	println("Endogenous Retirement (Energy) Module")

	### Expressions ###

	@expression(EP, eNewCapEnergy[t in 1:num_stages, y in RET_CAP_ENERGY],
		if y in NEW_CAP_ENERGY
			EP[:vCAPENERGY][t,y]
		else
			0*EP[:vZERO]
		end
	)

	# Construct and add the endogenous retirement constraint expressions
	@expression(EP, eRetCapTrackEnergy[t in 1:num_stages, y in RET_CAP_ENERGY], sum(EP[:vRETCAPENERGY][s,y] for s=1:t))

	@expression(EP, eNewCapTrackEnergy[t in 1:num_stages, y in RET_CAP_ENERGY], sum(eNewCapEnergy[s,y] for s=1:get_retirement_stage(t, inputs[t]["dfGen"][y,:Lifetime], stage_lens)))

	@expression(EP, eMinRetCapTrackEnergy[t in 1:num_stages, y in RET_CAP_ENERGY], sum((inputs[s]["dfGen"][y,:Min_Retired_Energy_Cap_MW]/inputs[s]["dfGen"][y,:Cap_Size]) for s=1:t))

	### Constratints ###

	@constraint(EP, cLifetimeRetEnergy[t in 1:num_stages, y in RET_CAP_ENERGY], eNewCapTrackEnergy[t,y] + eMinRetCapTrackEnergy[t,y]  <= eRetCapTrackEnergy[t,y])

end


function multistage_endogenous_retirement_discharge!(EP::Model, inputs::Dict, setup::Dict)

	println("Endogenous Retirement (Discharge) Module")

	NEW_CAP = inputs[1]["NEW_CAP"] # Set of all resources eligible for new capacity
	RET_CAP = inputs[1]["RET_CAP"] # Set of all resources eligible for capacity retirements

	multi_stage_settings = setup["MultiStageSettingsDict"]
	num_stages = num_stages = multi_stage_settings["NumStages"];
	stage_lens = multi_stage_settings["StageLengths"]

	### Expressions ###

	@expression(EP, eNewCap[t in 1:num_stages, y in RET_CAP],
		if y in NEW_CAP
			EP[:vCAP][t,y]
		else
			0*EP[:vZERO]
		end
	)

	# Construct and add the endogenous retirement constraint expressions
	@expression(EP, eRetCapTrack[t in 1:num_stages, y in RET_CAP], sum(EP[:vRETCAP][s,y] for s=1:t))

	@expression(EP, eNewCapTrack[t in 1:num_stages, y in RET_CAP], sum(eNewCap[s,y] for s=1:get_retirement_stage(t, inputs[t]["dfGen"][y,:Lifetime], stage_lens)))

	@expression(EP, eMinRetCapTrack[t in 1:num_stages, y in RET_CAP],
			sum((inputs[s]["dfGen"][y,:Min_Retired_Cap_MW]/inputs[s]["dfGen"][y,:Cap_Size]) for s=1:t)
	)

	### Constraints ###

	@constraint(EP, cLifetimeRet[t in 1:num_stages, y in RET_CAP], eNewCapTrack[t,y] + eMinRetCapTrack[t,y]  <= eRetCapTrack[t,y])

end



function multistage_energy_share_requirement_master!(EP::Model, inputs::Dict, setup::Dict)

	println("Multistage ESR policy module")

	num_stages = setup["MultiStageSettingsDict"]["NumStages"];

	nW,stage_subperiods,nco2cap,nesr = map_from_stage_to_subperiod(inputs,num_stages)

	@variable(EP,vESRbudget[w in 1:nW, ESR=1:nesr[w]])

	@constraint(EP,cESRShare_master[t in 1:num_stages, ESR=1:inputs[t]["nESR"]],sum(vESRbudget[w,ESR] for w in stage_subperiods[t]) ==0)

end

function multistage_co2_cap_master!(EP::Model,inputs::Dict,setup::Dict)

	println("Multistage CO2 Cap module")

	num_stages = setup["MultiStageSettingsDict"]["NumStages"];

	nW,stage_subperiods,nco2cap,nesr = map_from_stage_to_subperiod(inputs,num_stages)

	@variable(EP,vCO2budget[w=1:nW, cap=1:nco2cap[w]])

	if setup["CO2Cap"] == 1
		if setup["NumScenarios"]<=1

			@constraint(EP, cCO2Emissions_systemwide_master[t in 1:num_stages, cap=1:inputs[t]["NCO2Cap"]], sum(vCO2budget[w,cap] for w in stage_subperiods[t]) == sum(inputs[t]["dfMaxCO2"][z,cap] for z=findall(x->x==1, inputs[t]["dfCO2CapZones"][:,cap])))

		else

			@constraint(EP,cCO2Emissions_systemwide_master[t in 1:num_stages, cap=1:inputs[t]["NCO2Cap"], s=1:inputs[t]["NumScenarios"]], sum(vCO2budget[w,cap] for p in intersect(stage_subperiods[t],subperiods2scenarios(nW, inputs[t]["NumScenarios"])[s,1]:subperiods2scenarios(nW, inputs[t]["NumScenarios"])[s,2])) == sum(inputs[t]["dfMaxCO2"][z,cap] for z=findall(x->x==1, inputs[t]["dfCO2CapZones"][:,cap])))

		end
	elseif setup["CO2Cap"] == 2

			@constraint(EP, cCO2Emissions_systemwide_master[t in 1:num_stages, cap=1:inputs[t]["NCO2Cap"]],
				sum(vCO2budget[w,cap] for w in stage_subperiods[t]) == sum(inputs[t]["dfMaxCO2Rate"][z,cap] * sum(inputs[t]["omega"][h]*inputs[t]["pD"][h,z] for h=1:inputs[t]["T"]) for z = findall(x->x==1, inputs[t]["dfCO2CapZones"][:,cap])))

	elseif setup["CO2Cap"] == 3

		@constraint(EP, cCO2Emissions_systemwide_master[t in 1:num_stages, cap=1:inputs[t]["NCO2Cap"]],sum(vCO2budget[w,cap] for w in stage_subperiods[t])==0)

	end

end

function multistage_maximum_capacity_requirement!(EP::Model,inputs::Dict,setup::Dict)

	println("Maximum Capacity Requirement Module")
	num_stages = setup["MultiStageSettingsDict"]["NumStages"]
	# if input files are present, add maximum capacity requirement slack variables
	if haskey(inputs[1], "MaxCapPriceCap")
		@variable(EP, vMaxCap_slack[t in 1:num_stages, maxcap = 1:inputs[t]["NumberOfMaxCapReqs"]]>=0)
		EP[:eMaxCapRes] -= vMaxCap_slack

		@expression(EP, eCMaxCap_slack[t in 1:num_stages, maxcap = 1:inputs[t]["NumberOfMaxCapReqs"]], inputs[t]["MaxCapPriceCap"][maxcap] * EP[:vMaxCap_slack][t,maxcap])
		@expression(EP, eTotalCMaxCapSlack[t in 1:num_stages], sum(EP[:eCMaxCap_slack][t,maxcap] for maxcap = 1:inputs[t]["NumberOfMaxCapReqs"]))
		
	end
	
	@constraint(EP, cZoneMaxCapReq[t in 1:num_stages, maxcap = 1:inputs[t]["NumberOfMaxCapReqs"]], EP[:eMaxCapRes][t,maxcap] >= inputs[t]["MaxCapReq"][maxcap])

end


function multistage_minimum_capacity_requirement!(EP::Model,inputs::Dict,setup::Dict)

	println("Minimum Capacity Requirement Module")
	num_stages = setup["MultiStageSettingsDict"]["NumStages"]

	# if input files are present, add minimum capacity requirement slack variables
	if haskey(inputs[1], "MinCapPriceCap")
		@variable(EP, vMinCap_slack[t in 1:num_stages, mincap = 1:inputs[t]["NumberOfMinCapReqs"]]>=0)
		EP[:eMinCapRes] += vMinCap_slack
	
		@expression(EP, eCMinCap_slack[t in 1:num_stages, mincap = 1:inputs[t]["NumberOfMinCapReqs"]], inputs[t]["MinCapPriceCap"][mincap] * EP[:vMinCap_slack][t,mincap])
		@expression(EP, eTotalCMinCapSlack[t in 1:num_stages], sum(EP[:eCMinCap_slack][t,mincap] for mincap = 1:inputs[t]["NumberOfMinCapReqs"]) )
			
	end
	
	@constraint(EP, cZoneMinCapReq[t in 1:num_stages, mincap = 1:inputs[t]["NumberOfMinCapReqs"]], EP[:eMinCapRes][t,mincap] >= inputs[t]["MinCapReq"][mincap])

end

function multistage_retrofit!(EP::Model, inputs::Dict, setup::Dict)

	println("Multistage Retrofit Resources Module")

	num_stages = setup["MultiStageSettingsDict"]["NumStages"]
	G = inputs["G"]   # Number of resources (generators, storage, DR, and DERs)
	RESOURCES = inputs[1]["RESOURCES"] # Set of all resources by name
	RETRO = inputs[1]["RETRO"] # Set of all retrofit resources by ID
	NEW_CAP = inputs[1]["NEW_CAP"] # Set of all resources eligible for capacity expansion by ID
	RET_CAP = inputs[1]["RET_CAP"] # Set of all resources eligible for capacity retirements by ID
	RETRO_SOURCES = inputs[1]["RETROFIT_SOURCES"] # Source technologies by name for each retrofit [1:G]
	RETRO_SOURCE_IDS = inputs[1]["RETROFIT_SOURCE_IDS"] # Source technologies by ID for each retrofit [1:G]
	RETRO_EFFICIENCY = inputs[1]["RETROFIT_EFFICIENCIES"] # Ratio of installed retrofit capacity to source capacity [0:1] (indexed by retrofit tech r, source # i)
	NUM_RETRO_SOURCES = inputs[1]["NUM_RETROFIT_SOURCES"] # Number of possible sources for a given retrofit resource
	### Expressions ###

	# Retired capacity of all retirement-eligible resources
	@expression(EP, eRetroRetireCap[t in 1:num_stages, y in RET_CAP],
		EP[:vRETCAP][t,y]*inputs[t]["dfGen"][y, :Cap_Size]
	)

	# One-to-Many Retrofit Mapping: Sum of capacity being retrofitted from a resource to all of its possible destination retrofit technologies
	@expression(EP, eRetroRetireCapMap[t in 1:num_stages, y in RET_CAP],
		sum( EP[:vRETROFIT][t,y,r]*inputs[t]["dfGen"][y, :Cap_Size] for r in intersect( findall(x->in(RESOURCES[y],RETRO_SOURCES[x]),1:G), NEW_CAP ); init=0 )
	)

	# Many-to-One Retrofit Mapping: For a given retrofit technology, sum of retrofit capacity from all of its possible sources
	@expression(EP, eRetroInstallCapMap[t in 1:num_stages, r in intersect(RETRO, NEW_CAP)],
		sum( RETRO_SOURCE_IDS[r][i] in RET_CAP ? EP[:vRETROFIT][t,RETRO_SOURCE_IDS[r][i], r]*inputs[t]["dfGen"][RETRO_SOURCE_IDS[r][i], :Cap_Size]*RETRO_EFFICIENCY[r][i] : 0 for i in 1:NUM_RETRO_SOURCES[r]; init=0 )
	)

	
	# Installed capacity of all retrofit resources
	@expression(EP, eRetroInstallCap[t in 1:num_stages, r in intersect(RETRO, NEW_CAP)],
		EP[:vCAP][t,r]*inputs[t]["dfGen"][r, :Cap_Size]
	)
	
	### Constraints ###

	# (One-to-Many) Sum of retrofitted capacity from a given source technology must not exceed the retired capacity of that technology. (Retrofitting is included within retirement, not a distinct category)
	# TO DO: Add term for decommissioned capacity on RHS and make it an equality constraint
	@constraint(EP, cRetroSource[t in 1:num_stages, y in RET_CAP], eRetroRetireCap[t,y] >= eRetroRetireCapMap[t,y])

	# (Many-to-One) New installed capacity of retrofit technology r must be equal to the (efficiency-downscaled) sum of capacity retrofitted to technology r from source technologies yr
	@constraint(EP, cRetroDest[t in 1:num_stages, r in intersect(RETRO, NEW_CAP)], eRetroInstallCapMap[t,r] == eRetroInstallCap[t,r])

	return EP

end



function multistage_investment_charge!(EP::Model, inputs::Dict, setup::Dict)

	println("Multistage Charge Investment Module")
	num_stages = setup["MultiStageSettingsDict"]["NumStages"]

	STOR_ASYMMETRIC = inputs[1]["STOR_ASYMMETRIC"] # Set of storage resources with asymmetric (separte) charge/discharge capacity components

	NEW_CAP_CHARGE = inputs[1]["NEW_CAP_CHARGE"] # Set of asymmetric charge/discharge storage resources eligible for new charge capacity
	RET_CAP_CHARGE = inputs[1]["RET_CAP_CHARGE"] # Set of asymmetric charge/discharge storage resources eligible for charge capacity retirements

	### Variables ###

	## Storage capacity built and retired for storage resources with independent charge and discharge power capacities (STOR=2)

	# New installed charge capacity of resource "y"
	@variable(EP, vCAPCHARGE[t in 1:num_stages, y in NEW_CAP_CHARGE] >= 0)

	# Retired charge capacity of resource "y" from existing capacity
	@variable(EP, vRETCAPCHARGE[t in 1:num_stages, y in RET_CAP_CHARGE] >= 0)

	### Expressions ###
	@expression(EP, eExistingCapCharge[y in STOR_ASYMMETRIC], inputs[1]["dfGen"][y,:Existing_Charge_Cap_MW]*EP[:vZERO])

	@expression(EP, eTotalCapCharge[t in 1:num_stages, y in STOR_ASYMMETRIC],
	if (y in intersect(NEW_CAP_CHARGE, RET_CAP_CHARGE))
		eExistingCapCharge[y] + inputs[t]["dfGen"][y,:Cap_Size]*sum(EP[:vCAPCHARGE][s,y] - EP[:vRETCAPCHARGE][s,y] for s in 1:t)
	elseif (y in setdiff(NEW_CAP_CHARGE, RET_CAP_CHARGE))
		eExistingCapCharge[y] + inputs[t]["dfGen"][y,:Cap_Size]*sum(EP[:vCAPCHARGE][s,y] for s in 1:t)
	elseif (y in setdiff(RET_CAP_CHARGE, NEW_CAP_CHARGE))
		eExistingCapCharge[y] - inputs[t]["dfGen"][y,:Cap_Size]*sum(EP[:vRETCAPCHARGE][s,y] for s in 1:t)
	else
		eExistingCapCharge[y]
	end
	)
	#This expression is needed to bound the capacity that can be retired in each stage
	@expression(EP, eMaxRetCapCharge[t in 1:num_stages, y in RET_CAP_CHARGE], 
	if t== 1
		eExistingCapCharge[y]
	else
		eTotalCapCharge[t-1,y]
	end
	)

	## Objective Function Expressions ##

	# Fixed costs for resource "y" = annuitized investment cost plus fixed O&M costs
	# If resource is not eligible for new charge capacity, fixed costs are only O&M costs
	@expression(EP, eCFixCharge[t in 1:num_stages, y in STOR_ASYMMETRIC],
		if y in NEW_CAP_CHARGE # Resources eligible for new charge capacity
			inputs[t]["dfGen"][y,:Inv_Cost_Charge_per_MWyr]*inputs[t]["dfGen"][y,:Cap_Size]*vCAPCHARGE[t,y] + inputs[t]["dfGen"][y,:Fixed_OM_Cost_Charge_per_MWyr]*eTotalCapCharge[t,y]
		else
			inputs[t]["dfGen"][y,:Fixed_OM_Cost_Charge_per_MWyr]*eTotalCapCharge[t,y]
		end
	)

	# Sum individual resource contributions to fixed costs to get total fixed costs
	@expression(EP, eTotalCFixCharge[t in 1:num_stages], sum(EP[:eCFixCharge][t,y] for y in STOR_ASYMMETRIC))

	### Constratints ###

	## Constraints on retirements and capacity additions
	#Cannot retire more charge capacity than existing charge capacity
	# @constraint(EP,cAuxPositiveTotChargeCap[t in 1:num_stages,y in STOR_ASYMMETRIC], eTotalCapCharge[t,y]>=0)
	@constraint(EP, cMaxRetCharge[t in 1:num_stages, y in RET_CAP_CHARGE], inputs[t]["dfGen"][y,:Cap_Size]*vRETCAPCHARGE[t,y] <= eMaxRetCapCharge[t,y])

  	#Constraints on new built capacity

	# Constraint on maximum charge capacity (if applicable) [set input to -1 if no constraint on maximum charge capacity]
	# DEV NOTE: This constraint may be violated in some cases where Existing_Charge_Cap_MW is >= Max_Charge_Cap_MWh and lead to infeasabilty
    @constraint(EP, cMaxCapCharge[t in 1:num_stages, y in intersect(inputs[t]["dfGen"][inputs[t]["dfGen"].Max_Charge_Cap_MW.>0,:R_ID], STOR_ASYMMETRIC)], eTotalCapCharge[t,y] <= inputs[t]["dfGen"][y,:Max_Charge_Cap_MW])

	# Constraint on minimum charge capacity (if applicable) [set input to -1 if no constraint on minimum charge capacity]
	# DEV NOTE: This constraint may be violated in some cases where Existing_Charge_Cap_MW is <= Min_Charge_Cap_MWh and lead to infeasabilty
    @constraint(EP, cMinCapCharge[t in 1:num_stages, y in intersect(inputs[t]["dfGen"][inputs[t]["dfGen"].Min_Charge_Cap_MW.>0,:R_ID], STOR_ASYMMETRIC)], eTotalCapCharge[t,y] >= inputs[t]["dfGen"][y,:Min_Charge_Cap_MW])


end


function multistage_investment_energy!(EP::Model, inputs::Dict, setup::Dict)

	println("Multistage Storage Investment Module")

	num_stages = setup["MultiStageSettingsDict"]["NumStages"]

	STOR_ALL = inputs[1]["STOR_ALL"] # Set of all storage resources
	NEW_CAP_ENERGY = inputs[1]["NEW_CAP_ENERGY"] # Set of all storage resources eligible for new energy capacity
	RET_CAP_ENERGY = inputs[1]["RET_CAP_ENERGY"] # Set of all storage resources eligible for energy capacity 

	### Variables ###

	## Energy storage reservoir capacity (MWh capacity) built/retired for storage with variable power to energy ratio (STOR=1 or STOR=2)

	# New installed energy capacity of resource "y"
	@variable(EP, vCAPENERGY[t in 1:num_stages, y in NEW_CAP_ENERGY] >= 0)

	# Retired energy capacity of resource "y" from existing capacity
	@variable(EP, vRETCAPENERGY[t in 1:num_stages, y in RET_CAP_ENERGY] >= 0)

	### Expressions ###
	@expression(EP, eExistingCapEnergy[y in STOR_ALL], inputs[1]["dfGen"][y,:Existing_Cap_MWh]*EP[:vZERO])

	@expression(EP, eTotalCapEnergy[t in 1:num_stages, y in STOR_ALL],
	if (y in intersect(NEW_CAP_ENERGY, RET_CAP_ENERGY))
		eExistingCapEnergy[y] + inputs[t]["dfGen"][y,:Cap_Size]*sum(EP[:vCAPENERGY][s,y] - EP[:vRETCAPENERGY][s,y] for s in 1:t)
	elseif (y in setdiff(NEW_CAP_ENERGY, RET_CAP_ENERGY))
		eExistingCapEnergy[y] + inputs[t]["dfGen"][y,:Cap_Size]*sum(EP[:vCAPENERGY][s,y] for s in 1:t)
	elseif (y in setdiff(RET_CAP_ENERGY, NEW_CAP_ENERGY))
		eExistingCapEnergy[y] - inputs[t]["dfGen"][y,:Cap_Size]*sum(EP[:vRETCAPENERGY][s,y] for s in 1:t)
	else
		eExistingCapEnergy[y]
	end
	)

	#This expression is needed to bound the capacity that can be retired in each stage
	@expression(EP, eMaxRetCapEnergy[t in 1:num_stages, y in RET_CAP_ENERGY], 
	if t== 1
		eExistingCapEnergy[y]
	else
		eTotalCapEnergy[t-1,y]
	end
	)

	## Objective Function Expressions ##

	# Fixed costs for resource "y" = annuitized investment cost plus fixed O&M costs
	# If resource is not eligible for new energy capacity, fixed costs are only O&M costs
	@expression(EP, eCFixEnergy[t in 1:num_stages, y in STOR_ALL],
		if y in NEW_CAP_ENERGY # Resources eligible for new capacity
			inputs[t]["dfGen"][y,:Inv_Cost_per_MWhyr]*inputs[t]["dfGen"][y,:Cap_Size]*vCAPENERGY[t,y] + inputs[t]["dfGen"][y,:Fixed_OM_Cost_per_MWhyr]*eTotalCapEnergy[t,y]
		else
			inputs[t]["dfGen"][y,:Fixed_OM_Cost_per_MWhyr]*eTotalCapEnergy[t,y]
		end
	)

	# Sum individual resource contributions to fixed costs to get total fixed costs
	@expression(EP, eTotalCFixEnergy[t in 1:num_stages], sum(EP[:eCFixEnergy][t,y] for y in STOR_ALL))

	### Constraints ###

	## Constraints on retirements and capacity additions
	# Cannot retire more energy capacity than existing energy capacity
	# @constraint(EP,cAuxPositiveTotEnergyCap[t in 1:num_stages,y in STOR_ALL], eTotalCapEnergy[t,y]>=0)
	@constraint(EP, cMaxRetEnergy[t in 1:num_stages, y in RET_CAP_ENERGY], inputs[t]["dfGen"][y,:Cap_Size]*vRETCAPENERGY[t,y] <= eMaxRetCapEnergy[t,y])

	## Constraints on new built energy capacity
	# Constraint on maximum energy capacity (if applicable) [set input to -1 if no constraint on maximum energy capacity]
	# DEV NOTE: This constraint may be violated in some cases where Existing_Cap_MWh is >= Max_Cap_MWh and lead to infeasabilty
	@constraint(EP, cMaxCapEnergy[t in 1:num_stages, y in intersect(inputs[t]["dfGen"][inputs[t]["dfGen"].Max_Cap_MWh.>0,:R_ID], STOR_ALL)], eTotalCapEnergy[t,y] <= inputs[t]["dfGen"][y,:Max_Cap_MWh])

	# Constraint on minimum energy capacity (if applicable) [set input to -1 if no constraint on minimum energy apacity]
	# DEV NOTE: This constraint may be violated in some cases where Existing_Cap_MWh is <= Min_Cap_MWh and lead to infeasabilty
	@constraint(EP, cMinCapEnergy[t in 1:num_stages, y in intersect(inputs[t]["dfGen"][inputs[t]["dfGen"].Min_Cap_MWh.>0,:R_ID], STOR_ALL)], eTotalCapEnergy[t,y] >= inputs[t]["dfGen"][y,:Min_Cap_MWh])

	# Max and min constraints on energy storage capacity built (as proportion to discharge power capacity)
	@constraint(EP, cMinCapEnergyDuration[t in 1:num_stages, y in STOR_ALL], eTotalCapEnergy[t,y] >= inputs[t]["dfGen"][y,:Min_Duration] * EP[:eTotalCap][t,y])
	@constraint(EP, cMaxCapEnergyDuration[t in 1:num_stages, y in STOR_ALL], eTotalCapEnergy[t,y] <= inputs[t]["dfGen"][y,:Max_Duration] * EP[:eTotalCap][t,y])

end


function multistage_investment_transmission!(EP::Model, inputs::Dict, setup::Dict)

	println("Multistage transmission investment module")

	L = inputs[1]["L"]     # Number of transmission lines

	NetworkExpansion = setup["NetworkExpansion"]

	if NetworkExpansion == 1
		# Network lines and zones that are expandable have non-negative maximum reinforcement inputs
		EXPANSION_LINES = inputs[1]["EXPANSION_LINES"]
	end

	num_stages = setup["MultiStageSettingsDict"]["NumStages"]

	### Variables ###

	if NetworkExpansion == 1
		# Transmission network capacity reinforcements per line
		@variable(EP, vNEW_TRANS_CAP[t in 1:num_stages, l in EXPANSION_LINES] >= 0)
	end

	### Expressions ###

	@expression(EP, eExistingTransMax[l=1:L], inputs[1]["pTrans_Max"][l]*EP[:vZERO])
	
	## Transmission power flow and loss related expressions:
	# Total availabile maximum transmission capacity is the sum of existing maximum transmission capacity plus new transmission capacity
	if NetworkExpansion == 1
		@expression(EP, eAvail_Trans_Cap[t in 1:num_stages, l=1:L],
			if l in EXPANSION_LINES
				eExistingTransMax[l] + sum(vNEW_TRANS_CAP[s,l]*inputs[t]["pLine_Size"][l] for s in 1:t)
			else
				eExistingTransMax[l]
			end
		)
	else
		@expression(EP, eAvail_Trans_Cap[t in 1:num_stages, l=1:L], eExistingTransMax[l])
	end


	## Objective Function Expressions ##

	if NetworkExpansion == 1
		@expression(EP, eTotalCNetworkExp[t in 1:num_stages], sum(vNEW_TRANS_CAP[t,l]*inputs[t]["pLine_Size"][l]*inputs[t]["pC_Line_Reinforcement"][l] for l in EXPANSION_LINES))
	end

	## End Objective Function Expressions ##

	### Constraints ###

  	#@constraint(EP,cAuxPositiveTotTransCap[t in 1:num_stages,l=1:L],eAvail_Trans_Cap[t,l]>=0)

	# If network expansion is used:
	if NetworkExpansion == 1
		# Constrain maximum single-stage line capacity reinforcement for lines eligible for expansion
		@constraint(EP, cMaxLineReinforcement[t in 1:num_stages, l in EXPANSION_LINES], vNEW_TRANS_CAP[t,l]*inputs[t]["pLine_Size"][l] <= inputs[t]["pMax_Line_Reinforcement"][l])
		@constraint(EP, cMaxFlowPossible[t in 1:num_stages, l in EXPANSION_LINES], eAvail_Trans_Cap[t,l] <= inputs[t]["pTrans_Max_Possible"][l])
		
		if haskey(inputs[1], "dfMaxPathTrans")
			# Constrain maximum single-stage line capacity reinforcement for transmission paths
			@constraint(EP,cMaxPathReinforcement[t in 1:num_stages, i in 1:size(inputs[t]["dfMaxPathTrans"],1)], sum(vNEW_TRANS_CAP[t,l]*inputs[t]["pLine_Size"][l] for l in intersect(EXPANSION_LINES,findall(inputs[t]["TransmissionPath"].==inputs[t]["dfMaxPathTrans"][i,:transmission_path_name]))) <= inputs[t]["dfMaxPathTrans"][i,:Path_Max_Reinforcement])

			# Constrain maximum possible flow for transmission paths
			@constraint(EP,cMaxPathFlowPossible[t in 1:num_stages, i in 1:size(inputs[t]["dfMaxPathTrans"],1)], sum(eAvail_Trans_Cap[t,l] for l in intersect(EXPANSION_LINES,findall(inputs[t]["TransmissionPath"].==inputs[t]["dfMaxPathTrans"][i,:transmission_path_name]))) <= inputs[t]["dfMaxPathTrans"][i,:Path_Max_Flow_Possible])
		end
	end
	
	#END network expansion contraints



end

function multistage_investment_discharge!(EP::Model, inputs::Dict, setup::Dict)

	num_stages = setup["MultiStageSettingsDict"]["NumStages"]

	println("Multistage Investment Discharge Module")

	G = inputs[1]["G"] # Number of resources (generators, storage, DR, and DERs)

	NEW_CAP = inputs[1]["NEW_CAP"] # Set of all resources eligible for new capacity
	RET_CAP = inputs[1]["RET_CAP"] # Set of all resources eligible for capacity retirements
	RETRO = inputs[1]["RETRO"]     # Set of all retrofit resources
	# Additional retrofit information if necessary
	if !isempty(RETRO)
		NUM_RETRO_SOURCES = inputs[1]["NUM_RETROFIT_SOURCES"]       # The number of source resources for each retrofit resource
		RETRO_SOURCES = inputs[1]["RETROFIT_SOURCES"]               # Source technologies (Resource Name) for each retrofit [1:G]
		RETRO_SOURCE_IDS = inputs[1]["RETROFIT_SOURCE_IDS"]         # Source technologies (IDs) for each retrofit [1:G]
		RETRO_INV_CAP_COSTS = inputs[1]["RETROFIT_INV_CAP_COSTS"]   # The set of investment costs (capacity $/MWyr) of each retrofit by source
		RETRO_EFFICIENCY = inputs[1]["RETROFIT_EFFICIENCIES"]       # Ratio of installed retrofit capacity to retired source capacity [0:1]
	end

	### Variables ###

	# Retired capacity of resource "y" from existing capacity
	@variable(EP, vRETCAP[t in 1:num_stages, y in RET_CAP] >= 0);

    # New installed capacity of resource "y"
	@variable(EP, vCAP[t in 1:num_stages, y in NEW_CAP] >= 0);

	# Capacity from source resource "yr" that is being retrofitted into capacity of retrofit resource "r"
	if !isempty(RETRO)
		# Dependent iterators only allowed in forward sequence, so we reconstruct retrofit destinations from sources.
		ALL_SOURCES = intersect(collect(Set(collect(Iterators.flatten(RETRO_SOURCE_IDS)))),RET_CAP)
		DESTS_BY_SOURCE = [ y in ALL_SOURCES ? intersect(findall(x->in(inputs["RESOURCES"][y],RETRO_SOURCES[x]), 1:G), findall(x->x in NEW_CAP, 1:G)) : []  for y in 1:G]
		@variable(EP, vRETROFIT[t in 1:num_stages,yr in ALL_SOURCES, r in DESTS_BY_SOURCE[yr]] >= 0);     # Capacity retrofitted from source technology y to retrofit technology r
	end

	@expression(EP, eExistingCap[y in 1:G], inputs[1]["dfGen"][y,:Existing_Cap_MW]*EP[:vZERO])

	@expression(EP, eTotalCap[t in 1:num_stages, y in 1:G],
	if y in intersect(NEW_CAP, RET_CAP) # Resources eligible for new capacity and retirements
		eExistingCap[y] + inputs[t]["dfGen"][y,:Cap_Size]*sum(EP[:vCAP][s,y] - EP[:vRETCAP][s,y] for s in 1:t)
	elseif y in setdiff(NEW_CAP, RET_CAP) # Resources eligible for only new capacity
			eExistingCap[y] + inputs[t]["dfGen"][y,:Cap_Size]*sum(EP[:vCAP][s,y] for s in 1:t)
	elseif y in setdiff(RET_CAP, NEW_CAP) # Resources eligible for only capacity retirements
			eExistingCap[y] - inputs[t]["dfGen"][y,:Cap_Size]*sum(EP[:vRETCAP][s,y] for s in 1:t)
	else # Resources not eligible for new capacity or retirements
		eExistingCap[y]
	end
	)

	#This expression is needed to bound the capacity that can be retired in each stage
	@expression(EP, eMaxRetCap[t in 1:num_stages, y in RET_CAP], 
	if t==1
		eExistingCap[y]
	else
		eTotalCap[t-1,y];
	end)

	## Objective Function Expressions ##

	# Fixed costs for resource "y" = annuitized investment cost plus fixed O&M costs
	# If resource is not eligible for new capacity, fixed costs are only O&M costs
	@expression(EP, eCFix[t in 1:num_stages,y in 1:G],
		if y in setdiff(NEW_CAP, RETRO) # Resources eligible for new capacity (Non-Retrofit)
				inputs[t]["dfGen"][y,:Inv_Cost_per_MWyr]*inputs[t]["dfGen"][y,:Cap_Size]*vCAP[t,y] + inputs[t]["dfGen"][y,:Fixed_OM_Cost_per_MWyr]*eTotalCap[t,y]
		elseif y in intersect(NEW_CAP, RETRO) # Resources eligible for new capacity (Retrofit yr -> y)
				sum( RETRO_SOURCE_IDS[y][i] in RET_CAP ? RETRO_INV_CAP_COSTS[y][i]*inputs[t]["dfGen"][y,:Cap_Size]*vRETROFIT[t,RETRO_SOURCE_IDS[y][i],y]*RETRO_EFFICIENCY[y][i] : 0 for i in 1:NUM_RETRO_SOURCES[y]) + inputs[t]["dfGen"][y,:Fixed_OM_Cost_per_MWyr]*eTotalCap[t,y]
		else
			inputs[t]["dfGen"][y,:Fixed_OM_Cost_per_MWyr]*eTotalCap[t,y]
		end
	)


	# Sum individual resource contributions to fixed costs to get total fixed costs
	@expression(EP, eTotalCFix[t in 1:num_stages], sum(EP[:eCFix][t,y] for y in 1:G))

	### Constratints ###

	## Constraints on retirements and capacity additions
	# Cannot retire more capacity than existing capacity
	
	#@constraint(EP,cAuxPositiveTotCap[t in 1:num_stages,y in 1:G], eTotalCap[t,y]>=0)
	@constraint(EP, cMaxRet[t in 1:num_stages, y in RET_CAP], inputs[t]["dfGen"][y,:Cap_Size]*vRETCAP[t,y] <= eMaxRetCap[t,y])
	
	## Constraints on new built capacity
	# Constraint on maximum capacity (if applicable) [set input to -1 if no constraint on maximum capacity]
	# DEV NOTE: This constraint may be violated in some cases where Existing_Cap_MW is >= Max_Cap_MW and lead to infeasabilty
	@constraint(EP, cMaxCap[t in 1:num_stages, y in intersect(inputs[t]["dfGen"][inputs[t]["dfGen"].Max_Cap_MW.>0, :R_ID], 1:G)],  eTotalCap[t,y] <= inputs[t]["dfGen"][y, :Max_Cap_MW])

	# Constraint on minimum capacity (if applicable) [set input to -1 if no constraint on minimum capacity]
	# DEV NOTE: This constraint may be violated in some cases where Existing_Cap_MW is <= Min_Cap_MW and lead to infeasabilty
	@constraint(EP, cMinCap[t in 1:num_stages, y in intersect(inputs[t]["dfGen"][inputs[t]["dfGen"].Min_Cap_MW.>0, :R_ID], 1:G)],  eTotalCap[t,y] >= inputs[t]["dfGen"][y, :Min_Cap_MW])

	if setup["MinCapReq"] == 1
		@expression(EP, eMinCapResInvest[t in 1:num_stages, mincap in 1:inputs[t]["NumberOfMinCapReqs"]], sum(eTotalCap[t,y] for y in inputs[t]["dfGen"][inputs[t]["dfGen"][!, Symbol("MinCapTag_$mincap")] .== 1, :R_ID]))
		EP[:eMinCapRes] += eMinCapResInvest
	end

	if setup["MaxCapReq"] == 1
		@expression(EP, eMaxCapResInvest[t in 1:num_stages, maxcap = 1:inputs[t]["NumberOfMaxCapReqs"]], sum(eTotalCap[t,y] for y in inputs[t]["dfGen"][inputs[t]["dfGen"][!, Symbol("MaxCapTag_$maxcap")] .== 1, :R_ID]))
		EP[:eMaxCapRes] += eMaxCapResInvest
	end


end

function multistage_long_duration_storage_master!(EP::Model,inputs::Dict,setup::Dict)
	
	num_stages = setup["MultiStageSettingsDict"]["NumStages"]

	nW = sum(inputs[t]["REP_PERIOD"] for t in 1:num_stages);

	STOR_LONG_DURATION = inputs[1]["STOR_LONG_DURATION"]

	dfPeriodMap = Dict(t=>inputs[t]["Period_Map"] for t in 1:num_stages) # Dataframe that maps modeled periods to representative periods

	NPeriods = Dict(t=>size(inputs[t]["Period_Map"])[1] for t in 1:num_stages) # Number of modeled periods
	NPeriods[0] = 0;

	MODELED_PERIODS_INDEX = Dict(t=>1:NPeriods[t] for t in 1:num_stages);

	ALL_MODELED_PERIODS_INDEX = 1:sum(NPeriods[t] for t in 1:num_stages);

	REP_PERIODS_INDEX = Dict(t=>MODELED_PERIODS_INDEX[t][dfPeriodMap[t][!,:Rep_Period] .== MODELED_PERIODS_INDEX[t]] for t in 1:num_stages)

	rescale_all_periods(r,t) = sum(NPeriods[s] for s in 1:t-1;init=0) + r; #take an index from MODELED_PERIODS_INDEX[t] and maps it to an index in ALL_MODELED_PERIODS_INDEX

	rescale_rep_periods(r,t) = sum(inputs[s]["REP_PERIOD"] for s in 1:t-1;init=0) + dfPeriodMap[t][r,:Rep_Period_Index];

	### Variables ###

	# Variables to define inter-period energy transferred between modeled periods

	# State of charge of storage at beginning of each modeled period n
	@variable(EP, vSOCw[y in STOR_LONG_DURATION, n in ALL_MODELED_PERIODS_INDEX] >= 0)

	# Build up in storage inventory over each representative period w
	# Build up inventory can be positive or negative
	@variable(EP, vdSOC[y in STOR_LONG_DURATION, w=1:nW])

	### Constraints ###

	# Storage at beginning of period w = storage at beginning of period w-1 + storage built up in period w (after n representative periods)
	## Multiply storage build up term from prior period with corresponding weight
	@constraint(EP, cSoCBalLongDurationStorage[t in 1:num_stages, y in STOR_LONG_DURATION, n in MODELED_PERIODS_INDEX[t]],
					vSOCw[y, rescale_all_periods(mod1(n+1, NPeriods[t]),t)] == vSOCw[y,rescale_all_periods(n,t)] + vdSOC[y,rescale_rep_periods(n,t)])

	# Storage at beginning of each modeled period cannot exceed installed energy capacity
	@constraint(EP, cSoCBalLongDurationStorageUpper[t in 1:num_stages, y in STOR_LONG_DURATION, n in MODELED_PERIODS_INDEX[t]],
	vSOCw[y,rescale_all_periods(n,t)] <= EP[:eTotalCapEnergy][t,y])
	

end

function multistage_hydro_inter_period_linkage_master!(EP::Model,inputs::Dict,setup::Dict)

	STOR_HYDRO_LONG_DURATION = inputs[1]["STOR_HYDRO_LONG_DURATION"]

	num_stages = setup["MultiStageSettingsDict"]["NumStages"]

	nW = sum(inputs[t]["REP_PERIOD"] for t in 1:num_stages);

	dfPeriodMap = Dict(t=>inputs[t]["Period_Map"] for t in 1:num_stages) # Dataframe that maps modeled periods to representative periods

	NPeriods = Dict(t=>size(inputs[t]["Period_Map"])[1] for t in 1:num_stages) # Number of modeled periods
	NPeriods[0] = 0;

	MODELED_PERIODS_INDEX = Dict(t=>1:NPeriods[t] for t in 1:num_stages);

	ALL_MODELED_PERIODS_INDEX = 1:sum(NPeriods[t] for t in 1:num_stages);

	REP_PERIODS_INDEX = Dict(t=>MODELED_PERIODS_INDEX[t][dfPeriodMap[t][!,:Rep_Period] .== MODELED_PERIODS_INDEX[t]] for t in 1:num_stages)

	rescale_all_periods(r,t) = sum(NPeriods[s] for s in 1:t-1;init=0) + r; #take an index from MODELED_PERIODS_INDEX[t] and maps it to an index in ALL_MODELED_PERIODS_INDEX

	rescale_rep_periods(r,t) = sum(inputs[s]["REP_PERIOD"] for s in 1:t-1;init=0) + dfPeriodMap[t][r,:Rep_Period_Index];


	### Variables ###

	# Variables to define inter-period energy transferred between modeled periods

	# State of charge of storage at beginning of each modeled period n
	@variable(EP, vSOC_HYDROw[y in STOR_HYDRO_LONG_DURATION, n in ALL_MODELED_PERIODS_INDEX] >= 0)

	# Build up in storage inventory over each representative period w
	# Build up inventory can be positive or negative
	@variable(EP, vdSOC_HYDRO[y in STOR_HYDRO_LONG_DURATION, w=1:nW])

	### Constraints ###

	# Storage at beginning of period w = storage at beginning of period w-1 + storage built up in period w (after n representative periods)
	## Multiply storage build up term from prior period with corresponding weight
	@constraint(EP, cHydroReservoirLongDurationStorage[t in 1:num_stages, y in STOR_HYDRO_LONG_DURATION, n in MODELED_PERIODS_INDEX[t]],
					vSOC_HYDROw[y, rescale_all_periods(mod1(n+1, NPeriods[t]),t)] == vSOC_HYDROw[y,rescale_all_periods(n,t)] + vdSOC_HYDRO[y,rescale_rep_periods(n,t)])

	# Storage at beginning of each modeled period cannot exceed installed energy capacity
	@constraint(EP, cHydroReservoirLongDurationStorageUpper[t in 1:num_stages, y in STOR_HYDRO_LONG_DURATION, n in MODELED_PERIODS_INDEX[t]],vSOC_HYDROw[y,rescale_all_periods(n,t)] <= inputs[t]["dfGen"][y,:Hydro_Energy_to_Power_Ratio]*EP[:eTotalCap][t,y])

end
