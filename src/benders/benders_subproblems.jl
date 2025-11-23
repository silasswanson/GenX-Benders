
function generate_subproblem(setup::Dict,inputs::Dict,master_vars::Vector{String},master_cons::Vector{String},OPTIMIZER::MOI.OptimizerWithAttributes)

	if setup["MultiStage"]==1
		
		EP = Model(OPTIMIZER)

		@variable(EP, vZERO == 1);

		println("Built generator capacity")

		load_built_discharge_cap!(EP,setup,inputs);

		println("Built transmission capacity")

		load_built_transmission_cap!(EP,setup,inputs);
	
		if !isempty(inputs["STOR_ALL"])
			println("Built energy capacity")
			load_built_energy_cap!(EP,setup,inputs);
		end
	
		if !isempty(inputs["STOR_ASYMMETRIC"])
			println("Built charge capacity")
			load_built_charge_cap!(EP,setup,inputs);
		end

		setup["MultiStage"]=0;
	
		load_operational_model!(EP,setup,inputs);
		
		setup["MultiStage"]=1;

	else

		EP = generate_model(setup,inputs,OPTIMIZER);

	end

	master_cons_to_be_deleted = intersect(name.(all_constraints(EP,include_variable_in_set_constraints=false)),master_cons);

	for sc in master_cons_to_be_deleted
		delete(EP,constraint_by_name(EP,sc))
		unregister(EP,Symbol(splitfun(sc)))
	end
	
	if setup["BD_Mode"]=="full"
		# do nothing
	else
		if setup["CO2Cap"] >= 1
			co2_cap_decomp!(EP,inputs,setup)
		end

		if setup["EnergyShareRequirement"] >= 1
			energy_share_requirement_decomp!(EP, inputs, setup)
		end

		#Because sub-problems correspond to one sub-period, LDS modelling will not have been included, and LDS slack have been skipped. 
		#However, we need them. So we define these things here:
		if !isempty(inputs["STOR_LONG_DURATION"]) || !isempty(inputs["STOR_HYDRO_LONG_DURATION"])
			lds_slack!(EP,inputs,setup)
		end

		if !isempty(inputs["STOR_LONG_DURATION"])
			long_duration_storage_decomp!(EP,inputs,setup)
		end

		if !isempty(inputs["STOR_HYDRO_LONG_DURATION"])
			hydro_inter_period_linkage_decomp!(EP,inputs,setup)
		end


	end

	if setup["MultiStage"] == 1
		DF = 1 / (1 + setup["MultiStageSettingsDict"]["WACC"])^(setup["MultiStageSettingsDict"]["StageLengths"][inputs["CurStage"]] * (inputs["CurStage"] - 1));

		@objective(EP,Min, DF*inputs["OPEXMULT"]*EP[:eObj])
	else
		@objective(EP,Min, EP[:eObj])
	end


	master_vars_sub = intersect(name.(all_variables(EP)),master_vars);
	for sv in master_vars_sub
		if has_lower_bound(variable_by_name(EP,sv))
			delete_lower_bound(variable_by_name(EP,sv))
		end
		if has_upper_bound(variable_by_name(EP,sv))
			delete_upper_bound(variable_by_name(EP,sv))
		end
		set_objective_coefficient(EP,variable_by_name(EP,sv),0)
	end
	
	set_objective_coefficient(EP,EP[:vZERO],0)
	set_silent(EP)
	
	return EP, master_vars_sub
end



function load_operational_model!(EP::Model,setup::Dict,inputs::Dict)

	T = inputs["T"]     # Number of time steps (hours)
	Z = inputs["Z"]     # Number of zones

	@expression(EP, ePowerBalance[t=1:T, z=1:Z], 0)

	@expression(EP, eObj, 0)

	@expression(EP, eGenerationByZone[z=1:Z, t=1:T], 0)

	# Initialize Capacity Reserve Margin Expression
	if setup["CapacityReserveMargin"] > 0
		@expression(EP, eCapResMarBalance[res=1:inputs["NCapacityReserveMargin"], t=1:T], 0)
	end
	
	# Energy Share Requirement
	if setup["EnergyShareRequirement"] >= 1
		@expression(EP, eESR[ESR=1:inputs["nESR"]], 0)
	end
	
	discharge!(EP, inputs, setup)

	non_served_energy!(EP, inputs, setup)

	if setup["UCommit"] > 0
		ucommit!(EP, inputs, setup)
	end

	emissions!(EP, inputs)

	if setup["Reserves"] > 0
		reserves!(EP, inputs, setup)
	end

	if Z > 1
		operations_transmission!(EP, inputs, setup)
	end

	if inputs["REP_PERIOD"] > 1 && (!isempty(inputs["STOR_LONG_DURATION"]) || !isempty(inputs["STOR_HYDRO_LONG_DURATION"]))
		lds_slack!(EP,inputs,setup)
	end

	if !isempty(inputs["VRE"])
		curtailable_variable_renewable!(EP, inputs, setup)
	end

	if !isempty(inputs["MUST_RUN"])
		must_run!(EP, inputs, setup)
	end

	if !isempty(inputs["STOR_ALL"])
		println("Storage Resources Module")

		storage_all!(EP, inputs, setup)

		# Include Long Duration Storage only when modeling representative periods and long-duration storage
		if inputs["REP_PERIOD"] > 1 && !isempty(inputs["STOR_LONG_DURATION"])
			long_duration_storage!(EP, inputs,setup)
		end

		if !isempty(inputs["STOR_ASYMMETRIC"])
			storage_asymmetric!(EP, inputs, setup)
		end
		
		if !isempty(inputs["STOR_SYMMETRIC"])
			storage_symmetric!(EP, inputs, setup)
		end
		
		# ESR Lossses
		if setup["EnergyShareRequirement"] >= 1
			if setup["IncludeLossesInESR"] == 1
				@expression(EP, eESRStor[ESR=1:inputs["nESR"]], sum(inputs["dfESR"][z,ESR]*sum(EP[:eELOSS][y] for y in intersect(inputs["dfGen"][inputs["dfGen"].Zone.==z,:R_ID],inputs["STOR_ALL"])) for z=findall(x->x>0,inputs["dfESR"][:,ESR])))
				EP[:eESR] -= eESRStor
			end
		end

		# Capacity Reserves Margin policy
		if setup["CapacityReserveMargin"] > 0
			@expression(EP, eCapResMarBalanceStor[res=1:inputs["NCapacityReserveMargin"], t=1:T], sum(inputs["dfGen"][y,Symbol("CapRes_$res")] * (EP[:vP][y,t] - EP[:vCHARGE][y,t])  for y in inputs["STOR_ALL"]))
			EP[:eCapResMarBalance] += eCapResMarBalanceStor
		end

	end

	if !isempty(inputs["HYDRO_RES"])
		hydro_res!(EP, inputs, setup)
	end

	if inputs["REP_PERIOD"] > 1 && !isempty(inputs["STOR_HYDRO_LONG_DURATION"])
		hydro_inter_period_linkage!(EP, inputs,setup)
	end

	if !isempty(inputs["FLEX"])
		flexible_demand!(EP, inputs, setup)
	end

	if !isempty(inputs["THERM_ALL"])
		thermal!(EP, inputs, setup)
	end

	if setup["CO2Cap"] >= 1
		co2_cap!(EP, inputs, setup)
	end

	if setup["EnergyShareRequirement"] >= 1
		energy_share_requirement!(EP, inputs, setup)
	end
	
	if setup["CapacityReserveMargin"] > 0
		cap_reserve_margin!(EP, inputs, setup)
	end

	@constraint(EP, cPowerBalance[t=1:T, z=1:Z], EP[:ePowerBalance][t,z] == inputs["pD"][t,z])

	return EP
end



function co2_cap_decomp!(EP::Model, inputs::Dict, setup::Dict)
	dfGen = inputs["dfGen"]
	SEG = inputs["SEG"]  # Number of lines
	G = inputs["G"]     # Number of resources (generators, storage, DR, and DERs)
	T = inputs["T"]     # Number of time steps (hours)
	Z = inputs["Z"]     # Number of zones

	w = inputs["SubPeriod"];
	
	
	delete(EP,vec(EP[:cCO2Emissions_systemwide]))
	unregister(EP,:cCO2Emissions_systemwide)
	
	### Variables ###

	@variable(EP,vCO2budget[[w],1:inputs["NCO2Cap"]])

	### Constraints ###

	## Mass-based: Emissions constraint in absolute emissions limit (tons)
	if setup["CO2Cap"] == 1
		@constraint(EP, cCO2Emissions_systemwide[cap=1:inputs["NCO2Cap"]],
			sum(inputs["omega"][t]/inputs["scenario_probability"] * EP[:eEmissionsByZone][z,t] for z=findall(x->x==1, inputs["dfCO2CapZones"][:,cap]), t=1:T)- EP[:vCO2Cap_slack][cap] <= vCO2budget[w,cap]
		)

	## Load + Rate-based: Emissions constraint in terms of rate (tons/MWh)
	elseif setup["CO2Cap"] == 2 ##This part moved to non_served_energy.jl
		@constraint(EP, cCO2Emissions_systemwide[cap=1:inputs["NCO2Cap"]],
			sum(inputs["omega"][t] * EP[:eEmissionsByZone][z,t] for z=findall(x->x==1, inputs["dfCO2CapZones"][:,cap]), t=1:T) - EP[:vCO2Cap_slack][cap] <= vCO2budget[w,cap]
			- sum(inputs["dfMaxCO2Rate"][z,cap] * sum(inputs["omega"][t] * sum(EP[:vNSE][s,t,z] for s in 1:SEG) for t=1:T) for z = findall(x->x==1, inputs["dfCO2CapZones"][:,cap])) +
			sum(inputs["dfMaxCO2Rate"][z,cap] * setup["StorageLosses"] *  EP[:eELOSSByZone][z] for z=findall(x->x==1, inputs["dfCO2CapZones"][:,cap]))
		)

	## Generation + Rate-based: Emissions constraint in terms of rate (tons/MWh)
	elseif (setup["CO2Cap"]==3)
		@constraint(EP, cCO2Emissions_systemwide[cap=1:inputs["NCO2Cap"]],
			sum(inputs["omega"][t] * EP[:eEmissionsByZone][z,t] for z=findall(x->x==1, inputs["dfCO2CapZones"][:,cap]), t=1:T) - EP[:vCO2Cap_slack][cap] <= vCO2budget[w,cap] + 
			sum(inputs["dfMaxCO2Rate"][z,cap] * inputs["omega"][t] * EP[:eGenerationByZone][z,t] for t=1:T, z=findall(x->x==1, inputs["dfCO2CapZones"][:,cap]))
		)
	end 


end



function energy_share_requirement_decomp!(EP::Model, inputs::Dict, setup::Dict)

	delete(EP,EP[:cESRShare])
	unregister(EP,:cESRShare)
	
	w = inputs["SubPeriod"];

	T = inputs["T"]

	weighted_demand = [sum(inputs["omega"][t]*inputs["pD"][t,z] for t in 1:T) for z in 1:inputs["Z"]];

	dfGen = inputs["dfGen"]
	IncludeLossesInESR = setup["IncludeLossesInESR"]
	STOR_ALL = inputs["STOR_ALL"]

	### Variables ###

	@variable(EP,vESRbudget[[w],1:inputs["nESR"]])

	### Constraints ###

	# if IncludeLossesInESR==1
		
	# 	# @constraint(EP,cESRShare2[ESR=1:inputs["nESR"]], EP[:eESRDischarge][ESR]+ sum(inputs["dfESR"][z,ESR]*weighted_demand[z] for z=findall(x->x>0,inputs["dfESR"][:,ESR]))-EP[:eESRStor][ESR] -EP[:eESRTran][ESR] + EP[:vESR_slack][ESR] >= vESRbudget[ESR,w])

	# 	@constraint(EP, cESRShare[ESR=1:inputs["nESR"]], sum(inputs["omega"][t]*dfGen[y,Symbol("ESR_$ESR")]*EP[:vP][t,y] for t in 1:T, y=dfGen[findall(x->x>0,dfGen[!,Symbol("ESR_$ESR")]),:R_ID]) -EP[:eESRStor][ESR] -EP[:eESRTran][ESR] + EP[:vESR_slack][ESR] >= vESRbudget[ESR,w])

	# else
	# 	#@constraint(EP,cESRShare2[ESR=1:inputs["nESR"]], EP[:eESRDischarge][ESR]+ sum(inputs["dfESR"][z,ESR]*weighted_demand[z] for z=findall(x->x>0,inputs["dfESR"][:,ESR])) + EP[:vESR_slack][ESR] >= vESRbudget[ESR,w])

	# 	@constraint(EP, cESRShare[ESR=1:inputs["nESR"]], sum(inputs["omega"][t]*dfGen[y,Symbol("ESR_$ESR")]*EP[:vP][t,y] for t in 1:T, y=dfGen[findall(x->x>0,dfGen[!,Symbol("ESR_$ESR")]),:R_ID]) + EP[:vESR_slack][ESR] >= vESRbudget[ESR,w])
	# end

	
	@constraint(EP, cESRShare[ESR=1:inputs["nESR"]], EP[:eESR][ESR] >= vESRbudget[w,ESR])


end


function long_duration_storage_decomp!(EP::Model,inputs::Dict,setup::Dict)

	w = inputs["SubPeriod"];

	r = inputs["SubPeriod_Index"]

	dfGen = inputs["dfGen"]

	STOR_LONG_DURATION = inputs["STOR_LONG_DURATION"]
	STOR_SHORT_DURATION = inputs["STOR_SHORT_DURATION"]
	START_SUBPERIODS = inputs["START_SUBPERIODS"]

	hours_per_subperiod = inputs["hours_per_subperiod"] #total number of hours per subperiod

	if haskey(EP,:cSoCBalLongDurationStorageStart)
		delete.(EP,EP[:cSoCBalLongDurationStorageStart])
		unregister(EP,:cSoCBalLongDurationStorageStart)
	end

	if haskey(EP,:cSoCBalLongDurationStorageSub)
		delete.(EP,EP[:cSoCBalLongDurationStorageSub])
		unregister(EP,:cSoCBalLongDurationStorageSub)
	end

	if haskey(EP,:cSoCBalStart)
		delete.(EP,EP[:cSoCBalStart])
		unregister(EP,:cSoCBalStart)
		#### Even if there is only one period in the sub-problem, we are actually modeling LDS storage so this constraint has to be restricted to STOR_SHORT_DURATION (superseding the if-condition in storage_all.jl)
		@constraint(EP, cSoCBalStart[t in START_SUBPERIODS, y in STOR_SHORT_DURATION], EP[:vS][y,t] ==
			(1-dfGen[y,:Self_Disch])*EP[:vS][y,t+hours_per_subperiod-1] - (1/dfGen[y,:Eff_Down] * EP[:vP][y,t])
			+ (dfGen[y,:Eff_Up]*EP[:vCHARGE][y,t]))
	end


	### Variables ###

	# Variables to define inter-period energy transferred between modeled periods

	# State of charge of storage at beginning of each modeled period n
	@variable(EP, vSOCw[y in STOR_LONG_DURATION, [r]] >= 0)

	# Build up in storage inventory over each representative period w
	# Build up inventory can be positive or negative
	@variable(EP, vdSOC[y in STOR_LONG_DURATION, [w]])

	@variable(EP, vLDS_Start_slack[[w], y in STOR_LONG_DURATION])

	@variable(EP, vLDS_Sub_slack[y in STOR_LONG_DURATION,[r]])
	
	@constraint(EP,cSlackLDS_Start_Up[[w], y in STOR_LONG_DURATION],vLDS_Start_slack[w,y] <= EP[:vLDS_SLACK_MAX][1])
	@constraint(EP,cSlackLDS_Start_Lo[[w], y in STOR_LONG_DURATION],-vLDS_Start_slack[w,y]<= EP[:vLDS_SLACK_MAX][1])
	@constraint(EP,cSlackLDS_Sub_Up[y in STOR_LONG_DURATION, [r]],vLDS_Sub_slack[y,r]<= EP[:vLDS_SLACK_MAX][1])
	@constraint(EP,cSlackLDS_Sub_Lo[y in STOR_LONG_DURATION, [r]],-vLDS_Sub_slack[y,r]<= EP[:vLDS_SLACK_MAX][1])

	### Constraints ###

	# Links last time step with first time step, ensuring position in hour 1 is within eligible change from final hour position
	# Modified initial state of storage for long-duration storage - initialize wth value carried over from last period
	@constraint(EP, cSoCBalLongDurationStorageStart[[w], y in STOR_LONG_DURATION],
				    EP[:vS][y,1] == (1-dfGen[y,:Self_Disch])*( EP[:vS][y,hours_per_subperiod] - vdSOC[y,w])-(1/dfGen[y,:Eff_Down]*EP[:vP][y,1])+(dfGen[y,:Eff_Up]*EP[:vCHARGE][y,1]+ vLDS_Start_slack[w,y]))

	# Initial storage level for representative periods must also adhere to sub-period storage inventory balance
	# Initial storage = Final storage - change in storage inventory across representative period
	@constraint(EP, cSoCBalLongDurationStorageSub[y in STOR_LONG_DURATION, [r]],
					vSOCw[y,r] == EP[:vS][y,hours_per_subperiod] - vdSOC[y,w] +vLDS_Sub_slack[y,r])

end




function hydro_inter_period_linkage_decomp!(EP::Model,inputs::Dict,setup::Dict)
	
	w = inputs["SubPeriod"];

	r = inputs["SubPeriod_Index"]

	dfGen = inputs["dfGen"]

	STOR_HYDRO_LONG_DURATION = inputs["STOR_HYDRO_LONG_DURATION"]
	STOR_HYDRO_SHORT_DURATION = inputs["STOR_HYDRO_SHORT_DURATION"]

	START_SUBPERIODS = inputs["START_SUBPERIODS"]

	hours_per_subperiod = inputs["hours_per_subperiod"] #total number of hours per subperiod


	if haskey(EP,:cHydroReservoirLongDurationStorageStart)
		delete.(EP,EP[:cHydroReservoirLongDurationStorageStart])
		unregister(EP,:cHydroReservoirLongDurationStorageStart)
	end

	if haskey(EP,:cHydroReservoirLongDurationStorageSub)
		delete.(EP,EP[:cHydroReservoirLongDurationStorageSub])
		unregister(EP,:cHydroReservoirLongDurationStorageSub)
	end

	if haskey(EP,:cHydroReservoirStart)
		delete.(EP,EP[:cHydroReservoirStart])
		unregister(EP,:cHydroReservoirStart)
		#### Even if there is only one period in the sub-problem, we are actually modeling LDS storage so this constraint has to be restricted to STOR_HYDRO_SHORT_DURATION
		@constraint(EP, cHydroReservoirStart[y in STOR_HYDRO_SHORT_DURATION,t in START_SUBPERIODS], EP[:vS_HYDRO][y,t] == EP[:vS_HYDRO][y, hoursbefore(hours_per_subperiod,t,1)]- (1/dfGen[y,:Eff_Down]*EP[:vP][y,t]) - vSPILL[y,t] + inputs["pP_Max"][y,t]*EP[:eTotalCap][y])

	end

	### Variables ###

	# Variables to define inter-period energy transferred between modeled periods

	# State of charge of storage at beginning of each modeled period n
	@variable(EP, vSOC_HYDROw[y in STOR_HYDRO_LONG_DURATION, [r]] >= 0)

	# Build up in storage inventory over each representative period w
	# Build up inventory can be positive or negative
	@variable(EP, vdSOC_HYDRO[y in STOR_HYDRO_LONG_DURATION, [w]])

	@variable(EP, vHydro_Start_slack[[w], y in STOR_HYDRO_LONG_DURATION])

	@variable(EP, vHydro_Sub_slack[y in STOR_HYDRO_LONG_DURATION, [r]])

	@constraint(EP,cSlackHydro_Start_Up[[w], y in STOR_HYDRO_LONG_DURATION],vHydro_Start_slack[w,y] <= EP[:vLDS_SLACK_MAX][1])
	@constraint(EP,cSlackHydro_Start_Lo[[w], y in STOR_HYDRO_LONG_DURATION],-vHydro_Start_slack[w,y]<= EP[:vLDS_SLACK_MAX][1])
	@constraint(EP,cSlackHydro_Sub_Up[y in STOR_HYDRO_LONG_DURATION, [r]],vHydro_Sub_slack[y,r]<= EP[:vLDS_SLACK_MAX][1])
	@constraint(EP,cSlackHydro_Sub_Lo[y in STOR_HYDRO_LONG_DURATION, [r]],-vHydro_Sub_slack[y,r]<= EP[:vLDS_SLACK_MAX][1])

	### Constraints ###

	# Links last time step with first time step, ensuring position in hour 1 is within eligible change from final hour position
	# Modified initial state of storage for long-duration storage - initialize wth value carried over from last period
	@constraint(EP, cHydroReservoirLongDurationStorageStart[[w], y in STOR_HYDRO_LONG_DURATION],
				    EP[:vS_HYDRO][y,1] == (EP[:vS_HYDRO][y,hours_per_subperiod]-vdSOC_HYDRO[y,w])-(1/dfGen[y,:Eff_Down]*EP[:vP][y,1])-EP[:vSPILL][y,1]+inputs["pP_Max"][y,1]*EP[:eTotalCap][y]+ vHydro_Start_slack[w,y])

	# Initial storage level for representative periods must also adhere to sub-period storage inventory balance
	# Initial storage = Final storage - change in storage inventory across representative period
	@constraint(EP, cHydroReservoirLongDurationStorageSub[y in STOR_HYDRO_LONG_DURATION, [r]],
					vSOC_HYDROw[y,r] == EP[:vS_HYDRO][y,hours_per_subperiod] - vdSOC_HYDRO[y,w] + vHydro_Sub_slack[y,r])

end



function operations_transmission!(EP::Model,inputs::Dict,setup::Dict)

	println("Transmission Module")

	T = inputs["T"]     # Number of time steps (hours)
	Z = inputs["Z"]     # Number of zones
	L = inputs["L"]     # Number of transmission lines

	UCommit = setup["UCommit"]
	CapacityReserveMargin = setup["CapacityReserveMargin"]
	EnergyShareRequirement = setup["EnergyShareRequirement"]
	IncludeLossesInESR = setup["IncludeLossesInESR"]

	## sets and indices for transmission losses and expansion
	TRANS_LOSS_SEGS = inputs["TRANS_LOSS_SEGS"] # Number of segments used in piecewise linear approximations quadratic loss functions - can only take values of TRANS_LOSS_SEGS =1, 2
	LOSS_LINES = inputs["LOSS_LINES"] # Lines for which loss coefficients apply (are non-zero);

	### Variables ###

	# Power flow on each transmission line "l" at hour "t"
	@variable(EP, vFLOW[l=1:L,t=1:T]);

  	if (TRANS_LOSS_SEGS==1)  #loss is a constant times absolute value of power flow
		# Positive and negative flow variables
		@variable(EP, vTAUX_NEG[l in LOSS_LINES,t=1:T] >= 0)
		@variable(EP, vTAUX_POS[l in LOSS_LINES,t=1:T] >= 0)

		if UCommit == 1
			# Single binary variable to ensure positive or negative flows only
			@variable(EP, vTAUX_POS_ON[l in LOSS_LINES,t=1:T],Bin)
			# Continuous variable representing product of binary variable (vTAUX_POS_ON) and avail transmission capacity
			@variable(EP, vPROD_TRANSCAP_ON[l in LOSS_LINES,t=1:T]>=0)
		end
	else # TRANS_LOSS_SEGS>1
		# Auxiliary variables for linear piecewise interpolation of quadratic losses
		@variable(EP, vTAUX_NEG[l in LOSS_LINES, s=0:TRANS_LOSS_SEGS, t=1:T] >= 0)
		@variable(EP, vTAUX_POS[l in LOSS_LINES, s=0:TRANS_LOSS_SEGS, t=1:T] >= 0)
		if UCommit == 1
			# Binary auxilary variables for each segment >1 to ensure segments fill in order
			@variable(EP, vTAUX_POS_ON[l in LOSS_LINES, s=1:TRANS_LOSS_SEGS, t=1:T], Bin)
			@variable(EP, vTAUX_NEG_ON[l in LOSS_LINES, s=1:TRANS_LOSS_SEGS, t=1:T], Bin)
		end
    end

	# Transmission losses on each transmission line "l" at hour "t"
	@variable(EP, vTLOSS[l in LOSS_LINES,t=1:T] >= 0)

	### Expressions ###

	# Net power flow outgoing from zone "z" at hour "t" in MW
    @expression(EP, eNet_Export_Flows[z=1:Z,t=1:T], sum(inputs["pNet_Map"][l,z] * vFLOW[l,t] for l=1:L))

	# Losses from power flows into or out of zone "z" in MW
    	@expression(EP, eLosses_By_Zone[z=1:Z,t=1:T], sum(abs(inputs["pNet_Map"][l,z]) * (1/2) *vTLOSS[l,t] for l in LOSS_LINES))


	## Power Balance Expressions ##

	@expression(EP, ePowerBalanceNetExportFlows[t=1:T, z=1:Z],
		-eNet_Export_Flows[z,t])
	@expression(EP, ePowerBalanceLossesByZone[t=1:T, z=1:Z],
		-eLosses_By_Zone[z,t])

	EP[:ePowerBalance] += ePowerBalanceLossesByZone
	EP[:ePowerBalance] += ePowerBalanceNetExportFlows

	# Capacity Reserves Margin policy
	if CapacityReserveMargin > 0
		if Z > 1 
			@expression(EP, eCapResMarBalanceTrans[res=1:inputs["NCapacityReserveMargin"], t=1:T], sum(inputs["dfTransCapRes_excl"][l,res] * inputs["dfDerateTransCapRes"][l,res]* EP[:vFLOW][l,t] for l in 1:L))
			EP[:eCapResMarBalance] -= eCapResMarBalanceTrans
		end
	end

	### Constraints ###

  	## Power flow and transmission (between zone) loss related constraints

	# Maximum power flows, power flow on each transmission line cannot exceed maximum capacity of the line at any hour "t"
	# Allow expansion of transmission capacity for lines eligible for reinforcement
	@constraints(EP, begin
		cMaxFlow_out[l=1:L, t=1:T], vFLOW[l,t] <= EP[:eAvail_Trans_Cap][l]
		cMaxFlow_in[l=1:L, t=1:T], vFLOW[l,t] >= -EP[:eAvail_Trans_Cap][l]
	end)

	# Transmission loss related constraints - linear losses as a function of absolute value
	if TRANS_LOSS_SEGS == 1

		@constraints(EP, begin
			# Losses are alpha times absolute values
			cTLoss[l in LOSS_LINES, t=1:T], vTLOSS[l,t] == inputs["pPercent_Loss"][l]*(vTAUX_POS[l,t]+vTAUX_NEG[l,t])

			# Power flow is sum of positive and negative components
			cTAuxSum[l in LOSS_LINES, t=1:T], vTAUX_POS[l,t]-vTAUX_NEG[l,t] == vFLOW[l,t]

			# Sum of auxiliary flow variables in either direction cannot exceed maximum line flow capacity
			cTAuxLimit[l in LOSS_LINES, t=1:T], vTAUX_POS[l,t]+vTAUX_NEG[l,t] <= EP[:eAvail_Trans_Cap][l]
		end)

		if UCommit == 1
			# Constraints to limit phantom losses that can occur to avoid discrete cycling costs/opportunity costs due to min down
			@constraints(EP, begin
				cTAuxPosUB[l in LOSS_LINES, t=1:T], vTAUX_POS[l,t] <= vPROD_TRANSCAP_ON[l,t]

				# Either negative or positive flows are activated, not both
				cTAuxNegUB[l in LOSS_LINES, t=1:T], vTAUX_NEG[l,t] <= EP[:eAvail_Trans_Cap][l]-vPROD_TRANSCAP_ON[l,t]

				# McCormick representation of product of continuous and binary variable
				# (in this case, of: vPROD_TRANSCAP_ON[l,t] = EP[:eAvail_Trans_Cap][l] * vTAUX_POS_ON[l,t])
				# McCormick constraint 1
				[l in LOSS_LINES,t=1:T], vPROD_TRANSCAP_ON[l,t] <= inputs["pTrans_Max_Possible"][l]*vTAUX_POS_ON[l,t]

				# McCormick constraint 2
				[l in LOSS_LINES,t=1:T], vPROD_TRANSCAP_ON[l,t] <= EP[:eAvail_Trans_Cap][l]

				# McCormick constraint 3
				[l in LOSS_LINES,t=1:T], vPROD_TRANSCAP_ON[l,t] >= EP[:eAvail_Trans_Cap][l]-(1-vTAUX_POS_ON[l,t])*inputs["pTrans_Max_Possible"][l]
			end)
		end

	end # End if(TRANS_LOSS_SEGS == 1) block

	# When number of segments is greater than 1
	if (TRANS_LOSS_SEGS > 1)
		## between zone transmission loss constraints
		# Losses are expressed as a piecewise approximation of a quadratic function of power flows across each line
		# Eq 1: Total losses are function of loss coefficient times the sum of auxilary segment variables across all segments of piecewise approximation
		# (Includes both positive domain and negative domain segments)
		@constraint(EP, cTLoss[l in LOSS_LINES, t=1:T], vTLOSS[l,t] ==
							(inputs["pTrans_Loss_Coef"][l]*sum((2*s-1)*(inputs["pTrans_Max_Possible"][l]/TRANS_LOSS_SEGS)*vTAUX_POS[l,s,t] for s=1:TRANS_LOSS_SEGS)) +
							(inputs["pTrans_Loss_Coef"][l]*sum((2*s-1)*(inputs["pTrans_Max_Possible"][l]/TRANS_LOSS_SEGS)*vTAUX_NEG[l,s,t] for s=1:TRANS_LOSS_SEGS)) )
		# Eq 2: Sum of auxilary segment variables (s >= 1) minus the "zero" segment (which allows values to go negative)
		# from both positive and negative domains must total the actual power flow across the line
		@constraints(EP, begin
			cTAuxSumPos[l in LOSS_LINES, t=1:T], sum(vTAUX_POS[l,s,t] for s=1:TRANS_LOSS_SEGS)-vTAUX_POS[l,0,t]  == vFLOW[l,t]
			cTAuxSumNeg[l in LOSS_LINES, t=1:T], sum(vTAUX_NEG[l,s,t] for s=1:TRANS_LOSS_SEGS) - vTAUX_NEG[l,0,t]  == -vFLOW[l,t]
		end)
		if UCommit == 0 || UCommit == 2
			# Eq 3: Each auxilary segment variables (s >= 1) must be less than the maximum power flow in the zone / number of segments
			@constraints(EP, begin
				cTAuxMaxPos[l in LOSS_LINES, s=1:TRANS_LOSS_SEGS, t=1:T], vTAUX_POS[l,s,t] <= (inputs["pTrans_Max_Possible"][l]/TRANS_LOSS_SEGS)
				cTAuxMaxNeg[l in LOSS_LINES, s=1:TRANS_LOSS_SEGS, t=1:T], vTAUX_NEG[l,s,t] <= (inputs["pTrans_Max_Possible"][l]/TRANS_LOSS_SEGS)
			end)
		else # Constraints that can be ommitted if problem is convex (i.e. if not using MILP unit commitment constraints)
			# Eqs 3-4: Ensure that auxilary segment variables do not exceed maximum value per segment and that they
			# "fill" in order: i.e. one segment cannot be non-zero unless prior segment is at it's maximum value
			# (These constraints are necessary to prevents phantom losses in MILP problems)
			@constraints(EP, begin
				cTAuxOrderPos1[l in LOSS_LINES, s=1:TRANS_LOSS_SEGS, t=1:T], vTAUX_POS[l,s,t] <=  (inputs["pTrans_Max_Possible"][l]/TRANS_LOSS_SEGS)*vTAUX_POS_ON[l,s,t]
				cTAuxOrderNeg1[l in LOSS_LINES, s=1:TRANS_LOSS_SEGS, t=1:T], vTAUX_NEG[l,s,t] <=  (inputs["pTrans_Max_Possible"][l]/TRANS_LOSS_SEGS)*vTAUX_NEG_ON[l,s,t]
				cTAuxOrderPos2[l in LOSS_LINES, s=1:(TRANS_LOSS_SEGS-1), t=1:T], vTAUX_POS[l,s,t] >=  (inputs["pTrans_Max_Possible"][l]/TRANS_LOSS_SEGS)*vTAUX_POS_ON[l,s+1,t]
				cTAuxOrderNeg2[l in LOSS_LINES, s=1:(TRANS_LOSS_SEGS-1), t=1:T], vTAUX_NEG[l,s,t] >=  (inputs["pTrans_Max_Possible"][l]/TRANS_LOSS_SEGS)*vTAUX_NEG_ON[l,s+1,t]
			end)

			# Eq 5: Binary constraints to deal with absolute value of vFLOW.
			@constraints(EP, begin
				# If flow is positive, vTAUX_POS segment 0 must be zero; If flow is negative, vTAUX_POS segment 0 must be positive
				# (and takes on value of the full negative flow), forcing all vTAUX_POS other segments (s>=1) to be zero
				cTAuxSegmentZeroPos[l in LOSS_LINES, t=1:T], vTAUX_POS[l,0,t] <= inputs["pTrans_Max_Possible"][l]*(1-vTAUX_POS_ON[l,1,t])

				# If flow is negative, vTAUX_NEG segment 0 must be zero; If flow is positive, vTAUX_NEG segment 0 must be positive
				# (and takes on value of the full positive flow), forcing all other vTAUX_NEG segments (s>=1) to be zero
				cTAuxSegmentZeroNeg[l in LOSS_LINES, t=1:T], vTAUX_NEG[l,0,t] <= inputs["pTrans_Max_Possible"][l]*(1-vTAUX_NEG_ON[l,1,t])
			end)
		end
	end # End if(TRANS_LOSS_SEGS > 0) block

	# ESR Lossses
    if EnergyShareRequirement >= 1 && IncludeLossesInESR ==1
        @expression(EP, eESRTran[ESR=1:inputs["nESR"]],
                    sum(inputs["dfESR"][z,ESR]*sum(inputs["omega"][t]*EP[:eLosses_By_Zone][z,t] for t in 1:T) for z=findall(x->x>0,inputs["dfESR"][:,ESR])))
        EP[:eESR] -= eESRTran
	end


end




function generate_full_operational_subproblem(setup::Dict,inputs::Dict,OPTIMIZER::MOI.OptimizerWithAttributes,master_vars::Vector{String},master_cons::Vector{String})
	if setup["MultiStage"]==0
		EP,master_vars_sub = generate_subproblem(setup,inputs,master_vars,master_cons,OPTIMIZER);
	else
		EP, master_vars_sub = initialize_dist_helpers(setup,inputs,master_vars,master_cons);
	end
    return EP,master_vars_sub
end

function generate_decomp_operational_subproblems(setup::Dict,inputs_decomp::Dict,OPTIMIZER::MOI.OptimizerWithAttributes,master_vars::Vector{String},master_cons::Vector{String})
	## Start pre-solve timer
	subproblem_generation_time = time()
	nW = length(inputs_decomp);
	EP = Dict();
	master_vars_sub=Dict();

	for w in 1:nW
		EP[w],master_vars_sub[w] = generate_subproblem(setup,inputs_decomp[w],master_vars,master_cons,OPTIMIZER)
	end

	## Record pre-solver time
	subproblem_generation_time = time() - subproblem_generation_time
	println("Decomposed operational subproblems generation took $subproblem_generation_time seconds")
	return EP,master_vars_sub
end

function initialize_dist_helpers(setup::Dict,inputs_decomp::Dict,master_vars::Vector{String},master_cons::Vector{String})

    ##### Initialize a distributed arrays of JuMP models
	## Start pre-solve timer
	subproblem_generation_time = time()
    
    helpers_all = distribute([Dict() for i in 1:length(inputs_decomp)]);

    @sync for p in workers()
        @async @spawnat p begin
            W_local = localindices(helpers_all)[1];
            inputs_local = [inputs_decomp[k] for k in W_local];
			SUBPROB_OPTIMIZER =  configure_benders_subprob_solver(setup["Solver"], setup["settings_path"]);
            init_local_helper!(setup,inputs_local,localpart(helpers_all),master_vars,master_cons,SUBPROB_OPTIMIZER);
        end
    end

	# master_vars_sub = Dict();
	# for i in eachindex(helpers_all)
	# 	w = helpers_all[i]["SubPeriod"];
	# 	master_vars_sub[w] = helpers_all[i]["master_vars_sub"];
	# end

	p_id = workers();
    np_id = length(p_id);

    master_vars_sub = [Dict() for k in 1:np_id];

    @sync for k in 1:np_id
              @async master_vars_sub[k]= @fetchfrom p_id[k] get_local_master_vars(localpart(helpers_all))
    end

	master_vars_sub = merge(master_vars_sub...);

    ## Record pre-solver time
	subproblem_generation_time = time() - subproblem_generation_time
	println("Distributed operational subproblems generation took $subproblem_generation_time seconds")

    return helpers_all,master_vars_sub

end

function init_local_helper!(setup::Dict,inputs_local::Vector{Dict{Any,Any}},helper_local::Vector{Dict{Any,Any}},master_vars::Vector{String},master_cons::Vector{String},OPTIMIZER::MOI.OptimizerWithAttributes)

    nW = length(inputs_local)

    for i=1:nW
		EP, master_vars_sub = generate_subproblem(setup,inputs_local[i],master_vars,master_cons,OPTIMIZER);
        helper_local[i]["Model"] = EP;
        helper_local[i]["master_vars_sub"] = master_vars_sub
        helper_local[i]["SubPeriod"] = inputs_local[i]["SubPeriod"];
    end
end
