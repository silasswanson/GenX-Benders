@doc raw"""
	write_emissions(path::AbstractString, inputs::Dict, setup::Dict, EP::Model)

Function for reporting time-dependent CO$_2$ emissions by zone.

"""
function write_emissions(path::AbstractString, inputs::Dict, setup::Dict, EP::Model)
	dfGen = inputs["dfGen"]
	G = inputs["G"]     # Number of resources (generators, storage, DR, and DERs)
	T = inputs["T"]     # Number of time steps (hours)
	Z = inputs["Z"]     # Number of zones

	scale_factor = setup["ParameterScale"] == 1 ? ModelScalingFactor : 1


	if (setup["WriteShadowPrices"]==1 || setup["UCommit"]==0 || (setup["UCommit"]==2 && (setup["Reserves"]==0 || (setup["Reserves"]>0 && inputs["pDynamic_Contingency"]==0)))) # fully linear model
		# CO2 emissions by zone

		if setup["CO2Cap"]>=1
			if inputs["NumScenarios"]<=1
				# Dual variable of CO2 constraint = shadow price of CO2
				tempCO2Price = zeros(Z,inputs["NCO2Cap"])
				if has_duals(EP) == 1
					for cap in 1:inputs["NCO2Cap"]
						for z in findall(x->x==1, inputs["dfCO2CapZones"][:,cap])
							tempCO2Price[z,cap] = (-1) * dual.(EP[:cCO2Emissions_systemwide])[cap]
							# when scaled, The objective function is in unit of Million US$/kton, thus k$/ton, to get $/ton, multiply 1000
							tempCO2Price[z,cap] *= scale_factor
						end
					end
				end
				dfEmissions = hcat(DataFrame(Zone = 1:Z), DataFrame(tempCO2Price, :auto), DataFrame(AnnualSum = Array{Float64}(undef, Z)))
				auxNew_Names=[Symbol("Zone"); [Symbol("CO2_Price_$cap") for cap in 1:inputs["NCO2Cap"]]; Symbol("AnnualSum")]
				rename!(dfEmissions,auxNew_Names)
				for i in 1:Z
					dfEmissions[i,:AnnualSum] = sum(inputs["omega"].*value.(EP[:eEmissionsByZone][i,:]))*scale_factor
				end
			else
				# Dual variable of CO2 constraint = shadow price of CO2
				tempCO2Price = zeros(Z,inputs["NCO2Cap"]*inputs["NumScenarios"])
				if has_duals(EP) == 1
					for s in 1:inputs["NumScenarios"]
						for cap in 1:inputs["NCO2Cap"]
							for z in findall(x->x==1, inputs["dfCO2CapZones"][:,cap])
								tempCO2Price[z,(s-1)*inputs["NCO2Cap"]+cap] = (-1) * dual.(EP[:cCO2Emissions_systemwide])[cap,s]
								# when scaled, The objective function is in unit of Million US$/kton, thus k$/ton, to get $/ton, multiply 1000
								tempCO2Price[z,(s-1)*inputs["NCO2Cap"]+cap] *= scale_factor
							end
						end
					end
				end
				dfEmissions = hcat(DataFrame(Zone = 1:Z), DataFrame([tempCO2Price Array{Float64}(undef, Z,inputs["NumScenarios"])], :auto))
				auxNew_Names=[Symbol("Zone"); [Symbol("CO2_Price_$(cap)_Scenario_$s") for cap in 1:inputs["NCO2Cap"], s=1:inputs["NumScenarios"]][:]; [Symbol("AnnualSum_Scenario_$s") for s=1:inputs["NumScenarios"]]]
				rename!(dfEmissions,auxNew_Names)
				scenario_times = timeseries2scenarios(T, inputs["NumScenarios"])
				for i in 1:Z
					for s in 1:inputs["NumScenarios"]
						dfEmissions[i,Symbol("AnnualSum_Scenario_$s")] = sum(inputs["omega"][t]*value(EP[:eEmissionsByZone][i,t]) for t in scenario_times[s,1]:scenario_times[s,2])*scale_factor
					end
				end
			end

		else
			if inputs["NumScenarios"]<=1
				dfEmissions = DataFrame(Zone = 1:Z, AnnualSum = Array{Float64}(undef, Z))
				for i in 1:Z
					dfEmissions[i,:AnnualSum] = sum(inputs["omega"].*value.(EP[:eEmissionsByZone][i,:]))*scale_factor
				end
			else
				dfEmissions = DataFrame([1:Z Array{Float64}(undef,Z,inputs["NumScenarios"])],:auto);
				auxNew_Names = [Symbol("Zone");[Symbol("AnnualSum_Scenario_$s") for s=1:inputs["NumScenarios"]]];
				rename!(dfEmissions,auxNew_Names)
				scenario_times = timeseries2scenarios(T, inputs["NumScenarios"])
				for i in 1:Z
					for s in inputs["NumScenarios"]
						dfEmissions[i,Symbol("AnnualSum_Scenario_$s")] = sum(inputs["omega"][t]*value.(EP[:eEmissionsByZone][i,t]) for t in scenario_times[s,1]:scenario_times[s,2])*scale_factor
					end
				end
			end
		end

		dfEmissions = hcat(dfEmissions, DataFrame(value.(EP[:eEmissionsByZone])*scale_factor, :auto))


		if setup["CO2Cap"]>=1
			if inputs["NumScenarios"]<=1
				auxNew_Names=[Symbol("Zone");[Symbol("CO2_Price_$cap") for cap in 1:inputs["NCO2Cap"]];Symbol("AnnualSum");[Symbol("t$t") for t in 1:T]]
				rename!(dfEmissions,auxNew_Names)
				total = DataFrame(["Total" zeros(1,inputs["NCO2Cap"]) sum(dfEmissions[!,:AnnualSum]) fill(0.0, (1,T))], :auto)
				for t in 1:T
					total[:,t+inputs["NCO2Cap"]+2] .= sum(dfEmissions[:,Symbol("t$t")][1:Z])
				end
			else
				auxNew_Names=[Symbol("Zone");[Symbol("CO2_Price_$(cap)_Scenario_$s") for cap in 1:inputs["NCO2Cap"], s in 1:inputs["NumScenarios"]][:];[Symbol("AnnualSum_Scenario_$s") for s=1:inputs["NumScenarios"]];[Symbol("t$t") for t in 1:T]]
				rename!(dfEmissions,auxNew_Names)
				total = DataFrame(["Total" zeros(1,inputs["NCO2Cap"]*inputs["NumScenarios"]) [sum(dfEmissions[!,Symbol("AnnualSum_Scenario_$s")]) for s in 1:inputs["NumScenarios"]]' fill(0.0, (1,T))], :auto)
				for t in 1:T
					total[:,t+inputs["NCO2Cap"]*inputs["NumScenarios"]+inputs["NumScenarios"]+1] .= sum(dfEmissions[:,Symbol("t$t")][1:Z])
				end
			end

		else
			if inputs["NumScenarios"] <= 1
				auxNew_Names=[Symbol("Zone"); Symbol("AnnualSum"); [Symbol("t$t") for t in 1:T]]
				rename!(dfEmissions,auxNew_Names)
				total = DataFrame(["Total" sum(dfEmissions[!,:AnnualSum]) fill(0.0, (1,T))], :auto)
				for t in 1:T
					total[:,t+2] .= sum(dfEmissions[:,Symbol("t$t")][1:Z])
				end
			else
				auxNew_Names=[Symbol("Zone");[Symbol("AnnualSum_Scenario_$s") for s=1:inputs["NumScenarios"]];[Symbol("t$t") for t in 1:T]]
				rename!(dfEmissions,auxNew_Names)
				total = DataFrame(["Total" [sum(dfEmissions[!,Symbol("AnnualSum_Scenario_$s")]) for s in 1:inputs["NumScenarios"]]' fill(0.0, (1,T))], :auto)
				for t in 1:T
					total[:,t+inputs["NumScenarios"]+1] .= sum(dfEmissions[:,Symbol("t$t")][1:Z])
				end
			end
		end

        rename!(total,auxNew_Names)
        dfEmissions = vcat(dfEmissions, total)


## Aaron - Combined elseif setup["Dual_MIP"]==1 block with the first block since they were identical. Why do we have this third case? What is different about it?
	else
		if inputs["NumScenarios"]<= 1
			#CO2 emissions by zone
			dfEmissions = hcat(DataFrame(Zone = 1:Z), DataFrame(AnnualSum = Array{Float64}(undef, Z)))
			for i in 1:Z
				dfEmissions[i,:AnnualSum] = sum(inputs["omega"].*value.(EP[:eEmissionsByZone][i,:])) * scale_factor
			end
			dfEmissions = hcat(dfEmissions, DataFrame(value.(EP[:eEmissionsByZone]) * scale_factor, :auto))
			auxNew_Names=[Symbol("Zone");Symbol("AnnualSum");[Symbol("t$t") for t in 1:T]]
			rename!(dfEmissions,auxNew_Names)
			total = DataFrame(["Total" sum(dfEmissions[!,:AnnualSum]) fill(0.0, (1,T))], :auto)
			for t in 1:T
				total[:,t+2] .= sum(dfEmissions[:,Symbol("t$t")][1:Z])
			end
			rename!(total,auxNew_Names)
			dfEmissions = vcat(dfEmissions, total)
		else
			#CO2 emissions by zone
			dfEmissions = DataFrame([1:Z Array{Float64}(undef, Z,inputs["NumScenarios"])],:auto);
			auxNew_Names=[Symbol("Zone");[Symbol("AnnualSum_Scenario_$s") for s in 1: inputs["NumScenarios"]]];
			rename!(dfEmissions,auxNew_Names)
			scenario_times = timeseries2scenarios(T, inputs["NumScenarios"])
			for i in 1:Z
				for s in 1:inputs["NumScenarios"]
					dfEmissions[i,Symbol("AnnualSum_Scenario_$s")] = sum(inputs["omega"][t]*value(EP[:eEmissionsByZone][i,t]) for t in scenario_times[s,1]:scenario_times[s,2]) * scale_factor
				end
			end
			dfEmissions = hcat(dfEmissions, DataFrame(value.(EP[:eEmissionsByZone]) * scale_factor, :auto))
			auxNew_Names=[Symbol("Zone");[Symbol("AnnualSum_Scenario_$s") for s in 1: inputs["NumScenarios"]];[Symbol("t$t") for t in 1:T]]
			rename!(dfEmissions,auxNew_Names)
			total = DataFrame(["Total" [sum(dfEmissions[!,Symbol("AnnualSum_Scenario_$s")]) for s in 1:inputs["NumScenarios"]]' fill(0.0, (1,T))], :auto)
			for t in 1:T
				total[:,t+inputs["NumScenarios"]+1] .= sum(dfEmissions[:,Symbol("t$t")][1:Z])
			end
			rename!(total,auxNew_Names)
			dfEmissions = vcat(dfEmissions, total)
		end
	end
	CSV.write(joinpath(path, "emissions.csv"), dftranspose(dfEmissions, false), writeheader=false)
end
