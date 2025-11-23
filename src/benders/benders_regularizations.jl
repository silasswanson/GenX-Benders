function solve_l2_level_set_problem(EP::Model,master_vars::Vector{String},master_sol::NamedTuple,master_sol_best::NamedTuple,LB,UB,γ)

	m = length(master_vars);

	@constraint(EP,cLevel_set,EP[:eObj] + sum(EP[:vTHETA])<=LB+γ*(UB-LB))

	# eQuadObjFunL2 = QuadExpr(AffExpr(0.0));
	# for i in 1:m
	# 	add_to_expression!(eQuadObjFunL2,(variable_by_name(EP,master_vars[i])-master_sol_best.values[master_vars[i]]),(variable_by_name(EP,master_vars[i])-master_sol_best.values[master_vars[i]]))
	# end
	eQuadObjFunL2 = sum((variable_by_name(EP,s)-master_sol_best.values[s])^2 for s in master_vars)/m;

	@objective(EP,Min,eQuadObjFunL2)

    optimize!(EP)

	
	if has_values(EP)
		neg_cap_bool = check_negative_capacities(EP);
		
		if neg_cap_bool
			println("***Resolving the l2 problem because of negative capacities***")

			@constraint(EP,cPosTotalCap,EP[:eTotalCap].>=0)
			if haskey(EP,:eTotalCapEnergy)
				@constraint(EP,cPosTotalCapEnergy,EP[:eTotalCapEnergy].>=0)
			end
			if haskey(EP,:eTotalCapCharge)
				@constraint(EP,cPosTotalCapCharge,EP[:eTotalCapCharge].>=0)
			end
			if haskey(EP,:eAvail_Trans_Cap)
				@constraint(EP,cPosAvailTransCap,EP[:eAvail_Trans_Cap].>=0)
			end
			optimize!(EP)

			if has_values(EP)
				neg_cap_bool = check_negative_capacities(EP);
				if neg_cap_bool
					@warn  "Resolving did not work, skipping l2 regularization step"
					delete.(EP,EP[:cPosTotalCap])
					unregister(EP,:cPosTotalCap)
					if haskey(EP,:eTotalCapEnergy)
						delete.(EP,EP[:cPosTotalCapEnergy])
						unregister(EP,:cPosTotalCapEnergy)
					end
					if haskey(EP,:eTotalCapCharge)
						delete.(EP,EP[:cPosTotalCapCharge])
						unregister(EP,:cPosTotalCapCharge)
					end
					if haskey(EP,:eAvail_Trans_Cap)
						delete.(EP,EP[:cPosAvailTransCap])
						unregister(EP,:cPosAvailTransCap)
					end
				else
					master_sol = (;master_sol..., inv_cost=value(EP[:eObj]), values=Dict([s=>value(variable_by_name(EP,s)) for s in master_vars]), theta = value.(EP[:vTHETA]))
					delete.(EP,EP[:cPosTotalCap])
					unregister(EP,:cPosTotalCap)
					if haskey(EP,:eTotalCapEnergy)
						delete.(EP,EP[:cPosTotalCapEnergy])
						unregister(EP,:cPosTotalCapEnergy)
					end
					if haskey(EP,:eTotalCapCharge)
						delete.(EP,EP[:cPosTotalCapCharge])
						unregister(EP,:cPosTotalCapCharge)
					end
					if haskey(EP,:eAvail_Trans_Cap)
						delete.(EP,EP[:cPosAvailTransCap])
						unregister(EP,:cPosAvailTransCap)
					end
				end
			else
				@warn  "The l2 problem solution failed"
				delete.(EP,EP[:cPosTotalCap])
				unregister(EP,:cPosTotalCap)
				if haskey(EP,:eTotalCapEnergy)
					delete.(EP,EP[:cPosTotalCapEnergy])
					unregister(EP,:cPosTotalCapEnergy)
				end
				if haskey(EP,:eTotalCapCharge)
					delete.(EP,EP[:cPosTotalCapCharge])
					unregister(EP,:cPosTotalCapCharge)
				end
				if haskey(EP,:eAvail_Trans_Cap)
					delete.(EP,EP[:cPosAvailTransCap])
					unregister(EP,:cPosAvailTransCap)
				end
			end
		else
			master_sol = (;master_sol..., inv_cost=value(EP[:eObj]), values=Dict([s=>value(variable_by_name(EP,s)) for s in master_vars]), theta = value.(EP[:vTHETA]))
		end

	else
		@warn  "The l2 problem solution failed"
	end



	delete(EP,EP[:cLevel_set])
	unregister(EP,:cLevel_set)
	
	@objective(EP,Min, EP[:eObj] + sum(EP[:vTHETA]))

	return master_sol

end


function solve_int_level_set_problem(EP::Model,master_vars::Vector{String},master_sol::NamedTuple,LB,UB,γ)
	
	m = length(master_vars);

	@constraint(EP,cLevel_set,EP[:eObj] + sum(EP[:vTHETA])<=LB+γ*(UB-LB))

	@objective(EP,Min, 0*sum(EP[:vTHETA]))

    optimize!(EP)

	if has_values(EP)
		
		master_sol = (;master_sol..., inv_cost=value(EP[:eObj]), values=Dict([s=>value(variable_by_name(EP,s)) for s in master_vars]), theta = value.(EP[:vTHETA]))

	else

		if !has_values(EP)
			@warn  "the interior level set problem solution failed"
		else
			master_sol = (;master_sol..., inv_cost=value(EP[:eObj]), values=Dict([s=>value(variable_by_name(EP,s)) for s in master_vars]), theta = value.(EP[:vTHETA]))
		end
	end


	delete(EP,EP[:cLevel_set])
	unregister(EP,:cLevel_set)
	@objective(EP,Min, EP[:eObj] + sum(EP[:vTHETA]))
	
	return master_sol

end




function solve_trust_region_problem(EP::Model,master_vars::Vector{String},master_sol::NamedTuple,master_sol_best::NamedTuple,γ::Float64)
			
	m = length(master_vars);
	@constraint(EP,cTrustRegion1[i=1:m], variable_by_name(EP,master_vars[i]) - master_sol_best.values[master_vars[i]]  <= max(γ*abs(master_sol_best.values[master_vars[i]]),0.5)) ### γ*maximum(abs(master_sol.values[s]) for s in master_vars)
	@constraint(EP,cTrustRegion2[i=1:m], -variable_by_name(EP,master_vars[i]) + master_sol_best.values[master_vars[i]]  <= max(γ*abs(master_sol_best.values[master_vars[i]]),0.5))

	optimize!(EP)

	
	if has_values(EP)
		neg_cap_bool = check_negative_capacities(EP);
		
		if neg_cap_bool
			println("***Resolving the trust region problem with Crossover=1 because of negative capacities***")
			set_attribute(EP, "Crossover", 1)
			#set_attribute(EP, "BarHomogeneous", 1)
			optimize!(EP)
			if has_values(EP)
				neg_cap_bool = check_negative_capacities(EP);
				if neg_cap_bool
					@warn  "Resolving did not work, skipping the trust region step"
					set_attribute(EP, "Crossover", 0)
					sol = deepcopy(master_sol)	
				else
					sol = (inv_cost=value(EP[:eObj]), values=Dict([s=>value(variable_by_name(EP,s)) for s in master_vars]), theta = value.(EP[:vTHETA]))
					set_attribute(EP, "Crossover", 0)
					#set_attribute(EP, "BarHomogeneous", -1)
				end
			else
				@warn  "The trust region problem solution failed"
				set_attribute(EP, "Crossover", 0)
				sol = deepcopy(master_sol)	
			end
		else
			sol = (inv_cost=value(EP[:eObj]), values=Dict([s=>value(variable_by_name(EP,s)) for s in master_vars]), theta = value.(EP[:vTHETA]))
		end

	else
			@warn  "The trust region problem solution failed"
			sol = deepcopy(master_sol)	
	end

	delete(EP,EP[:cTrustRegion1])
	unregister(EP,:cTrustRegion1)
	delete(EP,EP[:cTrustRegion2])
	unregister(EP,:cTrustRegion2)

	return 	sol

end
