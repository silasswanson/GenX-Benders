

function solve_full_operational_subproblem(EP::Model,master_sol::NamedTuple,master_vars_sub::Vector{String})

    sub_results = Dict();
	
	sub_results[1] = solve_subproblem(EP, master_sol, master_vars_sub);

    return sub_results 

end

function solve_decomp_operational_subproblems(EP::Dict,master_sol::NamedTuple,master_vars_sub::Dict)

    sub_results = Dict();
    for w in keys(EP)
		sub_results[w] = solve_subproblem(EP[w],master_sol,master_vars_sub[w]);
    end
    
    return sub_results

end

function solve_dist_helpers(EP_helpers::DArray{Dict{Any, Any}, 1, Vector{Dict{Any, Any}}},master_sol::NamedTuple)

    p_id = workers();
    np_id = length(p_id);

    sub_results = [Dict() for k in 1:np_id];

    @sync for k in 1:np_id
              @async sub_results[k]= @fetchfrom p_id[k] solve_local_helper(localpart(EP_helpers),master_sol); ### This is equivalent to fetch(@spawnat p .....)
    end

	sub_results = merge(sub_results...);

    return sub_results
end
 
function solve_local_helper(helper_local::Vector{Dict{Any,Any}},master_sol::NamedTuple)

    local_sol=Dict();
    for m in helper_local
        EP = m["Model"];
        master_vars_sub = m["master_vars_sub"]
        w = m["SubPeriod"];
		local_sol[w] = solve_subproblem(EP,master_sol,master_vars_sub);
    end
    return local_sol
end


function solve_master_problem(EP::Model,master_vars::Vector{String})
	
	if any(is_integer.(all_variables(EP)))
		println("The master model is a MILP")
		optimize!(EP)
			if has_values(EP) #
				master_sol =  (LB = objective_value(EP), inv_cost =value(EP[:eObj]), values =Dict([s=>value.(variable_by_name(EP,s)) for s in master_vars]), theta = value.(EP[:vTHETA])) 
			else
				compute_conflict!(EP)
				list_of_conflicting_constraints = ConstraintRef[];
				for (F, S) in list_of_constraint_types(EP)
					for con in all_constraints(EP, F, S)
						if get_attribute(con, MOI.ConstraintConflictStatus()) == MOI.IN_CONFLICT
							push!(list_of_conflicting_constraints, con)
						end
					end
				end
				display(list_of_conflicting_constraints)
				@error "The master solution failed. This should not happen"
			end
	else 
		### The master model is an LP
		optimize!(EP)
		if has_values(EP)
			neg_cap_bool = check_negative_capacities(EP);
			
			if neg_cap_bool
				println("***Resolving the master problem with Crossover=1 because of negative capacities***")
				set_attribute(EP, "Crossover", 1)
				#set_attribute(EP, "BarHomogeneous", 1)
				optimize!(EP)
				if has_values(EP)
					master_sol =  (LB = objective_value(EP), inv_cost =value(EP[:eObj]), values =Dict([s=>value.(variable_by_name(EP,s)) for s in master_vars]), theta = value.(EP[:vTHETA])) 
					set_attribute(EP, "Crossover", 0)
					#set_attribute(EP, "BarHomogeneous", -1)
				else			
					println("The master problem solution failed, trying with BarHomogenous=1")
					set_attribute(EP, "BarHomogeneous", 1)
					optimize!(EP)
					if has_values(EP)
						master_sol =  (LB = objective_value(EP), inv_cost =value(EP[:eObj]), values =Dict([s=>value.(variable_by_name(EP,s)) for s in master_vars]), theta = value.(EP[:vTHETA])) 
						set_attribute(EP, "BarHomogeneous", -1)
					else
						@error "The master solution failed. This should not happen"
					end
				end
			else
				master_sol =  (LB = objective_value(EP), inv_cost =value(EP[:eObj]), values =Dict([s=>value.(variable_by_name(EP,s)) for s in master_vars]), theta = value.(EP[:vTHETA])) 
			end
		else
			println("The master problem solution failed, trying with BarHomogenous=1")
			set_attribute(EP, "BarHomogeneous", 1)
			optimize!(EP)
			if has_values(EP)
				master_sol =  (LB = objective_value(EP), inv_cost =value(EP[:eObj]), values =Dict([s=>value.(variable_by_name(EP,s)) for s in master_vars]), theta = value.(EP[:vTHETA])) 
				set_attribute(EP, "BarHomogeneous", -1)
			else
				@error "The master solution failed. This should not happen"
			end

		end
	end

	return master_sol

end

function solve_subproblem(EP::Model,master_sol::NamedTuple,master_vars_sub::Vector{String})

	
	fix_master_vars!(EP,master_sol,master_vars_sub)

	optimize!(EP)
	
	if has_values(EP)
		op_cost = objective_value(EP);
		lambda = [dual(FixRef(variable_by_name(EP,y))) for y in master_vars_sub];
		theta_coeff = 1;
		if haskey(EP,:eObjSlack)
			feasibility_slack = value(EP[:eObjSlack]);
		else
			feasibility_slack = 0.0;
		end
		
	else
		op_cost = 0;
		lambda = zeros(length(master_vars_sub));
		theta_coeff = 0;
		feasibility_slack = 0;
		@warn "The subproblem solution failed. This should not happen, double check the input files"
	end
    
	return (op_cost=op_cost,lambda = lambda,theta_coeff=theta_coeff,feasibility_slack=feasibility_slack)

end
