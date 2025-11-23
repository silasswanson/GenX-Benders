
function benders_decomposition(setup::Dict,benders_dict::Dict)

	EP_master = benders_dict["EP_master"];
	master_vars = benders_dict["master_vars"];

	EP_subprob = benders_dict["EP_subprob"];
	master_vars_sub = benders_dict["master_vars_sub"];

	if setup["MultiStage"]==1 && setup["BD_Mode"]=="full"
		setup["BD_Mode"] = "distributed";
	end

	
	EP_master, master_sol, LB_hist, UB_hist, cpu_time,feasibility_hist = benders_iterations(EP_master,master_vars,EP_subprob, master_vars_sub,setup);
		
	

	return (master_sol,EP_master,LB_hist=LB_hist,UB_hist=UB_hist,cpu_time=cpu_time,feasibility_hist = feasibility_hist)

end



function benders_iterations(EP_master::Model,master_vars::Vector{String},EP_subprob, master_vars_sub,setup::Dict)
	
	### This is a stabilised version of the Benders decomposition algorithm in:
	### A. Jacobson, F. Pecci, N. Sepulveda, Q. Xu, and J. Jenkins, “A computationally efficient Benders decomposition for energy systems planning problems with detailed operations and time-coupling constraints.” arXiv, Feb. 20, 2023. Accessed: Mar. 07, 2023. [Online]. Available: http://arxiv.org/abs/2302.10037

	## Start solver time
	solver_start_time = time()

	#### Algorithm parameters:
	
	MaxIter = setup["BD_MaxIter"]
    ConvTol = setup["BD_ConvTol"]
	MaxCpuTime = setup["BD_MaxCpuTime"]
	γ = setup["BD_StabParam"];
	stab_method = setup["BD_Stab_Method"];

	integer_routine_flag = false

	if setup["IntegerInvestments"] == 1 && setup["BD_Stab_Method"] != "off"
		all_master_vars = all_variables(EP_master);
		integer_vars = all_master_vars[is_integer.(all_master_vars)];
		binary_vars = all_master_vars[is_binary.(all_master_vars)];
		unset_integer.(integer_vars)
		unset_binary.(binary_vars)
		integer_routine_flag = true;
	end

    #### Initialize UB and LB
	master_sol = solve_master_problem(EP_master,master_vars);

    UB = Inf;
    LB = master_sol.LB;

    LB_hist = Float64[];
    UB_hist = Float64[];
    cpu_time = Float64[];
	feasibility_hist = Float64[];

	master_sol_best = deepcopy(master_sol);

    #### Run Benders iterations
    for k = 0:MaxIter
		
		start_subop_sol = time();
        if setup["BD_Mode"]=="full"
            subop_sol = solve_full_operational_subproblem(EP_subprob,master_sol,master_vars_sub);
        elseif setup["BD_Mode"]=="serial"
            subop_sol = solve_decomp_operational_subproblems(EP_subprob,master_sol,master_vars_sub);
        elseif setup["BD_Mode"]=="distributed"
            subop_sol = solve_dist_helpers(EP_subprob,master_sol);
        end
		cpu_subop_sol = time()-start_subop_sol;
		println("Solving the subproblems required $cpu_subop_sol seconds")

		UBnew = sum((subop_sol[w].theta_coeff==0 ? Inf : subop_sol[w].op_cost) for w in keys(subop_sol))+master_sol.inv_cost;
		if UBnew < UB
			master_sol_best = deepcopy(master_sol);
			UB = UBnew;
		end

		print("Updating the master problem....")
		time_start_update = time()
		if setup["BD_Mode"]=="full"
			update_master_problem_single_cut!(EP_master,subop_sol,master_sol,master_vars_sub)
		elseif setup["BD_Mode"]=="serial" || setup["BD_Mode"]=="distributed"
			update_master_problem_multi_cuts!(EP_master,subop_sol,master_sol,master_vars_sub)
		end
		time_master_update = time()-time_start_update
		println("done (it took $time_master_update s).")

		start_master_sol = time()
		unst_master_sol = solve_master_problem(EP_master,master_vars);
		cpu_master_sol = time()-start_master_sol;
		println("Solving the master problem required $cpu_master_sol seconds")

		LB = max(LB,unst_master_sol.LB);
		
		append!(LB_hist,LB)
        append!(UB_hist,UB)
		append!(feasibility_hist,sum(subop_sol[w].feasibility_slack for w in keys(subop_sol)))
        append!(cpu_time,time()-solver_start_time)

		if any(subop_sol[w].theta_coeff==0 for w in keys(subop_sol))
			println("***k = ", k,"      LB = ", LB,"     UB = ", UB,"       Gap = ", (UB-LB)/abs(LB),"       CPU Time = ",cpu_time[end])
		else
			println("k = ", k,"      LB = ", LB,"     UB = ", UB,"       Gap = ", (UB-LB)/abs(LB),"       CPU Time = ",cpu_time[end])
		end

        if (UB-LB)/abs(LB) <= ConvTol
			if integer_routine_flag
				println("*** Switching on integer constraints *** ")
				UB = Inf;
				set_integer.(integer_vars)
				set_binary.(binary_vars)
				master_sol = solve_master_problem(EP_master,master_vars);
				LB = master_sol.LB;
				master_sol_best = deepcopy(master_sol);
				integer_routine_flag = false;
			else
				break
			end
		elseif (cpu_time[end] >= MaxCpuTime)|| (k == MaxIter)
			break
		elseif UB==Inf
			master_sol = deepcopy(unst_master_sol);
		else
			if stab_method == "int_level_set"
				start_stab_method = time()
				println("Solving the interior level set problem with γ = $γ")
				if  setup["IntegerInvestments"] == 1 && integer_routine_flag==false
					unset_integer.(integer_vars)
					unset_binary.(binary_vars)
					for v in integer_vars
						fix(v,unst_master_sol.values[name(v)];force=true)
					end
					for v in binary_vars
						fix(v,unst_master_sol.values[name(v)];force=true)
					end
					master_sol = solve_int_level_set_problem(EP_master,master_vars,unst_master_sol,LB,UB,γ);
					unfix.(integer_vars)
					unfix.(binary_vars)
					set_integer.(integer_vars)
					set_binary.(binary_vars)
					set_lower_bound.(integer_vars,0.0)
					set_lower_bound.(binary_vars,0.0)
				else
					master_sol = solve_int_level_set_problem(EP_master,master_vars,unst_master_sol,LB,UB,γ);
				end
				cpu_stab_method = time()-start_stab_method;
				println("Solving the interior level set problem required $cpu_stab_method seconds")
			elseif stab_method == "l2_level_set"
				start_stab_method = time()
				println("Solving the l2 level set problem with γ = $γ")
				master_sol = solve_l2_level_set_problem(EP_master,master_vars,unst_master_sol,master_sol_best,LB,UB,γ);
				cpu_stab_method = time()-start_stab_method;
				println("Solving the l2 level set problem required $cpu_stab_method seconds")
			elseif occursin("trust_region",stab_method)
				start_stab_method = time()
				println("Solving the trust region problem with γ = $γ")
				master_sol = solve_trust_region_problem(EP_master,master_vars,unst_master_sol,master_sol_best,γ);
				cpu_stab_method = time()-start_stab_method;
				println("Solving the trust region problem required $cpu_stab_method seconds")
			else
				master_sol = deepcopy(unst_master_sol);
			end

		end

    end

	return (EP_master=EP_master,master_sol = master_sol_best,LB_hist = LB_hist,UB_hist = UB_hist,cpu_time = cpu_time,feasibility_hist = feasibility_hist)
end


function update_master_problem_single_cut!(EP::Model,subop_sol::Dict,master_sol::NamedTuple,master_vars_sub::Vector{String})
    
    @constraint(EP,subop_sol[1].theta_coeff*EP[:vTHETA][1] >= subop_sol[1].op_cost + sum(subop_sol[1].lambda[i]*(variable_by_name(EP,master_vars_sub[i]) - master_sol.values[master_vars_sub[i]]) for i in 1:length(master_vars_sub)));  

end

function update_master_problem_multi_cuts!(EP::Model,subop_sol::Dict,master_sol::NamedTuple,master_vars_sub::Dict)
    
	W = keys(subop_sol);
	
    @constraint(EP,[w in W],subop_sol[w].theta_coeff*EP[:vTHETA][w] >= subop_sol[w].op_cost + sum(subop_sol[w].lambda[i]*(variable_by_name(EP,master_vars_sub[w][i]) - master_sol.values[master_vars_sub[w][i]]) for i in 1:length(master_vars_sub[w])));

       
end
