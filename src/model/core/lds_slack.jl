function lds_slack!(EP::Model, inputs::Dict,setup::Dict)

    println("Including slacks for LDS constraints")

    @variable(EP,vLDS_SLACK_MAX[w in 1:inputs["REP_PERIOD"]]);

	@constraint(EP,cPosSlack[w in 1:inputs["REP_PERIOD"]],vLDS_SLACK_MAX[w]>=0)
    
    # PenaltyValue = (inputs["Weights"]/inputs["H"])*maximum(inputs["dfGen"].Inv_Cost_per_MWyr);

    PenaltyValue = 2*(inputs["Weights"]/inputs["H"])*inputs["Voll"][1] ;    

	@expression(EP,eObjSlack,sum(PenaltyValue[w]*vLDS_SLACK_MAX[w] for w in 1:inputs["REP_PERIOD"]))
    
    EP[:eObj] += eObjSlack

end