function integer_investments!(EP::Model, inputs::Dict, setup::Dict)

   
    if setup["NetworkExpansion"] == 1
        set_integer.(EP[:vNEW_TRANS_CAP])
	end
    
    set_integer.(EP[:vRETCAP])
    
    set_integer.(EP[:vCAP])

    if setup["MultiStage"]==1
        if !isempty(inputs[1]["STOR_ALL"])
            set_integer.(EP[:vCAPENERGY])
            set_integer.(EP[:vRETCAPENERGY])
        end
        
        if !isempty(inputs[1]["STOR_ASYMMETRIC"])
            set_integer.(EP[:vCAPCHARGE])
            set_integer.(EP[:vRETCAPCHARGE])
        end
    else
        if !isempty(inputs["STOR_ALL"])
            set_integer.(EP[:vCAPENERGY])
            set_integer.(EP[:vRETCAPENERGY])
        end
        
        if !isempty(inputs["STOR_ASYMMETRIC"])
            set_integer.(EP[:vCAPCHARGE])
            set_integer.(EP[:vRETCAPCHARGE])
        end
    end

end