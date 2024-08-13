@doc raw"""
	investment_charge_tes!(EP::Model, inputs::Dict)

This function defines the expressions and constraints keeping track of total available storage charge capacity across all resources as well as constraints on capacity retirements. The function also adds investment and fixed O\&M related costs related to charge capacity to the objective function.

The total capacity of each resource is defined as the sum of the existing capacity plus the newly invested capacity minus any retired capacity.

```math
\begin{aligned}
& \Delta^{total,charge}_{y,z} =(\overline{\Delta^{charge}_{y,z}}+\Omega^{charge}_{y,z}-\Delta^{charge}_{y,z}) \forall y \in \mathcal{O}^{asym}, z \in \mathcal{Z}
\end{aligned}
```

One cannot retire more capacity than existing capacity.
```math
\begin{aligned}
&\Delta^{charge}_{y,z} \leq \overline{\Delta^{charge}_{y,z}}
	\hspace{4 cm}  \forall y \in \mathcal{O}^{asym}, z \in \mathcal{Z}
\end{aligned}
```

For resources where $\overline{\Omega_{y,z}^{charge}}$ and $\underline{\Omega_{y,z}^{charge}}$ is defined, then we impose constraints on minimum and maximum power capacity.
```math
\begin{aligned}
& \Delta^{total,charge}_{y,z} \leq \overline{\Omega}^{charge}_{y,z}
	\hspace{4 cm}  \forall y \in \mathcal{O}^{asym}, z \in \mathcal{Z} \\
& \Delta^{total,charge}_{y,z}  \geq \underline{\Omega}^{charge}_{y,z}
	\hspace{4 cm}  \forall y \in \mathcal{O}^{asym}, z \in \mathcal{Z}
\end{aligned}
```

In addition, this function adds investment and fixed O&M related costs related to charge capacity to the objective function:
```math
\begin{aligned}
& 	\sum_{y \in \mathcal{O}^{asym} } \sum_{z \in \mathcal{Z}}
	\left( (\pi^{INVEST,charge}_{y,z} \times    \Omega^{charge}_{y,z})
	+ (\pi^{FOM,charge}_{y,z} \times  \Delta^{total,charge}_{y,z})\right)
\end{aligned}
```
"""
function investment_charge_tes!(EP::Model, inputs::Dict, setup::Dict)
    println("TES Charge Investment Module")

    gen = inputs["RESOURCES"]

    MultiStage = setup["MultiStage"]

    TES = inputs["TES"] # Set of storage resources with asymmetric (separte) charge/discharge capacity components

    NEW_CAP_CHARGE = inputs["NEW_CAP_CHARGE_TES"] # Set of asymmetric charge/discharge storage resources eligible for new charge capacity
    RET_CAP_CHARGE = inputs["RET_CAP_CHARGE_TES"] # Set of asymmetric charge/discharge storage resources eligible for charge capacity retirements

    ### Variables ###

    ## Storage capacity built and retired for storage resources with independent charge and discharge power capacities (STOR=2)

    # New installed charge capacity of resource "y"
    @variable(EP, vCAPCHARGE_TES[y in NEW_CAP_CHARGE]>=0)

    # Retired charge capacity of resource "y" from existing capacity
    @variable(EP, vRETCAPCHARGE_TES[y in RET_CAP_CHARGE]>=0)

    if MultiStage == 1
        @variable(EP, vEXISTINGCAPCHARGE_TES[y in TES]>=0)
    end

    ### Expressions ###

    if MultiStage == 1
        @expression(EP, eExistingCapCharge_TES[y in TES], vEXISTINGCAPCHARGE_TES[y])
    else
        @expression(EP,
            eExistingCapCharge_TES[y in TES],
            existing_charge_cap_mw(gen[y]))
    end

    @expression(EP, eTotalCapCharge_TES[y in TES],
        if (y in intersect(NEW_CAP_CHARGE, RET_CAP_CHARGE))
            eExistingCapCharge_TES[y] + EP[:vCAPCHARGE_TES][y] - EP[:vRETCAPCHARGE_TES][y]
        elseif (y in setdiff(NEW_CAP_CHARGE, RET_CAP_CHARGE))
            eExistingCapCharge_TES[y] + EP[:vCAPCHARGE_TES][y]
        elseif (y in setdiff(RET_CAP_CHARGE, NEW_CAP_CHARGE))
            eExistingCapCharge_TES[y] - EP[:vRETCAPCHARGE_TES][y]
        else
            eExistingCapCharge_TES[y] + EP[:vZERO]
        end)

    ## Objective Function Expressions ##

    # Fixed costs for resource "y" = annuitized investment cost plus fixed O&M costs
    # If resource is not eligible for new charge capacity, fixed costs are only O&M costs
    @expression(EP, eCFixCharge_TES[y in TES],
        if y in NEW_CAP_CHARGE # Resources eligible for new charge capacity
            inv_cost_charge_per_mwyr(gen[y]) * vCAPCHARGE_TES[y] +
            fixed_om_cost_charge_per_mwyr(gen[y]) * eTotalCapCharge_TES[y]
        else
            fixed_om_cost_charge_per_mwyr(gen[y]) * eTotalCapCharge_TES[y]
        end)

    # Sum individual resource contributions to fixed costs to get total fixed costs
    @expression(EP, eTotalCFixCharge_TES, sum(EP[:eCFixCharge_TES][y] for y in TES))

    # Add term to objective function expression
    if MultiStage == 1
        # OPEX multiplier scales fixed costs to account for multiple years between two model stages
        # We divide by OPEXMULT since we are going to multiply the entire objective function by this term later,
        # and we have already accounted for multiple years between stages for fixed costs.
        add_to_expression!(EP[:eObj], (1 / inputs["OPEXMULT"]), eTotalCFixCharge_TES)
    else
        add_to_expression!(EP[:eObj], eTotalCFixCharge_TES)
    end

    ### Constratints ###
    # Charge capacity for TES is 3x the discharge capacity
    @constraint(EP, cChargeCapRatio[y in TES],
    EP[:eTotalCapCharge_TES][y] == (3*EP[:eTotalCap][y]))
    
    if MultiStage == 1
        # Existing capacity variable is equal to existing capacity specified in the input file
        @constraint(EP,
            cExistingCapCharge_TES[y in TES],
            EP[:vEXISTINGCAPCHARGE_TES][y]==existing_charge_cap_mw(gen[y]))
    end

    ## Constraints on retirements and capacity additions
    #Cannot retire more charge capacity than existing charge capacity
    @constraint(EP,
        cMaxRetCharge_TES[y in RET_CAP_CHARGE],
        vRETCAPCHARGE_TES[y]<=eExistingCapCharge_TES[y])

    #Constraints on new built capacity

    # Constraint on maximum charge capacity (if applicable) [set input to -1 if no constraint on maximum charge capacity]
    # DEV NOTE: This constraint may be violated in some cases where Existing_Charge_Cap_MW is >= Max_Charge_Cap_MWh and lead to infeasabilty
    @constraint(EP,
        cMaxCapCharge_TES[y in intersect(ids_with_positive(gen, max_charge_cap_mw),
            TES)],
        eTotalCapCharge_TES[y]<=max_charge_cap_mw(gen[y]))

    # Constraint on minimum charge capacity (if applicable) [set input to -1 if no constraint on minimum charge capacity]
    # DEV NOTE: This constraint may be violated in some cases where Existing_Charge_Cap_MW is <= Min_Charge_Cap_MWh and lead to infeasabilty
    @constraint(EP,
        cMinCapCharge_TES[y in intersect(ids_with_positive(gen, min_charge_cap_mw),
            TES)],
        eTotalCapCharge_TES[y]>=min_charge_cap_mw(gen[y]))
end
