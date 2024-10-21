@doc raw"""
	investment_energy_tes!(EP::Model, inputs::Dict)

This function defines the expressions and constraints keeping track of total available TES storage charge capacity across all resources as well as constraints on capacity retirements. The function also adds investment and fixed O\&M related costs related to charge capacity to the objective function.

The total capacity of each resource is defined as the sum of the existing capacity plus the newly invested capacity minus any retired capacity.

```math
\begin{aligned}
& \Delta^{total,energy}_{y,z} =(\overline{\Delta^{energy}_{y,z}}+\Omega^{energy}_{y,z}-\Delta^{energy}_{y,z}) \forall y \in \mathcal{O}, z \in \mathcal{Z}
\end{aligned}
```

One cannot retire more capacity than existing capacity.

```math
\begin{aligned}
&\Delta^{energy}_{y,z} \leq \overline{\Delta^{energy}_{y,z}}
		\hspace{4 cm}  \forall y \in \mathcal{O}, z \in \mathcal{Z}
\end{aligned}
```

For resources where $\overline{\Omega_{y,z}^{energy}}$ and $\underline{\Omega_{y,z}^{energy}}$ is defined, then we impose constraints on minimum and maximum power capacity.

```math
\begin{aligned}
& \Delta^{total,energy}_{y,z} \leq \overline{\Omega}^{energy}_{y,z}
	\hspace{4 cm}  \forall y \in \mathcal{O}, z \in \mathcal{Z} \\
& \Delta^{total,energy}_{y,z}  \geq \underline{\Omega}^{energy}_{y,z}
	\hspace{4 cm}  \forall y \in \mathcal{O}, z \in \mathcal{Z}
\end{aligned}
```

In addition, this function adds investment and fixed O\&M related costs related to charge capacity to the objective function:

```math
\begin{aligned}
& 	\sum_{y \in \mathcal{O} } \sum_{z \in \mathcal{Z}}
	\left( (\pi^{INVEST,energy}_{y,z} \times    \Omega^{energy}_{y,z})
	+ (\pi^{FOM,energy}_{y,z} \times  \Delta^{total,energy}_{y,z})\right)
\end{aligned}
```
"""
function investment_energy_tes!(EP::Model, inputs::Dict, setup::Dict)
    println("TES Storage Investment Module")

    gen = inputs["RESOURCES"]

    MultiStage = setup["MultiStage"]

    TES = inputs["TES"] # Set of all storage resources
    NEW_CAP_ENERGY = inputs["NEW_CAP_ENERGY_TES"] # Set of all storage resources eligible for new energy capacity
    RET_CAP_ENERGY = inputs["RET_CAP_ENERGY_TES"] # Set of all storage resources eligible for energy capacity retirements

    ### Variables ###

    # New installed energy capacity of resource "y"
    @variable(EP, vCAPENERGY_TES[y in NEW_CAP_ENERGY]>=0)

    # Retired energy capacity of resource "y" from existing capacity
    @variable(EP, vRETCAPENERGY_TES[y in RET_CAP_ENERGY]>=0)

    if MultiStage == 1
        @variable(EP, vEXISTINGCAPENERGY_TES[y in TES]>=0)
    end

    ### Expressions ###

    if MultiStage == 1
        @expression(EP, eExistingCapEnergy_TES[y in TES], vEXISTINGCAPENERGY_TES[y])
    else
        @expression(EP, eExistingCapEnergy_TES[y in TES], existing_cap_mwh(gen[y]))
    end

    @expression(EP, eTotalCapEnergy_TES[y in TES],
        if (y in intersect(NEW_CAP_ENERGY, RET_CAP_ENERGY))
            eExistingCapEnergy_TES[y] + EP[:vCAPENERGY_TES][y] - EP[:vRETCAPENERGY_TES][y]
        elseif (y in setdiff(NEW_CAP_ENERGY, RET_CAP_ENERGY))
            eExistingCapEnergy_TES[y] + EP[:vCAPENERGY_TES][y]
        elseif (y in setdiff(RET_CAP_ENERGY, NEW_CAP_ENERGY))
            eExistingCapEnergy_TES[y] - EP[:vRETCAPENERGY_TES][y]
        else
            eExistingCapEnergy_TES[y] + EP[:vZERO]
        end)

    ## Objective Function Expressions ##

    # Fixed costs for resource "y" = annuitized investment cost plus fixed O&M costs
    # If resource is not eligible for new energy capacity, fixed costs are only O&M costs
    @expression(EP, eCFixEnergy_TES[y in TES],
        if y in NEW_CAP_ENERGY # Resources eligible for new capacity
            inv_cost_per_mwhyr(gen[y]) * vCAPENERGY_TES[y] +
            fixed_om_cost_per_mwhyr(gen[y]) * eTotalCapEnergy_TES[y]
        else
            fixed_om_cost_per_mwhyr(gen[y]) * eTotalCapEnergy_TES[y]
        end)

    # Sum individual resource contributions to fixed costs to get total fixed costs
    @expression(EP, eTotalCFixEnergy_TES, sum(EP[:eCFixEnergy_TES][y] for y in TES))

    # Add term to objective function expression
    if MultiStage == 1
        # OPEX multiplier scales fixed costs to account for multiple years between two model stages
        # We divide by OPEXMULT since we are going to multiply the entire objective function by this term later,
        # and we have already accounted for multiple years between stages for fixed costs.
        add_to_expression!(EP[:eObj], (1 / inputs["OPEXMULT"]), eTotalCFixEnergy_TES)
    else
        add_to_expression!(EP[:eObj], eTotalCFixEnergy_TES)
    end

    ### Constraints ###

    if MultiStage == 1
        @constraint(EP,
            cExistingCapEnergy_TES[y in TES],
            EP[:vEXISTINGCAPENERGY_TES][y]==existing_cap_mwh(gen[y]))
    end

    ## Constraints on retirements and capacity additions
    # Cannot retire more energy capacity than existing energy capacity
    @constraint(EP,
        cMaxRetEnergy_TES[y in RET_CAP_ENERGY],
        vRETCAPENERGY_TES[y]<=eExistingCapEnergy_TES[y])

    ## Constraints on new built energy capacity
    # Constraint on maximum energy capacity (if applicable) [set input to -1 if no constraint on maximum energy capacity]
    # DEV NOTE: This constraint may be violated in some cases where Existing_Cap_MWh is >= Max_Cap_MWh and lead to infeasabilty
    @constraint(EP,
        cMaxCapEnergy_TES[y in intersect(ids_with_positive(gen, max_cap_mwh), TES)],
        eTotalCapEnergy_TES[y]<=max_cap_mwh(gen[y]))

    # Constraint on minimum energy capacity (if applicable) [set input to -1 if no constraint on minimum energy apacity]
    # DEV NOTE: This constraint may be violated in some cases where Existing_Cap_MWh is <= Min_Cap_MWh and lead to infeasabilty
    @constraint(EP,
        cMinCapEnergy_TES[y in intersect(ids_with_positive(gen, min_cap_mwh), TES)],
        eTotalCapEnergy_TES[y]>=min_cap_mwh(gen[y]))

    # Max and min constraints on energy storage capacity built (as proportion to discharge power capacity)
    @constraint(EP,
        cMinCapEnergyDuration_TES[y in TES],
        EP[:eTotalCapEnergy_TES][y]>=min_duration(gen[y]) * EP[:eTotalCap][y])
    @constraint(EP,
        cMaxCapEnergyDuration_TES[y in TES],
        EP[:eTotalCapEnergy_TES][y]<=max_duration(gen[y]) * EP[:eTotalCap][y])
end
