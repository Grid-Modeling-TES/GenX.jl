@doc raw"""
	storage_all_tes!(EP::Model, inputs::Dict, setup::Dict)

Sets up variables and constraints common to all storage resources. See ```storage()``` in ```storage.jl``` for description of constraints.
"""
function storage_all_tes!(EP::Model, inputs::Dict, setup::Dict)
    # Setup variables, constraints, and expressions common to all storage resources
    println("TES Storage Core Resources Module")

    gen = inputs["RESOURCES"]
    OperationalReserves = setup["OperationalReserves"]
    CapacityReserveMargin = setup["CapacityReserveMargin"]

    virtual_discharge_cost = inputs["VirtualChargeDischargeCost"]

    T = inputs["T"]     # Number of time steps (hours)
    Z = inputs["Z"]     # Number of zones
    p = inputs["hours_per_subperiod"]

    TES = inputs["TES"]
    STOR_SHORT_DURATION = inputs["STOR_SHORT_DURATION_TES"]
    representative_periods = inputs["REP_PERIOD"]

    START_SUBPERIODS = inputs["START_SUBPERIODS"]
    INTERIOR_SUBPERIODS = inputs["INTERIOR_SUBPERIODS"]

    hours_per_subperiod = inputs["hours_per_subperiod"] #total number of hours per subperiod

    ### Variables ###

    # Storage level of resource "y" at hour "t" [MWh] on zone "z" - unbounded
    @variable(EP, vS_TES[y in TES, t = 1:T]>=0)

    # Energy withdrawn from grid by resource "y" at hour "t" [MWh] on zone "z"
    @variable(EP, vCHARGE_TES[y in TES, t = 1:T]>=0)

    # TES "discharge" in MWh
    @variable(EP, vUSE_TES[y = TES, t in 1:T]>=0)

    ### Expressions ###

    ## Objective Function Expressions ##

    #Variable costs of "charging" for technologies "y" during hour "t" in zone "z"
    @expression(EP,
        eCVar_in_TES[y in TES, t = 1:T],
        inputs["omega"][t]*var_om_cost_per_mwh_in(gen[y])*vCHARGE_TES[y, t])

    # Sum individual resource contributions to variable charging costs to get total variable charging costs
    @expression(EP, eTotalCVarInT_TES[t = 1:T], sum(eCVar_in_TES[y, t] for y in TES))
    @expression(EP, eTotalCVarIn_TES, sum(eTotalCVarInT_TES[t] for t in 1:T))
    EP[:eObj] += eTotalCVarIn_TES


    # Subtract TES revenue from objective function
    scale_factor = setup["ParameterScale"] == 1 ? 10^6 : 1  # If ParameterScale==1, costs are in millions of $
    @expression(EP,
        eHydrogenValue_TES[y in TES, t in 1:T],
        (inputs["omega"][t] * EP[:vUSE_TES][y, t] / tes_mwh_per_mmbtu(gen[y]) *
        heat_price_per_mmbtu(gen[y])/scale_factor))
    @expression(EP,
        eTotalHydrogenValueT_TES[t in 1:T],
        sum(eHydrogenValue_TES[y, t] for y in TES))
    @expression(EP, eTotalHydrogenValue_TES, sum(eTotalHydrogenValueT_TES[t] for t in 1:T))
    EP[:eObj] -= eTotalHydrogenValue_TES
    

    ## Power Balance Expressions ##
    @expression(EP, ePowerBalanceTES[t in 1:T, z in 1:Z],
        sum(EP[:vCHARGE_TES][y, t] for y in intersect(TES, resources_in_zone_by_rid(gen, z))))

    # TES consumes electricity so their vCHARGE is subtracted from power balance.
    EP[:ePowerBalance] -= ePowerBalanceTES

    
    ### Energy Share Requirement Policy ###
    # Since we're using vUSE to denote TES consumption, we subtract this from the eESR Energy Share Requirement balance to increase demand for clean resources if desired
    # TES demand is only accounted for in an ESR that the TES resources is tagged in in Generates_data.csv (e.g. ESR_N > 0) and
    # a share of TES demand equal to df[y,:ESR_N] must be met by resources qualifying for ESR_N for each TES resource y.
    # Using vCHARGE_TES instead of vUSE
    if setup["EnergyShareRequirement"] >= 1
        @expression(EP,
            eElectrolyzerESR_TES[ESR in 1:inputs["nESR"]],
            sum(inputs["omega"][t] * EP[:vCHARGE_TES][y, t]
            for y in intersect(TES, ids_with_policy(gen, esr, tag = ESR)),
            t in 1:T))
        EP[:eESR] -= eElectrolyzerESR_TES
    end


    ### Constraints ###
    
    # TES always discharges when there's a non-zero charge. This constraint sets the minimum heat production to the current charge divived by the maximum duration
    @constraints(EP, begin
        [y in TES, t in 1:T],
        EP[:vUSE_TES][y, t] >= (EP[:vS_TES][y ,t]/max_duration(gen[y]))
    end)

    # Capacity factor constraint
    @constraints(EP, begin
        [y in TES],
        sum(inputs["omega"][t] * EP[:vUSE_TES][y, t] for t in 1:T) >= sum(inputs["omega"][t]*min_power(gen[y])*EP[:eTotalCap][y] for t in 1:T)
    end)

    ### Remove vP (TES does not produce power so vP = 0 for all periods)
    @constraints(EP, begin
        [y in TES, t in 1:T], EP[:vP][y, t] == 0
    end)


    ## Storage energy capacity and state of charge related constraints:
    @constraint(EP, cCapRatio[y in TES],
        EP[:eTotalCapEnergy_TES][y] == (max_duration(gen[y])*EP[:eTotalCap][y]))


    # Links state of charge in first time step with decisions in last time step of each subperiod
    # We use a modified formulation of this constraint (cSoCBalLongDurationStorageStart_TES) when operations wrapping and long duration storage are being modeled
    if representative_periods > 1 && !isempty(inputs["STOR_LONG_DURATION_TES"])
        CONSTRAINTSET = STOR_SHORT_DURATION
    else
        CONSTRAINTSET = TES
    end
    @constraint(EP,
        cSoCBalStart_TES[t in START_SUBPERIODS, y in CONSTRAINTSET],
        EP[:vS_TES][y,
            t]==
        EP[:vS_TES][y, t + hours_per_subperiod - 1] -
        (1 / efficiency_down(gen[y]) * EP[:vUSE_TES][y, t])
        +
        (EP[:vCHARGE_TES][y, t]) -
        (self_discharge(gen[y]) * EP[:vS_TES][y, t + hours_per_subperiod - 1]))

    @constraints(EP,
        begin

            # Maximum energy stored must be less than energy capacity
            [y in TES, t in 1:T], EP[:vS_TES][y, t] <= EP[:eTotalCapEnergy_TES][y]

            # energy stored for the next hour
            cSoCBalInterior_TES[t in INTERIOR_SUBPERIODS, y in TES],
            EP[:vS_TES][y, t] ==
            EP[:vS_TES][y, t - 1] - (1 / efficiency_down(gen[y]) * EP[:vUSE_TES][y, t]) +
            (efficiency_up(gen[y]) * EP[:vCHARGE_TES][y, t]) -
            (self_discharge(gen[y]) * EP[:vS_TES][y, t - 1])
        end)

    # Storage discharge and charge power (and reserve contribution) related constraints:
    @constraints(EP,
        begin
            [y in TES, t = 1:T], 
            EP[:vUSE_TES][y, t] <= EP[:eTotalCap][y]
            [y in TES, t = 1:T],
            EP[:vUSE_TES][y, t] <=
            EP[:vS_TES][y, hoursbefore(hours_per_subperiod, t, 1)] *
            efficiency_down(gen[y])
        end)

    
    ### Maximum ramp up and down between consecutive hours - Not applicable as ramp rates are 100% between hours for TES

    ### Minimum heat production constraint (if any)
    #kmmbtu_to_mmbtu = 10^3
    #@constraint(EP,
    #    cHydrogenMin_TES[y in TES],
    #    sum(inputs["omega"][t] * EP[:eTesProduction][y, t] / tes_mwh_per_mmbtu(gen[y])
    #    for t in 1:T)>=tes_min_kmmbtu(gen[y]) * kmmbtu_to_mmbtu)


    ### Hydrogen Hourly Supply Matching Constraint (Constraint #6) ###
    # Requires generation from qualified resources (indicated by Qualified_Hydrogen_Supply==1 in the resource .csv files)
    # from within the same zone as the electrolyzers are located to be >= hourly consumption from electrolyzers in the zone
    # (and any charging by qualified storage within the zone used to help increase electrolyzer utilization).
    #if setup["HydrogenHourlyMatching"] == 1
    #    HYDROGEN_ZONES = unique(zone_id.(gen.Electrolyzer))
    #    QUALIFIED_SUPPLY = ids_with(gen, qualified_hydrogen_supply)
    #    @constraint(EP, cHourlyMatching[z in HYDROGEN_ZONES, t in 1:T],
    #        sum(EP[:vP][y, t]
    #        for y in intersect(resources_in_zone_by_rid(gen,z), QUALIFIED_SUPPLY))>=sum(EP[:vUSE][y,t]
    #        for y in intersect(resources_in_zone_by_rid(gen,z), ELECTROLYZERS)) + sum(EP[:vCHARGE][y,t]
    #        for y in intersect(resources_in_zone_by_rid(gen,z), QUALIFIED_SUPPLY, STORAGE)) + sum(EP[:vCHARGE_TES][y,t]
    #        for y in intersect(resources_in_zone_by_rid(gen,z), QUALIFIED_SUPPLY, TES)))
    #end


    # Attemps at constraining production in each hour are below - these tend to make the model non-linear.
    #=     @variable(EP, vUSE_TES[y in TES, t in 1:T]>=0)
    @variable(EP, vd1[y in TES, t in 1:T]>=0)
    @variable(EP, vd2[y in TES, t in 1:T]>=0)
    MMM = 15000

    @constraints(EP, begin
        [y in TES, t in 1:T],
        EP[:vUSE_TES][y, t] <= EP[:vS_TES][y, t]
        [y in TES, t in 1:T],
        EP[:vUSE_TES][y, t] <= (min_power(gen[y]) * EP[:eTotalCap][y])
        [y in TES, t in 1:T],
        EP[:vUSE_TES][y, t] >= (EP[:vS_TES][y, t] - (MMM*(1-vd1[y, t])))
        [y in TES, t in 1:T],
        EP[:vUSE_TES][y, t] >= ((min_power(gen[y]) * EP[:eTotalCap][y]) - (MMM*(1-vd2[y, t])))
        [y in TES, t in 1:T],
        (vd1[y, t] + vd2[y, t]) == 1
    end) =#

#=     @variable(EP, vBinary[y in TES, t in 1:T] in Semiinteger(0, 1))
    MMM = 15000
    @constraints(EP, begin
        [y in TES, t in 1:T],
        ((min_power(gen[y]) * EP[:eTotalCap][y]) -  EP[:vS_TES][y, t]) <= (MMM*vBinary[y, t])

        [y in TES, t in 1:T],
        (EP[:vS_TES][y, t] - (min_power(gen[y]) * EP[:eTotalCap][y])) <= (MMM*(1-vBinary[y, t]))
    end)

    @constraints(EP, begin
        [y in TES, t in 1:T],
        EP[:vUSE_TES][y, t] <= EP[:vS_TES][y, t]
        [y in TES, t in 1:T],
        EP[:vUSE_TES][y, t] <= (min_power(gen[y]) * EP[:eTotalCap][y])
        [y in TES, t in 1:T],
        EP[:vUSE_TES][y, t] >= (EP[:vS_TES][y, t] - (MMM*(1-vBinary[y, t])))
        [y in TES, t in 1:T],
        EP[:vUSE_TES][y, t] >= ((min_power(gen[y]) * EP[:eTotalCap][y]) - (MMM*vBinary[y, t]))
    end) =#

    #@expression(EP, eTesProduction[y in TES, t in 1:T],
    #    if EP[:vS_TES][y, t] <= (min_power(gen[y]) * EP[:eTotalCap][y])
    #        EP[:vS_TES][y, t]
    #    else
    #        min_power(gen[y]) * EP[:eTotalCap][y]
    #    end)
    #min(EP[:vS_TES][y, t],(min_power(gen[y]) * EP[:eTotalCap][y])))
    #@expression(EP, eTesProduction[y in TES, t in 1:T],
    #    ifelse(EP[:vS_TES][y, t] <= (min_power(gen[y]) * EP[:eTotalCap][y]), EP[:vS_TES][y, t], min_power(gen[y]) * EP[:eTotalCap][y]))
    


end
