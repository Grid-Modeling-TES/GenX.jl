@doc raw"""
	storage_asymmetric_tes!(EP::Model, inputs::Dict, setup::Dict)

Sets up variables and constraints specific to storage resources with asymmetric charge and discharge capacities. See ```storage()``` in ```storage.jl``` for description of constraints.
"""
function storage_asymmetric_tes!(EP::Model, inputs::Dict, setup::Dict)
    # Set up additional variables, constraints, and expressions associated with storage resources with asymmetric charge & discharge capacity
    # (e.g. most chemical, thermal, and mechanical storage options with distinct charge & discharge components/processes)
    # STOR = 2 corresponds to storage with distinct power and energy capacity decisions and distinct charge and discharge power capacity decisions/ratings

    println("TES Storage Resources with Asmymetric Charge/Discharge Capacity Module")

    OperationalReserves = setup["OperationalReserves"]
    CapacityReserveMargin = setup["CapacityReserveMargin"]

    T = inputs["T"]     # Number of time steps (hours)

    TES = inputs["TES"]

    ### Constraints ###

    # Storage discharge and charge power (and reserve contribution) related constraints for symmetric storage resources:
    if OperationalReserves == 1
        storage_asymmetric_operational_reserves_tes!(EP, inputs, setup)
    else
        if CapacityReserveMargin > 0
            # Maximum charging rate (including virtual charging to move energy held in reserve back to available storage) must be less than charge power rating
            @constraint(EP,
                [y in TES, t in 1:T],
                EP[:vCHARGE_TES][y, t] + EP[:vCAPRES_charge_TES][y, t]<=EP[:eTotalCapCharge_TES][y])
        else
            # Maximum charging rate (including virtual charging to move energy held in reserve back to available storage) must be less than charge power rating
            @constraint(EP,
                [y in TES, t in 1:T],
                EP[:vCHARGE_TES][y, t]<=EP[:eTotalCapCharge_TES][y])
        end
    end
end

@doc raw"""
	storage_asymmetric_operational_reserves_tes!(EP::Model, inputs::Dict)

Sets up variables and constraints specific to storage resources with asymmetric charge and discharge capacities when reserves are modeled. See ```storage()``` in ```storage.jl``` for description of constraints.
"""
function storage_asymmetric_operational_reserves_tes!(EP::Model, inputs::Dict, setup::Dict)
    T = inputs["T"]
    CapacityReserveMargin = setup["CapacityReserveMargin"] > 0

    TES = inputs["TES"]
    STOR_ASYM_REG = intersect(TES, inputs["REG"]) # Set of asymmetric storage resources with REG reserves

    vCHARGE_TES = EP[:vCHARGE_TES]
    vREG_charge = EP[:vREG_charge]
    eTotalCapCharge_TES = EP[:eTotalCapCharge_TES]

    expr = extract_time_series_to_expression(vCHARGE_TES, TES)
    add_similar_to_expression!(expr[STOR_ASYM_REG, :], vREG_charge[STOR_ASYM_REG, :])
    if CapacityReserveMargin
        vCAPRES_charge_TES = EP[:vCAPRES_charge_TES]
        add_similar_to_expression!(expr[TES, :],
            vCAPRES_charge_TES[TES, :])
    end
    @constraint(EP, [y in TES, t in 1:T], expr[y, t]<=eTotalCapCharge_TES[y])
end
