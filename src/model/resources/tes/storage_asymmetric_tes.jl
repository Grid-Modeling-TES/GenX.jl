@doc raw"""
	storage_asymmetric_tes!(EP::Model, inputs::Dict, setup::Dict)

Sets up variables and constraints specific to storage resources with asymmetric charge and discharge capacities. See ```storage()``` in ```storage.jl``` for description of constraints.
"""
function storage_asymmetric_tes!(EP::Model, inputs::Dict, setup::Dict)
    # Set up additional variables, constraints, and expressions associated with storage resources with asymmetric charge & discharge capacity
    println("TES Storage Resources with Asmymetric Charge/Discharge Capacity Module")

    OperationalReserves = setup["OperationalReserves"]
    CapacityReserveMargin = setup["CapacityReserveMargin"]

    T = inputs["T"]     # Number of time steps (hours)

    TES = inputs["TES"]

    ### Constraints ###

    # Maximum charging rate must be less than charge power rating
    @constraint(EP,
        [y in TES, t in 1:T],
        EP[:vCHARGE_TES][y, t]<=EP[:eTotalCapCharge_TES][y])
end

