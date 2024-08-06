@doc raw"""
	write_tes(path::AbstractString, inputs::Dict,setup::Dict, EP::Model)

Function for writing the capacities of different storage technologies, including hydro reservoir, flexible storage tech etc.
"""
function write_tes(path::AbstractString, inputs::Dict, setup::Dict, EP::Model)
    gen = inputs["RESOURCES"]
    zones = zone_id.(gen)

    T = inputs["T"]     # Number of time steps (hours)
    G = inputs["G"]
    TES = inputs["TES"]

    # Storage level (state of charge) of each resource in each time step
    dfStorage = DataFrame(Resource = inputs["RESOURCE_NAMES"], Zone = zones)
    storagevcapvalue = zeros(G, T)

    if !isempty(inputs["TES"])
        storagevcapvalue[TES, :] = value.(EP[:vS_TES][TES, :])
    end

    if setup["ParameterScale"] == 1
        storagevcapvalue *= ModelScalingFactor
    end

    dfStorage = hcat(dfStorage, DataFrame(storagevcapvalue, :auto))
    auxNew_Names = [Symbol("Resource"); Symbol("Zone"); [Symbol("t$t") for t in 1:T]]
    rename!(dfStorage, auxNew_Names)
    CSV.write(joinpath(path, "tes.csv"), dftranspose(dfStorage, false), header = false)
end
