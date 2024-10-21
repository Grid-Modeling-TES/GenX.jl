@doc raw"""
	write_tes(path::AbstractString, inputs::Dict,setup::Dict, EP::Model)

Function for writing the capacities of different storage technologies, including hydro reservoir, flexible storage tech etc.
"""
function write_tes(path::AbstractString, inputs::Dict, setup::Dict, EP::Model)
    gen = inputs["RESOURCES"]
    price = locational_marginal_price(EP, inputs, setup)

    T = inputs["T"]     # Number of time steps (hours)
    G = inputs["G"]
    TES = inputs["TES"]
    Z = inputs["Z"]

    # Storage level (state of charge) of each resource in each time step
    if !isempty(inputs["TES"])
        storagevcapvalue = zeros(length(TES)*5, T)
        column_counter = 1
        Resource_header = ["TES_Charge_Zone1", "TES_State_of_charge_Zone1", "TES_Production_MWh_Zone1", "TES_MMBtu_Production_Zone1", "TES_Cost_of_Electricity_Zone1"]
        Zones = [1,1,1,1,1]
        for z in 1:Z
            Y_ZONE = resources_in_zone_by_rid(gen, z)
            TES_ZONE = intersect(inputs["TES"], Y_ZONE)

            if !isempty(TES_ZONE)
                storagevcapvalue[column_counter, :] = value.(EP[:vCHARGE_TES][TES_ZONE, :])
                storagevcapvalue[column_counter+1, :] = value.(EP[:vS_TES][TES_ZONE, :])
                storagevcapvalue[column_counter+2, :] = value.(EP[:vUSE_TES][TES_ZONE, :])
                storagevcapvalue[column_counter+3, :] = 0.001 * (value.(EP[:vUSE_TES][TES_ZONE, :]) / (tes_mwh_per_mmbtu.(gen.Tes))[z])
                storagevcapvalue[column_counter+4, :] = ((value.(EP[:vCHARGE_TES][TES_ZONE, :]).data) .* price[z, :])[z, :]
                column_counter += 5
                if z > 1
                    push!(Resource_header, "TES_Charge_Zone$z")
                    push!(Resource_header, "TES_State_of_charge_Zone$z")
                    push!(Resource_header, "TES_Production_Zone$z")
                    push!(Resource_header, "TES_MMBtu_Production_Zone$z")
                    push!(Resource_header, "TES_Cost_of_Electricity_Zone$z")
                    push!(Zones, z)
                    push!(Zones, z)
                    push!(Zones, z)
                    push!(Zones, z)
                    push!(Zones, z)
                end
            end
        end  
        dfStorage = DataFrame(Resource = Resource_header, Zone = Zones, AnnualSum = Array{Union{Missing, Float64}}(undef, (length(TES)*5)))
        
        if setup["ParameterScale"] == 1
            storagevcapvalue *= ModelScalingFactor
        end
    
        dfStorage.AnnualSum .= storagevcapvalue * inputs["omega"]
    
        dfStorage = hcat(dfStorage, DataFrame(storagevcapvalue, :auto))
        auxNew_Names = [Symbol("Resource"); Symbol("Zone"); Symbol("Total"); [Symbol("t$t") for t in 1:T]]
        rename!(dfStorage, auxNew_Names)
        CSV.write(joinpath(path, "tes.csv"), dftranspose(dfStorage, false), header = false)
    end
end
