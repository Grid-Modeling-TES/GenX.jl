function write_heat_prices(path::AbstractString, inputs::Dict, setup::Dict, EP::Model)
    scale_factor = setup["ParameterScale"] == 1 ? 10^6 : 1  # If ParameterScale==1, costs are in millions of $
    dfHeatPrice = DataFrame(Heat_Price_Per_MMBTU = convert(Array{Float64},
        dual.(EP[:cHydrogenMin_TES]) * scale_factor))

    CSV.write(joinpath(path, "heat_prices.csv"), dfHeatPrice)
    return nothing
end
