# Three Zones with Thermal Energy Storage

This is a one-year example with hourly resolution which contains three zones representing Massachusetts, Connecticut, and Maine. It is designed to show thermal energy storage module in GenX. The sixteen represented resources include natural gas, solar PV, wind, TES and lithium-ion battery storage.

To run the model, first navigate to the example directory:

- Using a Julia REPL:

```bash
$ julia
julia> cd("example_systems/2b_three_zones_w_tes/")
```

- Using a terminal or command prompt:
```bash
$ cd example_systems/2b_three_zones_w_tes/
``` 
   
Next, ensure that your settings in `settings/GenX_settings.yml` are correct. The default settings use the solver `HiGHS`, time domain reduced input data (`TimeDomainReduction: 1`) and minimum capacity requirement policy (`MinCapReq: 1`) as specified in the `policies/Minimum_capacity_requirement.csv` file. Other optional policies include a capacity reserve margin, an energy share requirement (such as renewable portfolio standard (RPS) or clean electricity standard (CES) policies), a CO2 emissions cap, and a maximum capacity requirement policy (see the documentation for more details). You can also enforce the "three pillars" constraint (turned off by default) for TES by setting (`TESHourlyMatching: 1`).

Once the settings are confirmed, run the model with the `Run.jl` script in the example directory:

- Using a Julia REPL (recommended)
```julia
julia> include("Run.jl")
```
- Using a terminal or command prompt:
```bash
$ julia Run.jl
```

Once the model has completed, results will write to the `results` directory.

