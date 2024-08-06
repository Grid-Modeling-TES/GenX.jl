@doc raw"""
	long_duration_storage_tes!(EP::Model, inputs::Dict, setup::Dict)
This function creates variables and constraints enabling modeling of long duration storage resources when modeling representative time periods.\

**Storage inventory balance at beginning of each representative period**
The constraints in this section are used to approximate the behavior of long-duration energy storage technologies when approximating annual grid operations by modeling operations over representative periods. Previously, the state of charge balance for storage (as defined in ```storage_all_tes()```) assumed that state of charge at the beginning and end of each representative period has to be the same. In other words, the amount of energy built up or consumed by storage technology $o$ in zone $z$ over the representative period $m$, $\Delta Q_{o,z,m} = 0$. This assumption implicitly excludes the possibility of transferring energy from one representative period to the other which could be cost-optimal when the capital cost of energy storage capacity is relatively small. To model long-duration energy storage using representative periods, we replace the state of charge equation, such that the first term on the right hand side accounts for change in storage inventory associated with representative period $m$ ($\Delta Q_{o,z,m}$), which could be positive (net accumulation) or negative (net reduction).

```math
\begin{aligned}
& \Gamma_{o,z,(m-1)\times \tau^{period}+1 } =\left(1-\eta_{o,z}^{loss}\right)\times \left(\Gamma_{o,z,m\times \tau^{period}} -\Delta Q_{o,z,m}\right) -  \\
& \frac{1}{\eta_{o,z}^{discharge}}\Theta_{o,z,(m-1)\times \tau^{period}+1} + \eta_{o,z}^{charge}\Pi_{o,z,(m-1)\times \tau^{period}+1} \quad \forall o \in \mathcal{O}^{LDES}, z \in \mathcal{Z}, m \in \mathcal{M}
\end{aligned}
```
By definition $\mathcal{T}^{start}=\{\left(m-1\right) \times \tau^{period}+1 | m \in \mathcal{M}\}$, which implies that this constraint is defined for all values of $t \in T^{start}$. \
**Storage inventory change input periods**
We need additional variables and constraints to approximate energy exchange between representative periods, while accounting for their chronological occurence in the original input time series data and the possibility that two representative periods may not be adjacent to each other (see Figure below). To implement this, we introduce a new variable $Q_{o,z, n}$ that models inventory of storage technology $o \in O$ in zone $z$ in each input period $n \in \mathcal{N}$. Additionally we define a function mapping, $f: n \rightarrow m$, that uniquely maps each input period $n$ to its corresponding representative period $m$. This mapping is available as an output of the process used to identify representative periods (E.g. k-means clustering [Mallapragada et al., 2018](https://www.sciencedirect.com/science/article/pii/S0360544218315238?casa_token=I-6GVNMtAVIAAAAA:G8LFXFqXxRGrXHtrzmiIGm02BusIUmm83zKh8xf1BXY81-dTnA9p2YI1NnGuzlYBXsxK12by)).

![Modeling inter-period energy exchange via long-duration storage when using representative period temporal resolution to approximate annual grid operations](../../assets/LDES_approach.png)
*Figure. Modeling inter-period energy exchange via long-duration storage when using representative period temporal resolution to approximate annual grid operations*

The following two equations define the storage inventory at the beginning of each input period $n+1$ as
the sum of storage inventory at begining of previous input period $n$ plus change in storage inventory for that period.
The latter is approximated by the change in storage inventory in the corresponding representative period,
identified per the mapping $f(n)$.
If the input period is also a representative period,
then a second constraint enforces that initial storage level estimated by the intra-period storage balance constraint should equal the initial storage level estimated from the inter-period storage balance constraints.

```math
\begin{aligned}
& Q_{o,z,n+1} = Q_{o,z,n} + \Delta Q_{o,z,f(n)}
\quad \forall  o \in \mathcal{O}^{LDES}, z \in \mathcal{Z}, n \in \mathcal{N}
\end{aligned}
```

```math
\begin{aligned}
& Q_{o,z,n} =\Gamma_{o,z,f(n)\times \tau^{period}} - \Delta Q_{o,z,m}
\quad \forall  o \in \mathcal{O}^{LDES}, z \in \mathcal{Z}, n \in   \mathcal{N}^{rep},
\end{aligned}
```

Finally, the next constraint enforces that the initial storage level for each input period $n$ must be less than the installed energy capacity limit. This constraint ensures that installed energy storage capacity is consistent with the state of charge during both the operational time periods $t$ during each sample period $m$ as well as at the start of each chronologically ordered input period $n$ in the full annual time series.

```math
\begin{aligned}
    Q_{o,z,n} \leq \Delta^{total, energy}_{o,z}
\quad \forall n \in \mathcal{N}, o \in \mathcal{O}^{LDES}
\end{aligned}
```

If the capacity reserve margin constraint is enabled, a similar set of constraints is used to track the evolution of the energy held in reserve across representative periods. The main linking constraint is as follows:
```math
\begin{aligned}
& \Gamma^{CRM}_{o,z,(m-1)\times \tau^{period}+1 } =\left(1-\eta_{o,z}^{loss}\right)\times \left(\Gamma^{CRM}_{o,z,m\times \tau^{period}} -\Delta Q_{o,z,m}\right) +  \\
& \frac{1}{\eta_{o,z}^{discharge}}\Theta^{CRM}_{o,z,(m-1)\times \tau^{period}+1} - \eta_{o,z}^{charge}\Pi^{CRM}_{o,z,(m-1)\times \tau^{period}+1} \quad \forall o \in \mathcal{O}^{LDES}, z \in \mathcal{Z}, m \in \mathcal{M}
\end{aligned}
```
All other constraints are identical to those used to track the actual state of charge, except with the new variables $Q^{CRM}_{o,z,n}$ and $\Delta Q^{CRM}_{o,z,n}$ used in place of $Q_{o,z,n}$ and $\Delta Q_{o,z,n}$, respectively.
"""
function long_duration_storage_tes!(EP::Model, inputs::Dict, setup::Dict)
    println("TES Long Duration Storage Module")

    gen = inputs["RESOURCES"]

    CapacityReserveMargin = setup["CapacityReserveMargin"]

    REP_PERIOD = inputs["REP_PERIOD"]     # Number of representative periods

    STOR_LONG_DURATION = inputs["STOR_LONG_DURATION_TES"]

    hours_per_subperiod = inputs["hours_per_subperiod"] #total number of hours per subperiod

    dfPeriodMap = inputs["Period_Map"] # Dataframe that maps modeled periods to representative periods
    NPeriods = size(inputs["Period_Map"])[1] # Number of modeled periods

    MODELED_PERIODS_INDEX = 1:NPeriods
    REP_PERIODS_INDEX = MODELED_PERIODS_INDEX[dfPeriodMap[!, :Rep_Period] .== MODELED_PERIODS_INDEX]

    ### Variables ###

    # Variables to define inter-period energy transferred between modeled periods

    # State of charge of storage at beginning of each modeled period n
    @variable(EP, vSOCw_TES[y in STOR_LONG_DURATION, n in MODELED_PERIODS_INDEX]>=0)

    # Build up in storage inventory over each representative period w
    # Build up inventory can be positive or negative
    @variable(EP, vdSOC_TES[y in STOR_LONG_DURATION, w = 1:REP_PERIOD])

    if CapacityReserveMargin > 0
        # State of charge held in reserve for storage at beginning of each modeled period n
        @variable(EP, vCAPRES_socw_TES[y in STOR_LONG_DURATION, n in MODELED_PERIODS_INDEX]>=0)

        # Build up in storage inventory held in reserve over each representative period w
        # Build up inventory can be positive or negative
        @variable(EP, vCAPRES_dsoc_TES[y in STOR_LONG_DURATION, w = 1:REP_PERIOD])
    end

    ### Constraints ###

    # Links last time step with first time step, ensuring position in hour 1 is within eligible change from final hour position
    # Modified initial state of storage for long-duration storage - initialize wth value carried over from last period
    # Alternative to cSoCBalStart constraint which is included when not modeling operations wrapping and long duration storage
    # Note: tw_min = hours_per_subperiod*(w-1)+1; tw_max = hours_per_subperiod*w
    @constraint(EP,
        cSoCBalLongDurationStorageStart_TES[w = 1:REP_PERIOD, y in STOR_LONG_DURATION],
        EP[:vS_TES][y,
            hours_per_subperiod * (w - 1) + 1]==(1 - self_discharge(gen[y])) *
                                                (EP[:vS_TES][y, hours_per_subperiod * w] -
                                                 vdSOC_TES[y, w])
                                                -
                                                (1 / efficiency_down(gen[y]) * EP[:vP][
            y, hours_per_subperiod * (w - 1) + 1]) +
                                                (efficiency_up(gen[y]) * EP[:vCHARGE_TES][
            y, hours_per_subperiod * (w - 1) + 1]))

    # Storage at beginning of period w = storage at beginning of period w-1 + storage built up in period w (after n representative periods)
    ## Multiply storage build up term from prior period with corresponding weight
    @constraint(EP,
        cSoCBalLongDurationStorage_TES[y in STOR_LONG_DURATION, r in MODELED_PERIODS_INDEX],
        vSOCw_TES[y,
            mod1(r + 1, NPeriods)]==vSOCw_TES[y, r] +
                                    vdSOC_TES[y, dfPeriodMap[r, :Rep_Period_Index]])

    # Storage at beginning of each modeled period cannot exceed installed energy capacity
    @constraint(EP,
        cSoCBalLongDurationStorageUpper_TES[y in STOR_LONG_DURATION,
            r in MODELED_PERIODS_INDEX],
        vSOCw_TES[y, r]<=EP[:eTotalCapEnergy_TES][y])

    # Initial storage level for representative periods must also adhere to sub-period storage inventory balance
    # Initial storage = Final storage - change in storage inventory across representative period
    @constraint(EP,
        cSoCBalLongDurationStorageSub_TES[y in STOR_LONG_DURATION, r in REP_PERIODS_INDEX],
        vSOCw_TES[y,
            r]==EP[:vS_TES][y, hours_per_subperiod * dfPeriodMap[r, :Rep_Period_Index]] -
                vdSOC_TES[y, dfPeriodMap[r, :Rep_Period_Index]])

    # Capacity Reserve Margin policy
    if CapacityReserveMargin > 0
        # LDES Constraints for storage held in reserve

        # Links last time step with first time step, ensuring position in hour 1 is within eligible change from final hour position
        # Modified initial virtual state of storage for long-duration storage - initialize wth value carried over from last period
        # Alternative to cVSoCBalStart_TES constraint which is included when not modeling operations wrapping and long duration storage
        # Note: tw_min = hours_per_subperiod*(w-1)+1; tw_max = hours_per_subperiod*w
        @constraint(EP,
            cVSoCBalLongDurationStorageStart_TES[w = 1:REP_PERIOD, y in STOR_LONG_DURATION],
            EP[:vCAPRES_socinreserve_TES][y,
                hours_per_subperiod * (w - 1) + 1]==(1 - self_discharge(gen[y])) *
                                                    (EP[:vCAPRES_socinreserve_TES][
                y, hours_per_subperiod * w] - vCAPRES_dsoc_TES[y, w])
                                                    +
                                                    (1 / efficiency_down(gen[y]) *
                                                     EP[:vCAPRES_discharge_TES][
                y, hours_per_subperiod * (w - 1) + 1]) -
                                                    (efficiency_up(gen[y]) *
                                                     EP[:vCAPRES_charge_TES][
                y, hours_per_subperiod * (w - 1) + 1]))

        # Storage held in reserve at beginning of period w = storage at beginning of period w-1 + storage built up in period w (after n representative periods)
        ## Multiply storage build up term from prior period with corresponding weight
        @constraint(EP,
            cVSoCBalLongDurationStorage_TES[y in STOR_LONG_DURATION,
                r in MODELED_PERIODS_INDEX],
            vCAPRES_socw_TES[y,
                mod1(r + 1, NPeriods)]==vCAPRES_socw_TES[y, r] +
                                        vCAPRES_dsoc_TES[y, dfPeriodMap[r, :Rep_Period_Index]])

        # Initial reserve storage level for representative periods must also adhere to sub-period storage inventory balance
        # Initial storage = Final storage - change in storage inventory across representative period
        @constraint(EP,
            cVSoCBalLongDurationStorageSub_TES[y in STOR_LONG_DURATION, r in REP_PERIODS_INDEX],
            vCAPRES_socw_TES[y,r]==EP[:vCAPRES_socinreserve_TES][y,
                hours_per_subperiod * dfPeriodMap[r, :Rep_Period_Index]] -
                    vCAPRES_dsoc_TES[y, dfPeriodMap[r, :Rep_Period_Index]])

        # energy held in reserve at the beginning of each modeled period acts as a lower bound on the total energy held in storage
        @constraint(EP,
            cSOCMinCapResLongDurationStorage_TES[y in STOR_LONG_DURATION,
                r in MODELED_PERIODS_INDEX],
            vSOCw_TES[y, r]>=vCAPRES_socw_TES[y, r])
    end
end
