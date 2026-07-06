"""
    plot_cv_current(pnpresult, model;
                    species=nothing,
                    fig_size=(650, 400),
                    scale=cm^2/mA,
                    color_gradient=true)

Plot a CV current trace from a `pnpresult` and return only the figure.

- If `species` is provided, plots `currents(pnpresult, species) .* scale`.
- If `species === nothing`, defaults to the first entry of `model.cspecies`.
- `color_gradient=true` applies a red gradient along the voltage axis.
- Returns: `fig`.
"""
function plot_cv_current(pnpresult, model;
                         species=nothing,
                         fig_size=(650, 400),
                         scale=cm^2/mA,
                         color_gradient=true)

    ic = model.cspecies
    sp = (species === nothing) ? ic[6] : species

    fig = Figure(size = fig_size)
    ax  = Axis(fig[1, 1],
               ylabel = L"I (mA/cm^2)",
               xlabel = L"\phi (V \; \mathrm{vs}\; SHE)",
               # limits = ((-0.5, 0.9), (-2e-20, 2e-20)),
    )

    I = currents(pnpresult, sp) .* scale

    if color_gradient
        cols = RGBf.(range(0, 1, length(pnpresult.voltages)), 0.0, 0.0)
        lines!(ax, pnpresult.voltages, I; color=cols)
    else
        lines!(ax, pnpresult.voltages, I)
    end

    return fig
end



"""
    plot_conc_time_electrode(result, bulk;
                             nspecies=7,
                             fig_size=(650, 400),
                             scale=(mol/dm^3),
                             log_y=true,
                             legend=true)

Plot electrode-adjacent concentrations vs time for multiple species.

Assumptions
- `result.tsol.t` is the time vector.
- `result.tsol[i, 1, t]` accesses the concentration of species `i` at the electrode node (index 1)
  at time index `t`.
- `bulk` is an array-like object whose entries have `:name` and `:color` fields.

Behavior
- If `log_y=true`, negative/zero values are replaced by `NaN` before plotting.

Returns
- `fig`
"""
function plot_conc_time_electrode(result, bulk;
                                  nspecies=7,
                                  fig_size=(650, 400),
                                  scale=(mol/dm^3),
                                  log_y=true,
                                  legend=true)

    species = getproperty.(bulk, :name)
    colors  = getproperty.(bulk, :color)

    times = result.tsol.t
    nt    = length(times)

    conc_at_electrode = [result.tsol[i, 1, t] / scale for i in 1:nspecies, t in 1:nt]

    fig = Figure(size = fig_size)
    ax = Axis(fig[1, 1],
              xlabel = L"time / s",
              ylabel = L"c_{i,\,\mathrm{electrode}} / (\mathrm{mol/dm^3})",
              yscale = (log_y ? log10 : identity))

    for i in 1:nspecies
        y = conc_at_electrode[i, :]
        yplot = log_y ? (map(c -> (c > 0 ? c : NaN), y)) : y
        lines!(ax, times, yplot; color = colors[i], label = string(species[i]))
    end

    if legend
        Legend(fig[1, 2], ax; labelsize=10, backgroundcolor=RGBA(1.0, 1.0, 1.0, 0.5))
    end

    return fig
end

function plot_pressure_varied_sweep(P_recs;
                                   species=iohminus,
                                   fig_size=(1600, 900),
                                   scale=cm^2/mA,
                                   limits=nothing
)

    fig = Figure(size = fig_size)
    if limits !== nothing
        ax = Axis(fig[1, 1],
                  xlabel = L"\phi (V \; \mathrm{vs}\; SHE)",
                  ylabel = L"I (mA/cm^2)",
                  limits = limits,
                  )
    else
        ax = Axis(fig[1, 1],
                xlabel = L"\phi (V \; \mathrm{vs}\; SHE)",
                ylabel = L"I (mA/cm^2)",
                )
    end
    n = length(P_recs)
    cols = [RGB(1 - t, 0, t) for t in LinRange(0, 1, n)]

    plots  = Any[]
    labels = String[]

    for j in 1:n
        p, rec = P_recs[j]
        label = "$(p)\t pCO2(atm)"
        push!(labels, label)

        I = currents(rec, species) .* scale
        push!(plots, lines!(ax, rec.voltages, I; color = cols[j]))
    end

    Legend(fig[1, 2], plots, labels, "Theoretical"; framevisible=true)
    return fig
end

function capsplot_κ(vis, sys; κ_values = [0, 1, 5, 10, 20, 30, 40])
    n = length(κ_values)
    colors = [RGB(0, 0, i/n) for i in 1:n]
    dls = LiquidElectrolytes.DLCapSweepResult[]

    sys_local = deepcopy(sys)
    data = electrolytedata(sys_local)

    κ_original = copy(data.κ)
    cbulk_original = copy(data.c_bulk)

    try
        for (j, κval) in enumerate(κ_values)
            data.κ .= κval
            data.c_bulk .= [0.05, 0.05] * ufac"mol/dm^3"

            result = caps(sys_local)
            push!(dls, result)

            scalarplot!(
                vis,
                result.voltages,
                result.dlcaps / (μF / cm^2);
                linestyle = :solid,
                color = colors[j],
                clear = false,
                label = "κ = $(κval)"
            )
        end
    catch e
        @warn "capsplot_κ failed" exception=e
        rethrow()
    finally
        data.κ .= κ_original
        data.c_bulk .= cbulk_original
    end

    return dls
end
