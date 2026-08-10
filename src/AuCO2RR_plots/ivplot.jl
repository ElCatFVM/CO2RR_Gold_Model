using CSV, DataFrames, CairoMakie, GridVisualize




function iv_curve_axis(ivresult;
    cutoff = -0.4,
    showlegend = true,
    #title = "IV Curve",
    species = iohminus,
    data_dir = "../data/catmap_CO2R_data",
)

    # ---- data from simulation ----
    v_all = ivresult.voltages
    mask  = v_all .< cutoff
    volts = v_all[mask]

    # currents(ivresult, species) 
    I_sim = abs.(currents(ivresult, species))[mask] .* (cm^2/mA)
    #I_sim = max.(I_sim, eps(Float64))

    # ---- load csvs ----
    table  = readdlm(joinpath(data_dir, "IV-Ringe-digitized.csv"), ',', Float64, '\n')
    df_v   = table[:, 1]
    df_I   = abs.(table[:, 2])
    df_I   = max.(df_I, eps(Float64))

    table2 = readdlm(joinpath(data_dir, "Ringe-theorical.csv"), ',', Float64, '\n')
    df2_v  = table2[:, 1]
    df2_I  = abs.(table2[:, 2])
    df2_I  = max.(df2_I, eps(Float64))

    table3 = readdlm(joinpath(data_dir, "Ringe-experimental.csv"), ',', Float64, '\n')
    df3_v  = table3[:, 1]
    df3_I  = abs.(table3[:, 2])
    df3_I  = max.(df3_I, eps(Float64))

    # ---- FIGURE STYLE (match conc_vs_voltage_axis) ----
    fig = Figure(size=(960, 540))
    ax = Axis(fig[1, 1];
        xlabel = lab_voltage,
        # `|I|`, not `I`: `abs` above folds the sign away so the log axis can show it.
        ylabel = lab_current_abs,
        yscale = log10,
        limits = ((-1.25, -0.50), (1e-11, 1e2)),
    )

	xt = [-1.2, -1.0, -0.8, -0.6]
    ax.xticks = (xt, [@sprintf("%.1f", x) for x in xt])

    yt_vals = 10.0 .^ (0:-3:-9)
    yt_lbls = [powlab(0), powlab(-3), powlab(-6), powlab(-9)]
    ax.yticks = (yt_vals, yt_lbls)
    # Axis style (spinewidth, sizes, fonts, grid) now comes from the global theme.

    # ---- PLOTS ----
    # simulation line
    lines!(ax, volts, I_sim; color=:green, linewidth=LW_LINE, label="e⁻, we")


	"""	
    # digitized (cross)
    scatter!(ax, df_v, df_I;
        marker = :xcross,
        markersize = 12,
        color = :red,
        label = "Ringe et. al"
    )
	"""

    # theoretical (circle)
    scatter!(ax, df2_v, df2_I;
        marker = :circle,
        markersize = 8,
        color = :blue,
        label = "Ringe et. al : Theoretical"
    )

    # experimental (triangle up)
    scatter!(ax, df3_v, df3_I;
        marker = :utriangle,
        markersize = 12,
        color = :magenta,
        label = "Ringe et. al : Experimental"
    )

    showlegend && axislegend(ax, position=:rt)

    return fig
end


function conc_vs_voltage_axis(
    result, m;
    grid,
    useonly_pH::Bool = false,
    showlegend::Bool = false,
)
    bulk = m.bulk
    species  = getproperty.(bulk, :name)
    colors   = getproperty.(bulk, :color)
    nspecies = length(species)

    tsol  = LiquidElectrolytes.voltages_solutions(result)
    vgrid = result.voltages

    xcoords    = grid.components[XCoordinates]
    ielectrode = argmin(xcoords)

    scale = 1.0 / (mol / dm^3)
    nv = length(vgrid)
    conc_electrode = fill(NaN, nspecies, nv)

    for (j, v) in enumerate(vgrid)
        sol = tsol(v)
        sol === nothing && continue
        @inbounds for ia in 1:nspecies
            conc_electrode[ia, j] = sol[ia, ielectrode] * scale
        end
    end

    fig = Figure(size=(960, 540))
    ax = Axis(fig[1, 1];
        xlabel = lab_voltage,
        ylabel = lab_conc_surface,
        yscale  = log10,
        limits = ((-1.25, -0.50), (1e-11, 1e1)),
    )

    xt = [-1.2, -1.0, -0.8, -0.6]
    ax.xticks = (xt, [@sprintf("%.1f", x) for x in xt])

    yt_vals = 10.0 .^ (0:-3:-9)
    yt_lbls = [powlab(0), powlab(-3), powlab(-6), powlab(-9)]
    ax.yticks = (yt_vals, yt_lbls)
    # Axis style now comes from the global theme.

    if useonly_pH
        iH = findfirst(isequal("H⁺"), species)
        iH === nothing && error("H⁺ not found in species list.")
        y = max.(conc_electrode[iH, :], eps(Float64))
        lines!(ax, vgrid, y; color=colors[iH], linewidth=LW_LINE, label=species[iH])
    else
        for ia in 1:nspecies
            y = max.(conc_electrode[ia, :], eps(Float64))
            lines!(ax, vgrid, y; color=colors[ia], linewidth=LW_LINE, label=species[ia])
        end
    end

    showlegend && axislegend(ax, position=:rt)

    return (fig=fig, ax=ax, species=species, colors=colors, conc_electrode=conc_electrode, vgrid=vgrid)
end


function addplot_ax!(
    ax, sol, vshow, m;
    grid,
    useonly_pH::Bool=false,
    scale = 1.0 / (mol / dm^3),
    clear::Bool=true,
)
    bulk = m.bulk
    species = getproperty.(bulk, :name)
    colors  = getproperty.(bulk, :color)

    x = grid.components[XCoordinates] .+ 1.0e-14

    clear && empty!(ax)

    if useonly_pH
        iH = findfirst(isequal("H⁺"), species)
        iH === nothing && error("H⁺ not found in species list.")
        y = log10.(sol[iH, :] .* scale)
        lines!(ax, x, y; color=colors[iH], linewidth=LW_LINE, label=species[iH])
    else
        nc = min(size(sol, 1), length(species), length(colors))
        for ia in 1:nc
            y = log10.(sol[ia, :] .* scale)
            lines!(ax, x, y; color=colors[ia], linewidth=LW_LINE, label=species[ia])
        end
    end

    return ax
end

function plot1d_makie(
    result, vshow, m;
    grid,
    L,
    useonly_pH::Bool=false,
    df_compare=nothing,
    fig_size=(960, 540),
)
    bulk = m.bulk
    tsol = LiquidElectrolytes.voltages_solutions(result)

    fig = Figure(size=fig_size)
    ax  = Axis(fig[1, 1];
        xlabel = lab_distance,
        ylabel = rich(
            "log", subscript("10"), " ",
            rich("a", font = :bold_italic), subscript("i", font = :italic)
        ),
        xscale = log10,
        limits = ((1e-11, L*1.2), (-11, 1)),
    )

    # Axis style now comes from the global theme.

    sol = tsol(vshow)
    sol === nothing && error("No solution available at voltage $vshow")

    addplot_ax!(ax, sol, vshow, m; grid=grid, useonly_pH=useonly_pH, clear=true)

    if df_compare !== nothing
        addplot_ax!(ax, df_compare, vshow, m; grid=grid, useonly_pH=useonly_pH, clear=false)
    end

    return (fig=fig, ax=ax, species=getproperty.(bulk, :name), colors=getproperty.(bulk, :color))
end

function conc_vs_voltage_axis_compare(result, m;
    grid,
    useonly_pH=false,
    showlegend=true,
    compare=false
)
    bulk = m.bulk
    species  = getproperty.(bulk, :name)
    colors   = getproperty.(bulk, :color)
    nspecies = length(species)

    tsol  = LiquidElectrolytes.voltages_solutions(result)
    vgrid = result.voltages

    xcoords    = grid.components[XCoordinates]
    ielectrode = argmin(xcoords)

    scale = 1.0 / (mol / dm^3)
    nv = length(vgrid)
    conc_electrode = fill(NaN, nspecies, nv)

    for (j, v) in enumerate(vgrid)
        sol = tsol(v)
        sol === nothing && continue
        @inbounds for ia in 1:nspecies
            conc_electrode[ia, j] = sol[ia, ielectrode] * scale
        end
    end

    fig = Figure(size=(960, 540))
    ax = Axis(fig[1, 1];
        xlabel = lab_voltage,
        ylabel = lab_conc_surface,
        yscale = log10,
        limits = ((-1.25, -0.50), (1e-11, 1e1)),
    )

    xt = [-1.2, -1.0, -0.8, -0.6]
    ax.xticks = (xt, [@sprintf("%.1f", x) for x in xt])

    yt_vals = 10.0 .^ (0:-3:-9)
    yt_lbls = [powlab(0), powlab(-3), powlab(-6), powlab(-9)]
    ax.yticks = (yt_vals, yt_lbls)
    # Axis style now comes from the global theme.

    if useonly_pH
        iH = findfirst(isequal("H⁺"), species)
        iH === nothing && error("H⁺ not found in species list.")
        y = max.(conc_electrode[iH, :], eps(Float64))
        lines!(ax, vgrid, y; color=colors[iH], linewidth=LW_LINE, label=species[iH])
    else
        for ia in 1:nspecies
            y = max.(conc_electrode[ia, :], eps(Float64))
            lines!(ax, vgrid, y; color=colors[ia], linewidth=LW_LINE, label=species[ia])
        end
    end

    if compare
        path = "../data/catmap_CO2R_data/voltage-conc.csv"

        df = CSV.read(path, DataFrame; delim=',', ignorerepeated=true)

        # if parsed as single column like: "Index,Voltage,Concentration"
        if length(names(df)) == 1 && occursin(",", String(first(names(df))))
            df = CSV.read(path, DataFrame; delim=',', header=1, ignorerepeated=true)
        end

        # if actually TSV
        if length(names(df)) == 1 && occursin("\t", String(first(names(df))))
            df = CSV.read(path, DataFrame; delim='\t', ignorerepeated=true)
        end

        rename!(df, Dict(n => Symbol(replace(strip(String(n)), '\ufeff' => "")) for n in names(df)))

        for col in (:Index, :Voltage, :Concentration)
            hasproperty(df, col) || error("Catmap CSV header parse failed. Got $(names(df))")
        end

        df.Index = Int.(df.Index)
        df.Voltage = Float64.(df.Voltage)
        df.Concentration = Float64.(df.Concentration)

        for g in groupby(df, :Index)
            idx = first(g.Index)
            p = sortperm(g.Voltage)

            V = g.Voltage[p]
            y = max.(g.Concentration[p], eps(Float64))

            if 1 <= idx <= nspecies
                lines!(ax, V, y; color=colors[idx], linewidth=LW_EXP, linestyle=:dashdot)
            else
                lines!(ax, V, y; linewidth=LW_EXP, linestyle=:dashdot)
            end
        end
    end

    showlegend && axislegend(ax; position=:rt)
    return (fig=fig, ax=ax)
end




function plot1d(result, vshow, m; grid, L, df_compare=nothing)
    bulk = m.bulk
    tsol = LiquidElectrolytes.voltages_solutions(result)

    vis = GridVisualizer(;
        size    = (600, 300),
        clear   = true,
        legend  = :rt,
        limits  = (-11, 1),
        xlimits = (1e-11, L*1.2),
        # Plain strings, not the `lab_*` rich labels: this is a `GridVisualizer`, which
        # takes only `String`. Parentheses for the unit, to match the Makie figures.
        xlabel  = "Distance from Electrode x (m)",
        ylabel  = "log₁₀ aᵢ",
        xscale  = :log,
    )

    sol = tsol(vshow)
    sol === nothing && error("No solution available at voltage $vshow")

    addplot_solution!(vis, sol, vshow, m; grid=grid)

    if df_compare !== nothing
        addplot_df!(vis, df_compare, m)
    end

    return reveal(vis)
end

function addplot_solution!(vis, sol, vshow, m; grid)
    bulk = m.bulk
    species = getproperty.(bulk, :name)
    colors  = getproperty.(bulk, :color)

    scale = 1.0 / (mol / dm^3)
    title = @sprintf("Φ_we=%+1.2f [V vs. SHE]", vshow)

    x = grid.components[XCoordinates] .+ 1.0e-14

    scalarplot!(vis, x, log10.(sol[1, :] .* scale);
        color=colors[1], label=species[1], clear=true, title=title)

    nc = min(size(sol, 1), length(species), length(colors))
    for ia in 2:nc
        scalarplot!(vis, x, log10.(sol[ia, :] .* scale);
            color=colors[ia], label=species[ia], clear=false)
    end

    return vis
end


function addplot_df!(vis, df::DataFrame, m)
    bulk = m.bulk
    species  = getproperty.(bulk, :name)
    colors   = getproperty.(bulk, :color)
    nspecies = length(species)

    nms = names(df)
    @assert length(nms) ≥ 3
    rename!(df, Dict(nms[1]=>:Index, nms[2]=>:Distance, nms[3]=>:Concentration))

    idxs = unique(skipmissing(df.Index))

    for idx_raw in idxs
        idx = try
            Int(idx_raw)
        catch
            try
                parse(Int, String(idx_raw))
            catch
                continue
            end
        end

        (1 <= idx <= nspecies) || continue

        mask = (df.Index .== idx_raw)
        x = Float64.(df.Distance[mask])
        y = max.(Float64.(df.Concentration[mask]), eps(Float64))

        if maximum(x) > 1e-3
            x .*= 1e-6
        end
        x .+= 1e-14

        p = sortperm(x)
        scalarplot!(vis, x[p], log10.(y[p]);
            color=colors[idx], clear=false, label="")
    end

    return vis
end

function plot1d_movie(result, m; grid, L, step=5, file="concentrations.gif", framerate=2)
    bulk = m.bulk
    tsol = LiquidElectrolytes.voltages_solutions(result)

    vis = GridVisualizer(;
        size    = (650, 400),
        clear   = true,
        legend  = :rt,
        limits  = (-11, 1),
        xlimits = (1e-11, L*1.2),
        # Plain strings, not the `lab_*` rich labels: this is a `GridVisualizer`, which
        # takes only `String`. Parentheses for the unit, to match the Makie figures.
        xlabel  = "Distance from Electrode x (m)",
        ylabel  = "log₁₀ aᵢ",
        xscale  = :log,
    )

    vrange = result.voltages[end:-step:1]

    movie(vis, file=file, framerate=framerate) do vis
        for vshow in vrange
            sol = tsol(vshow)
            sol === nothing && continue
            addplot_solution!(vis, sol, vshow, m; grid=grid)
            reveal(vis)
        end
    end

    return file
end











function plotcurr_over_L(results::Dict{Int,Any};
    species=ico, cutoff=-0.4, title="IV vs L"
)
    vis = GridVisualizer(;
        size   = (800, 500),
        title  = title,
        xlabel = lab_voltage,
        ylabel = lab_current,
        legend = :rt,
        yscale = :log,
    )

    items = sort(collect(results); by=first)
    n = length(items)

    cols = Makie.resample_cmap(:cool, n)   
    for (k, (L, ivres)) in enumerate(items)
        volts = ivres.voltages
        mask  = volts .< cutoff

        scalarplot!(vis,
            volts[mask],
            abs.(currents(ivres, species))[mask] .* cm^2/mA;
            clear = false,
            label = "L = $(L) μm",
            color = cols[k],              
        )
    end

    reveal(vis)
end




"""
plot_iv_with_ringe_refs(ivresult; species=iohminus, cutoff=-0.4, paths=..., vis_kwargs...)

Plot an IV curve from `ivresult` and overlay three reference datasets from CSV files:
- IV-Ringe-digitized.csv (digitized)
- Ringe-theorical.csv (theoretical)
- Ringe-experimental.csv (experimental)

Assumes the reference CSVs have two numeric columns: voltage, current.
"""
function plot_iv_with_ringe_refs(
    ivresult;
    species = iohminus,
    cutoff  = -0.4,
    paths = (
        digitized   = "../data/catmap_CO2R_data/IV-Ringe-digitized.csv",
        theoretical = "../data/catmap_CO2R_data/Ringe-theorical.csv",
        experimental= "../data/catmap_CO2R_data/Ringe-experimental.csv",
    ),
    vis_kwargs...
)
    volts_all = ivresult.voltages
    mask = volts_all .< cutoff
    volts = volts_all[mask]
    I_sim = abs.(currents(ivresult, species))[mask] .* (cm^2/mA)

    vis = GridVisualizer(;
        size   = get(vis_kwargs, :size, (600, 400)),
        title  = get(vis_kwargs, :title, "IV Curve"),
        xlabel = get(vis_kwargs, :xlabel, lab_voltage),
        ylabel = get(vis_kwargs, :ylabel, lab_current),
        legend = get(vis_kwargs, :legend, :lb),
        yscale = get(vis_kwargs, :yscale, :log),
    )

    scalarplot!(
        vis,
        volts,
        I_sim;
        clear = true,
        linestyle = :solid,
        label = "sim",
    )

    function _read2col(path)
        tbl = readdlm(path, ',', Float64, '\n')
        return (voltage = tbl[:, 1], current = tbl[:, 2])
    end

    d1 = _read2col(paths.digitized)
    d2 = _read2col(paths.theoretical)
    d3 = _read2col(paths.experimental)

    scalarplot!(
        vis, d1.voltage, d1.current;
        clear=false, linewidth=0, markershape=:cross, markersize=8, markevery=1, color =:red,
        label="Ringe et. al (digitized)",
    )

    scalarplot!(
        vis, d2.voltage, d2.current;
        clear=false, linewidth=0, markershape=:circle, markersize=4, markevery=1, color =:blue,
        label="Ringe et. al (theoretical)",
    )

    scalarplot!(
        vis, d3.voltage, d3.current;
        clear=false, linewidth=0, markershape=:utriangle, markersize=8, markevery=1, color =:magenta,
        label="Ringe et. al (experimental)",
    )

    return reveal(vis)
end


# =====================================================================
# Added from scripts/row_interaction_script.jl  (only added, nothing removed)
# Batch 3: electrode activity vs voltage (self-contained; take electrolyte/grid).
# =====================================================================

# ── moved from cell f2a1829d (electrode_activity_vs_voltage) ──
function electrode_activity_vs_voltage(result, grid, m;
    c_ref = 1.0 * ufac"mol/dm^3",
    node_selector = :minx,
)
    electrolyte = m.elydata
    tsol = LiquidElectrolytes.voltages_solutions(result)
    vgrid = LiquidElectrolytes.voltages(result)

    xcoords = grid.components[XCoordinates]
    ielectrode = node_selector === :minx ? argmin(xcoords) : argmax(xcoords)

    cspecies = electrolyte.cspecies
    ip       = LiquidElectrolytes.pressure_index(electrolyte)

    nspecies = length(cspecies)
    nv       = length(vgrid)

    γ_e = fill(NaN, nspecies, nv)
    a_e = fill(NaN, nspecies, nv)
    c_e = fill(NaN, nspecies, nv)

    for (j, U) in enumerate(vgrid)
        sol = tsol(U)
        sol === nothing && continue

        unode = view(sol, :, ielectrode)
        pnode = unode[ip]

        γ = zeros(eltype(sol), size(sol, 1))

        electrolyte.actcoeff!(γ, unode, pnode, electrolyte)

        for (k, ic) in enumerate(cspecies)
            cval = unode[ic]
            γval = γ[ic]
            aval = γval * (cval / c_ref)

            c_e[k, j] = cval / c_ref
            γ_e[k, j] = γval
            a_e[k, j] = aval
        end
    end

    return (
        voltages = vgrid,
        ielectrode = ielectrode,
        cspecies = cspecies,
        gamma_electrode = γ_e,
        activity_electrode = a_e,
        concentration_scaled = c_e,
    )
end

# ── moved from cell 11d6598c (activity_vs_voltage_axis) ──
function activity_vs_voltage_axis(
    result, m;
    grid,
    ipressure,
    model_type::String, # "Stefan(MPB)" "DGML"
    useonly_pH::Bool = false,
    showlegend::Bool = false,
)
    bulk = m.bulk
    electrolyte = m.elydata
    species  = getproperty.(bulk, :name)
    colors   = getproperty.(bulk, :color)
    nspecies = length(species)

    tsol  = LiquidElectrolytes.voltages_solutions(result)
    vgrid = result.voltages

    xcoords    = grid.components[XCoordinates]
    ielectrode = argmin(xcoords)

    scale = 1.0 / (mol / dm^3)
    nv = length(vgrid)

    activity_electrode = fill(NaN, nspecies, nv)
    gamma_electrode    = fill(NaN, nspecies, nv)
    conc_electrode     = fill(NaN, nspecies, nv)

    cspecies = electrolyte.cspecies

    for (j, v) in enumerate(vgrid)
        sol = tsol(v)
        sol === nothing && continue

        # Extract spatial unit (x=0) our unit vector and magnitude
        cnode = zeros(eltype(sol), maximum(cspecies))
        for ic in cspecies
            cnode[ic] = sol[ic, ielectrode]
        end
        pnode = sol[ipressure, ielectrode]

        v0 = electrolyte.v0
        bar_c = 1.0 / v0
        RT = electrolyte.RT

        Phi = sum(cnode[ic] * electrolyte.v[ic] for ic in cspecies)
        solvent_frac = max(1.0 - Phi, eps(Float64))

        for ia in 1:nspecies
            c = sol[ia, ielectrode]
            v_a = electrolyte.v[ia]
            size_ratio = v_a / v0

            term_conc = c / bar_c

            if model_type == "DMGL_γ"
				# DGML Model: Treats the electrolyte as an ideal incompressible mixture.
                # term_press: Captures the mechanical pressure penalty scaled by the specific volume difference.
                # Species with v_a > 0 are physically repelled by local pressure gradients (Barodiffusion).
                # Dimensionless species (v_a = 0) feel pure pressure-correction without steric linkage.
                term_press  = exp((1.0 - size_ratio) * pnode / (bar_c * RT))
                term_steric = solvent_frac^(-size_ratio)
                a_eff_thermo = term_conc * term_press * term_steric

            elseif model_type == "Stefan_γ"
				# Stefan's Model (Bikerman-Freise): Applies a global lattice-based steric penalty.
                # All species, regardless of their actual size (even v_a = 0), are subjected to the exact same penalty (1 - Φ)^-1.
                # This causes the unphysical coupled depletion of point-charge species when supporting cations overcrowd.
                a_eff_thermo = term_conc * (solvent_frac^(-1.0))

            else
                error("Invalid model_type. Use 'DGML' or 'Stefan'.")
            end

            c_scale = c * scale
            a_eff_scaled = a_eff_thermo * (bar_c * scale)

            conc_electrode[ia, j]     = c
            gamma_electrode[ia, j]    = a_eff_scaled / max(c_scale, eps(Float64))
            activity_electrode[ia, j] = a_eff_scaled
        end
    end

    fig = Figure(size=(1200, 800))

    ax = Axis(fig[1,1];
        xlabel = lab_voltage,
        ylabel = lab_activity_surface,
        limits = ((-1.25, -0.50), nothing),
        yscale = log10
    )

    if useonly_pH
        iH = findfirst(isequal("H⁺"), species)
        iH === nothing && error("H⁺ not found in species list.")

        valid_mask = .!isnan.(activity_electrode[iH, :])
        x_data = vgrid[valid_mask]
        y_data = max.(activity_electrode[iH, valid_mask], eps(Float64))

        lines!(ax, x_data, y_data; color=colors[iH], linewidth=LW_LINE, label=species[iH])
    else
        for ia in 1:nspecies
            valid_mask = .!isnan.(activity_electrode[ia, :])
            x_data = vgrid[valid_mask]
            y_data = max.(activity_electrode[ia, valid_mask], eps(Float64))

            if sum(valid_mask) > 0
                lines!(ax, x_data, y_data; color=colors[ia], linewidth=LW_LINE, label=species[ia])
            end
        end
    end

    showlegend && axislegend(ax, position=:rt)

    return (
        fig=fig,
        species=species,
        colors=colors,
        conc_electrode=conc_electrode,
        gamma_electrode=gamma_electrode,
        activity_electrode=activity_electrode,
        vgrid=vgrid,
    )
end

# ==========================================================================
# Panel builders for the polarization figure.
#
# These read the CSVs written by `sweeps_csv.jl`'s publication exporters rather
# than a live result: that is the pipeline the published figure stands on, and
# it means a figure can be redrawn without re-running a sweep. The column
# contract is the exporters' — `Voltage` plus `Current`, `Capacitance`, or one
# column per species.
#
# Each draws into a `fig[i, j]` and returns its `Axis`, so panels can be
# composed by `plot_iv_summary` or used one at a time.
# ==========================================================================

"""
    species_columns(df)

Species names of a wide species-vs-voltage CSV: every column except `Voltage`, in file
order.

The exporters write the species names as the column headers, so the plot reads them from
the file instead of from a list repeated at the call site — a list that would have to be
kept in the same order as the matrix rows and silently mislabels every curve when it is not.
"""
species_columns(df) = [n for n in names(df) if n != "Voltage"]

"""
    panel_polarization!(fig, panel_pos, df; refs, kwargs...)

Polarization curve `|I|` against potential, log scale.

`df` needs `Voltage` and `Current`, as [`export_polarization_csv`](@ref) writes them. `refs`
is a vector of `(df, label, marker, colour)` tuples drawn as scatter behind the simulation —
digitised literature points, typically.

`abs` because the axis is logarithmic; the label says `|I|` for the same reason.
"""
function panel_polarization!(
        fig, panel_pos, df;
        ref = nothing,
        ref_label = "CatINT",
        ref_scale = 1.0,
        ref_color = colorant"#999999",
        ref_lw = LW_DASH,
        points = nothing,
        points_label = "Experiment",
        points_color = colorant"#00916a",
        markersize = 12,
        scale = cm^2 / mA,
        color = colorant"#2980B9",
        label = "MPNP@LiquidElectrolytes.jl",
        lw = LW_SOLID,
        xlabel = lab_voltage,
        ylabel = lab_current_co_abs,
        yaxisposition = :left,
        limits = ((-1.5, -0.4), (1.0e-10, 1.0e2)),
        yticks = (10.0 .^ (0:-5:-10), [powlab(0), powlab(-5), powlab(-10)]),
        xticks = LinearTicks(5),
        legend_position = :lb,
        showlegend = true,
    )
    ax = Axis(
        panel_pos;
        xlabel = xlabel, ylabel = ylabel, yaxisposition = yaxisposition,
        yscale = log10, limits = limits, xticks = xticks, yticks = yticks,
    )
    lines!(ax, df.Voltage, abs.(df.Current) .* scale; color, linewidth = lw)
    # `ref_scale`, separate from `scale`: the simulated column is stored in SI and needs
    # `cm^2 / mA`, while a digitised literature curve is already in mA cm⁻². One shared
    # factor would silently move one of the two by five orders of magnitude.
    ref === nothing || lines!(
        ax, ref.Voltage, abs.(ref.Current) .* ref_scale;
        color = ref_color, linewidth = ref_lw, linestyle = :dash,
    )
    points === nothing || scatter!(
        ax, points.Voltage, abs.(points.Current) .* ref_scale;
        marker = :circle, markersize = markersize, color = points_color,
    )

    if showlegend
        elems = Any[LineElement(color = color, linewidth = lw)]
        labels = Any[label]
        if ref !== nothing
            push!(elems, LineElement(color = ref_color, linewidth = ref_lw, linestyle = :dash))
            push!(labels, ref_label)
        end
        if points !== nothing
            push!(elems, MarkerElement(color = points_color, marker = :circle,
                                       markersize = markersize + 2))
            push!(labels, points_label)
        end
        axislegend(
            ax, elems, labels;
            position = legend_position, labelsize = 28,
            framevisible = true, backgroundcolor = (:white, 0.5),
            framecolor = (:black, 0.5), patchsize = (45, 22),
        )
    end
    return ax
end

"""
    panel_species_vs_voltage!(fig, panel_pos, df; ylabel, kwargs...)

One curve per species against potential, log scale — the shared body of the concentration
and activity panels, which differ only in their label and limits.

`ref` overlays a second table with the same columns, dashed, for a comparison model.
Species missing from `ref` are skipped rather than erroring, because a reference table
often carries fewer of them.
"""
function panel_species_vs_voltage!(
        fig, panel_pos, df;
        ylabel,
        species = nothing,
        ref = nothing,
        ref_value_col = :Concentration,
        colors = nothing,
        ref_colors = nothing,
        annotate = nothing,
        annotate_fontsize = 32,
        lw = LW_SOLID,
        ref_lw = LW_DASH,
        xlabel = lab_voltage,
        yaxisposition = :left,
        limits = ((-1.5, -0.5), (1.0e-11, 1.0e1)),
        yticks = (10.0 .^ (0:-3:-9), [powlab(0), powlab(-3), powlab(-6), powlab(-9)]),
        xticks = LinearTicks(5),
        showlegend = false,
    )
    # Explicit `species` pins the order, which matters when `ref` is long-format and keyed
    # by a numeric index into that same list.
    sp = species === nothing ? species_columns(df) : species
    # `get`, not indexing: a CSV may carry a column this map has never seen, and a missing
    # colour should cost that one curve its identity, not the whole figure.
    cols = colors === nothing ? [get(SPECIES_COLORS, s, colorant"#555555") for s in sp] : colors
    rcols = ref_colors === nothing ?
        [get(SPECIES_PASTEL, s, colorant"#BBBBBB") for s in sp] : ref_colors

    ax = Axis(
        panel_pos;
        xlabel = xlabel, ylabel = ylabel, yaxisposition = yaxisposition,
        yscale = log10, limits = limits, xticks = xticks, yticks = yticks,
    )
    for (i, s) in enumerate(sp)
        y = max.(df[!, s], eps(Float64))     # the log axis cannot take an exact zero
        lines!(ax, df.Voltage, y; color = cols[i], linewidth = lw, label = s)
    end
    if ref !== nothing
        for (i, s) in enumerate(sp)
            xr, yr = _ref_series(ref, i, s, ref_value_col)
            xr === nothing && continue
            lines!(
                ax, xr, max.(yr, eps(Float64));
                color = rcols[i], linewidth = ref_lw, linestyle = :dash,
            )
        end
    end

    # Direct labels beat a legend for seven overlapping curves, but where each one fits is a
    # property of the data, not something to derive: pass the positions as
    # `Dict("K⁺" => (-1.35, 1e-1), …)` in data coordinates. Species absent from the dict are
    # simply not labelled, so a partial dict is fine.
    if annotate !== nothing
        for (i, s) in enumerate(sp)
            haskey(annotate, s) || continue
            x, y = annotate[s]
            text!(
                ax, x, y;
                text = species_rich(s), color = cols[i],
                fontsize = annotate_fontsize, font = :bold, align = (:center, :center),
            )
        end
    end

    showlegend && axislegend(ax; position = :rb)
    return ax
end

"""
    _ref_series(ref, i, name, value_col)

One species' reference curve, from a table in either layout.

**Wide**: one column per species, named as in the main table — what the exporters here
write. **Long**: an `Index` column numbering species in the order they were passed, plus
`Voltage` and `value_col` — what the digitised CatINT tables use.

Returns `(nothing, nothing)` when the species is absent, because a reference table routinely
carries fewer species than the model does and a missing one should drop its dashed curve
rather than the figure.
"""
function _ref_series(ref, i, name, value_col)
    cols = names(ref)
    if String(name) in cols
        return ref.Voltage, ref[!, String(name)]
    elseif "Index" in cols && String(value_col) in cols
        m = ref[!, :Index] .== i
        any(m) || return (nothing, nothing)
        return ref[m, :Voltage], ref[m, value_col]
    end
    return (nothing, nothing)
end

"""
    species_rich(name)

Species name as rich text, with digits subscripted and charges superscripted.

The CSV headers already carry Unicode sub/superscripts (`HCO₃⁻`), which render but sit on
the baseline at the wrong size. This rebuilds the few known names properly and passes
anything else through unchanged.
"""
function species_rich(name)
    name == "K⁺" && return rich("K", superscript("+"))
    name == "H⁺" && return rich("H", superscript("+"))
    name == "OH⁻" && return rich("OH", superscript("−"))
    name == "HCO₃⁻" && return rich("HCO", subscript("3"), superscript("−"))
    name == "CO₃²⁻" && return rich("CO", subscript("3"), superscript("2−"))
    name == "CO₂" && return rich("CO", subscript("2"))
    return rich(name)
end

"""
    panel_dlcap!(fig, panel_pos, df; kwargs...)

Differential capacitance against potential referred to the pzc, on a linear axis.

`shift` is subtracted from the voltage column so the curve is read against `U_pzc`; the
published figure shifts by 0.16 V.

`scale` converts the stored value. `capscalc` and the CSV exporter work in SI, i.e. F m⁻²,
while capacitance is quoted in μF cm⁻² — a factor of 100. Getting this wrong is not obvious
from the figure: 0.17 F m⁻² and 17 μF cm⁻² are the same number in different clothes, and
only the axis label says which one is on the page.

`label` is drawn inside the axis at `label_pos` (relative coordinates), which is how the
published panel names its model rather than spending a legend on one curve.
"""
function panel_dlcap!(
        fig, panel_pos, df;
        shift = 0.0,
        scale = cm^2 / μF,
        color = colorant"#2980B9",
        lw = LW_SOLID,
        xlabel = lab_voltage_pzc,
        ylabel = lab_capacitance,
        yaxisposition = :left,
        limits = ((-0.9, 0.9), (0, 25)),
        xticks = LinearTicks(5),
        yticks = Makie.automatic,
        label = "",
        label_pos = (0.8, 0.6),
        label_fontsize = 28,
    )
    ax = Axis(
        panel_pos;
        xlabel = xlabel, ylabel = ylabel,
        yaxisposition = yaxisposition, limits = limits,
        xticks = xticks, yticks = yticks,
    )
    lines!(ax, df.Voltage .- shift, df.Capacitance .* scale; color, linewidth = lw)
    isempty(label) || text!(
        ax, Point2f(label_pos...);
        text = label, space = :relative, color = color,
        fontsize = label_fontsize, font = :bold, align = (:center, :center),
    )
    return ax
end

"""
    species_table(vgrid, M, species)

A `species × voltage` matrix as the wide table the panels read: `Voltage` plus one column
per species.

Same layout the publication CSV exporters write, so a figure can be fed from a live result
or from a file without the plotting code knowing which.
"""
function species_table(vgrid, M, species)
    size(M, 1) == length(species) || error(
        "matrix has $(size(M, 1)) rows but $(length(species)) species names"
    )
    df = DataFrame(Voltage = collect(vgrid))
    for (i, s) in enumerate(species)
        df[!, String(s)] = vec(M[i, :])
    end
    return df
end

"""
    polarization_table(ivresult; species, scale)

`Voltage`/`Current` table from a live IV result, in the layout
[`panel_polarization!`](@ref) reads.

`species` defaults to OH⁻, matching the CSV exporter and valid on a cathodic sweep — see the
note on [`export_polarization_csv`](@ref).
"""
polarization_table(ivresult; species = iohminus, scale = 1.0) =
    DataFrame(
        Voltage = collect(ivresult.voltages),
        Current = currents(ivresult, species) .* scale,
    )

"""
    surface_conc_table(ivresult, m; grid)

Surface concentrations against potential, as a wide table.

Reads the same computation [`conc_vs_voltage_axis`](@ref) plots, so the table and the
standalone figure cannot disagree; the figure it builds on the way is discarded.
"""
function surface_conc_table(ivresult, m; grid)
    out = conc_vs_voltage_axis(ivresult, m; grid = grid)
    return species_table(out.vgrid, out.conc_electrode, out.species)
end

"""
    surface_activity_table(ivresult, m; grid, model_type, ipressure)

Surface activities against potential, as a wide table.

`model_type` selects the activity model and must match the one the electrolyte was built
with — `"Stefan_γ"` or `"DMGL_γ"`. It is not derived from the electrolyte because the
activity expression lives in the plotting layer, so passing the wrong one silently produces
a plausible curve for the wrong model.
"""
function surface_activity_table(
        ivresult, m;
        grid,
        model_type::String,
        ipressure = LiquidElectrolytes.pressure_index(m.elydata),
    )
    out = activity_vs_voltage_axis(
        ivresult, m; grid = grid, ipressure = ipressure, model_type = model_type,
    )
    return species_table(out.vgrid, out.activity_electrode, out.species)
end

"""
    plot_iv_summary_from_result(ivresult, m; grid, model_type, df_cdl, kwargs...)

[`plot_iv_summary`](@ref) with panels (a)–(c) built from a live IV result instead of from
CSV.

Panel (d) still takes `df_cdl`: differential capacitance comes from a separate
`dlcapsweep`, not from an IV sweep, so there is nothing in `ivresult` to build it from.
The reference overlays — CatINT and the experimental points — stay file-driven, since they
are digitised literature data.

All of [`plot_iv_summary`](@ref)'s keywords apply.
"""
function plot_iv_summary_from_result(
        ivresult, m;
        grid,
        model_type::String,
        df_cdl,
        pol_species = iohminus,
        ipressure = LiquidElectrolytes.pressure_index(m.elydata),
        kwargs...,
    )
    return plot_iv_summary(;
        df_pol = polarization_table(ivresult; species = pol_species),
        df_act = surface_activity_table(ivresult, m; grid, model_type, ipressure),
        df_conc = surface_conc_table(ivresult, m; grid),
        df_cdl = df_cdl,
        kwargs...,
    )
end

"""
    plot_iv_summary(; df_pol, df_act, df_conc, df_cdl, kwargs...)

The four-panel polarization figure: current, surface activity, surface concentration and
differential capacitance, all against potential.

Composed from [`panel_polarization!`](@ref), [`panel_species_vs_voltage!`](@ref) and
[`panel_dlcap!`](@ref), each of which is usable on its own. Was a `let` block in
`plot_publication_notebook.jl` with its own axis styling and its own species list; the
panels now share the package style and read their species from the CSV headers.

Rows 1 and 2 share the potential axis and only the bottom row is labelled, so the four read
as one figure rather than four plots. Capacitance is the exception: it is normally measured
over a wider window, so its axis is not linked — pass `link_cdl = true` if it should be.

Every `df_*` is a table written by the matching exporter in `sweeps_csv.jl`; the `*_ref`
tables overlay a second model.
"""
function plot_iv_summary(;
        df_pol,
        df_act,
        df_conc,
        df_cdl,
        pol_ref = nothing,
        pol_points = nothing,
        act_ref = nothing,
        conc_ref = nothing,
        act_labels = nothing,
        conc_labels = nothing,
        species = ["K⁺", "H⁺", "HCO₃⁻", "CO₃²⁻", "CO₂", "OH⁻", "CO"],
        cdl_shift = 0.16,
        cdl_label = "MPNP",
        pol_scale = cm^2 / mA,
        act_limits = ((-1.5, -0.5), (1.0e-11, 1.0e4)),
        act_yticks = (10.0 .^ (3:-3:-9),
                      [powlab(3), powlab(0), powlab(-3), powlab(-6), powlab(-9)]),
        conc_limits = ((-1.5, -0.5), (1.0e-11, 1.0e1)),
        fig_size = (1400, 1100),
        figure_padding = (40, 60, 30, 60),
        panel_labels = ("(a)", "(b)", "(c)", "(d)"),
        showlegend = true,
    )
    fig = Figure(size = fig_size, figure_padding = figure_padding)

    ax_pol = panel_polarization!(
        fig, fig[1, 1], df_pol;
        ref = pol_ref, points = pol_points, scale = pol_scale, xlabel = "", showlegend,
    )
    # Right-hand y axes in column 2 so each row's two labels sit on the outside of the
    # figure rather than back to back down the middle.
    ax_act = panel_species_vs_voltage!(
        fig, fig[1, 2], df_act;
        ylabel = lab_activity_surface, species, ref = act_ref, annotate = act_labels,
        xlabel = "", yaxisposition = :right,
        limits = act_limits, yticks = act_yticks,
    )
    ax_con = panel_species_vs_voltage!(
        fig, fig[2, 1], df_conc;
        ylabel = lab_conc_surface, species, ref = conc_ref, annotate = conc_labels,
        limits = conc_limits,
    )
    ax_cdl = panel_dlcap!(
        fig, fig[2, 2], df_cdl;
        shift = cdl_shift, label = cdl_label, yaxisposition = :right,
    )

    # Each panel keeps its own voltage window — (a) reaches further negative than (b)/(c),
    # and (d) is against U − U_pzc entirely — so the axes are deliberately not linked. Only
    # the x *label* is dropped from the top row, since the tick values still differ.
    ax_pol.xlabelvisible = false
    ax_act.xlabelvisible = false

    if panel_labels !== nothing
        for (pos, lbl, pad) in zip(
                (fig[1, 1, TopLeft()], fig[1, 2, TopLeft()],
                 fig[2, 1, TopLeft()], fig[2, 2, TopLeft()]),
                panel_labels, (130, 0, 130, 0),
            )
            Label(
                pos, lbl;
                fontsize = 30, font = :bold, padding = (0, pad, 8, 0),
                halign = :right, valign = :bottom,
            )
        end
    end

    rowgap!(fig.layout, 20)
    colgap!(fig.layout, 30)
    resize_to_layout!(fig)
    return (fig = fig, pol = ax_pol, act = ax_act, conc = ax_con, cdl = ax_cdl)
end






