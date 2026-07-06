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
        ylabel = lab_current,
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
    lines!(ax, volts, I_sim; color=:green, linewidth=5, label="e⁻, we")


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
        ylabel = rich(rich("c", font = :italic), subscript("i"), superscript("+"), "  (M)"),
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
        lines!(ax, vgrid, y; color=colors[iH], linewidth=5, label=species[iH])
    else
        for ia in 1:nspecies
            y = max.(conc_electrode[ia, :], eps(Float64))
            lines!(ax, vgrid, y; color=colors[ia], linewidth=5, label=species[ia])
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
        lines!(ax, x, y; color=colors[iH], linewidth=5, label=species[iH])
    else
        nc = min(size(sol, 1), length(species), length(colors))
        for ia in 1:nc
            y = log10.(sol[ia, :] .* scale)
            lines!(ax, x, y; color=colors[ia], linewidth=5, label=species[ia])
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
        xlabel = "Distance from electrode [m]",
        ylabel = rich("log", subscript("10"), " ", rich("c", font = :italic), "(", rich("a", font = :italic), subscript("i"), ")"),
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
        ylabel = rich(rich("c", font = :italic), subscript("i"), superscript("+"), "  (M)"),
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
        lines!(ax, vgrid, y; color=colors[iH], linewidth=3, label=species[iH])
    else
        for ia in 1:nspecies
            y = max.(conc_electrode[ia, :], eps(Float64))
            lines!(ax, vgrid, y; color=colors[ia], linewidth=3, label=species[ia])
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
                lines!(ax, V, y; color=colors[idx], linewidth=2, linestyle=:dashdot)
            else
                lines!(ax, V, y; linewidth=2, linestyle=:dashdot)
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
        xlabel  = "Distance from electrode [m]",
        ylabel  = "log c(aᵢ)",
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
        xlabel  = "Distance from electrode [m]",
        ylabel  = "log c(aᵢ)",
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
        ylabel = rich(rich("ã", font = :italic), subscript("i")),
        limits = ((-1.25, -0.50), nothing),
        yscale = log10
    )

    if useonly_pH
        iH = findfirst(isequal("H⁺"), species)
        iH === nothing && error("H⁺ not found in species list.")

        valid_mask = .!isnan.(activity_electrode[iH, :])
        x_data = vgrid[valid_mask]
        y_data = max.(activity_electrode[iH, valid_mask], eps(Float64))

        lines!(ax, x_data, y_data; color=colors[iH], linewidth=5, label=species[iH])
    else
        for ia in 1:nspecies
            valid_mask = .!isnan.(activity_electrode[ia, :])
            x_data = vgrid[valid_mask]
            y_data = max.(activity_electrode[ia, valid_mask], eps(Float64))

            if sum(valid_mask) > 0
                lines!(ax, x_data, y_data; color=colors[ia], linewidth=5, label=species[ia])
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






