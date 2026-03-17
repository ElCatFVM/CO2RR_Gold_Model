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
        xlabel = L"\mathbf{\text{U}\ \mathrm{vs.}\ \text{SHE}\ \mathrm{(V)}}",
        ylabel = L"\mathbf{I}\;(\mathrm{mA/cm^2})",
        yscale = log10,
        limits = ((-1.25, -0.50), (1e-11, 1e2)),
    )

	xt = [-1.2, -1.0, -0.8, -0.6]
    ax.xticks = (xt, [@sprintf("%.1f", x) for x in xt])

    yt_vals = 10.0 .^ (0:-3:-9)
    yt_lbls = [L"10^{0}", L"10^{-3}", L"10^{-6}", L"10^{-9}"]
    ax.yticks = (yt_vals, yt_lbls)

    ax.spinewidth = 5.5
    ax.xtickwidth = 2.0
    ax.ytickwidth = 2.0
    ax.xticksize  = 8
    ax.yticksize  = 8
    ax.xlabelsize = 25
    ax.ylabelsize = 25
    ax.xticklabelsize = 25
    ax.yticklabelsize = 25

    ax.xgridvisible = false
    ax.ygridvisible = false
    ax.xlabelpadding = 10
    ax.ylabelpadding = 10
    ax.xlabelfont = :bold

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
    result;
    bulk,
    grid,
    useonly_pH::Bool = false,
    showlegend::Bool = false,
)
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
        xlabel = L"\mathbf{\text{U}\ \mathrm{vs.}\ \text{SHE}\ (V)}",
        ylabel = L"\mathbf{c_i^{+}}\;(\mathrm{M})",
        yscale  = log10,
        limits = ((-1.25, -0.50), (1e-11, 1e1)),
    )

    xt = [-1.2, -1.0, -0.8, -0.6]
    ax.xticks = (xt, [@sprintf("%.1f", x) for x in xt])

    yt_vals = 10.0 .^ (0:-3:-9)
    yt_lbls = [L"10^{0}", L"10^{-3}", L"10^{-6}", L"10^{-9}"]
    ax.yticks = (yt_vals, yt_lbls)

    ax.spinewidth = 5.5
    ax.xtickwidth = 2.0
    ax.ytickwidth = 2.0
    ax.xticksize  = 8
    ax.yticksize  = 8
    ax.xlabelsize = 25
    ax.ylabelsize = 25
    ax.xticklabelsize = 25
    ax.yticklabelsize = 25
    ax.xgridvisible = false
    ax.ygridvisible = false
    ax.xlabelpadding = 10
    ax.ylabelpadding = 10
    ax.xlabelfont = :bold

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
    ax, sol, vshow;
    bulk,
    grid,
    useonly_pH::Bool=false,
    scale = 1.0 / (mol / dm^3),
    clear::Bool=true,
)
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
    result, vshow;
    bulk,
    grid,
    L,
    useonly_pH::Bool=false,
    df_compare=nothing,
    fig_size=(960, 540),
)
    tsol = LiquidElectrolytes.voltages_solutions(result)

    fig = Figure(size=fig_size)
    ax  = Axis(fig[1, 1];
        xlabel = "Distance from electrode [m]",
        ylabel = L"\log_{10} c(a_i)",
        xscale = log10,
        limits = ((1e-11, L*1.2), (-11, 1)),
    )

    ax.spinewidth = 5.5
    ax.xtickwidth = 2.0
    ax.ytickwidth = 2.0
    ax.xticksize  = 8
    ax.yticksize  = 8
    ax.xlabelsize = 25
    ax.ylabelsize = 25
    ax.xticklabelsize = 25
    ax.yticklabelsize = 25
    ax.xgridvisible = false
    ax.ygridvisible = false
    ax.xlabelpadding = 10
    ax.ylabelpadding = 10

    sol = tsol(vshow)
    sol === nothing && error("No solution available at voltage $vshow")

    addplot_ax!(ax, sol, vshow; bulk=bulk, grid=grid, useonly_pH=useonly_pH, clear=true)

    if df_compare !== nothing
        addplot_ax!(ax, df_compare, vshow; bulk=bulk, grid=grid, useonly_pH=useonly_pH, clear=false)
    end

    return (fig=fig, ax=ax, species=getproperty.(bulk, :name), colors=getproperty.(bulk, :color))
end

function conc_vs_voltage_axis_compare(result; 
    bulk, 
    grid, 
    useonly_pH=false, 
    showlegend=true, 
    compare=false
)
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
        xlabel = L"\text{Voltage}\ U\ \mathrm{vs.}\ \text{SHE}\ (V)",
        ylabel = L"c_i^{+}\;(\mathrm{M})",
        yscale = log10,
        limits = ((-1.25, -0.50), (1e-11, 1e1)),
    )

    xt = [-1.2, -1.0, -0.8, -0.6]
    ax.xticks = (xt, [@sprintf("%.1f", x) for x in xt])

    yt_vals = 10.0 .^ (0:-3:-9)
    yt_lbls = [L"10^{0}", L"10^{-3}", L"10^{-6}", L"10^{-9}"]
    ax.yticks = (yt_vals, yt_lbls)

    ax.spinewidth = 2.5
    ax.xtickwidth = 2.0
    ax.ytickwidth = 2.0
    ax.xticksize  = 8
    ax.yticksize  = 8
    ax.xlabelsize = 30
    ax.ylabelsize = 30
    ax.xticklabelsize = 20
    ax.yticklabelsize = 20
    ax.xgridvisible = false
    ax.ygridvisible = false

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




function plot1d(result, vshow; bulk, grid, L, df_compare=nothing)
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

    addplot_solution!(vis, sol, vshow; bulk=bulk, grid=grid)

    if df_compare !== nothing
        addplot_df!(vis, df_compare; bulk=bulk)
    end

    return reveal(vis)
end

function addplot_solution!(vis, sol, vshow; bulk, grid)
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


function addplot_df!(vis, df::DataFrame; bulk)
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

function plot1d_movie(result; bulk, grid, L, step=5, file="concentrations.gif", framerate=2)
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
            addplot_solution!(vis, sol, vshow; bulk=bulk, grid=grid)
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
        xlabel = L"\phi_{we}\;(\mathrm{V\;vs\;SHE})",
        ylabel = L"I\;(\mathrm{mA/cm^2})",
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
        xlabel = get(vis_kwargs, :xlabel, L"\phi_{we} \, (\mathrm{V \; vs \; SHE})"),
        ylabel = get(vis_kwargs, :ylabel, L"I / (\mathrm{mA/cm^2})"),
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






