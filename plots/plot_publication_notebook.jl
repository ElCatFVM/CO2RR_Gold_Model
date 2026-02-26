### A Pluto.jl notebook ###
# v0.20.8

using Markdown
using InteractiveUtils

# ╔═╡ 1472eb23-8b3b-453b-9a91-a550a9988c54
begin
	using Pkg
 	Pkg.activate(joinpath(@__DIR__, ".."))	
	using CSV, DataFrames, Colors, LessUnitful
	using CairoMakie
	using Printf
end

# ╔═╡ f81c8540-591a-4a56-843f-f51085195e2e
@unitfactors mol dm m s K μm bar Pa eV μF V cm μA mA Å nm mm;

# ╔═╡ 1abfb78d-7291-4950-97df-aeb789378bfc
GoldModel = CSV.read("../data/dataplotfiles/iv_GoldModel.csv", DataFrame);

# ╔═╡ 9e7ebac4-b131-4362-9c2e-a07070715df6
LandStorfer = CSV.read("../data/dataplotfiles/iv_LandStorfer.csv", DataFrame);

# ╔═╡ 270509a2-433d-42af-886b-983f226f3229
function plot_IV_from_folder(dir="iv_csv";
                             xlim=(-1.3, 0.9),
                             pressures=[0.1, 0.3, 0.6, 1.0],
                             yscale_log=false)

    files = sort(filter(f -> endswith(f, ".csv"), readdir(dir)))
    @assert !isempty(files) "No CSV files found in '$(dir)'."

    # ------------------------------
    # Figure & Axis styling
    # ------------------------------
    fig = Figure(size=(1050, 500))
    ax = Axis(fig[1, 1];
	          xlabel = L"φ~(\mathrm{V~vs~SHE})",
	          ylabel = L"j~(\mathrm{mA\,cm^{-2}})",
	          limits = ((-1.3, 0.9), (-5.5, 1.8)),
	          xlabelsize = 25, ylabelsize = 25,
	          xgridvisible = false, 
			  ygridvisible = false,
			  spinewidth = 4.5,
			  xtickwidth = 4.5, 
			  ytickwidth = 4.5, 
			  xticklabelsize = 25, yticklabelsize = 25,
    )

    # ------------------------------
    # Color map and data plotting
    # ------------------------------
    n = length(files)
    cols = Makie.resample_cmap(:winter, n)
    plots, labels = Makie.AbstractPlot[], String[]

    for (i, fname) in enumerate(files)
        df = CSV.read(joinpath(dir, fname), DataFrame)
        x = hasproperty(df, :φ) ? df.φ : df[:, 1]
        y = hasproperty(df, :I) ? df.I : df[:, 2]
        mask = (xlim[1] .<= x) .& (x .<= xlim[2])

        label = @sprintf("%.1f pCO₂ atm", pressures[i])
        line = lines!(ax, x[mask], y[mask]; color=cols[i], linewidth=4)
        push!(plots, line)
        push!(labels, label)
    end

    # ------------------------------
    # Region highlighting (same as experimental plot)
    # ------------------------------
    vspan!(ax, -1.30, -0.70, color=(colorant"#87CEFA", 0.30))  # blue region
    text!(ax, -1.00, 1.55, text="CO₂ reduction reaction",
          font="sans-bold", align=(:center, :top),
          fontsize=16, color="#1E90FF")

    vspan!(ax, -0.20, 0.10, color=(colorant"#F7DC6F", 0.30))   # yellow region
    text!(ax, -0.05, 1.55, text="CO oxidation\nreaction\nwith CO₃²⁻",
          font="sans-bold", align=(:center, :top),
          fontsize=16, color=:orange)

    vspan!(ax, 0.10, 0.90, color=(colorant"#F1948A", 0.30))    # red region
    text!(ax, 0.50, 1.55, text="CO oxidation reaction with H₂O",
          font="sans-bold", align=(:center, :top),
          fontsize=16, color="#FF6F61")

    # ------------------------------
    # Axis ticks and limits
    # ------------------------------
    ax.xticks = -1.5:0.3:1.0
    ax.yticks = 1:-1:-5
    # ------------------------------
    # Legend inside axis (same style)
    # ------------------------------
    leg = Legend(fig[1, 1], plots, labels, "Theoretical";
        framevisible = false,
        halign = :right, valign = :bottom,
        padding = (0, 0, 0, 0),
        tellwidth = false, tellheight = false, labelsize = 20, titlesize = 23,
    )
    translate!(leg.blockscene, -40, 40, 0)

    fig
end

# ╔═╡ a9bb3083-d499-4462-845b-c41963cae1e5
fig_2 = plot_IV_from_folder("iv_csv"; yscale_log = false)

# ╔═╡ 40b4e182-7aa2-4518-842a-e70dd9dced0d
let
    fig = Figure(size = (1050, 500))
    ax = Axis(fig[1, 1];
	          xlabel = L"φ~(\mathrm{V~vs~SHE})",
	          ylabel = L"j~(\mathrm{mA\,cm^{-2}})",
	          limits = ((-1.3, 0.9), (-5.5, 1.8)),
	          xlabelsize = 25, ylabelsize = 25,
	          xgridvisible = false, 
			  ygridvisible = false,
			  spinewidth = 4.5,
			  xtickwidth = 4.5, 
			  ytickwidth = 4.5, 
			  xticklabelsize = 25, yticklabelsize = 25,
    )

    # Experimental Data Plotting based on M.T.M Koper
    raw = CSV.read("Langmuir 2021, 37, 5707−5716/Figure_3.csv", DataFrame; header=false)
    pres = vec(Matrix(raw[1:1, :]))
	pressures=[0.1, 0.3, 0.6, 1.0]
    sub = Matrix(raw[4:end, :])
    num = map(x -> x === missing ? NaN : parse(Float64, x), sub)
    num_df = DataFrame(num, :auto)
    npairs = size(num_df, 2) ÷ 2
    cols1 = resample_cmap(:winter, npairs)

    plot_objs1 = []
    labels1 = String[]
    for j in 1:npairs
        xcol, ycol = 2j - 1, 2j
        label = @sprintf("%.1f pCO₂ atm", pressures[j])        
		push!(labels1, label)
        line = lines!(ax, num_df[!, xcol], ((num_df[!, ycol])); color = cols1[j], linewidth = 4)
        push!(plot_objs1, line)
    end

	#---
	vspan!(ax, -1.30, -0.70, color=(colorant"#87CEFA", 0.30))  # lightskyblue, alpha=0.3
	text!(ax, -1.00, 1.55, text="CO₂ reduction reaction", font = "sans-bold", align=
		  (:center, :top), fontsize=16, color="#1E90FF")

	vspan!(ax, -0.20,  0.10, color=(colorant"#F7DC6F", 0.30))  # light yellow
	text!(ax,  -0.05,  1.55, text="CO oxidation\nreaction\nwith CO₃²⁻", font = 
		  "sans-bold", align=(:center, :top), fontsize=16, color=:orange)
	
	vspan!(ax,  0.10,  0.90, color=(colorant"#F1948A", 0.30))  # light red
	text!(ax,  0.50,  1.55, text="CO oxidation reaction with H₂O", font = "sans-bold", 
		  align=(:center, :top), fontsize=16, color="#FF6347")


	#---
	ax.xticks = -1.5:0.3:1.0
	ax.yticks = 1:-1:-5

	
	leg = Legend(fig[1, 1], plot_objs1, labels1, "Experimental";
	    framevisible = false,
	    halign = :right, valign = :bottom,  labelsize = 20, titlesize = 23,
	    padding = (0, 0, 0, 0),
	    tellwidth = false, tellheight = false
	)
	translate!(leg.blockscene,  -40, 40, 0)
	fig
end

# ╔═╡ 4911ad1e-e75a-4d41-83dc-675b3f26ce39
begin
fig = Figure(size=(800, 400))
ax = Axis(fig[1, 1], xlabel = L"φ \text{(V vs SHE)}",ylabel = L"\log|j|\;(\mathrm{mA\,cm^{-2}})",
	          xlabelsize = 18, ylabelsize = 18,
	          xgridvisible = false, 
			  ygridvisible = false,
			  spinewidth = 4.5,
			  xtickwidth = 4.5, 
			  ytickwidth = 4.5, 
     	  	  xticklabelsize = 17, yticklabelsize = 17,
	        )

	lines!(ax, GoldModel.x, log.(abs.(GoldModel.y)), color="#dc143c", label="File 1",  linewidth=4)
	lines!(ax, LandStorfer.x, log.(abs.(LandStorfer.y)), color="#5559fa", label="File 2",  linewidth=4)
	text!(ax, -1.37, -4.0, text="Volume-fraction Model", color="#5559fa", fontsize=19, font = "sans-bold")
	vlines!(ax, -1.5:0.5:0.0; color=(:gray, 0.3), linestyle=:solid, linewidth=2)
	hlines!(ax, -20:10:0; color=(:gray, 0.3), linestyle=:solid, linewidth=2)
	text!(ax, -0.85, -1.00, text="Pressure-dependent Model", color="#dc143c", fontsize=19, font = "sans-bold")
	ax.xticks = -2:0.25:0
	ax.yticks = 0:-5:-25
#axislegend(ax)
fig
end

# ╔═╡ dcb38fd2-23b6-481e-83b7-4f3ee332c4ee
begin
	#curr(J, ix) = [F * abs(j[ix]) for j in J]
	
	function plotcurr(result; df = nothing)
	    scale = 1 / (mol / dm^3)
	    volts = result.voltages[result.voltages .< -0.4]
	    vis = GridVisualizer(;
	                         size = (600, 400),
	                         tilte = "IV Curve",
	                         xlabel = L"\phi_{we} \, (\mathrm{V \; vs \; SHE})",
	                         ylabel = L"I / (\mathrm{mA/cm^2})",       
	                         legend = :lb,
							 yscale = :log,
		)
							 
	    scalarplot!(vis,
	                volts,
	                abs.(currents(ivresult, iohminus))[result.voltages .< -0.4] .* cm^2/mA;
	                color = :green,
	                clear = false,
	                linestyle = :solid,
	                label = "e⁻, we")
		if !isnothing(df)
			scalarplot!(vis,
						df[:voltage],
						df[:current],
						clear = false,
						linewidth = 0,
						markershape = :cross,
						markersize = 8,
						markevery = 1,
						color = :red,
						label = "Ringe et. al")
		end
		
	    reveal(vis)
	end
end

# ╔═╡ ee091126-671e-46d8-8ed5-9457aeefd8fb
function project_root()
    return normpath(joinpath(@__DIR__, ".."))
end

# ╔═╡ 4b02a8c0-a854-436a-9fdc-99ecf637e858
function output_dir()
    return joinpath(project_root(), "data", "output")
end

# ╔═╡ 59b13ba5-c85d-4e31-bd06-e1541b426205
function capsplot_fixed(
    result, title;
    is_Landstorfer::Bool = false,
    nshow::Int = 201,
    xlimits_L=(-1.0, 1.0),
    ylimits_L=(0, 100),
    show_cdl0::Bool=true,
)

    fig = Figure(size = (1050, 500))
    ax = Axis(fig[1, 1];
        xlabel = L"\phi~(\mathrm{V~vs~}\phi_{pzc})",
        ylabel = L"C_{dl}~(\mu \mathrm{F\,cm^{-2}})",
        xlabelsize = 25,
        ylabelsize = 25,
        xticklabelsize = 25,
        yticklabelsize = 25,
        xgridvisible = false,
        ygridvisible = false,
        spinewidth = 4.5,
        xtickwidth = 4.5,
        ytickwidth = 4.5,
        limits = (xlimits_L, ylimits_L)
    )

    nres = length(result)
    hmol = 1 / max(nres, 1)

    plot_objs = Any[]
    labels = String[]

    for i in 1:nres
        c = RGB(i * hmol, 0, 1 - i * hmol)

        v   = result[i].voltage_range
        cap = vec(result[i].dlcaps)

        n = min(nshow, length(v), length(cap))
        v   = v[1:n]
        cap = cap[1:n] / (μF / cm^2)

        # --- label ---
        lbl = if hasproperty(result[i], :molarity)
            "$(result[i].molarity) M"
        elseif hasproperty(result[i], :comb)
            "×$(result[i].comb)"
        else
            "Run $i"
        end

        line = lines!(ax, v, cap;
                      color=c,
                      linewidth=4)

        push!(plot_objs, line)
        push!(labels, lbl)

        # --- PZC marker ---
        if show_cdl0
            scatter!(ax, [0.0],
                     [result[i].cdl0] / (μF / cm^2);
                     color=c,
                     markersize=10)
        end
    end

    # 강조 영역 (앞 스타일 유지)
    vspan!(ax, -1.30, -0.70, color=(colorant"#87CEFA", 0.25))
    vspan!(ax, -0.20,  0.10, color=(colorant"#F7DC6F", 0.25))
    vspan!(ax,  0.10,  0.90, color=(colorant"#F1948A", 0.25))

    leg = Legend(fig[1, 1], plot_objs, labels, title;
        framevisible = false,
        halign = :right,
        valign = :bottom,
        labelsize = 20,
        titlesize = 23,
        tellwidth = false,
        tellheight = false
    )

    translate!(leg.blockscene, -40, 40, 0)

    return fig
end

# ╔═╡ 12a083c3-dde1-4a99-8d66-c7c43afc6b68
let


    csv_paths = [
        raw"../data/output/DLCap_Dirichlet_DMGL_γ_pb_model_1.csv",
        raw"../data/output/DLCap_Dirichlet_DMGL_γ_pb_model_2.csv",
        raw"../data/output/DLCap_Dirichlet_DMGL_γ_pb_model_3.csv",
    ]

    xcol = :Voltage
    ycol = :Capacitance
	n = length(csv_paths)
    colors = Makie.resample_cmap(:cool, n)
    labels = ["model 1", "model 2", "model 3"]

    fig = Figure(size = (480, 660))  # not too tall
    axs = Axis[]

    xlims = (-0.8, 0.4)
    ylims = (0, 200)   
    for i in 1:3
        ax = Axis(fig[i, 1];
            ylabel = (i == 2 ? "Cdl (μF cm⁻²)" : ""), 
            xlabel = (i == 3 ? "Voltage (V)" : ""),   
            title  = labels[i],

            limits = (xlims[1], xlims[2], ylims[1], ylims[2]),

            xgridvisible = false,
            ygridvisible = false,

            xlabelsize = 18, ylabelsize = 18,
            titlesize = 16,

            spinewidth = 3.5,
            xtickwidth = 3.5,
            ytickwidth = 3.5,
            xticklabelsize = 15,
            yticklabelsize = 15,
        )

        # Hide x tick labels for upper panels (shared x-axis look)
        if i < 3
            ax.xticklabelsvisible = false
            ax.xlabelvisible = false
        end

        push!(axs, ax)
    end

    # Link x-axes so zoom/pan stays consistent (true shared x)
    linkxaxes!(axs...)

    for (i, p) in enumerate(csv_paths)
        df = CSV.read(p, DataFrame; header=1)

        lines!(
            axs[i],
            df[!, xcol],
            df[!, ycol] / (μF / cm^2),
            color = colors[i],
            linewidth = 3,
        )
    end

    # Tight-ish spacing
    rowgap!(fig.layout, 8)
    colgap!(fig.layout, 8)

    display(fig)
end

# ╔═╡ Cell order:
# ╠═1472eb23-8b3b-453b-9a91-a550a9988c54
# ╠═f81c8540-591a-4a56-843f-f51085195e2e
# ╠═1abfb78d-7291-4950-97df-aeb789378bfc
# ╠═9e7ebac4-b131-4362-9c2e-a07070715df6
# ╠═270509a2-433d-42af-886b-983f226f3229
# ╠═a9bb3083-d499-4462-845b-c41963cae1e5
# ╠═40b4e182-7aa2-4518-842a-e70dd9dced0d
# ╠═4911ad1e-e75a-4d41-83dc-675b3f26ce39
# ╠═dcb38fd2-23b6-481e-83b7-4f3ee332c4ee
# ╠═ee091126-671e-46d8-8ed5-9457aeefd8fb
# ╠═4b02a8c0-a854-436a-9fdc-99ecf637e858
# ╠═59b13ba5-c85d-4e31-bd06-e1541b426205
# ╠═12a083c3-dde1-4a99-8d66-c7c43afc6b68
