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
    raw = CSV.read("../data/Langmuir_CV_data/Figure_3.csv", DataFrame; header=false)
    pres = vec(Matrix(raw[1:1, :]))
	pressures = ["Ar sat", "0.1", "0.2", "0.3", "0.5", "0.6", "1.0"]
    sub = Matrix(raw[4:end, :])
    num = map(x -> x === missing ? NaN : parse(Float64, x), sub)
    num_df = DataFrame(num, :auto)
    npairs = size(num_df, 2) ÷ 2
    cols1 = resample_cmap(:winter, npairs)

    plot_objs1 = []
    labels1 = String[]
    for j in 1:npairs
        xcol, ycol = 2j - 1, 2j
        label = pressures[j] * " pCO₂ atm"       
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

# ╔═╡ 8545d818-d255-4e8a-af8f-72fdaf9d4bdd


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
        raw"../data/output/DLCap_Robin_Stefan_γ_pb_All_species_1.csv",
        raw"../data/output/DLCap_Robin_Stefan_γ_pnp_All_species_1.csv",
        raw"../data/output/DLCap_Robin_Stefan_γ_pb_Potassium_only_1.csv",
	    raw"../data/output/DLCap_Robin_Stefan_γ_pnp_Potassium_only_1.csv",

    ]

    xcol = :Voltage
    ycol = :Capacitance
	n = length(csv_paths)
    colors = Makie.resample_cmap(:cool, n)
	#for i in length(csv_paths)
    labels = ["model 1", "model 2", "model 3", "model 4"]

    fig = Figure(size = (480, 660))  # not too tall
    axs = Axis[]

    xlims = (-0.8, 0.8)
    ylims = (0, 50)   
    for i in 1:n
        ax = Axis(fig[i, 1];
            ylabel = (i == n/2 ? "Cdl (μF cm⁻²)" : ""), 
            xlabel = (i == n ? "Voltage (V)" : ""),   
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

# ╔═╡ 95844847-6f46-4840-979d-c54282295f41
let


    csv_paths = [
       # raw"../data/output/DLCap_Robin_Stefan_γ_pb_All_species_1.csv",
       # raw"../data/output/DLCap_Robin_Stefan_γ_pnp_All_species_1.csv",
       # raw"../data/output/DLCap_Robin_Stefan_γ_pb_Potassium_only_1.csv",
	    raw"../data/output/DLCap_Robin_Stefan_γ_pnp_Potassium_only_1.csv",

    ]

    xcol = :Voltage
    ycol = :Capacitance
	n = length(csv_paths)
    colors = :lightgray
	#for i in length(csv_paths)
    #labels = ["model 1", "model 2", "model 3", "model 4"]

    fig = Figure(size = (480, 220))  # not too tall
    axs = Axis[]

    xlims = (-0.8, 0.8)
    ylims = (0, 50)   
    for i in 1:n
        ax = Axis(fig[i, 1];
            ylabel = (i == n ? "Cdl (μF cm⁻²)" : ""), 
            xlabel = (i == n ? "Voltage (V)" : ""),   
            #title  = labels[i],

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
            color = colors,#[i],
			
            linewidth = 3,
        )
    end

    # Tight-ish spacing
    rowgap!(fig.layout, 8)
    colgap!(fig.layout, 8)

    display(fig)
end

# ╔═╡ 54d0654f-9add-421b-ba54-90f276e3a10e
let


    csv_paths = [
        raw"../data/output/IV_Robin_DMGL_γ_pnp_ohminus.csv",
        raw"../data/output/IV_Robin_DMGL_γ_pnp_ohminus_koh_nonsol.csv",
        raw"../data/output/Polarization_Curve_Robin_Stefan_γ_pnp_Potassium_only.csv",
    ]

    xcol = :Voltage
    ycol = :Current

	n = length(csv_paths)
    colors = Makie.resample_cmap(:isoluminant_cm_70_c39_n256, n)
    labels = ["model 1", "model 2", "model 3"]

    fig = Figure(size = (480, 660))  # not too tall
    axs = Axis[]

    xlims = (-1.3, -0.4)
    ylims = (-20, 4)   
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
            log.(abs.(df[!, ycol].* (cm^2/mA))),
            color = colors[i],
            linewidth = 3,
        )
    end

    # Tight-ish spacing
    rowgap!(fig.layout, 8)
    colgap!(fig.layout, 8)

    display(fig)
end

# ╔═╡ d7ad452c-d3c1-485c-92ac-0b006dfb5acd
let


    csv_paths = [
       # raw"../data/output/IV_Robin_DMGL_γ_pnp_ohminus.csv",
       # raw"../data/output/IV_Robin_DMGL_γ_pnp_ohminus_koh_nonsol.csv",
        raw"../data/output/Polarization_Curve_Robin_Stefan_γ_pnp_Potassium_only.csv",
    ]

    xcol = :Voltage
    ycol = :Current

	n = length(csv_paths)
    colors = :skyblue#Makie.resample_cmap(:isoluminant_cm_70_c39_n256, n)
    #labels = ["model 1", "model 2", "model 3"]

    fig = Figure(size = (420, 250))  # not too tall
    axs = Axis[]

    xlims = (-1.3, -0.4)
    ylims = (-20, 4)   
    for i in 1:n
        ax = Axis(fig[i, 1];
            ylabel ="Current" ,
            xlabel ="Voltage" ,
           # title  = labels[i],

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



        push!(axs, ax)
    end

    # Link x-axes so zoom/pan stays consistent (true shared x)
    linkxaxes!(axs...)

    for (i, p) in enumerate(csv_paths)
        df = CSV.read(p, DataFrame; header=1)

        lines!(
            axs[i],
            df[!, xcol],
            log.(abs.(df[!, ycol].* (cm^2/mA))),
            color = colors,#[i],
            linewidth = 3,
        )
    end

    # Tight-ish spacing
    rowgap!(fig.layout, 8)
    colgap!(fig.layout, 8)

    display(fig)
end

# ╔═╡ 08bb3a34-88d9-4934-a9c5-a5d61fa77e7e
let
    csv_paths = [
        raw"../data/output/concentration_Robin_Stefan_γ.csv",
        raw"../data/output/concentration_Robin_DMGL_γ.csv",
    ]

    species = ["K⁺","H⁺","HCO₃⁻","CO₃²⁻","CO₂","OH⁻","CO"]
    colors  = [:orange, :gray, :brown, :violet, :red, :green, :blue]

    fig = Figure(size=(700, 800))

    for i in 1:length(csv_paths)
        df = CSV.read(csv_paths[i], DataFrame)

        vgrid = df[!, :Voltage]  # <- CSV 전압 컬럼명이 다르면 여기만 바꾸세요
        conc_electrode = permutedims(Matrix(df[!, Symbol.(species)]))  # (7, N)

        ax = Axis(fig[i, 1];
            xlabel = L"\mathbf{\text{U}\ \mathrm{vs.}\ \text{SHE}\ (V)}",
            ylabel = L"\mathbf{c_i^{+}}\;(\mathrm{M})",
            yscale = log10,
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

        if i < 2
            ax.xticklabelsvisible = false
            ax.xlabelvisible = false
        end
		
        for ia in 1:7
            y = max.(conc_electrode[ia, :], eps(Float64))
            lines!(ax, vgrid, y; color=colors[ia], linewidth=5, label=species[ia])
        end

		if i == 1
			text!(ax, -1.15, 0.3, text=L"\mathrm{K^+}",
				  color=colors[1], fontsize=18, font="sans-bold")
			text!(ax, -0.90, 0.0000002, text=L"\mathrm{H^+}", 
				  color=colors[2], fontsize=18, font = "sans-bold") 
			text!(ax, -1.05, 5e-11, text=L"\mathrm{CO_3^{2-}}",
				  color=colors[4], fontsize=18, font = "sans-bold") 
			text!(ax, -1.0, 2.5e-7, text=L"\mathrm{HCO_3^-}", 
				  color=colors[3], fontsize=18, font = "sans-bold") 
			text!(ax, -1.2, 5.2e-6, text=L"\mathrm{CO_2}", 
				  color=colors[5], fontsize=18, font = "sans-bold")
			text!(ax, -0.90, 8e-10, text=L"\mathrm{OH^-}",
				  color=colors[6], fontsize=18, font = "sans-bold") 
			text!(ax, -1.15, 0.00024, text=L"\mathrm{CO}", 
			 	  color=colors[7], fontsize=18, font = "sans-bold")
		else 
			text!(ax, -1.2, 0.3, text=L"\mathrm{K^+}",
			  	  color=colors[1], fontsize=18, font="sans-bold")
			text!(ax, -0.90, 0.0000004, text=L"\mathrm{H^+}", 
				  color=colors[2], fontsize=18, font = "sans-bold") 
			text!(ax, -1.05, 2e-9, text=L"\mathrm{CO_3^{2-}}",
				  color=colors[4], fontsize=18, font = "sans-bold") 
			text!(ax, -1.0, 7e-6, text=L"\mathrm{HCO_3^-}", 
				  color=colors[3], fontsize=18, font = "sans-bold") 
			text!(ax, -1.21, 5e-4, text=L"\mathrm{CO_2}", 
				  color=colors[5], fontsize=18, font = "sans-bold")
			text!(ax, -0.9, 7e-9, text=L"\mathrm{OH^-}",
				  color=colors[6], fontsize=18, font = "sans-bold") 
			text!(ax, -1.15, 0.003, text=L"\mathrm{CO}", 
			 	  color=colors[7], fontsize=18, font = "sans-bold")
		end
    end

    display(fig)
end

# ╔═╡ 89112dd8-08c8-4887-90db-e2e9c5ce0d88
let
    csv_paths = [
        raw"../data/output/Concentration_Robin_Stefan_γ_Potassium_only.csv",
    ]

    species = ["K⁺","H⁺","HCO₃⁻","CO₃²⁻","CO₂","OH⁻","CO"]
    colors  = [:orange, :gray, :brown, :violet, :red, :green, :blue]

    fig = Figure(size=(960, 400))

    for i in 1:length(csv_paths)
        df = CSV.read(csv_paths[i], DataFrame)

        vgrid = df[!, :Voltage] 
        conc_electrode = permutedims(Matrix(df[!, Symbol.(species)]))  

        ax = Axis(fig[i, 1];
            xlabel = L"\mathbf{\text{U}\ \mathrm{vs.}\ \text{SHE}\ (V)}",
            ylabel = L"\mathbf{c_i^{+}}\;(\mathrm{M})",
            yscale = log10,
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

        #ax.xticklabelsvisible = false
        #ax.xlabelvisible = false
		
        for ia in 1:7
            y = max.(conc_electrode[ia, :], eps(Float64))
            lines!(ax, vgrid, y; color=colors[ia], linewidth=5, label=species[ia])
        end
			text!(ax, -1.15, 0.3, text=L"\mathrm{K^+}",
				  color=colors[1], fontsize=18, font="sans-bold")
			text!(ax, -0.90, 0.0000002, text=L"\mathrm{H^+}", 
				  color=colors[2], fontsize=18, font = "sans-bold") 
			text!(ax, -1.05, 3e-11, text=L"\mathrm{CO_3^{2-}}",
				  color=colors[4], fontsize=18, font = "sans-bold") 
			text!(ax, -1.0, 1.5e-7, text=L"\mathrm{HCO_3^-}", 
				  color=colors[3], fontsize=18, font = "sans-bold") 
			text!(ax, -1.2, 5.2e-6, text=L"\mathrm{CO_2}", 
				  color=colors[5], fontsize=18, font = "sans-bold")
			text!(ax, -0.90, 10e-10, text=L"\mathrm{OH^-}",
				  color=colors[6], fontsize=18, font = "sans-bold") 
			text!(ax, -1.17, 0.00024, text=L"\mathrm{CO}", 
			 	  color=colors[7], fontsize=18, font = "sans-bold")
    end

    display(fig)
end

# ╔═╡ 17b95990-f6e8-4da5-909a-1a2d46d7494d
let
     FS_BIG   = 16   
     FS_SMALL = 16      
     SP_BIG   = 3.5
     SP_SMALL = 3.5
     TICKW_BIG   = 2.0
     TICKW_SMALL = 2.0
     TICKL_BIG   = 8
     TICKL_SMALL = 8

    function style_axis!(ax; big::Bool)
        if big
            ax.spinewidth = SP_BIG
            ax.xtickwidth = TICKW_BIG
            ax.ytickwidth = TICKW_BIG
            ax.xticksize  = TICKL_BIG
            ax.yticksize  = TICKL_BIG
            ax.xlabelsize = FS_BIG
            ax.ylabelsize = FS_BIG
            ax.xticklabelsize = FS_BIG
            ax.yticklabelsize = FS_BIG
            ax.xlabelpadding = 10
            ax.ylabelpadding = 10
            ax.xlabelfont = :bold
        else
            ax.spinewidth = SP_SMALL
            ax.xtickwidth = TICKW_SMALL
            ax.ytickwidth = TICKW_SMALL
            ax.xticksize  = TICKL_SMALL
            ax.yticksize  = TICKL_SMALL
            ax.xlabelsize = FS_SMALL
            ax.ylabelsize = FS_SMALL
            ax.xticklabelsize = FS_SMALL
            ax.yticklabelsize = FS_SMALL
        end
        ax.xgridvisible = false
        ax.ygridvisible = false
        return ax
    end


    df_conc = CSV.read(raw"../data/output/Concentration_Dirichlet_Stefan_γ_Potassium_only.csv", DataFrame)
    df_pol  = CSV.read(raw"../data/output/Polarization_Curve_Dirichlet_Stefan_γ_pnp_Potassium_only.csv", DataFrame)
    df_cdl  = CSV.read(raw"../data/output/DLCap_Dirichlet_Stefan_γ_pnp_Potassium_only_1.csv", DataFrame)


    fig = Figure(size = (1200, 720),figure_padding = (300, 300, 30, 45))

    ax_cdl = Axis(fig[1, 1],
        xlabel = L"\text{Voltage} \mathrm{(V)}",
		ylabel = L"\text{C_{dl}} \textrm{(μF\,cm^{-2})}",
        limits = (-0.8, 0.8, 0, 300),
    )
    style_axis!(ax_cdl; big=false)

    ax_pol = Axis(fig[1, 2],
        xlabel = L"\text{Voltage} \mathrm{(V)}",
		ylabel = L"|j| \; \textrm{(mA\,cm^{-2})}",
        limits = (-1.3, -0.4, 1e-10, 100),
		yscale = log10
				  
    )
    style_axis!(ax_pol; big=false)

    ax_con = Axis(fig[2, 1:2],
        xlabel = L"\text{Voltage} \mathrm{(V)}",
        ylabel = L"\mathbf{c_i^{+}}\;(\mathrm{M})",
        yscale = log10,
        limits = (-1.25, 0.00, 1e-11, 1e1),
    )
    style_axis!(ax_con; big=true)

	
    rowsize!(fig.layout, 1, Relative(0.33))
    rowsize!(fig.layout, 2, Relative(0.67))
    colsize!(fig.layout, 1, Relative(1))
    colsize!(fig.layout, 2, Relative(1))
    rowgap!(fig.layout, 14)
    colgap!(fig.layout, 22)


    species = ["K⁺","H⁺","HCO₃⁻","CO₃²⁻","CO₂","OH⁻","CO"]
    colors  = [:orange, :gray, :brown, :violet, :red, :green, :blue]

    vgrid = df_conc[!, :Voltage]
    conc_electrode = permutedims(Matrix(df_conc[!, Symbol.(species)]))  # (7, N)

    xt = [-1.2, -1.0, -0.8, -0.6]
    ax_con.xticks = (xt, [@sprintf("%.1f", x) for x in xt])
    ax_pol.xticks = (xt, [@sprintf("%.1f", x) for x in xt])


    yt_vals_con = 10.0 .^ (0:-3:-9)
    yt_labs_con = [L"10^{0}", L"10^{-3}", L"10^{-6}", L"10^{-9}"]
    ax_con.yticks = (yt_vals_con, yt_labs_con)

	yt_vals_pol = 10.0 .^ (0:-4:-15)
	yt_labs_pol = [L"10^{0}", L"10^{-5}", L"10^{-10}", L"10^{-15}"]
	ax_pol.yticks = (yt_vals_pol, yt_labs_pol)

	yt_vals_cdl = [10, 150, 300]
	yt_labs_cdl = [L"10", L"150", L"300"]
	ax_cdl.yticks = (yt_vals_cdl, yt_labs_cdl)
	xt = [-0.8, -0.4, 0.0, 0.4, 0.8]		
	ax_cdl.xticks = (xt, [@sprintf("%.1f", x) for x in xt])

	
    for ia in 1:7
        y = max.(conc_electrode[ia, :], eps(Float64))
        lines!(ax_con, vgrid, y; color=colors[ia], linewidth=3)
    end

	
	# Stefan_Robin_Potassium
    #text!(ax_con, -1.15, 0.3,      text=L"\mathrm{K^+}",        color=colors[1], fontsize=18)
    #text!(ax_con, -0.90, 2e-7,     text=L"\mathrm{H^+}",        color=colors[2], fontsize=18)
    #text!(ax_con, -1.05, 3e-11,    text=L"\mathrm{CO_3^{2-}}",  color=colors[4], fontsize=18)
    #text!(ax_con, -1.00, 1.5e-7,   text=L"\mathrm{HCO_3^-}",    color=colors[3], fontsize=18)
    #text!(ax_con, -1.20, 5.2e-6,   text=L"\mathrm{CO_2}",       color=colors[5], fontsize=18)
    #text!(ax_con, -0.90, 1e-9,     text=L"\mathrm{OH^-}",       color=colors[6], fontsize=18)
    #text!(ax_con, -1.17, 2.4e-4,   text=L"\mathrm{CO}",         color=colors[7], fontsize=18)
	
	text!(ax_con, -1.2, 0.3, text=L"\mathrm{K^+}",
			  color=colors[1], fontsize=18, font="sans-bold")
	text!(ax_con, -0.90, 0.0000004, text=L"\mathrm{H^+}", 
			  color=colors[2], fontsize=18, font = "sans-bold") 
	text!(ax_con, -1.05, 2e-9, text=L"\mathrm{CO_3^{2-}}",
			  color=colors[4], fontsize=18, font = "sans-bold") 
	text!(ax_con, -1.0, 7e-6, text=L"\mathrm{HCO_3^-}", 
			  color=colors[3], fontsize=18, font = "sans-bold") 
	text!(ax_con, -1.21, 5e-4, text=L"\mathrm{CO_2}", 
			  color=colors[5], fontsize=18, font = "sans-bold")
	text!(ax_con, -0.9, 7e-9, text=L"\mathrm{OH^-}",
			  color=colors[6], fontsize=18, font = "sans-bold") 
	text!(ax_con, -1.15, 0.003, text=L"\mathrm{CO}", 
		 	  color=colors[7], fontsize=18, font = "sans-bold")


	
    lines!(ax_pol,
        df_pol[!, :Voltage],
        (abs.(df_pol[!, :Current] .* (cm^2/mA))),
        color = :skyblue,
        linewidth = 3,
    )


    lines!(ax_cdl,
        df_cdl[!, :Voltage],
        df_cdl[!, :Capacitance] / (μF / cm^2),
        color = RGBAf(0.5, 0.2, 0.8, 0.5),
        linewidth = 3,
    )

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
# ╠═8545d818-d255-4e8a-af8f-72fdaf9d4bdd
# ╠═4911ad1e-e75a-4d41-83dc-675b3f26ce39
# ╠═dcb38fd2-23b6-481e-83b7-4f3ee332c4ee
# ╠═ee091126-671e-46d8-8ed5-9457aeefd8fb
# ╠═4b02a8c0-a854-436a-9fdc-99ecf637e858
# ╠═59b13ba5-c85d-4e31-bd06-e1541b426205
# ╠═12a083c3-dde1-4a99-8d66-c7c43afc6b68
# ╠═95844847-6f46-4840-979d-c54282295f41
# ╠═54d0654f-9add-421b-ba54-90f276e3a10e
# ╠═d7ad452c-d3c1-485c-92ac-0b006dfb5acd
# ╠═08bb3a34-88d9-4934-a9c5-a5d61fa77e7e
# ╠═89112dd8-08c8-4887-90db-e2e9c5ce0d88
# ╠═17b95990-f6e8-4da5-909a-1a2d46d7494d
