### A Pluto.jl notebook ###
# v0.20.25

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

# ╔═╡ b3f44506-52eb-491c-bc44-76c9c49a43cd
let
    fig = Figure(size = (1050, 500))
    ax = Axis(fig[1, 1];
        xlabel = L"\phi~(\mathrm{V~vs~SHE})",
        ylabel = L"j~(\mathrm{mA\,cm^{-2}})",
		#limits = ((-1.3, 0.9), (-7.5, 2.5)),
        xlabelsize = 25, ylabelsize = 25,
        xgridvisible = false,
        ygridvisible = false,
        spinewidth = 4.5,
        xtickwidth = 4.5,
        ytickwidth = 4.5,
        xticklabelsize = 25, yticklabelsize = 25,
    )


    df = CSV.read("../data/output/cv_profile_comp__σ_1000.0_Robin_DMGL_γ_pnp_All_species_Scanrate_0.05_Periods_1_sweep_range_-1.2-1.2_cv.csv", DataFrame)

    cols1 = :red

    plot_objs1 = Any[]
    labels1 = String[]


        line = lines!(
            ax,
            df.Voltage,
            df.Current./2;
            color = cols1,
            linewidth = 4
        )
        push!(plot_objs1, line)
    ax.xticks = -1.5:0.3:0.9
    #ax.yticks = 1:-3:-15
"""
    leg = Legend(fig[1, 1], plot_objs1, labels1, "Extracted";
        framevisible = false,
        halign = :right, valign = :bottom,
        labelsize = 20, titlesize = 23,
        padding = (0, 0, 0, 0),
        tellwidth = false, tellheight = false
    )
    translate!(leg.blockscene, -40, 40, 0)
"""
    fig
end

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
"""
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

"""
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

# ╔═╡ ccb13e5a-5ceb-4a56-9f3f-14307fc78115
let
    wanted    = ["0.1", "0.5", "1.0"]
    pressures = ["Ar sat", "0.1", "0.2", "0.3", "0.5", "0.6", "1.0"]

    raw    = CSV.read("../data/Langmuir_CV_data/Figure_3.csv", DataFrame; header=false)
    sub    = Matrix(raw[4:end, :])
    num    = map(x -> x === missing ? NaN : parse(Float64, x), sub)
    num_df = DataFrame(num, :auto)
    npairs = size(num_df, 2) ÷ 2

    # 원하는 압력에 해당하는 pair 인덱스만 추출 → [2, 5, 7]
    keep = findall(in(wanted), pressures[1:npairs])

    # 3개뿐이니 불연속(distinct) 색상
    cols1 = resample_cmap(:winter, npairs)

    fig, ax = with_theme(electrochemistry_theme()) do
        f = Figure(size = (1050, 500))
        a = Axis(f[1, 1];
                 xlabel = L"φ~(\mathrm{V~vs~SHE})",
                 ylabel = L"j~(\mathrm{mA\,cm^{-2}})",
                 limits = ((-1.3, 0.9), (-5.5, 1.8)),
                 xticks = -1.5:0.3:1.0,
                 yticks = 1:-1:-5)
        return f, a
    end

    plot_objs1 = []
    labels1    = String[]
    for (k, j) in enumerate(keep)
        xcol, ycol = 2j - 1, 2j
        push!(labels1, pressures[j] * " pCO₂ atm")
        line = lines!(ax, num_df[!, xcol], num_df[!, ycol];
                      color = cols1[mod1(k, length(cols1))], linewidth = 4)
        push!(plot_objs1, line)
    end

    leg = Legend(fig[1, 1], plot_objs1, labels1, "Experimental";
        framevisible = false,
        halign = :right, valign = :bottom, labelsize = 20, titlesize = 23,
        padding = (0, 0, 0, 0),
        tellwidth = false, tellheight = false)
    translate!(leg.blockscene, -40, 40, 0)
    fig
end

# ╔═╡ b0a4b942-5654-4b22-8849-90bf6c7f957f
let
    # 1. Figure & Axis Setup
    fig = Figure(size = (1100, 700))
    
    # Primary Axis: Concentrations (Log Scale)
    ax1 = Axis(fig[1, 1];
        yscale = log10, # Log scale is essential for multiple species
        xlabel = L"L~(\mu\mathrm{m})",
        ylabel = L"c_{i}^{\mathrm{surface}}~(\mathrm{M})",
        xlabelsize = 26, ylabelsize = 26,
        xticklabelsize = 22, yticklabelsize = 22,
        xgridvisible = false, ygridvisible = false,
        spinewidth = 3
    )

    # Secondary Axis: Current Density (Linear Scale)
    ax2 = Axis(fig[1, 1];
        ylabel = L"|j_{\mathrm{peak}}|~(\mathrm{mA\,cm^{-2}})",
        ylabelsize = 26, yticklabelsize = 22,
        ygridvisible = false,
        yaxisposition = :right,
        yticklabelcolor = :black,
        ylabelcolor = :black,
        spinewidth = 3
    )
    
    hidespines!(ax2, :l, :t, :b)
    hidexdecorations!(ax2)

    # 2. Data Loading & Extraction
    csv_filename = "../data/output/L_varied_uncompensated_.csv"
    df = CSV.read(csv_filename, DataFrame)
    L_list = sort(unique(df.L_um))
    
    # Identify all concentration columns
    conc_cols = filter(name -> occursin("_Surface_M", String(name)), names(df))
    
    # Prepare colormap for species
    species_colors = resample_cmap(:tab10, length(conc_cols))

    # 3. Processing and Plotting
    # A. Plot Each Species Concentration
    for (idx, col) in enumerate(conc_cols)
        clean_name = replace(String(col), "_Surface_M" => "")
        
        # Extract peak/representative value for each L
        vals = [maximum(df[df.L_um .== l, col]) for l in L_list]
        
        # Apply clamping for log scale safety
        vals_safe = map(v -> (v > 1e-15 ? v : 1e-15), vals)
        
        scatterlines!(ax1, L_list, vals_safe; 
            color = species_colors[idx], 
            linewidth = 3, 
            markersize = 10,
            label = clean_name)
    end

    # B. Plot Peak Current (as a reference)
    peak_currents = [maximum(abs.(df[df.L_um .== l, :Current_mAcm2])) for l in L_list]
    
    # Plot current with a distinct bold black dashed line
    p_current = scatterlines!(ax2, L_list, peak_currents; 
        color = (:black, 0.5), 
        linewidth = 5, 
        markersize = 14, 
        marker = :diamond,
        linestyle = :dash, 
        label = "Current")

    # 4. Final Styling & Legend
    # Adjust limits for log scale visibility
    ylims!(ax1, 1e-12, 10.0) 
    
    # Combined Legend
    Legend(fig[1, 2], ax1, "Species"; 
        framevisible = true, 
        labelsize = 18, 
        titlesize = 20)
    
    # Add a separate label for current if needed or include in legend
    Label(fig[0, 1], "Surface Concentration & Peak Current vs. Boundary Layer Thickness", 
          fontsize = 24, font = :bold)

    fig
end

# ╔═╡ fc8096e4-01ca-451a-87cd-7e2e0171a531
let
    # 1. Data Loading
    csv_filename = "../data/output/L_varied_compensated_timestep__σ_80.0_Robin_DMGL_γ_pnp_All_species_Scanrate_0.05_Periods_1_sweep_range_-1.2-0.8_cv.csv"
    df = CSV.read(csv_filename, DataFrame)
    L_list = sort(unique(df.L_um))
    
    # Identify concentration columns
    conc_cols = filter(name -> occursin("_Surface_M", String(name)), names(df))
    num_species = length(conc_cols)
    
    # 2. Figure Setup
    fig = Figure(size = (1000, 300 * (1 + num_species)))
    cols = resample_cmap(:viridis, length(L_list))
    
    axes = []

    # --- TOP PANEL: Current Density (Linear Scale) ---
    ax_current = Axis(fig[1, 1];
        ylabel = L"j~(\mathrm{mA\,cm^{-2}})",
        xticklabelsvisible = false,
        xgridvisible = false, ygridvisible = false,
        spinewidth = 4, xtickwidth = 4, ytickwidth = 4,
        ylabelsize = 24, yticklabelsize = 20
    )
    push!(axes, ax_current)

    # --- DYNAMIC PANELS: Concentrations (Safe Log Scale) ---
    for (i, col) in enumerate(conc_cols)
        clean_name = replace(String(col), "_Surface_M" => "")
        
        ax_conc = Axis(fig[i + 1, 1];
            yscale = log10,        
            ylabel = "[$clean_name] (M)", 
			limits = ((-1.3, 0.9), (1e-12, 1e4)),
			xlabel = i == num_species ? L"V_{\mathrm{eff}}~(\mathrm{V})" : "",
            xticklabelsvisible = i == num_species,
            xgridvisible = false, ygridvisible = false,
            spinewidth = 4, xtickwidth = 4, ytickwidth = 4,
            ylabelsize = 22, yticklabelsize = 20
        )
        # Set limits to avoid log(0) errors if data is empty or all clamped
        # limits!(ax_conc, nothing, (1e-12, 1.0)) 
        push!(axes, ax_conc)
    end

    # 3. Plotting Loop
    plot_elements = []
    for (j, L_val) in enumerate(L_list)
        subdf = df[df.L_um .== L_val, :]
        
        # Plot Current (Linear)
        ln = lines!(axes[1], subdf.Voltage_V, subdf.Current_mAcm2; 
            color = cols[j], linewidth = 4)
        push!(plot_elements, ln)

        # Plot Concentrations (Log with Clamping)
        for (i, col) in enumerate(conc_cols)
            y_raw = subdf[!, col]
            # Apply the clamping method to handle negative/zero values for log scale
            y_safe = map(c -> (c > 1e-12 ? c : 1e-12), y_raw)
            
            lines!(axes[i + 1], subdf.Voltage_V, y_safe; 
                color = cols[j], linewidth = 4)
        end
    end

    # 4. Global Styling
    linkxaxes!(axes...)
    
    Legend(fig[1, 1], plot_elements, ["L = $(round(l)) μm" for l in L_list], "Thickness";
        framevisible = false, halign = :right, valign = :top,
        tellwidth = false, tellheight = false, labelsize = 18, titlesize = 20)

    rowgap!(fig.layout, 10) 
    
    fig
end

# ╔═╡ 05eb8a6f-d7b5-4f37-a905-2209323afe2d
let
    # 1. Data Loading
    csv_filename = "../data/output/L_varied_compensated_timestep__σ_80.0_Robin_DMGL_γ_pnp_All_species_Scanrate_0.05_Periods_1_sweep_range_-1.2-0.8_cv.csv"
    df = CSV.read(csv_filename, DataFrame)
    L_list = sort(unique(df.L_um))
    
    # Identify concentration columns
    conc_cols = filter(name -> occursin("_Surface_M", String(name)), names(df))
    num_species = length(conc_cols)
    
    # 2. Figure Setup
    fig = Figure(size = (1000, 300 * (1 + num_species)))
    cols = resample_cmap(:viridis, length(L_list))
    
    axes = []

    # --- TOP PANEL: Current Density (Linear Scale) ---
    ax_current = Axis(fig[1, 1];
        ylabel = L"j~(\mathrm{mA\,cm^{-2}})",
        xticklabelsvisible = false,
        xgridvisible = false, ygridvisible = false,
        spinewidth = 4, xtickwidth = 4, ytickwidth = 4,
        ylabelsize = 24, yticklabelsize = 20
    )
    push!(axes, ax_current)

    # --- DYNAMIC PANELS: Concentrations (Safe Log Scale) ---
    for (i, col) in enumerate(conc_cols)
        clean_name = replace(String(col), "_Surface_M" => "")
        
        ax_conc = Axis(fig[i + 1, 1];
            yscale = log10,        
            ylabel = "[$clean_name] (M)", 
			limits = ((-1, 81), (1e-12, 1e4)),
            xlabel = i == num_species ? L"\text{Time}\;(s)" : "",
            xticklabelsvisible = i == num_species,
            xgridvisible = false, ygridvisible = false,
            spinewidth = 4, xtickwidth = 4, ytickwidth = 4,
            ylabelsize = 22, yticklabelsize = 20
        )
        # Set limits to avoid log(0) errors if data is empty or all clamped
        # limits!(ax_conc, nothing, (1e-12, 1.0)) 
        push!(axes, ax_conc)
    end

    # 3. Plotting Loop
    plot_elements = []
    for (j, L_val) in enumerate(L_list)
        subdf = df[df.L_um .== L_val, :]
        
        # Plot Current (Linear)
        ln = lines!(axes[1], subdf.Times_s, subdf.Current_mAcm2; 
            color = cols[j], linewidth = 4)
        push!(plot_elements, ln)

        # Plot Concentrations (Log with Clamping)
        for (i, col) in enumerate(conc_cols)
            y_raw = subdf[!, col]
            # Apply the clamping method to handle negative/zero values for log scale
            y_safe = map(c -> (c > 1e-128 ? c : 1e-128), y_raw)
            
            lines!(axes[i + 1], subdf.Times_s, y_safe; 
                color = cols[j], linewidth = 4)
        end
    end

    # 4. Global Styling
    linkxaxes!(axes...)
    
    Legend(fig[1, 1], plot_elements, ["L = $(round(l)) μm" for l in L_list], "Thickness";
        framevisible = false, halign = :right, valign = :top,
        tellwidth = false, tellheight = false, labelsize = 18, titlesize = 20)

    rowgap!(fig.layout, 10) 
    
    fig
end

# ╔═╡ 5bb2da63-cc7d-4fcf-98b0-d2f4fa1dc1ef
let
    # 1. Data Loading
    csv_filename = "../data/output/L_varied_uncompensated__σ_80.0_Robin_DMGL_γ_pnp_All_species_Scanrate_0.05_Periods_1_sweep_range_-1.2-0.8_cv.csv"
    df = CSV.read(csv_filename, DataFrame)
    L_list = sort(unique(df.L_um))
    
    # Identify concentration columns
    conc_cols = filter(name -> occursin("_Surface_M", String(name)), names(df))
    num_species = length(conc_cols)
    
    # 2. Figure Setup
    fig = Figure(size = (1000, 300 * (1 + num_species)))
    cols = resample_cmap(:viridis, length(L_list))
    
    axes = []

    # --- TOP PANEL: Current Density (Linear Scale) ---
    ax_current = Axis(fig[1, 1];
        ylabel = L"j~(\mathrm{mA\,cm^{-2}})",
        xticklabelsvisible = false,
        xgridvisible = false, ygridvisible = false,
        spinewidth = 4, xtickwidth = 4, ytickwidth = 4,
        ylabelsize = 24, yticklabelsize = 20
    )
    push!(axes, ax_current)

    # --- DYNAMIC PANELS: Concentrations (Safe Log Scale) ---
    for (i, col) in enumerate(conc_cols)
        clean_name = replace(String(col), "_Surface_M" => "")
        
        ax_conc = Axis(fig[i + 1, 1];
            yscale = log10,        
            ylabel = "[$clean_name] (M)", 
			limits = ((-1.3, 0.9), (1e-12, 1e4)),
            xlabel = i == num_species ? L"V_{\mathrm{eff}}~(\mathrm{V})" : "",
            xticklabelsvisible = i == num_species,
            xgridvisible = false, ygridvisible = false,
            spinewidth = 4, xtickwidth = 4, ytickwidth = 4,
            ylabelsize = 22, yticklabelsize = 20
        )
        # Set limits to avoid log(0) errors if data is empty or all clamped
        # limits!(ax_conc, nothing, (1e-12, 1.0)) 
        push!(axes, ax_conc)
    end

    # 3. Plotting Loop
    plot_elements = []
    for (j, L_val) in enumerate(L_list)
        subdf = df[df.L_um .== L_val, :]
        
        # Plot Current (Linear)
        ln = lines!(axes[1], subdf.Voltage_V, subdf.Current_mAcm2; 
            color = cols[j], linewidth = 4)
        push!(plot_elements, ln)

        # Plot Concentrations (Log with Clamping)
        for (i, col) in enumerate(conc_cols)
            y_raw = subdf[!, col]
            # Apply the clamping method to handle negative/zero values for log scale
            y_safe = map(c -> (c > 1e-128 ? c : 1e-128), y_raw)
            
            lines!(axes[i + 1], subdf.Voltage_V, y_safe; 
                color = cols[j], linewidth = 4)
        end
    end

    # 4. Global Styling
    linkxaxes!(axes...)
    
    Legend(fig[1, 1], plot_elements, ["L = $(round(l)) μm" for l in L_list], "Thickness";
        framevisible = false, halign = :right, valign = :top,
        tellwidth = false, tellheight = false, labelsize = 18, titlesize = 20)

    rowgap!(fig.layout, 10) 
    
    fig
end

# ╔═╡ 2ddeace6-663d-4fd1-8494-f1c88c19e628
let
    # 1. Data Loading
    csv_filename = "../data/output/L_varied_uncompensated__σ_80.0_Robin_DMGL_γ_pnp_All_species_Scanrate_0.05_Periods_1_sweep_range_-1.2-0.8_cv.csv"
    df = CSV.read(csv_filename, DataFrame)
    L_list = sort(unique(df.L_um))
    
    # Identify concentration columns
    conc_cols = filter(name -> occursin("_Surface_M", String(name)), names(df))
    num_species = length(conc_cols)
    
    # 2. Figure Setup
    fig = Figure(size = (1000, 300 * (1 + num_species)))
    cols = resample_cmap(:viridis, length(L_list))
    
    axes = []

    # --- TOP PANEL: Current Density (Linear Scale) ---
    ax_current = Axis(fig[1, 1];
        ylabel = L"j~(\mathrm{mA\,cm^{-2}})",
        xticklabelsvisible = false,
        xgridvisible = false, ygridvisible = false,
        spinewidth = 4, xtickwidth = 4, ytickwidth = 4,
        ylabelsize = 24, yticklabelsize = 20
    )
    push!(axes, ax_current)

    # --- DYNAMIC PANELS: Concentrations (Safe Log Scale) ---
    for (i, col) in enumerate(conc_cols)
        clean_name = replace(String(col), "_Surface_M" => "")
        
        ax_conc = Axis(fig[i + 1, 1];
            yscale = log10,        
            ylabel = "[$clean_name] (M)", 
			limits = ((-1, 81), (1e-12, 1e4)),
            xlabel = i == num_species ? L"\text{Time}\;(s)" : "",
            xticklabelsvisible = i == num_species,
            xgridvisible = false, ygridvisible = false,
            spinewidth = 4, xtickwidth = 4, ytickwidth = 4,
            ylabelsize = 22, yticklabelsize = 20, xlabelsize = 28
        )
        # Set limits to avoid log(0) errors if data is empty or all clamped
        # limits!(ax_conc, nothing, (1e-12, 1.0)) 
        push!(axes, ax_conc)
    end

    # 3. Plotting Loop
    plot_elements = []
    for (j, L_val) in enumerate(L_list)
        subdf = df[df.L_um .== L_val, :]
        
        # Plot Current (Linear)
        ln = lines!(axes[1], subdf.Times_s, subdf.Current_mAcm2; 
            color = cols[j], linewidth = 4)
        push!(plot_elements, ln)

        # Plot Concentrations (Log with Clamping)
        for (i, col) in enumerate(conc_cols)
            y_raw = subdf[!, col]
            # Apply the clamping method to handle negative/zero values for log scale
            y_safe = map(c -> (c > 1e-128 ? c : 1e-128), y_raw)
            
            lines!(axes[i + 1], subdf.Times_s, y_safe; 
                color = cols[j], linewidth = 4)
        end
    end

    # 4. Global Styling
    linkxaxes!(axes...)
    
    Legend(fig[1, 1], plot_elements, ["L = $(round(l)) μm" for l in L_list], "Thickness";
        framevisible = false, halign = :right, valign = :top,
        tellwidth = false, tellheight = false, labelsize = 18, titlesize = 20)

    rowgap!(fig.layout, 10) 
    
    fig
end

# ╔═╡ caa05490-0c7d-44ec-9be8-73f7a4473d8a
let
    # 1. Data Loading
    csv_filename = "../data/output/L_varied_uncompensated_.csv"
    df = CSV.read(csv_filename, DataFrame)
    L_list = sort(unique(df.L_um))
    
    # Identify concentration columns
    conc_cols = filter(name -> occursin("_Surface_M", String(name)), names(df))
    num_species = length(conc_cols)
    
    # 2. Figure Setup
    fig = Figure(size = (1000, 300 * (1 + num_species)))
    cols = resample_cmap(:viridis, length(L_list))
    
    axes = []

    # --- TOP PANEL: Current Density (Linear Scale) ---
    ax_current = Axis(fig[1, 1];
        ylabel = L"j~(\mathrm{mA\,cm^{-2}})",
        xticklabelsvisible = false,
        xgridvisible = false, ygridvisible = false,
        spinewidth = 4, xtickwidth = 4, ytickwidth = 4,
        ylabelsize = 24, yticklabelsize = 20
    )
    push!(axes, ax_current)

    # --- DYNAMIC PANELS: Concentrations (Safe Log Scale) ---
    for (i, col) in enumerate(conc_cols)
        clean_name = replace(String(col), "_Surface_M" => "")
        
        ax_conc = Axis(fig[i + 1, 1];
            yscale = log10,        
            ylabel = "[$clean_name] (M)", 
			limits = ((-1, 81), (1e-12, 1e4)),
            xlabel = i == num_species ? L"V_{\mathrm{eff}}~(\mathrm{V})" : "",
            xticklabelsvisible = i == num_species,
            xgridvisible = false, ygridvisible = false,
            spinewidth = 4, xtickwidth = 4, ytickwidth = 4,
            ylabelsize = 22, yticklabelsize = 20
        )
        # Set limits to avoid log(0) errors if data is empty or all clamped
        # limits!(ax_conc, nothing, (1e-12, 1.0)) 
        push!(axes, ax_conc)
    end

    # 3. Plotting Loop
    plot_elements = []
    for (j, L_val) in enumerate(L_list)
        subdf = df[df.L_um .== L_val, :]
        
        # Plot Current (Linear)
        ln = lines!(axes[1], subdf.Times_s, subdf.Current_mAcm2; 
            color = cols[j], linewidth = 4)
        push!(plot_elements, ln)

        # Plot Concentrations (Log with Clamping)
        for (i, col) in enumerate(conc_cols)
            y_raw = subdf[!, col]
            # Apply the clamping method to handle negative/zero values for log scale
            y_safe = map(c -> (c > 1e-128 ? c : 1e-128), y_raw)
            
            lines!(axes[i + 1], subdf.Times_s, y_safe; 
                color = cols[j], linewidth = 4)
        end
    end

    # 4. Global Styling
    linkxaxes!(axes...)
    
    Legend(fig[1, 1], plot_elements, ["L = $(round(l)) μm" for l in L_list], "Thickness";
        framevisible = false, halign = :right, valign = :top,
        tellwidth = false, tellheight = false, labelsize = 18, titlesize = 20)

    rowgap!(fig.layout, 10) 
    
    fig
end

# ╔═╡ 2b10e19b-c099-4e4e-9cc1-5be1255d03bd
1023

# ╔═╡ 4a68d318-ab6f-45b3-ad63-33d4a77c534d
let
    fig = Figure(size = (1050, 500))
    ax = Axis(fig[1, 1];
        xlabel = L"\phi~(\mathrm{V~vs~SHE})",
        ylabel = L"j~(\mathrm{mA\,cm^{-2}})",
        limits = ((-1.3, 0.9), (-7.5, 1.8)),
        xlabelsize = 25, ylabelsize = 25,
        xgridvisible = false,
        ygridvisible = false,
        spinewidth = 4.5,
        xtickwidth = 4.5,
        ytickwidth = 4.5,
        xticklabelsize = 25, yticklabelsize = 25,
    )

    # --------------------------------------------------
    # extracted CSV 읽기
    # --------------------------------------------------
    df = CSV.read("../data/output/pressure_varied_σ_80.0Robin_DMGL_γ_pnp_All_species_Scanrate_0.05_Periods_1.csv", DataFrame)

    # pressure
    plist = sort(unique(df.Pressure))


    cols1 = resample_cmap(:winter, length(plist))

    plot_objs1 = Any[]
    labels1 = String[]

    for (j, p) in enumerate(plist)
        subdf = df[df.Pressure .== p, :]


        label = @sprintf("%.1f pCO₂ atm", p)
        push!(labels1, label)

        line = lines!(
            ax,
            subdf.Voltage,
            subdf.Value./2;
            color = cols1[j],
            linewidth = 4
        )
        push!(plot_objs1, line)
    end
"""
    # --------------------------------------------------
    # reaction region shading
    # --------------------------------------------------
    vspan!(ax, -1.30, -0.70, color=(colorant"#87CEFA", 0.30))
    text!(ax, -1.00, 1.55,
        text = "CO₂ reduction reaction",
        font = "sans-bold",
        align = (:center, :top),
        fontsize = 16,
        color = "#1E90FF"
    )

    vspan!(ax, -0.20, 0.10, color=(colorant"#F7DC6F", 0.30))
    text!(ax, -0.05, 1.55,
        text = "CO oxidation\nreaction\nwith CO₃²⁻",
        font = "sans-bold",
        align = (:center, :top),
        fontsize = 16,
        color = :orange
    )

    vspan!(ax, 0.10, 0.90, color=(colorant"#F1948A", 0.30))
    text!(ax, 0.50, 1.55,
        text = "CO oxidation reaction with H₂O",
        font = "sans-bold",
        align = (:center, :top),
        fontsize = 16,
        color = "#FF6347"
    )
"""
    ax.xticks = -1.2:0.3:0.9
    ax.yticks = 1:-2:-7

    leg = Legend(fig[1, 1], plot_objs1, labels1, "Extracted";
        framevisible = false,
        halign = :right, valign = :bottom,
        labelsize = 20, titlesize = 23,
        padding = (0, 0, 0, 0),
        tellwidth = false, tellheight = false
    )
    translate!(leg.blockscene, -40, 40, 0)

    fig
end

# ╔═╡ 8545d818-d255-4e8a-af8f-72fdaf9d4bdd
let
    fig = Figure(size = (1050, 500))

    ax = Axis(fig[1, 1];
        xlabel = L"\phi~(\mathrm{V~vs~SHE})",
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

    # --------------------------------------------------
    # experimental CSV 읽기
    # --------------------------------------------------
    raw = CSV.read("../data/Langmuir_CV_data/Figure_3.csv", DataFrame; header=false)

    sub = Matrix(raw[4:end, :])
    num = map(x -> x === missing ? NaN : parse(Float64, x), sub)
    num_df = DataFrame(num, :auto)

    # Ar sat 제외
    pressures_exp = ["0.1", "0.2", "0.3", "0.5", "0.6", "1.0"]
    exp_indices = 2:7   # 1번째 pair(Ar sat) 제외

    cols1 = resample_cmap(:winter, length(exp_indices))

    plot_objs1 = Any[]
    labels1 = String[]

    for (k, j) in enumerate(exp_indices)
        xcol, ycol = 2j - 1, 2j
        label = pressures_exp[k] * " pCO₂ atm"
        push!(labels1, label)

        line = lines!(
            ax,
            num_df[!, xcol],
            num_df[!, ycol];
            color = cols1[k],
            linewidth = 4
        )
        push!(plot_objs1, line)
    end

    ax_ex = Axis(fig[1, 2];
        xlabel = L"\phi~(\mathrm{V~vs~SHE})",
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

    # --------------------------------------------------
    # simulation CSV 읽기
    # --------------------------------------------------
    df = CSV.read("../data/output/pressure_varied_σ_1000.0Robin_DMGL_γ_pnp_Potassium_only.csv", DataFrame)

    plist = sort(unique(df.Pressure))
    cols2 = resample_cmap(:winter, length(plist))

    plot_objs2 = Any[]
    labels2 = String[]

    for (j, p) in enumerate(plist)
        subdf = df[df.Pressure .== p, :]

        label = @sprintf("%.1f pCO₂ atm", p)
        push!(labels2, label)

        line = lines!(
            ax_ex,
            subdf.Voltage,
            subdf.Value ./ 2;
            color = cols2[j],
            linewidth = 4
        )
        push!(plot_objs2, line)
    end

    leg1 = Legend(fig[1, 1], plot_objs1, labels1, "Experiment";
        framevisible = false,
        halign = :right, valign = :bottom,
        labelsize = 20, titlesize = 23,
        tellwidth = false, tellheight = false
    )

    leg2 = Legend(fig[1, 2], plot_objs2, labels2, "Simulation";
        framevisible = false,
        halign = :right, valign = :bottom,
        labelsize = 20, titlesize = 23,
        tellwidth = false, tellheight = false
    )

    translate!(leg1.blockscene, -40, 40, 0)
    translate!(leg2.blockscene, -40, 40, 0)
	ax.xticks = -1.5:0.6:1.0 	; 	ax_ex.xticks = -1.5:0.6:1.0
	ax.yticks = 1:-1:-5 		;	ax_ex.yticks = 1:-1:-5

	
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

# ╔═╡ a09d4a6a-1017-4792-b684-7f9ddfeba83b
begin
		const FS_BIG   = 18
	    const FS_SMALL = 18
	    const SP_BIG   = 3.5
	    const SP_SMALL = 3.5
	    const TICKW_BIG   = 2.0
	    const TICKW_SMALL = 2.0
	    const TICKL_BIG   = 8
	    const TICKL_SMALL = 8
end

# ╔═╡ 5a090c2f-eb24-4453-b650-f4bc6c32c952
let
    df_prev = CSV.read(raw"../data/output/Concentration_Robin_Stefan_γ_Same_Size.csv", DataFrame)
    df_curr = CSV.read(raw"../data/output/Concentration_Robin_DMGL_γ_Same_Size.csv", DataFrame)

    species = ["K⁺","H⁺","HCO₃⁻","CO₃²⁻","CO₂","OH⁻","CO"]
    colors  = [:orange, :gray, :brown, :violet, :red, :green, :blue]

    V = df_curr[!, "Voltage"]

    floors = Dict{String, Float64}()
    for s in species
        vals = vcat(df_prev[!, s], df_curr[!, s])
        pos  = vals[vals .> 0]
        floors[s] = isempty(pos) ? 1e-30 : minimum(pos) * 0.1
    end

    fig = Figure(size = (950, 750), figure_padding = (25, 25, 30, 25))

    

    function format_axis!(ax; big=true)
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
            ax.ylabelfont = :bold
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

    ax_main = Axis(
        fig[1, 1],
        #xlabel = L"\text{Voltage} \mathrm{(V)}",
        ylabel = L"\textbf{Interfacial\ concentration}\qquad \mathbf{c_i^{\ddagger}}\;(\mathrm{M})",
        yscale = log10,
    )
		#xlabel = false
	ax_main.xticklabelsvisible = false
    ax_aux = Axis(
        fig[2, 1],
        xlabel = L"\text{Voltage} \mathrm{(V)}",
        ylabel = L"\textbf{log10(current / previous)}",
    )

    format_axis!(ax_main, big=true)
    format_axis!(ax_aux,  big=true)

    for (i, s) in enumerate(species)
        prev = df_prev[!, s]
        curr = df_curr[!, s]
        f    = floors[s]

        lines!(ax_main, V, curr;
            color = colors[i],
            linestyle = :solid,
            linewidth = 3.5
        )
		
        lines!(ax_main, V, prev;
            color = (:gray60, 0.8),
            linestyle = :dash,
            linewidth = 2
        )

        delta = log10.((curr .+ f) ./ (prev .+ f))
        lines!(ax_aux, V, delta;
            color = colors[i],
            linewidth = 3
        )
    end

    hlines!(ax_aux, [0.0]; color = :black, linestyle = :dot, linewidth = 2)

    linkxaxes!(ax_main, ax_aux)
    ylims!(ax_aux, -0.075, 0.075)
    ylims!(ax_main, 1e-12, 1e2)
    xlims!(ax_aux, -1.2, -0.6)

    text!(ax_main, -1.15, 3e-1, text=L"\mathrm{K^+}",
        color=colors[1], fontsize=18, font=:bold)
    text!(ax_main, -0.90, 2e-8, text=L"\mathrm{H^+}",
        color=colors[2], fontsize=18, font=:bold)
    text!(ax_main, -1.05, 5e-11, text=L"\mathrm{CO_3^{2-}}",
        color=colors[4], fontsize=18, font=:bold)
    text!(ax_main, -1.05, 1e-7, text=L"\mathrm{HCO_3^-}",
        color=colors[3], fontsize=18, font=:bold)
    text!(ax_main, -1.06, 9.2e-6, text=L"\mathrm{CO_2}",
        color=colors[5], fontsize=18, font=:bold)
    text!(ax_main, -0.90, 1e-9, text=L"\mathrm{OH^-}",
        color=colors[6], fontsize=18, font=:bold)
    text!(ax_main, -1.17, 2.4e-4, text=L"\mathrm{CO}",
        color=colors[7], fontsize=18, font=:bold)
	ax_main.yticks = ([1e-12, 1e-8, 1e-4, 1], [L"10^{-12}", L"10^{-8}", L"10^{-4}", L"10^{0}"])

    rowsize!(fig.layout, 1, Relative(0.72))
    rowsize!(fig.layout, 2, Relative(0.28))
    rowgap!(fig.layout, 14)
    colgap!(fig.layout, 22)

    resize_to_layout!(fig)
    fig
end

# ╔═╡ 224cf209-4aaa-49ca-a577-d77a597d256c
let
    df_prev = CSV.read(raw"../data/output/Activity_Curve_Robin_Stefan_γ_pnp_Same_Size.csv", DataFrame)
    df_curr = CSV.read(raw"../data/output/Activity_Curve_Robin_DMGL_γ_pnp_Same_Size.csv", DataFrame)

    species = ["K⁺","H⁺","HCO₃⁻","CO₃²⁻","CO₂","OH⁻","CO"]
    colors  = [:orange, :gray, :brown, :violet, :red, :green, :blue]

    V = df_curr[!, "Voltage"]

    floors = Dict{String, Float64}()
    for s in species
        vals = vcat(df_prev[!, s], df_curr[!, s])
        pos  = vals[vals .> 0]
        floors[s] = isempty(pos) ? 1e-30 : minimum(pos) * 0.1
    end

    fig = Figure(size = (950, 750), figure_padding = (25, 25, 30, 25))

    

    function format_axis!(ax; big=true)
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
            ax.ylabelfont = :bold
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

    ax_main = Axis(
        fig[1, 1],
        #xlabel = L"\text{Voltage} \mathrm{(V)}",
        ylabel = L"\textbf{Interfacial\ concentration}\qquad \mathbf{c_i^{\ddagger}}\;(\mathrm{M})",
        yscale = log10,
    )
		#xlabel = false
	ax_main.xticklabelsvisible = false
    ax_aux = Axis(
        fig[2, 1],
        xlabel = L"\text{Voltage} \mathrm{(V)}",
        ylabel = L"\textbf{log10(current / previous)}",
    )

    format_axis!(ax_main, big=true)
    format_axis!(ax_aux,  big=true)

    for (i, s) in enumerate(species)
        prev = df_prev[!, s]
        curr = df_curr[!, s]
        f    = floors[s]

        lines!(ax_main, V, curr;
            color = colors[i],
            linestyle = :solid,
            linewidth = 3.5
        )
		
        lines!(ax_main, V, prev;
            color = (:gray60, 0.8),
            linestyle = :dash,
            linewidth = 2
        )

        delta = log10.((curr .+ f) ./ (prev .+ f))
        lines!(ax_aux, V, delta;
            color = colors[i],
            linewidth = 3
        )
    end

    hlines!(ax_aux, [0.0]; color = :black, linestyle = :dot, linewidth = 2)

    linkxaxes!(ax_main, ax_aux)
    ylims!(ax_aux, -0.075, 0.075)
    ylims!(ax_main, 1e-10, 1e4)
    xlims!(ax_aux, -1.2, -0.6)

    text!(ax_main, -1.15, 50, text=L"\mathrm{K^+}",
        color=colors[1], fontsize=18, font=:bold)
    text!(ax_main, -0.875, 7e-7, text=L"\mathrm{H^+}",
        color=colors[2], fontsize=18, font=:bold)
    text!(ax_main, -1.05, 4e-9, text=L"\mathrm{CO_3^{2-}}",
        color=colors[4], fontsize=18, font=:bold)
    text!(ax_main, -1.05, 7e-5, text=L"\mathrm{HCO_3^-}",
        color=colors[3], fontsize=18, font=:bold)
    text!(ax_main, -1.10, 4.2e-4, text=L"\mathrm{CO_2}",
        color=colors[5], fontsize=18, font=:bold)
    text!(ax_main, -0.875, 1e-8, text=L"\mathrm{OH^-}",
        color=colors[6], fontsize=18, font=:bold)
    text!(ax_main, -1.17, 3.0e-2, text=L"\mathrm{CO}",
        color=colors[7], fontsize=18, font=:bold)

	ax_main.yticks = ([1e-10, 1e-6, 1e-2, 1e2], [L"10^{-10}", L"10^{-6}", L"10^{-2}", L"10^{2}"])


	
    rowsize!(fig.layout, 1, Relative(0.72))
    rowsize!(fig.layout, 2, Relative(0.28))
    rowgap!(fig.layout, 14)
    colgap!(fig.layout, 22)

    resize_to_layout!(fig)
    fig
end

# ╔═╡ 400709f7-2a7c-431a-96ce-99a678b840f7
let
    df_prev = CSV.read(raw"../data/catmap_CO2R_data/voltage-conc.csv", DataFrame)
	    df_curr = CSV.read(raw"../data/output/Concentration_Robin_Stefan_γ_Potassium_only.csv", DataFrame)

    species = ["K⁺","H⁺","HCO₃⁻","CO₃²⁻","CO₂","OH⁻","CO"]
    colors  = [:orange, :gray, :brown, :violet, :red, :green, :blue]

    V = df_curr[!, "Voltage"]

    fig = Figure(size = (950, 600), figure_padding = (25, 25, 30, 25))

    function format_axis!(ax; big=true)
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
            ax.ylabelfont = :bold
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

    ax_main = Axis(
        fig[1, 1],
        xlabel = L"\text{Voltage} \mathrm{(V)}",
        ylabel = L"\textbf{Interfacial\ concentration}\qquad \mathbf{c_i^{\ddagger}}\;(\mathrm{M})",
        yscale = log10,
		limits = (-1.25, -0.55, 1e-11, 1e1),

    )

    format_axis!(ax_main, big=true)

    # df_curr: species별 곡선
    for (i, s) in enumerate(species)
        curr = df_curr[!, s]
        lines!(ax_main, V, curr;
            color = colors[i],
            linewidth = 3.5
        )
    end

    # df_prev: Index=1:7 을 species 순서와 매칭해서 점으로 표시
    for i in 1:length(species)
        row = df_prev[df_prev.Index .== i, :]
        if nrow(row) > 0
            scatter!(ax_main,
                row[!, "Voltage"],
                row[!, "Concentration"] ;
                color = colors[i],
                markersize = 5,
                #strokecolor = :black,
                #strokewidth = 1.5
            )
        end
    end

    #text!(ax_main, -1.15, 50, text=L"\mathrm{K^+}",
    #    color=colors[1], fontsize=18, font=:bold)
    #text!(ax_main, -0.875, 7e-7, text=L"\mathrm{H^+}",
    #    color=colors[2], fontsize=18, font=:bold)
    #text!(ax_main, -1.05, 4e-9, text=L"\mathrm{CO_3^{2-}}",
    #    color=colors[4], fontsize=18, font=:bold)
    #text!(ax_main, -1.05, 7e-5, text=L"\mathrm{HCO_3^-}",
    #    color=colors[3], fontsize=18, font=:bold)
    #text!(ax_main, -1.10, 4.2e-4, text=L"\mathrm{CO_2}",
    #    color=colors[5], fontsize=18, font=:bold)
    #text!(ax_main, -0.875, 1e-8, text=L"\mathrm{OH^-}",
    #    color=colors[6], fontsize=18, font=:bold)
    #text!(ax_main, -1.17, 3.0e-2, text=L"\mathrm{CO}",
    #    color=colors[7], fontsize=18, font=:bold)

    ax_main.yticks = (
        [1e-10, 1e-6, 1e-2, 1e2],
        [L"10^{-10}", L"10^{-6}", L"10^{-2}", L"10^{2}"]
    )

    #xlims!(ax_main, -1.2, -0.6)
    #ylims!(ax_main, 1e-10, 1e4)

    resize_to_layout!(fig)
    fig
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
    FS_BIG   = 24
    FS_SMALL = 24
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
            ax.ylabelfont = :bold
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

    priv_color = (:gray60, 0.8)

    df_conc = CSV.read(raw"../data/output/Concentration_Robin_Stefan_γ_Same_Size.csv", DataFrame)
    df_pol  = CSV.read(raw"../data/output/Polarization_Curve_Robin_Stefan_γ_pnp_Same_Size.csv", DataFrame)
    df_cdl  = CSV.read(raw"../data/output/DLCap_Robin_Stefan_γ_pb_Same_Size_1.csv", DataFrame)
    df_act  = CSV.read(raw"../data/output/Activity_Curve_Robin_Stefan_γ_pnp_Same_Size.csv", DataFrame)

    df_conc_pr = CSV.read(raw"../data/output/Concentration_Robin_DMGL_γ_Same_Size.csv", DataFrame)
    df_pol_pr  = CSV.read(raw"../data/output/Polarization_Curve_Robin_DMGL_γ_pnp_Same_Size.csv", DataFrame)
    df_cdl_pr  = CSV.read(raw"../data/output/DLCap_Robin_DMGL_γ_pb_Same_Size_1.csv", DataFrame)
    df_act_pr  = CSV.read(raw"../data/output/Activity_Curve_Robin_DMGL_γ_pnp_Same_Size.csv", DataFrame)


	
    fig = Figure(size = (1200, 820), figure_padding = (40, 40, 30, 45))

    ax_cdl = Axis(fig[1, 1],
        xlabel = L"\quad φ\; \mathrm{(V vs SHE)}",
        ylabel = L"\text{C_{dl}} \textrm{(μF\,cm^{-2})}",
        limits = (-0.8, 0.8, 0, 25),
    )
    style_axis!(ax_cdl; big=false)

    ax_pol = Axis(fig[1, 2],
        xlabel = L"\quad φ\; \mathrm{(V vs SHE)}",
        ylabel = L"|j| \; \textrm{(mA\,cm^{-2})}",
        limits = (-1.3, -0.4, 1e-10, 100),
		yaxisposition = :right,
        yscale = log10
    )
    style_axis!(ax_pol; big=false)

    ax_con = Axis(fig[2, 1],
        xlabel = L"\quad φ\; \mathrm{(V vs SHE)}",
        ylabel = L"\text{Interfacial\ concentration}\; \mathbf{c_i^{\ddagger}}\;(\mathrm{M})",
        yscale = log10,
        limits = (-1.25, -0.5, 1e-11, 1e1),
    )
    style_axis!(ax_con; big=true)

    ax_act = Axis(fig[2, 2],
        xlabel = L"\quad φ\; \mathrm{(V vs SHE)}",
        ylabel = L"\text{Interfacial\ activity}\qquad \mathbf{a_i^{\ddagger}}",
        yscale = log10,
        limits = (-1.25, -0.5, 1e-11, 1e4),
    )
    style_axis!(ax_act; big=true)

    rowsize!(fig.layout, 1, Relative(0.28))
    rowsize!(fig.layout, 2, Relative(0.72))
    rowgap!(fig.layout, 14)
    colgap!(fig.layout, 22)

    species = ["K⁺","H⁺","HCO₃⁻","CO₃²⁻","CO₂","OH⁻","CO"]
    colors  = [:orange, :gray, :brown, :violet, :red, :green, :blue]
	colors_pastel = [
	    "#F7C97F",  
	    "#D3D3D3",  
	    "#D9C2A7", 
	    "#FFC5D3",  
	    "#FF746C",  
	    "#80EF80",  
	    "#AFCBFF",  
	]

    vgrid_con = df_conc[!, :Voltage]
    conc_electrode = permutedims(Matrix(df_conc[!, Symbol.(species)]))

    vgrid_act = df_act[!, :Voltage]
    act_electrode = permutedims(Matrix(df_act[!, Symbol.(species)]))

    vgrid_con_pr = df_conc_pr[!, :Voltage]
    conc_electrode_pr = permutedims(Matrix(df_conc_pr[!, Symbol.(species)]))

    vgrid_act_pr = df_act_pr[!, :Voltage]
    act_electrode_pr = permutedims(Matrix(df_act_pr[!, Symbol.(species)]))

    xt_bottom = [-1.2, -1.0, -0.8, -0.6]
    ax_con.xticks = (xt_bottom, [@sprintf("%.1f", x) for x in xt_bottom])
    ax_act.xticks = (xt_bottom, [@sprintf("%.1f", x) for x in xt_bottom])
    ax_pol.xticks = (xt_bottom, [@sprintf("%.1f", x) for x in xt_bottom])

    yt_vals_con = 10.0 .^ (0:-3:-9)
    yt_labs_con = [L"10^{0}", L"10^{-3}", L"10^{-6}", L"10^{-9}"]
    ax_con.yticks = (yt_vals_con, yt_labs_con)

    yt_vals_act = 10.0 .^ (3:-3:-9)
    yt_labs_act = [L"10^{3}", L"10^{0}", L"10^{-3}", L"10^{-6}", L"10^{-9}"]
    ax_act.yticks = (yt_vals_act, yt_labs_act)

    yt_vals_pol = 10.0 .^ (0:-5:-10)
    yt_labs_pol = [L"10^{0}", L"10^{-5}", L"10^{-10}"]
    ax_pol.yticks = (yt_vals_pol, yt_labs_pol)

    yt_vals_cdl = [0, 10, 20, 30]
    yt_labs_cdl = [L"0", L"10", L"20", L"30"]
    ax_cdl.yticks = (yt_vals_cdl, yt_labs_cdl)

    xt_top = [-0.8, -0.4, 0.0, 0.4, 0.8]
    ax_cdl.xticks = (xt_top, [@sprintf("%.1f", x) for x in xt_top])

	solid_linewidth = 3.5
	dash_linewidth = 2
	
    for ia in 1:7
        lines!(ax_con, vgrid_con, max.(conc_electrode[ia, :], eps(Float64));
            color=colors[ia], linewidth=solid_linewidth)
        lines!(ax_act, vgrid_act, max.(act_electrode[ia, :], eps(Float64));
            color=colors[ia], linewidth=solid_linewidth)
        lines!(ax_con, vgrid_con_pr, max.(conc_electrode_pr[ia, :], eps(Float64));
            color=colors_pastel[ia], linewidth=dash_linewidth, linestyle=:dash)
        lines!(ax_act, vgrid_act_pr, max.(act_electrode_pr[ia, :], eps(Float64));
            color=colors_pastel[ia], linewidth=dash_linewidth, linestyle=:dash)
    end

    lines!(ax_pol, df_pol[!, :Voltage], abs.(df_pol[!, :Current] .* (cm^2/mA));
        color="#8ED1C6", linewidth=solid_linewidth)
    lines!(ax_pol, df_pol_pr[!, :Voltage], abs.(df_pol_pr[!, :Current] .* (cm^2/mA));
        color=priv_color, linewidth=dash_linewidth, linestyle=:dash)

    lines!(ax_cdl, df_cdl[!, :Voltage], df_cdl[!, :Capacitance] / (μF / cm^2);
        color="#8ED1C6", linewidth=solid_linewidth)
    lines!(ax_cdl, df_cdl_pr[!, :Voltage], df_cdl_pr[!, :Capacitance] / (μF / cm^2);
        color=priv_color, linewidth=dash_linewidth, linestyle=:dash)

    model_elems = [
        LineElement(color = priv_color, linewidth = 1.8, linestyle = :dash),
        LineElement(color = :black, linewidth = 3, linestyle = :solid),
    ]
    model_labels = [
        "MPB",
        "DGML",
    ]

    Legend(
        fig[1, 1],
        [model_elems],
        [model_labels],
        ["Model"];
        framevisible = false,
        nbanks = 1,
        tellwidth = false,
        tellheight = false,
        halign = :right,
        valign = :bottom,
        margin = (10, 10, 10, 10),
    )
	text!(fig[2, 1], -1.15, 0.3, text=L"\mathrm{K^+}",
				  color=colors[1], fontsize=24, font="sans-bold")
	text!(fig[2, 1], -0.90, 0.000000009, text=L"\mathrm{H^+}", 
				  color=colors[2], fontsize=24, font = "sans-bold") 
	text!(fig[2, 1], -1.05, 2e-11, text=L"\mathrm{CO_3^{2-}}",
				  color=colors[4], fontsize=24, font = "sans-bold") 
	text!(fig[2, 1], -1.055, 4.00e-6, text=L"\mathrm{HCO_3^-}", 
				  color=colors[3], fontsize=24, font = "sans-bold") 
	text!(fig[2, 1], -1.23, 7.0e-6, text=L"\mathrm{CO_2}", 
				  color=colors[5], fontsize=24, font = "sans-bold")
	text!(fig[2, 1], -0.90, 9e-10, text=L"\mathrm{OH^-}",
				  color=colors[6], fontsize=24, font = "sans-bold") 
	text!(fig[2, 1], -1.17, 0.00024, text=L"\mathrm{CO}", 
			 	  color=colors[7], fontsize=24, font = "sans-bold")
    resize_to_layout!(fig)

	Label(fig[1, 1, TopLeft()], "(a)",
    fontsize = 25,
    font = :bold,
    padding = (0, 45, 8, 0),
    halign = :right,
    valign = :bottom)

Label(fig[1, 2, TopLeft()], "(b)",
    fontsize = 25,
    font = :bold,
    padding = (0, 45, 8, 0),
    halign = :right,
    valign = :bottom)

Label(fig[2, 1, TopLeft()], "(c)",
    fontsize = 25,
    font = :bold,
    padding = (0, 45, 8, 0),
    halign = :right,
    valign = :bottom)

Label(fig[2, 2, TopLeft()], "(d)",
    fontsize = 25,
    font = :bold,
    padding = (0, 45, 8, 0),
    halign = :right,
    valign = :bottom)


	fig
	
    #display(fig)
end

# ╔═╡ 0290138e-4839-49a0-9a77-050bb4ab2cc0
let
    FS_BIG   = 24
    FS_SMALL = 24
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

    df_conc = 
		CSV.read(raw"../data/output/Concentration_Robin_Stefan_γ_Potassium_only.csv", DataFrame)
    df_pol  = 
		CSV.read(raw"../data/output/Polarization_Curve_Robin_Stefan_γ_pnp_Potassium_only.csv", DataFrame)
    df_cdl  = 
		CSV.read(raw"../data/output/DLCap_Robin_Stefan_γ_pnp_Potassium_only_1.csv", DataFrame)
    df_act  = 
		CSV.read(raw"../data/output/Activity_Curve_Robin_Stefan_γ_pnp_Potassium_only.csv", DataFrame)


	df_conc_pr = 
		CSV.read(raw"../data/output/Concentration_Robin_DMGL_γ_Potassium_only.csv", DataFrame)
    df_pol_pr  = 
		CSV.read(raw"../data/output/Polarization_Curve_Robin_DMGL_γ_pnp_Potassium_only.csv", DataFrame)
    df_cdl_pr  = 
		CSV.read(raw"../data/output/DLCap_Robin_DMGL_γ_pnp_Potassium_only_1.csv", DataFrame)
    df_act_pr  = 
		CSV.read(raw"../data/output/Activity_Curve_Robin_DMGL_γ_pnp_Potassium_only.csv", DataFrame)
	

      fig = Figure(size = (1200, 820), figure_padding = (40, 40, 30, 45))
    priv_color = (:gray60, 0.8)

    ax_cdl = Axis(fig[1, 1],
        xlabel = L"\quad φ\; \mathrm{(V vs SHE)}",
        ylabel = L"\text{C_{dl}} \textrm{(μF\,cm^{-2})}",
        limits = (-0.8, 0.8, 0, 25),
    )
    style_axis!(ax_cdl; big=false)

    ax_pol = Axis(fig[1, 2],
        xlabel = L"\quad φ\; \mathrm{(V vs SHE)}",
        ylabel = L"|j| \; \textrm{(mA\,cm^{-2})}",
        limits = (-1.3, -0.4, 1e-10, 100),
		#yaxisposition = :right,
        yscale = log10
    )
    style_axis!(ax_pol; big=false)

    ax_con = Axis(fig[2, 1],
        xlabel = L"\quad φ\; \mathrm{(V vs SHE)}",
        ylabel = L"\text{Interfacial\ concentration}\quad \mathbf{c_i^{\ddagger}}\;(\mathrm{M})",
        yscale = log10,
        limits = (-1.25, -0.5, 1e-11, 1e1),
    )
    style_axis!(ax_con; big=true)

    ax_act = Axis(fig[2, 2],
        xlabel = L"\quad φ\; \mathrm{(V vs SHE)}",
        ylabel = L"\text{Interfacial\ activity}\qquad \mathbf{a_i^{\ddagger}}\; ",
        yscale = log10,
        limits = (-1.25, -0.5, 1e-11, 1e4),
		#yaxisposition = :right
    )
    style_axis!(ax_act; big=true)

    rowsize!(fig.layout, 1, Relative(0.28))
    rowsize!(fig.layout, 2, Relative(0.72))
    rowgap!(fig.layout, 14)
    colgap!(fig.layout, 22)

    species = ["K⁺","H⁺","HCO₃⁻","CO₃²⁻","CO₂","OH⁻","CO"]
    colors  = [:orange, :gray, :brown, :violet, :red, :green, :blue]
	colors_pastel = [
	    "#F7C97F",  
	    "#D3D3D3",  
	    "#D9C2A7", 
	    "#FFC5D3",  
	    "#FF746C",  
	    "#80EF80",  
	    "#AFCBFF",  
	]
    vgrid_con = df_conc[!, :Voltage]
    conc_electrode = permutedims(Matrix(df_conc[!, Symbol.(species)]))

    vgrid_act = df_act[!, :Voltage]
    act_electrode = permutedims(Matrix(df_act[!, Symbol.(species)]))

    vgrid_con_pr = df_conc_pr[!, :Voltage]
    conc_electrode_pr = permutedims(Matrix(df_conc_pr[!, Symbol.(species)]))

    vgrid_act_pr = df_act_pr[!, :Voltage]
    act_electrode_pr = permutedims(Matrix(df_act_pr[!, Symbol.(species)]))

    xt_bottom = [-1.2, -1.0, -0.8, -0.6]
    ax_con.xticks = (xt_bottom, [@sprintf("%.1f", x) for x in xt_bottom])
    ax_act.xticks = (xt_bottom, [@sprintf("%.1f", x) for x in xt_bottom])
    ax_pol.xticks = (xt_bottom, [@sprintf("%.1f", x) for x in xt_bottom])

    yt_vals_con = 10.0 .^ (0:-3:-9)
    yt_labs_con = [L"10^{0}", L"10^{-3}", L"10^{-6}", L"10^{-9}"]
    ax_con.yticks = (yt_vals_con, yt_labs_con)

    yt_vals_act = 10.0 .^ (3:-3:-9)
    yt_labs_act = [L"10^{3}", L"10^{0}", L"10^{-3}", L"10^{-6}", L"10^{-9}"]
    ax_act.yticks = (yt_vals_act, yt_labs_act)

    yt_vals_pol = 10.0 .^ (0:-5:-10)
    yt_labs_pol = [L"10^{0}", L"10^{-5}", L"10^{-10}"]
    ax_pol.yticks = (yt_vals_pol, yt_labs_pol)

    yt_vals_cdl = [0, 10, 20]
    yt_labs_cdl = [L"0", L"10", L"20"]
    ax_cdl.yticks = (yt_vals_cdl, yt_labs_cdl)

    xt_top = [-0.8, -0.4, 0.0, 0.4, 0.8]
    ax_cdl.xticks = (xt_top, [@sprintf("%.1f", x) for x in xt_top])

	solid_linewidth = 3.5
	dash_linewidth = 2
	
    for ia in 1:7
        lines!(ax_con, vgrid_con, max.(conc_electrode[ia, :], eps(Float64));
            color=colors[ia], linewidth=solid_linewidth)
        lines!(ax_act, vgrid_act, max.(act_electrode[ia, :], eps(Float64));
            color=colors[ia], linewidth=solid_linewidth)
        lines!(ax_con, vgrid_con_pr, max.(conc_electrode_pr[ia, :], eps(Float64));
            color=colors_pastel[ia], linewidth=dash_linewidth, linestyle=:dash)
        lines!(ax_act, vgrid_act_pr, max.(act_electrode_pr[ia, :], eps(Float64));
            color=colors_pastel[ia], linewidth=dash_linewidth, linestyle=:dash)
    end

    lines!(ax_pol, df_pol[!, :Voltage], abs.(df_pol[!, :Current] .* (cm^2/mA));
        color="#8ED1C6", linewidth=solid_linewidth)
    lines!(ax_pol, df_pol_pr[!, :Voltage], abs.(df_pol_pr[!, :Current] .* (cm^2/mA));
        color=priv_color, linewidth=dash_linewidth, linestyle=:dash)

    lines!(ax_cdl, df_cdl[!, :Voltage], df_cdl[!, :Capacitance] / (μF / cm^2);
        color="#8ED1C6", linewidth=solid_linewidth)
    lines!(ax_cdl, df_cdl_pr[!, :Voltage], df_cdl_pr[!, :Capacitance] / (μF / cm^2);
        color=priv_color, linewidth=dash_linewidth, linestyle=:dash)

    model_elems = [
        LineElement(color = priv_color, linewidth = 1.8, linestyle = :dash),
        LineElement(color = :black, linewidth = 3, linestyle = :solid),
    ]
    model_labels = [
        "MPB",
        "DGML",
    ]

    Legend(
        fig[1, 1],
        [model_elems],
        [model_labels],
        ["Model"];
        framevisible = false,
        nbanks = 1,
        tellwidth = false,
        tellheight = false,
        halign = :right,
        valign = :bottom,
        margin = (10, 10, 10, 10),
    )
	text!(fig[2, 1], -1.15, 0.3, text=L"\mathrm{K^+}",
				  color=colors[1], fontsize=24, font="sans-bold")
	text!(fig[2, 1], -0.90, 0.000000009, text=L"\mathrm{H^+}", 
				  color=colors[2], fontsize=24, font = "sans-bold") 
	text!(fig[2, 1], -1.05, 2e-11, text=L"\mathrm{CO_3^{2-}}",
				  color=colors[4], fontsize=24, font = "sans-bold") 
	text!(fig[2, 1], -1.055, 4.00e-6, text=L"\mathrm{HCO_3^-}", 
				  color=colors[3], fontsize=24, font = "sans-bold") 
	text!(fig[2, 1], -1.23, 7.0e-6, text=L"\mathrm{CO_2}", 
				  color=colors[5], fontsize=24, font = "sans-bold")
	text!(fig[2, 1], -0.90, 9e-10, text=L"\mathrm{OH^-}",
				  color=colors[6], fontsize=24, font = "sans-bold") 
	text!(fig[2, 1], -1.17, 0.00024, text=L"\mathrm{CO}", 
			 	  color=colors[7], fontsize=24, font = "sans-bold")
    resize_to_layout!(fig)

	Label(fig[1, 1, TopLeft()], "(a)",
    fontsize = 25,
    font = :bold,
    padding = (0, 30, 8, 0),
    halign = :right,
    valign = :bottom)

Label(fig[1, 2, TopLeft()], "(b)",
    fontsize = 25,
    font = :bold,
    padding = (0, 30, 8, 0),
    halign = :right,
    valign = :bottom)

Label(fig[2, 1, TopLeft()], "(c)",
    fontsize = 25,
    font = :bold,
    padding = (0, 30, 8, 0),
    halign = :right,
    valign = :bottom)

Label(fig[2, 2, TopLeft()], "(d)",
    fontsize = 25,
    font = :bold,
    padding = (0, 30, 8, 0),
    halign = :right,
    valign = :bottom)


	fig
	
    #display(fig)
end

# ╔═╡ Cell order:
# ╠═1472eb23-8b3b-453b-9a91-a550a9988c54
# ╠═f81c8540-591a-4a56-843f-f51085195e2e
# ╠═1abfb78d-7291-4950-97df-aeb789378bfc
# ╠═9e7ebac4-b131-4362-9c2e-a07070715df6
# ╠═270509a2-433d-42af-886b-983f226f3229
# ╠═a9bb3083-d499-4462-845b-c41963cae1e5
# ╠═b3f44506-52eb-491c-bc44-76c9c49a43cd
# ╠═40b4e182-7aa2-4518-842a-e70dd9dced0d
# ╠═ccb13e5a-5ceb-4a56-9f3f-14307fc78115
# ╠═b0a4b942-5654-4b22-8849-90bf6c7f957f
# ╠═fc8096e4-01ca-451a-87cd-7e2e0171a531
# ╠═05eb8a6f-d7b5-4f37-a905-2209323afe2d
# ╠═5bb2da63-cc7d-4fcf-98b0-d2f4fa1dc1ef
# ╠═2ddeace6-663d-4fd1-8494-f1c88c19e628
# ╠═caa05490-0c7d-44ec-9be8-73f7a4473d8a
# ╠═2b10e19b-c099-4e4e-9cc1-5be1255d03bd
# ╟─4a68d318-ab6f-45b3-ad63-33d4a77c534d
# ╟─8545d818-d255-4e8a-af8f-72fdaf9d4bdd
# ╠═dcb38fd2-23b6-481e-83b7-4f3ee332c4ee
# ╠═ee091126-671e-46d8-8ed5-9457aeefd8fb
# ╟─4b02a8c0-a854-436a-9fdc-99ecf637e858
# ╟─59b13ba5-c85d-4e31-bd06-e1541b426205
# ╟─95844847-6f46-4840-979d-c54282295f41
# ╠═a09d4a6a-1017-4792-b684-7f9ddfeba83b
# ╠═5a090c2f-eb24-4453-b650-f4bc6c32c952
# ╠═224cf209-4aaa-49ca-a577-d77a597d256c
# ╟─400709f7-2a7c-431a-96ce-99a678b840f7
# ╟─d7ad452c-d3c1-485c-92ac-0b006dfb5acd
# ╟─89112dd8-08c8-4887-90db-e2e9c5ce0d88
# ╠═17b95990-f6e8-4da5-909a-1a2d46d7494d
# ╠═0290138e-4839-49a0-9a77-050bb4ab2cc0
