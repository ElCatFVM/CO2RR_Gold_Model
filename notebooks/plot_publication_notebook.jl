### A Pluto.jl notebook ###
# v1.0.3

using Markdown
using InteractiveUtils

# ╔═╡ 1472eb23-8b3b-453b-9a91-a550a9988c54
begin
	using Pkg
 	Pkg.activate(joinpath(@__DIR__, ".."))	
	using CSV, DataFrames, Colors, LessUnitful
	using CairoMakie
	using Printf
	using Makie: rich, subscript, superscript, italic

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

# ╔═╡ 7ba7df52-eee1-4f36-a5d0-379f31c6c5ef
function electrochemistry_theme()
    Theme(
        size = (960, 540),
        fonts = Attributes(
            regular = Makie.to_font("DejaVu Sans"),
            bold    = Makie.to_font("DejaVu Sans Bold"),
        ),
        Axis = (
            spinewidth        = 5.5,
            xtickwidth        = 2.0,
            ytickwidth        = 2.0,
            xticksize         = 8,
            yticksize         = 8,
            xlabelsize        = 25,
            ylabelsize        = 25,
            xticklabelsize    = 25,
            yticklabelsize    = 25,
            xgridvisible      = false,
            ygridvisible      = false,
            xlabelpadding     = 10,
            ylabelpadding     = 10,
            xlabelfont        = :bold,
            ylabelfont        = :bold,
            yticks            = LinearTicks(5),
            xticks            = LinearTicks(4),
        ),
        Lines = (linewidth = 3,),
    )
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

    # extract only the pair indices for the desired pressures → [2, 5, 7]
    keep = findall(in(wanted), pressures[1:npairs])

    # only 3, so use distinct (discrete) colors
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
    # read extracted CSV
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
    # read experimental CSV
    # --------------------------------------------------
    raw = CSV.read("../data/Langmuir_CV_data/Figure_3.csv", DataFrame; header=false)

    sub = Matrix(raw[4:end, :])
    num = map(x -> x === missing ? NaN : parse(Float64, x), sub)
    num_df = DataFrame(num, :auto)

    # exclude Ar sat
    pressures_exp = ["0.1", "0.2", "0.3", "0.5", "0.6", "1.0"]
    exp_indices = 2:7   # exclude the 1st pair (Ar sat)

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
    # read simulation CSV
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

    # highlight region (keep the earlier style)
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
    const FS_BIG   = 22
    const FS_SMALL = 18
    const SP_BIG   = 3.5
    const SP_SMALL = 3.5
    const TICKW_BIG   = 2.0
    const TICKW_SMALL = 2.0
    const TICKL_BIG   = 8
    const TICKL_SMALL = 8

    function style_axis!(ax; big::Bool)
        fs = big ? FS_BIG   : FS_SMALL
        sp = big ? SP_BIG   : SP_SMALL
        tw = big ? TICKW_BIG : TICKW_SMALL
        tl = big ? TICKL_BIG : TICKL_SMALL

        ax.spinewidth     = sp
        ax.xtickwidth     = tw
        ax.ytickwidth     = tw
        ax.xticksize      = tl
        ax.yticksize      = tl
        ax.xlabelsize     = fs
        ax.ylabelsize     = fs
        ax.xticklabelsize = fs
        ax.yticklabelsize = fs
        ax.xlabelpadding  = 10
        ax.ylabelpadding  = 10
        if big
            ax.xlabelfont = :bold
            ax.ylabelfont = :bold
        end
        ax.xgridvisible = false
        ax.ygridvisible = false
        return ax
    end
end

# ╔═╡ 4d250254-b8a8-4ff2-95c6-81a02ddfd982
let
    FS_BIG   = 28
    FS_SMALL = 24
    SP_BIG   = 3.5
    SP_SMALL = 3.5
    TICKW_BIG   = 2.0
    TICKW_SMALL = 2.0
    TICKL_BIG   = 8
    TICKL_SMALL = 8

    solid_lw = 7.0
    dash_lw  = 4.5

    species   = ["K⁺","H⁺","HCO₃⁻","CO₃²⁻","CO₂","OH⁻","CO"]
    sp_colors = ["#E07B39","#888888","#7B5C3E","#222222",
                 "#C0392B","#27AE60","#2980B9"]
    sp_pastel = ["#F7C97F","#D3D3D3","#D9C2A7","#666666",
                 "#FF746C","#80EF80","#AFCBFF"]
    sp_rich = [
        rich("K", superscript("+")),
        rich("H", superscript("+")),
        rich("HCO", subscript("3"), superscript("−")),
        rich("CO", subscript("3"), superscript("2−")),
        rich("CO", subscript("2")),
        rich("OH", superscript("−")),
        rich("CO"),
    ]

    powlab(n) = rich("10", superscript(string(n)))

    fig = Figure(size = (1250, 750), figure_padding = (25, 30, 30, 25))

    function format_axis!(ax; big=true)
        fs = big ? FS_BIG : FS_SMALL
        sp = big ? SP_BIG : SP_SMALL
        tw = big ? TICKW_BIG : TICKW_SMALL
        tl = big ? TICKL_BIG : TICKL_SMALL
        ax.spinewidth     = sp
        ax.xtickwidth     = tw
        ax.ytickwidth     = tw
        ax.xticksize      = tl
        ax.yticksize      = tl
        ax.xlabelsize     = fs
        ax.ylabelsize     = fs
        ax.xticklabelsize = fs
        ax.yticklabelsize = fs
        ax.xlabelpadding  = 10
        ax.ylabelpadding  = 10
        ax.xgridvisible   = false
        ax.ygridvisible   = false
        return ax
    end

    # labels: physical quantity italic / abbreviations roman
    lab_phi   = rich(rich("U", font=:italic), "  (V vs. SHE)")
    lab_con   = rich("Interfacial Concentration\n\n\n",
                     rich("c", font=:italic),
                     subscript(rich("i", font=:italic)), superscript("‡"), "  (M)")
    lab_act   = rich("Interfacial Activity\n\n\n ",
                     rich("a", font=:italic),
                     subscript(rich("i", font=:italic)), superscript("‡"))
    # ratio: log(CatINT c / MPNP c)
    lab_ratio_c = rich("log(", superscript("DGML"), rich("c", font=:italic),
                       subscript(rich("i", font=:italic)), superscript("‡"),
                       " / ", superscript("MPNP"), rich("c", font=:italic),
                       subscript(rich("i", font=:italic)), superscript("‡"), ")")
    lab_ratio_a = rich("log(", superscript("DGML"), rich("a", font=:italic),
                       subscript(rich("i", font=:italic)), superscript("‡"),
                       " / ", superscript("MPNP"), rich("a", font=:italic),
                       subscript(rich("i", font=:italic)), superscript("‡"), ")")

    xt = [-1.5, -1.2, -0.9, -0.6]
    xtlab = [@sprintf("%.1f", x) for x in xt]

    # ── (a) Concentration ──────────────────────────────────────
    df_prev_a = CSV.read(raw"../data/output/Concentration_Robin_Stefan_γ_Same_Size.csv", DataFrame)  # MPNP (solid)
    df_curr_a = CSV.read(raw"../data/output/Concentration_Robin_DMGL_γ_Same_Size.csv", DataFrame)     # CatINT (dash)
    V_a = df_curr_a[!, "Voltage"]

    floors_a = Dict{String, Float64}()
    for s in species
        vals = vcat(df_prev_a[!, s], df_curr_a[!, s])
        pos  = vals[vals .> 0]
        floors_a[s] = isempty(pos) ? 1e-30 : minimum(pos) * 0.1
    end

    ax_main_a = Axis(fig[1, 1], ylabel = lab_con, yscale = log10)
    ax_main_a.xticklabelsvisible = false
    ax_aux_a  = Axis(fig[2, 1], xlabel = lab_phi, ylabel = lab_ratio_c)
    format_axis!(ax_main_a, big=true)
    format_axis!(ax_aux_a,  big=true)

    for (i, s) in enumerate(species)
        mpnp   = df_prev_a[!, s]   # solid
        catint = df_curr_a[!, s]   # dash
        f      = floors_a[s]
        lines!(ax_main_a, V_a, mpnp;   color = sp_colors[i], linestyle = :solid, linewidth = solid_lw)
        lines!(ax_main_a, V_a, catint; color = sp_pastel[i], linestyle = :dash,  linewidth = dash_lw)
        delta = log10.((catint .+ f) ./ (mpnp .+ f))   # log(CatINT/MPNP)
        lines!(ax_aux_a, V_a, delta; color = sp_colors[i], linewidth = solid_lw - 2)
    end
    hlines!(ax_aux_a, [0.0]; color = :black, linestyle = :dot, linewidth = 2)

    linkxaxes!(ax_main_a, ax_aux_a)
    ylims!(ax_aux_a, -0.075, 0.075)
    ylims!(ax_main_a, 1e-12, 1e2)
    xlims!(ax_aux_a, -1.5, -0.6)

    sp_pos_a = [
        (-1.35, 3e-1),  (-0.90, 1e-8),  (-0.85, 1e-4),  (-1.1, 5e-12),
        (-1.35,  2e-6),  (-0.90, 1e-9),  (-1.37, 10e-5),
    ]
    order_a = [1,2,3,4,5,6,7]
    for i in order_a
        x,y = sp_pos_a[i]
        text!(ax_main_a, x, y; text=sp_rich[i], color=sp_colors[i], fontsize=24, font=:bold)
    end
    ax_main_a.yticks = ([1e-12, 1e-8, 1e-4, 1], [powlab(-12), powlab(-8), powlab(-4), powlab(0)])
    ax_main_a.xticks = (xt, xtlab)
    ax_aux_a.xticks  = (xt, xtlab)

    # ── (b) Activity ───────────────────────────────────────────
    df_prev_b = CSV.read(raw"../data/output/Activity_Curve_Robin_Stefan_γ_pnp_Same_Size.csv", DataFrame)  # MPNP
    df_curr_b = CSV.read(raw"../data/output/Activity_Curve_Robin_DMGL_γ_pnp_Same_Size.csv", DataFrame)     # CatINT
    V_b = df_curr_b[!, "Voltage"]

    floors_b = Dict{String, Float64}()
    for s in species
        vals = vcat(df_prev_b[!, s], df_curr_b[!, s])
        pos  = vals[vals .> 0]
        floors_b[s] = isempty(pos) ? 1e-30 : minimum(pos) * 0.1
    end

    ax_main_b = Axis(fig[1, 2], ylabel = lab_act, yscale = log10)
    ax_main_b.xticklabelsvisible = false
    ax_aux_b  = Axis(fig[2, 2], xlabel = lab_phi, ylabel = lab_ratio_a)
    format_axis!(ax_main_b, big=true)
    format_axis!(ax_aux_b,  big=true)

    for (i, s) in enumerate(species)
        mpnp   = df_prev_b[!, s]
        catint = df_curr_b[!, s]
        f      = floors_b[s]
        lines!(ax_main_b, V_b, mpnp;   color = sp_colors[i], linestyle = :solid, linewidth = solid_lw)
        lines!(ax_main_b, V_b, catint; color = sp_pastel[i], linestyle = :dash,  linewidth = dash_lw)
        delta = log10.((catint .+ f) ./ (mpnp .+ f))
        lines!(ax_aux_b, V_b, delta; color = sp_colors[i], linewidth = solid_lw - 2)
    end
    hlines!(ax_aux_b, [0.0]; color = :black, linestyle = :dot, linewidth = 2)

    linkxaxes!(ax_main_b, ax_aux_b)
    ylims!(ax_aux_b, -0.075, 0.075)
    ylims!(ax_main_b, 1e-10, 1e4)
    xlims!(ax_aux_b, -1.5, -0.6)

    sp_pos_b = [
        (-1.35, 50),    (-0.875, 4e-7),  (-0.85, 70e-5),  (-1.15, 5e-10),
        (-1.33, 5e-4), (-0.875, 1e-8),  (-1.4, 3.0e-2),
    ]
    for i in 1:7
        x,y = sp_pos_b[i]
        text!(ax_main_b, x, y; text=sp_rich[i], color=sp_colors[i], fontsize=24, font=:bold)
    end
    ax_main_b.yticks = ([1e-10, 1e-6, 1e-2, 1e2], [powlab(-10), powlab(-6), powlab(-2), powlab(2)])
    ax_main_b.xticks = (xt, xtlab)
    ax_aux_b.xticks  = (xt, xtlab)

    Label(fig[1, 1, TopLeft()], "(a)", fontsize = 24, font = :bold, padding = (0, 5, 5, 0))
    Label(fig[1, 2, TopLeft()], "(b)", fontsize = 24, font = :bold, padding = (0, 5, 5, 0))

    rowsize!(fig.layout, 1, Relative(0.72))
    rowsize!(fig.layout, 2, Relative(0.28))
    rowgap!(fig.layout, 14)
    colgap!(fig.layout, 40)

    resize_to_layout!(fig)
    fig
end

# ╔═╡ ac4deeb9-812a-429c-afe9-389b22814ee6
let
    FS_BIG   = 32
    FS_SMALL = 32
    SP_BIG   = 3.5
    SP_SMALL = 3.5
    TICKW_BIG   = 2.0
    TICKW_SMALL = 2.0
    TICKL_BIG   = 8
    TICKL_SMALL = 8

    solid_lw = 7.0
    dash_lw  = 4.5


    df_prev = CSV.read(raw"../data/catmap_CO2R_data/voltage-conc.csv", DataFrame)              # CatINT (dash, long-format Index)
    df_curr = CSV.read(raw"../data/output/Concentration_Robin_DMGL_γ_All_species.csv", DataFrame)  # MPNP (solid)

    species   = ["K⁺","H⁺","HCO₃⁻","CO₃²⁻","CO₂","OH⁻","CO"]
    sp_colors = ["#E07B39","#888888","#7B5C3E","#222222",
                 "#C0392B","#27AE60","#2980B9"]
    sp_pastel = ["#F7C97F","#D3D3D3","#D9C2A7","#666666",
                 "#FF746C","#80EF80","#AFCBFF"]
    sp_rich = [
        rich("K", superscript("+")),
        rich("H", superscript("+")),
        rich("HCO", subscript("3"), superscript("−")),
        rich("CO", subscript("3"), superscript("2−")),
        rich("CO", subscript("2")),
        rich("OH", superscript("−")),
        rich("CO"),
    ]

    powlab(n) = rich("10", superscript(string(n)))

    V = df_curr[!, "Voltage"]

    fig = Figure(size = (800, 700), figure_padding = (25, 25, 30, 25))

    function format_axis!(ax; big=true)
        fs = big ? FS_BIG : FS_SMALL
        sp = big ? SP_BIG : SP_SMALL
        tw = big ? TICKW_BIG : TICKW_SMALL
        tl = big ? TICKL_BIG : TICKL_SMALL
        ax.spinewidth     = sp
        ax.xtickwidth     = tw
        ax.ytickwidth     = tw
        ax.xticksize      = tl
        ax.yticksize      = tl
        ax.xlabelsize     = fs
        ax.ylabelsize     = fs
        ax.xticklabelsize = fs
        ax.yticklabelsize = fs
        ax.xlabelpadding  = 10
        ax.ylabelpadding  = 10
        ax.xgridvisible   = false
        ax.ygridvisible   = false
        return ax
    end

    # labels: physical quantity italic / abbreviations·units roman  (concentration)
    lab_con = rich("Interfacial Concentration\n\n",
                   rich("c", font=:italic),
                   subscript(rich("i", font=:italic)), superscript("‡"), "  (M)")
    lab_phi = rich(rich("U", font=:italic), "  (V vs. SHE)")
    lab_ratio = rich("log(", superscript("CatINT"), rich("c", font=:italic),
                     subscript(rich("i", font=:italic)), superscript("‡"),
                     " / ", superscript("MPNP"), rich("c", font=:italic),
                     subscript(rich("i", font=:italic)), superscript("‡"), ")")

    ax_main = Axis(fig[1, 1];
        ylabel = lab_con,
        yscale = log10,
        limits = (-1.20, -0.55, 1e-11, 1e3),
    )
    ax_main.xticklabelsvisible = false

    # ── MPNP (solid) + CatINT (dash) ───────────────────────────
    for (i, s) in enumerate(species)
        curr = df_curr[!, s]
        lines!(ax_main, V, curr; color = sp_colors[i], linewidth = solid_lw)
        row = df_prev[df_prev.Index .== i, :]
        if nrow(row) > 0
            lines!(ax_main, row[!, "Voltage"], row[!, "Concentration"];
                color = sp_pastel[i], linewidth = dash_lw, linestyle = :dash)
        end
    end

    # species labels (rich)
    sp_pos = [
        (-1.15, 5),     (-0.875, 1e-8),  (-1.05, 30e-7),  (-1.05, 1e-11),
        (-1.2, 72e-7),    (-0.875, 7e-10), (-1.17, 3.0e-4),
    ]
    for i in 1:7
        x,y = sp_pos[i]
        text!(ax_main, x, y; text=sp_rich[i], color=sp_colors[i], fontsize=24, font=:bold)
    end

    # ── ratio aux: log(CatINT / MPNP) ──────────────────────────
    ax_aux = Axis(fig[2, 1];
        xlabel = lab_phi,
        ylabel = lab_ratio,
    )
    for (i, s) in enumerate(species)
        row_prev = df_prev[df_prev.Index .== i, :]
        if nrow(row_prev) > 0
            matched_curr = [df_curr[findmin(abs.(V .- v))[2], s] for v in row_prev[!, "Voltage"]]
            delta_y = log10.(row_prev[!, "Concentration"] ./ matched_curr)  # log(CatINT/MPNP)
            scatter!(ax_aux, row_prev[!, "Voltage"], delta_y;
                color = sp_colors[i], markersize = 10)
        end
    end
    hlines!(ax_aux, [0.0]; color = :black, linestyle = :dot, linewidth = 2)

    format_axis!(ax_main, big=true)
    format_axis!(ax_aux, big=true)

    linkxaxes!(ax_main, ax_aux)
    xlims!(ax_aux, -1.250, -0.55)
    ylims!(ax_aux, -0.3, 0.3)

    xt = [-1.4, -1.2, -1.0, -0.8, -0.6]
    xtlab = [@sprintf("%.1f", x) for x in xt]
    ax_main.yticks = ([1e-10, 1e-6, 1e-2, 1e2], [powlab(-10), powlab(-6), powlab(-2), powlab(2)])
    ax_main.xticks = (xt, xtlab)
    ax_aux.xticks  = (xt, xtlab)

    # ── legend: MPNP(solid) / CatINT(dash) ─────────────────────
    legend_entries = [
        LineElement(color = :black,  linestyle = :solid, linewidth = solid_lw),
        LineElement(color = :gray60, linestyle = :dash,  linewidth = dash_lw),
    ]
    legend_labels = ["MPNP@LiquidElectrolytes.jl", "CatINT"]

    axislegend(ax_main, legend_entries, legend_labels;
        position = :rt,
        orientation = :vertical,
        labelsize = 20,
        font = :bold,
        framevisible = true,
        backgroundcolor = (:white, 0.5),
        framecolor = (:black, 0.5),
    )

    rowsize!(fig.layout, 1, Relative(0.68))
    rowsize!(fig.layout, 2, Relative(0.32))
    rowgap!(fig.layout, 15)

    resize_to_layout!(fig)
    fig
end

# ╔═╡ 9d660147-670a-4728-b2ef-8bb978f6981b
let
    FS_BIG   = 32
    FS_SMALL = 32
    SP_BIG   = 3.5
    SP_SMALL = 3.5
    TICKW_BIG   = 2.0
    TICKW_SMALL = 2.0
    TICKL_BIG   = 8
    TICKL_SMALL = 8

    solid_lw = 7.0
    dash_lw  = 4.5

    df_sim = CSV.read(raw"../data/E.Acta_Cap_data/Landstorfer_NaClO4_0.005M.csv", DataFrame)  # dash
    df_exp = CSV.read(raw"../data/Valette_Cap_data/NaClO4_0.005M.csv", DataFrame)             # solid

    # row 1 = Voltage, row 2 = dlcaps (row-major)
    M_sim = Matrix(df_sim)
    M_exp = Matrix(df_exp)
    V_sim    = Float64.(M_sim[:, 1]);  caps_sim = Float64.(M_sim[:, 2])
    V_exp    = Float64.(M_exp[:, 1]);  caps_exp = Float64.(M_exp[:, 2])

    fig = Figure(size = (800, 700), figure_padding = (25, 25, 30, 25))

    function format_axis!(ax; big=true)
        fs = big ? FS_BIG : FS_SMALL
        sp = big ? SP_BIG : SP_SMALL
        tw = big ? TICKW_BIG : TICKW_SMALL
        tl = big ? TICKL_BIG : TICKL_SMALL
        ax.spinewidth     = sp
        ax.xtickwidth     = tw
        ax.ytickwidth     = tw
        ax.xticksize      = tl
        ax.yticksize      = tl
        ax.xlabelsize     = fs
        ax.ylabelsize     = fs
        ax.xticklabelsize = fs
        ax.yticklabelsize = fs
        ax.xlabelpadding  = 10
        ax.ylabelpadding  = 10
        ax.xgridvisible   = false
        ax.ygridvisible   = false
        return ax
    end

    # labels: physical quantity italic / abbreviations·units roman
    lab_phi = rich(rich("U", font=:italic), "  (V vs. SHE)")
    lab_cap = rich("Differential Capacitance\n\n",
                   rich("C", font=:italic),
                   subscript(rich("dl", font=:italic)),
                   "  (μF cm", superscript("−2"), ")")

    ax = Axis(fig[1, 1];
        xlabel = lab_phi,
        ylabel = lab_cap,
    )

    lines!(ax, V_exp, caps_exp; color = "#222222", linewidth = solid_lw)
    lines!(ax, V_sim, caps_sim; color = "#888888", linewidth = dash_lw, linestyle = :dash)

    format_axis!(ax, big=true)

    # ── legend: Experiment(solid) / Theory(dash) ───────────────
    legend_entries = [
        LineElement(color = :black,  linestyle = :solid, linewidth = solid_lw),
        LineElement(color = :gray60, linestyle = :dash,  linewidth = dash_lw),
    ]
    legend_labels = ["Experiment", "Theory"]

    axislegend(ax, legend_entries, legend_labels;
        position = :rt,
        orientation = :vertical,
        labelsize = 32,
        font = :bold,
        framevisible = true,
        backgroundcolor = (:white, 0.5),
        framecolor = (:black, 0.5),
    )

    resize_to_layout!(fig)
    fig
end

# ╔═╡ 400709f7-2a7c-431a-96ce-99a678b840f7
let
    FS_BIG   = 32
    FS_SMALL = 32
    SP_BIG   = 3.5
    SP_SMALL = 3.5
    TICKW_BIG   = 2.0
    TICKW_SMALL = 2.0
    TICKL_BIG   = 8
    TICKL_SMALL = 8

    solid_lw = 7.0
    dash_lw  = 4.5

    df_prev = CSV.read(raw"../data/catmap_CO2R_data/voltage-activ.csv", DataFrame)              # CatINT (dash, long-format Index)
    df_curr = CSV.read(raw"../data/output/Activity_Curve_Robin_Stefan_γ_pnp_All_species.csv", DataFrame)  # MPNP (solid)

    species   = ["K⁺","H⁺","HCO₃⁻","CO₃²⁻","CO₂","OH⁻","CO"]
    sp_colors = ["#E07B39","#888888","#7B5C3E","#222222",
                 "#C0392B","#27AE60","#2980B9"]
    sp_pastel = ["#F7C97F","#D3D3D3","#D9C2A7","#666666",
                 "#FF746C","#80EF80","#AFCBFF"]
    sp_rich = [
        rich("K", superscript("+")),
        rich("H", superscript("+")),
        rich("HCO", subscript("3"), superscript("−")),
        rich("CO", subscript("3"), superscript("2−")),
        rich("CO", subscript("2")),
        rich("OH", superscript("−")),
        rich("CO"),
    ]

    powlab(n) = rich("10", superscript(string(n)))

    V = df_curr[!, "Voltage"]

    fig = Figure(size = (800, 700), figure_padding = (25, 25, 30, 25))

    function format_axis!(ax; big=true)
        fs = big ? FS_BIG : FS_SMALL
        sp = big ? SP_BIG : SP_SMALL
        tw = big ? TICKW_BIG : TICKW_SMALL
        tl = big ? TICKL_BIG : TICKL_SMALL
        ax.spinewidth     = sp
        ax.xtickwidth     = tw
        ax.ytickwidth     = tw
        ax.xticksize      = tl
        ax.yticksize      = tl
        ax.xlabelsize     = fs
        ax.ylabelsize     = fs
        ax.xticklabelsize = fs
        ax.yticklabelsize = fs
        ax.xlabelpadding  = 10
        ax.ylabelpadding  = 10
        ax.xgridvisible   = false
        ax.ygridvisible   = false
        return ax
    end

    # labels: physical quantity italic / abbreviations·units roman
    lab_act = rich("Interfacial Activity\n\n",
                   rich("a", font=:italic),
                   subscript(rich("i", font=:italic)), superscript("‡"))
    lab_phi = rich(rich("U", font=:italic), "  (V vs. SHE)")
    lab_ratio = rich("log(", superscript("CatINT"), rich("a", font=:italic),
                     subscript(rich("i", font=:italic)), superscript("‡"),
                     " / ", superscript("MPNP"), rich("a", font=:italic),
                     subscript(rich("i", font=:italic)), superscript("‡"), ")")

    ax_main = Axis(fig[1, 1];
        ylabel = lab_act,
        yscale = log10,
        limits = (-1.25, -0.55, 1e-11, 1e3),
    )
    ax_main.xticklabelsvisible = false

    # ── MPNP (solid) + CatINT (dash) ───────────────────────────
    for (i, s) in enumerate(species)
        curr = df_curr[!, s]
        lines!(ax_main, V, curr; color = sp_colors[i], linewidth = solid_lw)
        row = df_prev[df_prev.Index .== i, :]
        if nrow(row) > 0
            lines!(ax_main, row[!, "Voltage"], row[!, "Concentration"];
                color = sp_pastel[i], linewidth = dash_lw, linestyle = :dash)
        end
    end

    # species labels (rich)
    sp_pos = [
        (-1.15, 30),     (-0.875, 2e-7),  (-1.05, 20e-7),  (-1.1, 0.5e-9),
        (-1.2, 9e-4),    (-0.875, 16e-9), (-1.17, 5000e-4),
    ]
    for i in 1:7
        x,y = sp_pos[i]
        text!(ax_main, x, y; text=sp_rich[i], color=sp_colors[i], fontsize=24, font=:bold)
    end

    # ── ratio aux: log(CatINT / MPNP) ──────────────────────────
    ax_aux = Axis(fig[2, 1];
        xlabel = lab_phi,
        ylabel = lab_ratio,
    )
    for (i, s) in enumerate(species)
        row_prev = df_prev[df_prev.Index .== i, :]
        if nrow(row_prev) > 0
            matched_curr = [df_curr[findmin(abs.(V .- v))[2], s] for v in row_prev[!, "Voltage"]]
            delta_y = log10.(row_prev[!, "Concentration"] ./ matched_curr)  # log(CatINT/MPNP)
            scatter!(ax_aux, row_prev[!, "Voltage"], delta_y;
                color = sp_colors[i], markersize = 10)
        end
    end
    hlines!(ax_aux, [0.0]; color = :black, linestyle = :dot, linewidth = 2)

    format_axis!(ax_main, big=true)
    format_axis!(ax_aux, big=true)

    linkxaxes!(ax_main, ax_aux)
    xlims!(ax_aux, -1.25, -0.55)
    ylims!(ax_aux, -0.3, 0.3)

    xt = [-1.2, -1.0, -0.8, -0.6]
    xtlab = [@sprintf("%.1f", x) for x in xt]
    ax_main.yticks = ([1e-10, 1e-6, 1e-2, 1e2], [powlab(-10), powlab(-6), powlab(-2), powlab(2)])
    ax_main.xticks = (xt, xtlab)
    ax_aux.xticks  = (xt, xtlab)

    # ── legend: MPNP(solid) / CatINT(dash) ─────────────────────
    legend_entries = [
        LineElement(color = :black,  linestyle = :solid, linewidth = solid_lw),
        LineElement(color = :gray60, linestyle = :dash,  linewidth = dash_lw),
    ]
    legend_labels = ["MPNP@LiquidElectrolytes.jl", "CatINT"]

    axislegend(ax_main, legend_entries, legend_labels;
        position = :rt,
        orientation = :vertical,
        labelsize = 20,
        font = :bold,
        framevisible = true,
        backgroundcolor = (:white, 0.5),
        framecolor = (:black, 0.5),
    )

    rowsize!(fig.layout, 1, Relative(0.68))
    rowsize!(fig.layout, 2, Relative(0.32))
    rowgap!(fig.layout, 15)

    resize_to_layout!(fig)
    fig
end

# ╔═╡ d7ad452c-d3c1-485c-92ac-0b006dfb5acd
let


    csv_paths = [
        raw"../data/output/IV_Robin_DMGL_γ_pnp_ohminus.csv",
        raw"../data/output/IV_Robin_DMGL_γ_pnp_ohminus_koh_nonsol.csv",
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
    FS_BIG   = 32
    FS_SMALL = 32
    SP_BIG   = 3.5
    SP_SMALL = 3.5
    TICKW_BIG   = 2.0
    TICKW_SMALL = 2.0
    TICKL_BIG   = 8
    TICKL_SMALL = 8

    function style_axis!(ax; big::Bool)
        fs = big ? FS_BIG : FS_SMALL
        sp = big ? SP_BIG : SP_SMALL
        tw = big ? TICKW_BIG : TICKW_SMALL
        tl = big ? TICKL_BIG : TICKL_SMALL
        ax.spinewidth     = sp
        ax.xtickwidth     = tw
        ax.ytickwidth     = tw
        ax.xticksize      = tl
        ax.yticksize      = tl
        ax.xlabelsize     = fs
        ax.ylabelsize     = fs
        ax.xticklabelsize = fs
        ax.yticklabelsize = fs
        ax.xlabelpadding  = 10
        ax.ylabelpadding  = 10
        ax.xgridvisible   = false
        ax.ygridvisible   = false
        return ax
    end

    col_mpnp   = "#2980B9"          # MPNP (solid, blue)
    col_catint = "#999999"          # CatINT (dash, gray) ── (a)(b)(c)
    col_dgml   = ("#999999", 0.9)   # (d) gray dotted
    col_exp    = "#222222"          # Experiment (dot)
    solid_lw = 7.0
    dash_lw  = 4.5

    # ── rich-text tick-label helper (avoid serif: rich instead of L"") ───────
    powlab(n) = rich("10", superscript(string(n)))

    # ── rich-text axis labels: quantity=italic, abbreviations/units=roman ──────
    lab_pol  = rich("Partial CO Current\n\n",
                    "|", rich("I", font=:italic), subscript("CO"),
                    "|  (mA cm", superscript("−2"), ")")
    lab_act  = rich("Interfacial Activity\n\n",
                    rich("a", font=:italic),
                    subscript(rich("i", font=:italic)), superscript("‡"))
    lab_con  = rich("Interfacial Concentration\n\n",
                    rich("c", font=:italic),
                    subscript(rich("i", font=:italic)), superscript("‡"),
                    "  (M)")
    lab_cdl  = rich("Differential Capacitance\n\n",
                    rich("C", font=:italic), subscript("dl"),
                    "  (μF cm", superscript("−2"), ")")
    lab_xshe = rich("Voltage ", rich("U", font=:italic), "  (V vs. SHE)")
    lab_xpzc = rich("Voltage [", rich("U", font=:italic), " − ",
                    rich("U", font=:italic), subscript("pzc"), "]  (V)")

    df_conc    = CSV.read(raw"../data/output/Concentration_Robin_Stefan_γ_Same_Size.csv", DataFrame)
    df_pol     = CSV.read(raw"../data/output/Polarization_Curve_Robin_Stefan_γ_pnp_Same_Size.csv", DataFrame)
    df_cdl     = CSV.read(raw"../data/output/DLCap_Robin_Stefan_γ_pb_Same_Size_1.csv", DataFrame)
    df_act     = CSV.read(raw"../data/output/Activity_Curve_Robin_Stefan_γ_pnp_Same_Size.csv", DataFrame)


	
    df_cap_exp_dgml = CSV.read(raw"../data/output/DLCap_Robin_DMGL_γ_pnp_All_species_1.csv", DataFrame)


	
    df_conc_pr = CSV.read(raw"../data/catmap_CO2R_data/voltage-conc.csv",    DataFrame)
    df_act_pr  = CSV.read(raw"../data/catmap_CO2R_data/voltage-activ.csv",   DataFrame)
    df_pol_pr  = CSV.read(raw"../data/catmap_CO2R_data/Ringe-theorical.csv", DataFrame; header=[:Voltage, :Current])
    # (a) experimental data (dot) — TODO: replace with the real path
    df_pol_exp = CSV.read(raw"../data/catmap_CO2R_data/Ringe-experimental.csv", DataFrame; header=[:Voltage, :Current])

    species   = ["K⁺","H⁺","HCO₃⁻","CO₃²⁻","CO₂","OH⁻","CO"]
    sp_colors = ["#E07B39","#888888","#7B5C3E","#222222",
                 "#C0392B","#27AE60","#2980B9"]
    sp_pastel = ["#F7C97F","#D3D3D3","#D9C2A7","#666666",
                 "#FF746C","#80EF80","#AFCBFF"]
    # rich species labels (defined once, used in both (b) and (c))
    sp_rich = [
        rich("K", superscript("+")),
        rich("H", superscript("+")),
        rich("HCO", subscript("3"), superscript("−")),
        rich("CO", subscript("3"), superscript("2−")),
        rich("CO", subscript("2")),
        rich("OH", superscript("−")),
        rich("CO"),
    ]

    vgrid_con = df_conc[!, :Voltage]
    vgrid_act = df_act[!, :Voltage]
    conc_electrode = permutedims(Matrix(df_conc[!, Symbol.(species)]))
    act_electrode  = permutedims(Matrix(df_act[!,  Symbol.(species)]))

    xt_bottom = [-1.4, -1.2, -1.0, -0.8, -0.6]
    xt_cdl    = [-0.8, -0.4,  0.0,  0.4,  0.8]

    fig = Figure(size = (1400, 1100), figure_padding = (40, 60, 30, 60))

    # (a) Partial CO Current — no x-label
    ax_pol = Axis(fig[1, 1];
        ylabel = lab_pol,
        limits = (-1.5, -0.4, 1e-10, 100),
        yscale = log10,
        yaxisposition = :left,
    )
    style_axis!(ax_pol; big=false)
    ax_pol.xticks = (xt_bottom, [@sprintf("%.1f", x) for x in xt_bottom])
    ax_pol.yticks = (10.0 .^ (0:-5:-10), [powlab(0), powlab(-5), powlab(-10)])
    ax_pol.xticklabelsvisible = true
    ax_pol.xlabelvisible      = false

    # (b) Interfacial Activity — right y-axis
    ax_act = Axis(fig[1, 2];
        ylabel = lab_act,
        yscale = log10,
        limits = (-1.5, -0.5, 1e-11, 1e4),
        yaxisposition = :right,
    )
    style_axis!(ax_act; big=false)
    ax_act.xticks = (xt_bottom, [@sprintf("%.1f", x) for x in xt_bottom])
    ax_act.yticks = (10.0 .^ (3:-3:-9),
                     [powlab(3), powlab(0), powlab(-3), powlab(-6), powlab(-9)])

    # (c) Interfacial Concentration
    ax_con = Axis(fig[2, 1];
        xlabel = lab_xshe,
        ylabel = lab_con,
        yscale = log10,
        limits = (-1.5, -0.5, 1e-11, 1e1),
        yaxisposition = :left,
    )
    style_axis!(ax_con; big=true)
    ax_con.xticks = (xt_bottom, [@sprintf("%.1f", x) for x in xt_bottom])
    ax_con.yticks = (10.0 .^ (0:-3:-9),
                     [powlab(0), powlab(-3), powlab(-6), powlab(-9)])
    ax_con.xticklabelsvisible = true
    ax_con.xlabelvisible      = true

    # (d) Differential Capacitance — right y-axis
    ax_cdl = Axis(fig[2, 2];
        xlabel = lab_xpzc,
        ylabel = lab_cdl,
        limits = (-0.9, 0.9, 0, 25),
        yaxisposition = :right,
    )
    style_axis!(ax_cdl; big=true)
    ax_cdl.xticks = (xt_cdl, [@sprintf("%.1f", x) for x in xt_cdl])
    ax_cdl.yticks = ([0, 10, 20, 30, 40, 50], ["0","10","20","30","40","50"])

    rowgap!(fig.layout, 20)
    colgap!(fig.layout, 30)

    # ── plots ───────────────────────────────────────────────────
    # (a) MPNP(solid) vs CatINT(dash) + Experiment(dot)
    lines!(ax_pol, df_pol[!, :Voltage],    abs.(df_pol[!, :Current]    .* (cm^2/mA));
        color=col_mpnp,   linewidth=solid_lw)
    lines!(ax_pol, df_pol_pr[!, :Voltage], abs.(df_pol_pr[!, :Current]);
        color=col_catint, linewidth=dash_lw, linestyle=:dash)
    scatter!(ax_pol, df_pol_exp[!, :Voltage], abs.(df_pol_exp[!, :Current]);
        color=col_exp, markersize=12)

    # (b)(c) species: MPNP(solid, species color) vs CatINT(dash, pastel)
    for ia in 1:7
        lines!(ax_act, vgrid_act, max.(act_electrode[ia, :], eps(Float64));
            color=sp_colors[ia], linewidth=solid_lw)
        ma = df_act_pr[!, :Index] .== ia
        any(ma) && lines!(ax_act, df_act_pr[ma, :Voltage],
            max.(df_act_pr[ma, :Concentration], eps(Float64));
            color=sp_pastel[ia], linewidth=dash_lw, linestyle=:dash)

        lines!(ax_con, vgrid_con, max.(conc_electrode[ia, :], eps(Float64));
            color=sp_colors[ia], linewidth=solid_lw)
        mc = df_conc_pr[!, :Index] .== ia
        any(mc) && lines!(ax_con, df_conc_pr[mc, :Voltage],
            max.(df_conc_pr[mc, :Concentration], eps(Float64));
            color=sp_pastel[ia], linewidth=dash_lw, linestyle=:dash)
    end

    # (d) MPNP (solid, blue) only
    lines!(ax_cdl, df_cdl[!, :Voltage] .- 0.16,
           df_cdl[!, :Capacitance] / (μF/cm^2);
        color=col_mpnp, linewidth=solid_lw)

    # ── legends ─────────────────────────────────────────────────
    # (a): MPNP / CatINT / Experiment
    pol_elems  = [LineElement(color=col_mpnp,   linewidth=solid_lw),
                  LineElement(color=col_catint, linewidth=dash_lw, linestyle=:dash),
                  MarkerElement(color=col_exp, marker=:circle, markersize=14)]
    pol_labels = ["MPNP@LiquidElectrolytes.jl", "CatINT", "Experiment"]

    axislegend(ax_pol, pol_elems, pol_labels;
        position        = :lb,
        orientation     = :vertical,
        labelsize       = 28,
        font            = :bold,
        framevisible    = true,
        backgroundcolor = (:white, 0.5),
        framecolor      = (:black, 0.5),
        patchsize       = (45, 22),
    )
    text!(ax_cdl, 0.4, 10; text="MPNP",
          color=col_mpnp, fontsize=28, font=:bold,
          align=(:center, :center))
    # ── species text labels — both (b) activity and (c) concentration ─
    sp_pos_con = [          # (coordinates are relative to (c) concentration)
        (-1.35,  1e-1),     # K⁺
        (-1.05,  9e-9),     # H⁺
        (-0.8,   2e-4),     # HCO₃⁻
        (-1.40,  2e-11),    # CO₃²⁻
        (-1.35,  7e-7),     # CO₂
        (-0.90,  9e-10),    # OH⁻
        (-1.35,  1e-4),     # CO
    ]
    for ia in 1:7
        x, y = sp_pos_con[ia]
        text!(ax_con, x, y; text=sp_rich[ia], color=sp_colors[ia],
              fontsize=32, font=:bold, offset=(0, 0))
    end
    # same labels on (b) activity too (only y positions adjusted to the activity scale)
    sp_pos_act = [
        (-1.35,  3e1),      # K⁺
        (-1.00,  9e-7),     # H⁺
        (-0.8,   6e-4),     # HCO₃⁻
        (-1.40,  2e-9),     # CO₃²⁻
        (-1.35,  7e-5),     # CO₂
        (-0.90,  9e-8),     # OH⁻
        (-1.35,  4e-2),     # CO
    ]
    for ia in 1:7
        x, y = sp_pos_act[ia]
        text!(ax_act, x, y; text=sp_rich[ia], color=sp_colors[ia],
              fontsize=32, font=:bold, offset=(0, 0))
    end

    # ── panel labels ────────────────────────────────────────────
    for (pos, lbl, pad) in [
        (fig[1, 1, TopLeft()], "(a)", 130),
        (fig[1, 2, TopLeft()], "(b)", 0),
        (fig[2, 1, TopLeft()], "(c)", 130),
        (fig[2, 2, TopLeft()], "(d)", 0),
    ]
        Label(pos, lbl;
            fontsize = 30, font = :bold,
            padding  = (0, pad, 8, 0),
            halign   = :right, valign = :bottom,
        )
    end

    resize_to_layout!(fig)
    fig
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

    priv_color = (:gray60, 0.6)

    cdl_mpb_color  = "#A7C7E7"
    cdl_dgml_color = "#C3A6E0"
    exp_color      = "#555555"

    # ── MPB (solid) — replaced with Potassium_only ───────────────────
    df_conc    = CSV.read(raw"../data/output/Concentration_Robin_Stefan_γ_Potassium_only.csv", DataFrame)
    df_pol     = CSV.read(raw"../data/output/Polarization_Curve_Robin_Stefan_γ_pnp_Potassium_only.csv", DataFrame)
    df_cdl     = CSV.read(raw"../data/output/DLCap_Robin_Stefan_γ_pnp_Potassium_only_1.csv", DataFrame)
    df_act     = CSV.read(raw"../data/output/Activity_Curve_Robin_Stefan_γ_pnp_Potassium_only.csv", DataFrame)
    df_cap_exp      = CSV.read(raw"../data/Valette_Cap_data/NaClO4_0.005M.csv", DataFrame)
    df_cap_exp_dgml = CSV.read(raw"../data/output/DLCap_Robin_DMGL_γ_pnp_All_species_1.csv", DataFrame)
    # ── CatINT (dotted) — same as the first, long-format ────────────
    # Index mapping: 1=K⁺ 2=H⁺ 3=HCO₃⁻ 4=CO₃²⁻ 5=CO₂ 6=OH⁻ 7=CO
    df_conc_pr  = CSV.read(raw"../data/catmap_CO2R_data/voltage-conc.csv",   DataFrame)  # Index,Voltage,Concentration
    df_act_pr   = CSV.read(raw"../data/catmap_CO2R_data/voltage-activ.csv",  DataFrame)  # Index,Voltage,Concentration
    df_pol_pr   = CSV.read(raw"../data/catmap_CO2R_data/Ringe-theorical.csv", DataFrame; header=[:Voltage, :Current])

    species       = ["K⁺","H⁺","HCO₃⁻","CO₃²⁻","CO₂","OH⁻","CO"]
    colors        = [:orange, :gray, :brown, :violet, :red, :green, :blue]
    colors_pastel = ["#F7C97F","#D3D3D3","#D9C2A7","#FFC5D3","#FF746C","#80EF80","#AFCBFF"]

    vgrid_con    = df_conc[!, :Voltage]
    vgrid_act    = df_act[!, :Voltage]
    conc_electrode    = permutedims(Matrix(df_conc[!, Symbol.(species)]))
    act_electrode     = permutedims(Matrix(df_act[!, Symbol.(species)]))

    # (a)(b)(c) shared model legend: CatINT(dash) / MPB@LiquidElectrolytes.jl(solid)
    model_elems  = [
        LineElement(color = priv_color, linewidth = 1.8, linestyle = :dash),
        LineElement(color = :black,     linewidth = 3,   linestyle = :solid),
    ]
    model_labels = ["CatINT", "MPB@LiquidElectrolyte.jl"]

    # (d) capacitance legend: distinguished by alpha+linewidth (no dash), 3 pastel colors
    cdl_mpb_lw  = 7    # previous model: thick and light
    cdl_dgml_lw = 3.5  # current model: thin and dark
    cdl_elems  = [
        LineElement(color = (cdl_mpb_color, 0.55), linewidth = cdl_mpb_lw),
        LineElement(color = (cdl_dgml_color, 1.0), linewidth = cdl_dgml_lw),
        MarkerElement(color = exp_color, marker = :circle, markersize = 10),
    ]
    cdl_labels = ["MPB@LiquidElectrolyte.jl", "DGML@LiquidElectrolyte.jl", "Experiment@Valette"]

    xt_bottom = [-1.4, -1.2, -1.0, -0.8, -0.6]
    xt_top    = [-0.8, -0.4,  0.0,  0.4,  0.8]

    # figure for the vertical stack
    fig = Figure(size = (900, 1500), figure_padding = (40, 40, 30, 45))

    xlabel_str = L"\textbf{Voltage}\ U \; \mathrm{(V \; vs. \; SHE)}"

    # ── (a) Partial CO Current  →  fig[1,1] ────────────────────
    ax_pol = Axis(fig[1, 1];
        ylabel = L"|I_{CO}| \; \mathrm{(mA \; cm^{-2})}",
        limits = (-1.5, -0.4, 1e-10, 100),
        yscale = log10,
    )
    style_axis!(ax_pol; big=false)
    ax_pol.xticks = (xt_bottom, [@sprintf("%.1f", x) for x in xt_bottom])
    ax_pol.yticks = (10.0 .^ (0:-5:-10),
                     [L"10^{0}", L"10^{-5}", L"10^{-10}"])
    ax_pol.xticklabelsvisible = false

    # ── (b) Interfacial Concentration  →  fig[2,1] ─────────────
    ax_con = Axis(fig[2, 1];
        ylabel = L"c_\alpha^{\ddagger} \; \mathrm{(M)}",
        yscale = log10,
        limits = (-1.5, -0.5, 1e-11, 1e1),
    )
    style_axis!(ax_con; big=true)
    ax_con.xticks = (xt_bottom, [@sprintf("%.1f", x) for x in xt_bottom])
    ax_con.yticks = (10.0 .^ (0:-3:-9),
                     [L"10^{0}", L"10^{-3}", L"10^{-6}", L"10^{-9}"])
    ax_con.xticklabelsvisible = false

    # ── (c) Interfacial Activity  →  fig[3,1]  (x-axis label here) ─
    ax_act = Axis(fig[3, 1];
        xlabel = xlabel_str,
        ylabel = L"a_\alpha^{\ddagger}",
        yscale = log10,
        limits = (-1.5, -0.5, 1e-11, 1e4),
    )
    style_axis!(ax_act; big=true)
    ax_act.xticks = (xt_bottom, [@sprintf("%.1f", x) for x in xt_bottom])
    ax_act.yticks = (10.0 .^ (3:-3:-9),
                     [L"10^{3}", L"10^{0}", L"10^{-3}", L"10^{-6}", L"10^{-9}"])

    # ── (d) Differential Capacitance  →  fig[4,1] ──────────────
    # x-axis is potential relative to pzc (U - U_pzc); remove y-axis limits (auto)
    xlabel_pzc = L"\textbf{Voltage}\ U - U_\mathrm{pzc} \; \mathrm{(V)}"
    ax_cdl = Axis(fig[4, 1];
        xlabel = xlabel_pzc,
        ylabel = L"C_\mathrm{dl} \; (\mu\mathrm{F \; cm^{-2}})",
    )
    style_axis!(ax_cdl; big=false)
    xt_cdl = [-0.8, -0.4, 0.0, 0.4, 0.8]
    ax_cdl.xticks = (xt_cdl .- 0.16, [@sprintf("%.1f", x) for x in xt_cdl])

    # row height ratio 2:3:3:2
    rowsize!(fig.layout, 1, Relative(2/10))
    rowsize!(fig.layout, 2, Relative(3/10))
    rowsize!(fig.layout, 3, Relative(3/10))
    rowsize!(fig.layout, 4, Relative(2/10))
    rowgap!(fig.layout, 14)

    solid_linewidth = 3.5
    dash_linewidth  = 2.0

    # ── MPB (solid) ────────────────────────────────────────────
    for ia in 1:7
        lines!(ax_con, vgrid_con, max.(conc_electrode[ia, :], eps(Float64));
            color=colors[ia], linewidth=solid_linewidth)
        lines!(ax_act, vgrid_act, max.(act_electrode[ia, :],  eps(Float64));
            color=colors[ia], linewidth=solid_linewidth)
    end

    # ── CatINT (dashed) — extract per species from long-format ───────────
    for ia in 1:7
        # concentration
        mc = df_conc_pr[!, :Index] .== ia
        if any(mc)
            lines!(ax_con, df_conc_pr[mc, :Voltage],
                max.(df_conc_pr[mc, :Concentration], eps(Float64));
                color=colors_pastel[ia], linewidth=dash_linewidth, linestyle=:dash)
        end
        # activity
        ma = df_act_pr[!, :Index] .== ia
        if any(ma)
            lines!(ax_act, df_act_pr[ma, :Voltage],
                max.(df_act_pr[ma, :Concentration], eps(Float64));
                color=colors_pastel[ia], linewidth=dash_linewidth, linestyle=:dash)
        end
    end

    # ── Partial CO current ─────────────────────────────────────
    lines!(ax_pol, df_pol[!, :Voltage],    abs.(df_pol[!, :Current]    .* (cm^2/mA));
        color="#8ED1C6", linewidth=solid_linewidth)
    lines!(ax_pol, df_pol_pr[!, :Voltage], abs.(df_pol_pr[!, :Current]);
        color=priv_color, linewidth=dash_linewidth, linestyle=:dash)

    # ── Differential capacitance — distinguished by alpha+linewidth (no dash) ──
    # MPB (previous model): thick and light, drawn first so it sits behind
    lines!(ax_cdl, df_cdl[!, :Voltage] .- 0.16, df_cdl[!, :Capacitance] / (μF/cm^2);
        color=(cdl_mpb_color, 0.55), linewidth=cdl_mpb_lw)
    # DGML (current model): thin and dark
    lines!(ax_cdl, df_cap_exp_dgml[!, :Voltage] .- 0.16, df_cap_exp_dgml[!, :Capacitance] / (μF/cm^2);
        color=(cdl_dgml_color, 1.0), linewidth=cdl_dgml_lw)
    # Experiment
    scatter!(ax_cdl, df_cap_exp[!, :Voltage] .+ 0.972, df_cap_exp[!, :Cdl];
        color=exp_color, markersize=7)

    # ── Legends ────────────────────────────────────────────────
    # (a) panel: shared model (CatINT / MPB)
    Legend(fig[1, 1], [model_elems], [model_labels], ["Model"];
        framevisible = false,
        nbanks       = 1,
        tellwidth    = false,
        tellheight   = false,
        halign       = :left,
        valign       = :bottom,
        margin       = (10, 10, 10, 10),
        labelsize    = 22,
        titlesize    = 24,
        patchsize    = (40, 20),
    )

    # (d) panel: MPB / DGML / Experiment
    Legend(fig[4, 1], [cdl_elems], [cdl_labels], ["Capacitance"];
        framevisible = false,
        nbanks       = 1,
        tellwidth    = false,
        tellheight   = false,
        halign       = :left,
        valign       = :top,
        margin       = (10, 10, 10, 10),
        labelsize    = 22,
        titlesize    = 24,
        patchsize    = (40, 20),
    )

    # species labels — (b) concentration panel (ax_con)
    text!(ax_con, -1.35,  3e-1;   text=L"\mathrm{K^+}",       color=colors[1], fontsize=24, font=:bold)
    text!(ax_con, -0.90,  9e-9;   text=L"\mathrm{H^+}",       color=colors[2], fontsize=24, font=:bold)
    text!(ax_con, -1.30,  2e-11;  text=L"\mathrm{CO_3^{2-}}", color=colors[4], fontsize=24, font=:bold)
    text!(ax_con, -1.055, 4e-6;   text=L"\mathrm{HCO_3^-}",   color=colors[3], fontsize=24, font=:bold)
    text!(ax_con, -1.3,  7e-7;   text=L"\mathrm{CO_2}",      color=colors[5], fontsize=24, font=:bold)
    text!(ax_con, -0.90,  9e-10;  text=L"\mathrm{OH^-}",      color=colors[6], fontsize=24, font=:bold)
    text!(ax_con, -1.35,  1e-4; text=L"\mathrm{CO}",        color=colors[7], fontsize=24, font=:bold)

    # outer labels — left of each panel (col 0), rotated vertically
    Label(fig[1, 0], "Partial CO Current";
        rotation = π/2, fontsize = 24, font = :bold, tellheight = false)
    Label(fig[2, 0], "Interfacial Concentration";
        rotation = π/2, fontsize = 24, font = :bold, tellheight = false)
    Label(fig[3, 0], "Interfacial Activity";
        rotation = π/2, fontsize = 24, font = :bold, tellheight = false)
    Label(fig[4, 0], "Differential Capacitance";
        rotation = π/2, fontsize = 24, font = :bold, tellheight = false)

    # panel labels
    for (pos, lbl) in [(fig[1,1,TopLeft()], "(a)"), (fig[2,1,TopLeft()], "(b)"),
                       (fig[3,1,TopLeft()], "(c)"), (fig[4,1,TopLeft()], "(d)")]
        Label(pos, lbl; fontsize=25, font=:bold,
              padding=(0,45,8,0), halign=:right, valign=:bottom)
    end

    resize_to_layout!(fig)
    fig
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
# ╠═7ba7df52-eee1-4f36-a5d0-379f31c6c5ef
# ╠═ccb13e5a-5ceb-4a56-9f3f-14307fc78115
# ╠═b0a4b942-5654-4b22-8849-90bf6c7f957f
# ╠═fc8096e4-01ca-451a-87cd-7e2e0171a531
# ╠═05eb8a6f-d7b5-4f37-a905-2209323afe2d
# ╠═5bb2da63-cc7d-4fcf-98b0-d2f4fa1dc1ef
# ╠═2ddeace6-663d-4fd1-8494-f1c88c19e628
# ╠═caa05490-0c7d-44ec-9be8-73f7a4473d8a
# ╠═2b10e19b-c099-4e4e-9cc1-5be1255d03bd
# ╠═4a68d318-ab6f-45b3-ad63-33d4a77c534d
# ╟─8545d818-d255-4e8a-af8f-72fdaf9d4bdd
# ╠═dcb38fd2-23b6-481e-83b7-4f3ee332c4ee
# ╠═ee091126-671e-46d8-8ed5-9457aeefd8fb
# ╟─4b02a8c0-a854-436a-9fdc-99ecf637e858
# ╟─59b13ba5-c85d-4e31-bd06-e1541b426205
# ╟─95844847-6f46-4840-979d-c54282295f41
# ╠═a09d4a6a-1017-4792-b684-7f9ddfeba83b
# ╠═4d250254-b8a8-4ff2-95c6-81a02ddfd982
# ╠═ac4deeb9-812a-429c-afe9-389b22814ee6
# ╠═9d660147-670a-4728-b2ef-8bb978f6981b
# ╠═400709f7-2a7c-431a-96ce-99a678b840f7
# ╠═d7ad452c-d3c1-485c-92ac-0b006dfb5acd
# ╟─89112dd8-08c8-4887-90db-e2e9c5ce0d88
# ╠═17b95990-f6e8-4da5-909a-1a2d46d7494d
# ╠═0290138e-4839-49a0-9a77-050bb4ab2cc0
