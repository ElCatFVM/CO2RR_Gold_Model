### A Pluto.jl notebook ###
# v1.0.0

using Markdown
using InteractiveUtils

# ╔═╡ a70cef7d-2a2f-4155-bdf3-fec9df94c63f
begin
    using Pkg
    Pkg.activate(joinpath(@__DIR__, ".."))
    using CSV
    using CairoMakie
    using Printf
    using DataFrames
	using DrWatson
end

# ╔═╡ baa3e08e-5d64-4c8f-9f6d-5fdb40e97bc5
begin
    using HypertextLiteral: @htl_str, @htl
    using UUIDs: uuid1
    using PlutoUI
end

# ╔═╡ 2623e9ed-252d-43a1-bc77-4c5bd9051653
sweepcomparedir(args...)=datadir("sweepcompare", args...)

# ╔═╡ 5691afd8-6583-41de-b33f-c752ef14815c
species = ["K⁺", "H⁺", "HCO₃⁻", "CO₃²⁻", "CO₂", "OH⁻", "CO"]


# ╔═╡ 492dcba1-44cd-424e-b890-532e3d33d058
colors = [:orange, :gray, :brown, :violet, :red, :green, :blue]


# ╔═╡ 55d73c4c-782a-49d6-b182-424356e2f379
function style!(ax)

    ax.spinewidth = 5.5
    ax.xtickwidth = 2.0
    ax.ytickwidth = 2.0
    ax.xticksize = 8
    ax.yticksize = 8
    ax.xlabelsize = 25
    ax.ylabelsize = 25
    ax.xticklabelsize = 25
    ax.yticklabelsize = 25
    ax.xgridvisible = false
    ax.ygridvisible = false
    ax.xlabelpadding = 10
    ax.ylabelpadding = 10
    return ax.xlabelfont = :bold


end

# ╔═╡ bbe2fa26-1b19-42c9-a2f7-fdd7cec85e65
function plot_iv_conc(csv_path)
    df = CSV.read(csv_path, DataFrame)

    fig = Figure(size = (700, 500))


    vgrid = df[!, :Voltage]  # <- if the CSV voltage column name differs, change only this

    ax = Axis(
        fig[1, 1];
        xlabel = L"\mathbf{\text{U}\ \mathrm{vs.}\ \text{SHE}\ (V)}",
        ylabel = L"\mathbf{c_i^{+}}\;(\mathrm{M})",
        yscale = log10,
        limits = ((-1.25, -0.5), (1.0e-11, 1.0e1)),
    )
    style!(ax)
    xt = [-1.2, -1.0, -0.8, -0.6]
    ax.xticks = (xt, [@sprintf("%.1f", x) for x in xt])

    yt_vals = 10.0 .^ (0:-3:-9)
    yt_lbls = [L"10^{0}", L"10^{-3}", L"10^{-6}", L"10^{-9}"]
    #  ax.yticks = (yt_vals, yt_lbls)

    for ia in 1:7
        y = max.(df[!, "c_" * species[ia]], eps(Float64))
        lines!(
            ax, vgrid, y; color = colors[ia], linewidth = 5,
            label = species[ia]
        )
    end
    axislegend(position = :rc, backgroundcolor = (:white, 0.5))
    return fig
end

# ╔═╡ 6886a522-fbe1-4251-8269-91cb8b8c5243
function plot_iv_conc(csv_path1, csv_path2)
    df1 = CSV.read(csv_path1, DataFrame)
    df2 = CSV.read(csv_path2, DataFrame)

    fig = Figure(size = (700, 500))


    vgrid1 = df1[!, :Voltage]  
    vgrid2 = df2[!, :Voltage] 

    ax = Axis(
        fig[1, 1];
        xlabel = L"\mathbf{\text{U}\ \mathrm{vs.}\ \text{SHE}\ (V)}",
        ylabel = L"\mathbf{c_i^{+}}\;(\mathrm{M})",
        yscale = log10,
        limits = ((-1.25, -0.5), (1.0e-11, 1.0e1)),
    )
    style!(ax)
    xt = [-1.2, -1.0, -0.8, -0.6]
    ax.xticks = (xt, [@sprintf("%.1f", x) for x in xt])

    yt_vals = 10.0 .^ (0:-3:-9)
    yt_lbls = [L"10^{0}", L"10^{-3}", L"10^{-6}", L"10^{-9}"]
    #  ax.yticks = (yt_vals, yt_lbls)

    for ia in 1:7
        y = max.(df1[!, "c_" * species[ia]], eps(Float64))
        lines!(
            ax, vgrid1, y; color = colors[ia], linewidth = 5,
            label = species[ia]
        )
        y = max.(df2[!, "c_" * species[ia]], eps(Float64))
        lines!(
            ax, vgrid2, y; color = colors[ia], linewidth = 5, linestyle = :dash
        )

    end
    axislegend(position = :rc, backgroundcolor = (:white, 0.5))
    return fig
end


# ╔═╡ f8929bed-15ef-483f-bd70-e4722460fa67
function plot_iv_act(csv_path)
    df = CSV.read(csv_path, DataFrame)

    fig = Figure(size = (700, 500))


    vgrid = df[!, :Voltage]  

    ax = Axis(
        fig[1, 1];
        xlabel = L"\mathbf{\text{U}\ \mathrm{vs.}\ \text{SHE}\ (V)}",
        ylabel = L"\mathbf{c_i^{+}}\;(\mathrm{M})",
        yscale = log10,
        limits = ((-1.25, -0.5), (1.0e-11, 1.0e6)),
    )
    style!(ax)
    xt = [-1.2, -1.0, -0.8, -0.6]
    ax.xticks = (xt, [@sprintf("%.1f", x) for x in xt])

    yt_vals = 10.0 .^ (0:-3:-9)
    yt_lbls = [L"10^{0}", L"10^{-3}", L"10^{-6}", L"10^{-9}"]
    #  ax.yticks = (yt_vals, yt_lbls)

    for ia in 1:7
        y = max.(df[!, "a_" * species[ia]], eps(Float64))
        lines!(
            ax, vgrid, y; color = colors[ia], linewidth = 5,
            label = species[ia]
        )
    end
    axislegend(position = :rc, backgroundcolor = (:white, 0.5))
    return fig
end


# ╔═╡ 7173dc96-f960-43e2-b7fc-4fddfe730ae1
function plot_iv_act(csv_path1, csv_path2)
    df1 = CSV.read(csv_path1, DataFrame)
    df2 = CSV.read(csv_path2, DataFrame)

    fig = Figure(size = (700, 500))


    vgrid1 = df1[!, :Voltage] 
    vgrid2 = df2[!, :Voltage]  
	
    ax = Axis(
        fig[1, 1];
        xlabel = L"\mathbf{\text{U}\ \mathrm{vs.}\ \text{SHE}\ (V)}",
        ylabel = L"\mathbf{a_i^{+}}\;(\mathrm{M})",
        yscale = log10,
        limits = ((-1.25, -0.5), (1.0e-11, 1.0e6)),
    )
    style!(ax)
    xt = [-1.2, -1.0, -0.8, -0.6]
    ax.xticks = (xt, [@sprintf("%.1f", x) for x in xt])

    yt_vals = 10.0 .^ (0:-3:-9)
    yt_lbls = [L"10^{0}", L"10^{-3}", L"10^{-6}", L"10^{-9}"]
    #  ax.yticks = (yt_vals, yt_lbls)

    for ia in 1:7
        y = max.(df1[!, "a_" * species[ia]], eps(Float64))
        lines!(
            ax, vgrid1, y; color = colors[ia], linewidth = 5,
            label = species[ia]
        )
        y = max.(df2[!, "a_" * species[ia]], eps(Float64))
        lines!(
            ax, vgrid2, y; color = colors[ia], linewidth = 5, linestyle = :dash
        )

    end
    axislegend(position = :rc, backgroundcolor = (:white, 0.5))
    return fig
end


# ╔═╡ 21d89357-65f2-4ed3-96ad-9075635ffdbb
function plot_iv_gamma(csv_path1, csv_path2)
    df1 = CSV.read(csv_path1, DataFrame)
    df2 = CSV.read(csv_path2, DataFrame)

    fig = Figure(size = (700, 500))


    vgrid1 = df1[!, :Voltage]  
    vgrid2 = df2[!, :Voltage]
	
    ax = Axis(
        fig[1, 1];
        xlabel = L"\mathbf{\text{U}\ \mathrm{vs.}\ \text{SHE}\ (V)}",
        ylabel = L"\mathbf{γ_i^{+}}\;(\mathrm{M})",
        #        yscale = log10,
    )
    xlims!(ax,(-1.25, -0.5))
	style!(ax)
    xt = [-1.2, -1.0, -0.8, -0.6]
    ax.xticks = (xt, [@sprintf("%.1f", x) for x in xt])

    yt_vals = 10.0 .^ (0:-3:-9)
    yt_lbls = [L"10^{0}", L"10^{-3}", L"10^{-6}", L"10^{-9}"]
    #  ax.yticks = (yt_vals, yt_lbls)

    for ia in 2:7
        y = max.(df1[!, "γ_" * species[ia]], eps(Float64))
        lines!(
            ax, vgrid1, y; color = colors[ia], linewidth = 5,
            label = species[ia]
        )
        y = max.(df2[!, "γ_" * species[ia]], eps(Float64))
        lines!(
            ax, vgrid2, y; color = colors[ia], linewidth = 5, linestyle = :dash
        )

    end
    axislegend(position = :rc, backgroundcolor = (:white, 0.5))
    return fig
end


# ╔═╡ 2b0ce48e-2fdb-4717-ada3-6a997380a1ed
md"""
## Concentrations: Stefan\_γ (full) vs DGML\_γ (dashed)
"""

# ╔═╡ 4dd33442-0641-42b8-9ee7-b5b9de38da71
plot_iv_conc(
    sweepcomparedir("sweep_iv_actcoeff=Stefan_γ!_bcmodel=Robin_hydrated=true_model=Gold.csv"),
    sweepcomparedir("sweep_iv_actcoeff=DGML_γ!_bcmodel=Robin_hydrated=true_model=Gold.csv"),
)

# ╔═╡ 7df90be5-a0f2-4e3b-ba75-d86196fa3e6c
md"""
## Activities: Stefan\_γ (full) vs DGML\_γ (dashed)
"""

# ╔═╡ 57569741-f178-4196-8b7f-2efc2224261d
plot_iv_act(
    sweepcomparedir("sweep_iv_actcoeff=Stefan_γ!_bcmodel=Robin_hydrated=true_model=Gold.csv"),
    sweepcomparedir("sweep_iv_actcoeff=DGML_γ!_bcmodel=Robin_hydrated=true_model=Gold.csv"),
)

# ╔═╡ 44a69c22-05da-4478-893e-05d51e0f1ba1
md"""
## Concentrations: DGML\_γ (full) vs Potassium\_γ (dashed)
"""

# ╔═╡ 9e2e5296-0ed5-4128-81ad-b546e17809f8
plot_iv_conc(
    sweepcomparedir("sweep_iv_actcoeff=DGML_γ!_bcmodel=Robin_hydrated=true_model=Gold.csv"),
    sweepcomparedir("sweep_iv_actcoeff=Potassium_γ!_bcmodel=Robin_hydrated=true_model=Gold.csv"),
)

# ╔═╡ 6b16b555-da01-44d5-b54b-52c5a0dd60d8
md"""
## Activities: DGML\_γ (full) vs Potassium\_γ (dashed)
"""

# ╔═╡ 92680105-811e-49d2-8d29-d772226777d2
plot_iv_act(
    sweepcomparedir("sweep_iv_actcoeff=DGML_γ!_bcmodel=Robin_hydrated=true_model=Gold.csv"),
    sweepcomparedir("sweep_iv_actcoeff=Potassium_γ!_bcmodel=Robin_hydrated=true_model=Gold.csv"),
)

# ╔═╡ 2f8e0e31-0123-4423-a465-4be05cfbac5f
plot_iv_conc(
    sweepcomparedir("sweep_iv_actcoeff=Stefan_γ!_bcmodel=Robin_hydrated=true_model=Gold.csv"),
    sweepcomparedir("sweep_iv_actcoeff=Potassium_γ!_bcmodel=Robin_hydrated=true_model=Gold.csv"),
)

# ╔═╡ 24b0e0f8-67e4-4cce-b712-15a645c5f86e
plot_iv_act(
    sweepcomparedir("sweep_iv_actcoeff=Stefan_γ!_bcmodel=Robin_hydrated=true_model=Gold.csv"),
    sweepcomparedir("sweep_iv_actcoeff=Potassium_γ!_bcmodel=Robin_hydrated=true_model=Gold.csv"),
)

# ╔═╡ daacc29e-20ac-4520-806d-5169356d8892
md"""
## IV curves
"""

# ╔═╡ 09bfabe4-ee58-4e84-af37-eedd7ad0d2b4
let
	df0=CSV.read(datadir("dataplotfiles","iv_GoldModel.csv"), DataFrame)
	df1=CSV.read(sweepcomparedir("sweep_iv_actcoeff=DGML_γ!_bcmodel=Robin_hydrated=true_model=Gold.csv"), DataFrame)
	df2=CSV.read(sweepcomparedir("sweep_iv_actcoeff=Stefan_γ!_bcmodel=Robin_hydrated=true_model=Gold.csv"), DataFrame)
	df3=CSV.read(sweepcomparedir("sweep_iv_actcoeff=Potassium_γ!_bcmodel=Robin_hydrated=true_model=Gold.csv"), DataFrame)

	fig=Figure(size=(700,300))
	ax=Axis(fig[1,1])


	lines!(ax,df0[!,"x"], -df0[!,"y"], label="GoldModel.csv")
	lines!(ax,df1[!,"Voltage"], -df1[!,"Current"]/10, label="DGML_γ")
	lines!(ax,df2[!,"Voltage"], -df2[!,"Current"]/10, label="Stefan_γ")
	lines!(ax,df3[!,"Voltage"], -df3[!,"Current"]/10, label="Potassium_γ")
	axislegend()
	fig
end

# ╔═╡ 8af12f1c-d35b-4cc9-8185-1bb5adbb69e8
html"""<hr>"""

# ╔═╡ 784b4c3e-bb2a-4940-a83a-ed5e5898dfd4
html"""<style>.dont-panic{ display: none }</style>"""

# ╔═╡ afe4745f-f9f1-4e23-8735-cbec6fb79c41
begin
    function floataside(text::Markdown.MD; top = 1)
        uuid = uuid1()
        return @htl(
            """
            		<style>


            		@media (min-width: calc(700px + 30px + 300px)) {
            			aside.plutoui-aside-wrapper-$(uuid) {

            	color: var(--pluto-output-color);
            	position:fixed;
            	right: 1rem;
            	top: $(top)px;
            	width: 400px;
            	padding: 10px;
            	border: 3px solid rgba(0, 0, 0, 0.15);
            	border-radius: 10px;
            	box-shadow: 0 0 11px 0px #00000010;
            	/* That is, viewport minus top minus Live Docs */
            	max-height: calc(100vh - 5rem - 56px);
            	overflow: auto;
            	z-index: 40;
            	background-color: var(--main-bg-color);
            	transition: transform 300ms cubic-bezier(0.18, 0.89, 0.45, 1.12);

            			}
            			aside.plutoui-aside-wrapper > div {
            #				width: 300px;
            			}
            		}
            		</style>

            		<aside class="plutoui-aside-wrapper-$(uuid)">
            		<div>
            		$(text)
            		</div>
            		</aside>

            		"""
        )
    end
    floataside(stuff; kwargs...) = floataside(md"""$(stuff)"""; kwargs...)
end;


# ╔═╡ b8fd36a7-d8d1-45f7-b66e-df9132168bfc
# https://discourse.julialang.org/t/adding-a-restart-process-button-in-pluto/76812/5
restart_button() = html"""
<script>
	const button = document.createElement("button")

	button.addEventListener("click", () => {
		editor_state_set(old_state => ({
			notebook: {
				...old_state.notebook,
				process_status: "no_process",
			},
		})).then(() => {
			window.requestAnimationFrame(() => {
				document.querySelector("#process_status a").click()
			})
		})
	})
	button.innerText = "Restart notebook"

	return button
</script>
""";

# ╔═╡ Cell order:
# ╠═a70cef7d-2a2f-4155-bdf3-fec9df94c63f
# ╠═2623e9ed-252d-43a1-bc77-4c5bd9051653
# ╠═5691afd8-6583-41de-b33f-c752ef14815c
# ╠═492dcba1-44cd-424e-b890-532e3d33d058
# ╠═55d73c4c-782a-49d6-b182-424356e2f379
# ╠═bbe2fa26-1b19-42c9-a2f7-fdd7cec85e65
# ╠═6886a522-fbe1-4251-8269-91cb8b8c5243
# ╠═f8929bed-15ef-483f-bd70-e4722460fa67
# ╠═7173dc96-f960-43e2-b7fc-4fddfe730ae1
# ╠═21d89357-65f2-4ed3-96ad-9075635ffdbb
# ╟─2b0ce48e-2fdb-4717-ada3-6a997380a1ed
# ╠═4dd33442-0641-42b8-9ee7-b5b9de38da71
# ╟─7df90be5-a0f2-4e3b-ba75-d86196fa3e6c
# ╠═57569741-f178-4196-8b7f-2efc2224261d
# ╠═44a69c22-05da-4478-893e-05d51e0f1ba1
# ╠═9e2e5296-0ed5-4128-81ad-b546e17809f8
# ╠═6b16b555-da01-44d5-b54b-52c5a0dd60d8
# ╠═92680105-811e-49d2-8d29-d772226777d2
# ╠═2f8e0e31-0123-4423-a465-4be05cfbac5f
# ╠═24b0e0f8-67e4-4cce-b712-15a645c5f86e
# ╟─daacc29e-20ac-4520-806d-5169356d8892
# ╠═09bfabe4-ee58-4e84-af37-eedd7ad0d2b4
# ╟─8af12f1c-d35b-4cc9-8185-1bb5adbb69e8
# ╟─baa3e08e-5d64-4c8f-9f6d-5fdb40e97bc5
# ╟─784b4c3e-bb2a-4940-a83a-ed5e5898dfd4
# ╟─afe4745f-f9f1-4e23-8735-cbec6fb79c41
# ╟─b8fd36a7-d8d1-45f7-b66e-df9132168bfc
