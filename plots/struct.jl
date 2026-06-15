function electrochemistry_theme()
    Theme(
        size = (960, 540),
        Axis = (
            spinewidth = 5.5,
            xtickwidth = 2.0,
            ytickwidth = 2.0,
            xticksize = 8,
            yticksize = 8,
            xlabelsize = 25,
            ylabelsize = 25,
            xticklabelsize = 25,
            yticklabelsize = 25,
            xgridvisible = false,
            ygridvisible = false,
            xlabelpadding = 10,
            ylabelpadding = 10,
            xlabelfont = :bold,
            ylabelfont = :bold,
			yticks = LinearTicks(5),
			xticks = LinearTicks(4)
        ),
        Lines = (
            linewidth = 3, 
        )
    )
end