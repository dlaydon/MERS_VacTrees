


MakeLayoutMatrix 		= function(NoCols, 	NoRows, 														AddTitle = TRUE, AddColNames = TRUE, 	AddRowNames = TRUE, BY_ROW = FALSE, RemoveCorner = TRUE)
{
	## Start by making a matrix. Then add space for colnames, rownames and full titles as necessary. 
	LayoutMatrix = matrix(1:(NoRows*NoCols), ncol = NoCols, byrow = BY_ROW)
	if (AddRowNames)
	{
		LayoutMatrix = LayoutMatrix + 1
		LayoutMatrix = cbind(rep(1, dim(LayoutMatrix)[1]), 		LayoutMatrix)
	}
	if (AddColNames) 
	{
		LayoutMatrix = LayoutMatrix + 1
		LayoutMatrix = rbind(rep(1, dim(LayoutMatrix)[2]), 		LayoutMatrix)
	}
	if (AddTitle) 
	{
		LayoutMatrix = LayoutMatrix + 1
		LayoutMatrix = rbind(rep(1, dim(LayoutMatrix)[2]), 		LayoutMatrix)
	}
	
	if (RemoveCorner)
		if (AddRowNames & AddColNames) ## remove weird "corner that ruins everyting in rest of plotting code (specifically position of text/rownames/colnames etc.) 
		{
			if (AddTitle) 	LayoutMatrix[2,1] = max(LayoutMatrix) + 1 	else  	### remove the corner (that is below the title) hence 2,1
							LayoutMatrix[1,1] = max(LayoutMatrix) + 1 			### remove the corner hence 1,1
		}
	return(LayoutMatrix)
}
MakeHeightsVecForLayout	= function(			NoRows, Ratio_Plots_To_Title = 4, Ratio_Plots_To_ColNames = 2, 	AddTitle = TRUE, AddColNames = TRUE)
{
	Heights = rep(1, NoRows)	
	if (AddColNames) 	Heights = c(1/Ratio_Plots_To_ColNames	, Heights)
	if (AddTitle)  		Heights = c(1/Ratio_Plots_To_Title		, Heights)
	return(Heights)	
}
MakeWidthsVecForLayout 	= function(NoCols, 	Ratio_Plots_To_RowNames = 2, 																			AddRowNames = TRUE)
{
	Widths = rep(1, NoCols)
	if (AddRowNames) Widths = c(1/Ratio_Plots_To_RowNames, Widths)
	return(Widths)	
}

SetUpMultiPlot 			= function(NoCols = NULL, NoRows = NULL, ColNames = NULL, RowNames = NULL, MultiPlot_Title = NULL, BY_ROW = FALSE, RemoveCorner = TRUE, 
		Ratio_Plots_To_Title = 4, Ratio_Plots_To_ColNames = 2, Ratio_Plots_To_RowNames = 2, ColName_Cex = 1.4, RowName_Cex = 1.4, MultiPlot_Title_Cex = 2, OG_Mar = OrigMAR, modelrun = ModelRun)
{
	if (is.null(NoCols)) NoCols = length(ColNames)
	if (is.null(NoRows)) NoRows = length(RowNames)
	if (is.null(MultiPlot_Title	)) AddTitle 	= FALSE	else  AddTitle 		= TRUE	
	if (is.null(ColNames		)) AddColNames 	= FALSE else  AddColNames 	= TRUE  
	if (is.null(RowNames		)) AddRowNames 	= FALSE else  AddRowNames 	= TRUE  
	
	LayoutMatrix 	= MakeLayoutMatrix			(NoCols = NoCols, NoRows = NoRows, AddTitle = AddTitle, AddColNames = AddColNames, AddRowNames = AddRowNames, BY_ROW = BY_ROW, RemoveCorner = RemoveCorner)
	HEIGHTS 		= MakeHeightsVecForLayout	(NoRows = NoRows, Ratio_Plots_To_Title = Ratio_Plots_To_Title, Ratio_Plots_To_ColNames = Ratio_Plots_To_ColNames, AddTitle = AddTitle, AddColNames = AddColNames)
	WIDTHS			= MakeWidthsVecForLayout	(NoCols = NoCols, Ratio_Plots_To_RowNames = Ratio_Plots_To_RowNames, AddRowNames = AddRowNames)
	layout(LayoutMatrix, heights = HEIGHTS, widths = WIDTHS)
	
	if (AddTitle)
	{	
		par(mar = rep(0, 4)) ### without this you often get the annoying "Error in plot.new() : figure margins too large" error. Set margins to zero for the title, then reset them to OrigMAR (global variable defined at start of script). 
		plot.new()
		text(0.5, 0.5, MultiPlot_Title , cex = MultiPlot_Title_Cex, font = 2)
	}
	if (AddColNames)
	{
#		par(mar = rep(0, 4))
		plot(NA, ylim = 0:1, xlim = c(0,length(ColNames)) , yaxs='i', xaxs='i', xaxt = "n", yaxt = "n", bty="n")
		text(1:length(ColNames) - 0.5 , 0.5, ColNames, cex = ColName_Cex, font = 2)
	}
	if (AddRowNames)
	{
#		par(mar = rep(0, 4))
		plot(NA, xlim = 0:1, ylim = c(0,length(RowNames)), yaxs='i', xaxs='i', xaxt = "n", yaxt = "n", bty="n")
		text(0.5, (1:length(RowNames) - 0.5), rev(RowNames), cex = RowName_Cex, font = 2, srt = 90) ### don't know why you have to reverse the rownames, but it otherwise the names come out in the wrong order. 
	}
	
	# reset par ("mar")
	par(mar = OG_Mar) ### 
}
GetTransparentColours 	= function(ColourString, Alpha = 0.5)
{
	Cols 			= col2rgb(ColourString, alpha = FALSE) ### colours to be set.
	TransparentCols = rgb(red = Cols[1,]/255, green = Cols[2,]/255, blue = Cols[3,]/255, alpha = Alpha)
	
	return (TransparentCols)
}
CloseOpenPlotDevices 	= function()
{
	if (!is.null(dev.list()))	for (openplot in 1:length(dev.list())) dev.off()
}



