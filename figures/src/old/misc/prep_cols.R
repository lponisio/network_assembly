

prep.cols <- function(dats){
	cInt <- classIntervals(dats, style="fisher")
	cPal <- tim.colors(24)
	cols  <- findColours(cInt,cPal)
	cols   <- ifelse(cols == "#00008F", "white", cols ) 
	return(cols)
}