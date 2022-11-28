plotUmatrix2<-function (Matrix = NULL, BestMatches = NULL, Cls = NULL, ClsColors = NULL, 
          BmSize = 2, DrawLegend = F,Title = NULL) 
{
  if (is.null(ClsColors)) 
    ClsColors = rainbow(length(unique(Cls)))
  Colormap = c("#3C6DF0", "#3C6DF0", "#3C6DF0", "#006602", 
                 "#006A02", "#006D01", "#007101", "#007501", "#007901", 
                 "#007C00", "#008000", "#068103", "#118408", "#0B8305", 
                 "#17860A", "#1D870D", "#228810", "#288A12", "#2E8B15", 
                 "#348D18", "#398E1A", "#3F8F1D", "#45911F", "#4A9222", 
                 "#509325", "#569527", "#5C962A", "#61982C", "#67992F", 
                 "#6D9A32", "#729C34", "#789D37", "#7E9F39", "#84A03C", 
                 "#89A13F", "#8FA341", "#95A444", "#9AA547", "#A0A749", 
                 "#A6A84C", "#ACAA4E", "#B1AB51", "#B7AC54", "#BDAE56", 
                 "#C3AF59", "#C8B15B", "#CEB25E", "#CBAF5C", "#C8AC59", 
                 "#C5A957", "#C3A654", "#C0A352", "#BDA050", "#BA9D4D", 
                 "#B7994B", "#B49648", "#B29346", "#AF9044", "#AC8D41", 
                 "#A98A3F", "#A6873C", "#A3843A", "#A08138", "#9E7E35", 
                 "#9B7B33", "#987830", "#95752E", "#92722B", "#8F6E29", 
                 "#8C6B27", "#8A6824", "#876522", "#84621F", "#815F1D", 
                 "#7E5C1B", "#7B5918", "#795616", "#765313", "#714E0F", 
                 "#6C480B", "#674307", "#6F4D15", "#785822", "#806230", 
                 "#896D3E", "#91774C", "#998159", "#A28C67", "#AA9675", 
                 "#B3A183", "#BBAB90", "#C3B59E", "#CCC0AC", "#D4CABA", 
                 "#DDD5C7", "#E5DFD5", "#E7E1D8", "#E9E4DB", "#EBE6DE", 
                 "#ECE8E1", "#EEEAE4", "#F0EDE7", "#F2EFEA", "#F4F1ED", 
                 "#F6F4F0", "#F8F6F3", "#F9F8F6", "#FBFAF9", "#FDFDFC", 
                 "#FFFFFF", "#FFFFFF", "#FEFEFE", "#FEFEFE", "#FEFEFE", 
                 "#FDFDFD", "#FDFDFD", "#FDFDFD", "#FCFCFC", "#FCFCFC", 
                 "#FCFCFC", "#FBFBFB", "#FBFBFB", "#FBFBFB", "#FAFAFA", 
                 "#FAFAFA", "#FAFAFA", "#F9F9F9", "#F9F9F9", "#FFFFFF", 
                 "#FFFFFF")
  
  BmSize = BmSize/2
  if (is.null(Matrix)) 
    stop("Matrix needs to be given")
  if (!is.matrix(Matrix)) {
    stop("Matrix has to be of type matrix")
  }
  if (length(dim(Matrix)) != 2) 
    stop("Matrix has to be a of type matrix, not an array")
  if (!is.null(Cls)) {
    if (length(ClsColors) < length(unique(Cls))) 
      ClsColors = c(ClsColors,rainbow((length(unique(Cls))-length(ClsColors))))
    if (length(unique(Cls)) > length(ClsColors)) 
      stop(paste("The amount of given Colors (ClsColors) is not enough for the Cls. The highest Cls Value is", 
                 length(unique(Cls)), "while the number of ClsColors only reaches to", 
                 length(ClsColors)))
  }
  
  quants = quantile(as.vector(Matrix), c(0.01, 0.5, 0.99),na.rm = T)
  minU = quants[1]
  maxU = quants[3]
  Matrix = (Matrix - quants[1])/(quants[3] - quants[1])
  indMax = which(Matrix > 1, arr.ind = T)
  indMin = which(Matrix < 0, arr.ind = T)
  if (length(indMax) > 0) 
    Matrix[indMax] = 1
  if (length(indMin) > 0) 
    Matrix[indMin] = 0

  if (!is.null(BestMatches)) {
    if (is.null(rownames(BestMatches))) 
      rownames(BestMatches) <- 1:nrow(BestMatches)
    BestMatchesFilter = rep(T, nrow(BestMatches))
  }
  nrows = nrow(Matrix)
  ncols = ncol(Matrix)
  colorx = Colormap
  for (i in 1:Nrlevels2) {
    Matrix[(Matrix >= levelBreaks[i]) & (Matrix <= levelBreaks[i +1])] = levelBreaks[i]
  }
  if (!is.null(BestMatches)) {
    dup = duplicated(BestMatches[, 1], fromLast = T) | duplicated(BestMatches[,1])
    duplicatedBestMatches = BestMatches[dup, , drop = F]
  }
  unnamedMatrix <- Matrix
  colnames(unnamedMatrix) <- NULL
  rownames(unnamedMatrix) <- NULL
  dfMatrix <- melt(unnamedMatrix)
  colnames(dfMatrix) <- c("y", "x", "z")
  dfMatrix2 = dfMatrix
  dfMatrix2[is.na(dfMatrix$z), 3] = 0
  MatrixPlot <- ggplot(dfMatrix, aes_string("x", "y")) + geom_tile(aes_string(fill = "z")) + 
    scale_fill_gradientn(colours = colorx, space = "Lab", 
                         na.value = "transparent") + stat_contour(aes_string(z = "z"), 
                                                                  data = dfMatrix2, bins = Nrlevels, size = 0.25, color = "black", 
                                                                  alpha = alpha, na.rm = F) + ylab("Lines (y)") + xlab("Columns (x)")
  if (!DrawLegend) 
    MatrixPlot <- MatrixPlot + theme(legend.position = "none")
  if (FixedRatio) 
    MatrixPlot <- MatrixPlot + coord_fixed(1, xlim = c(0.5, 
                                                       ncols + 0.5), ylim = c(0.5, nrows + 0.5), expand = 0) + 
    scale_y_reverse()
  else MatrixPlot <- MatrixPlot + coord_cartesian(xlim = c(0.5, 
                                                           ncols + 0.5), ylim = c(nrows + 0.5, 0.5), expand = 0) + 
    scale_y_reverse()
  
  if (!is.null(BestMatches)) {
    if (is.null(Cls)) 
      Cls <- rep(1, nrow(BestMatches))
    Cls = factor(Cls)
    uniqueClasses <- sort(na.last = T, unique(Cls))
    ClsColors = c("darkgreen", ClsColors)
    names(ClsColors) = 0:(length(ClsColors) - 1)
    MatrixPlot = MatrixPlot + scale_color_manual(values = ClsColors, 
                                                 name = "Clusters")
    if (!is.null(BestMatchesShape)) 
    d = data.frame(y = BestMatches[, 2], x = BestMatches[,3], class = Cls)
    MatrixPlot <- MatrixPlot + geom_point(aes_string(y = "y",x = "x", col = "factor(class)"), 
                                            data = d, stroke = BmSize)
  }
  return(MatrixPlot + ggtitle(Title) + scale_size_area())
}
