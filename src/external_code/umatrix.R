library(ggplot2)
library(reshape2)
library(dplyr)
library(plotly)
library(lattice)
library(rgl)

plot_umatrix<-function(res_esomtrain,single_cell=FALSE,emphasis_edges=FALSE,plot_ly=FALSE){
  quants = quantile(as.vector(res_esomtrain$Umatrix), c(0.15,0.975),na.rm = T)
  Matrix = (res_esomtrain$Umatrix - quants[1])/(quants[2] - quants[1])
  indMax = which(Matrix > 1, arr.ind = T)
  indMin = which(Matrix < 0, arr.ind = T)
  if (length(indMax) > 0) 
    Matrix[indMax] = 1
  if (length(indMin) > 0) 
    Matrix[indMin] = 0
  dt.umatrix<-melt(Matrix)
  matrix.bestm<-process_bestmatchs(res_esomtrain$BestMatches)
  
  if(emphasis_edges & single_cell){
    keept_nodes <- dt.umatrix[which(dt.umatrix$value<0.9),1:2]
    colnames(keept_nodes)<-c("V2","V3")
    matrix.bestm <- semi_join(matrix.bestm,keept_nodes,by=c("V2","V3"))
  }
  
  colfunc <- colorRampPalette(c("white","black"))
  g1<-ggplot(dt.umatrix,aes(x=Var2,y=Var1))+
    geom_tile(aes(fill=value))+
    scale_fill_gradientn(colours = colfunc(50))+
    stat_contour(aes(z = value),size = 0.25, color = "black")+
    theme(line = element_blank(),axis.text = element_blank(),title = element_blank(),panel.background=element_blank())
  if(single_cell){
    g1 <- g1+geom_point(data=matrix.bestm,aes(x=V3,y=V2,col=Key),stroke = 0.001,position = "jitter")
  }
  if(plot_ly){
    g1<-ggplotly(g1)
  }
  print(g1)
}

plot_umatrix_3D<- function(res_esomtrain,emphasis_edges=FALSE,gif_file=FALSE){
  # Pre-process grid
  umatplot <- process_umatrix(res_esomtrain$Umatrix)
  umatplot<-10*umatplot
  # Pre-process sphere
  bestm <- unique(process_bestmatchs(res_esomtrain$BestMatches))
  colnames(bestm) <- c("Group","Var1","Var2")
  sphere_df <- melt(umatplot)
  
  if(emphasis_edges){
    sphere_df <- sphere_df[sphere_df[,3]<(10*0.8),]
  }
  
  t1 <- inner_join(bestm,sphere_df,by=c("Var1","Var2"))
  # Color process
  zlim<-range(umatplot)
  zlen<-zlim[2]-zlim[1]+1
  colfunc <- colorRampPalette(c("white","black"))
  colorlut<-colfunc(zlen)
  col<-colorlut[umatplot-zlim[1]+1]
  colclass <- rainbow(length(unique(t1[,1])))
  colbestm <- colclass[as.factor(t1[,1])]
  t1$col<-colbestm
  # Plot
  open3d()
  rgl.surface(x = 1:nrow(umatplot), z = 1:ncol(umatplot), y = umatplot,color=col,lit=F)
  rgl.spheres(x=t1[,2],z=t1[,3],y=t1[,4],radius=0.2,color=t1$col)
  rgl.viewpoint(theta=15)
  
  if(gif_file){
    movie3d(spin3d(axis = c(0, 1, 0),rpm=2), duration = 30,
            dir = "~/Documents/VIB_work/vib_som/movie3d/")
  }
}

process_bestmatchs<-function(res_bestmatchs){
  matrix.bestm<-as.data.frame(res_bestmatchs)
  matrix.bestm[,1]<-as.character(matrix.bestm[,1])
  matrix.bestm[,2]<-as.numeric(as.character(matrix.bestm[,2]))
  matrix.bestm[,3]<-as.numeric(as.character(matrix.bestm[,3]))
  return(matrix.bestm)
}

process_umatrix <- function(res_umatrix){
  quants = quantile(as.vector(res_umatrix), c(0.15,0.975),na.rm = T)
  Matrix = (res_umatrix - quants[1])/(quants[2] - quants[1])
  indMax = which(Matrix > 1, arr.ind = T)
  indMin = which(Matrix < 0, arr.ind = T)
  if (length(indMax) > 0) 
    Matrix[indMax] = 1
  if (length(indMin) > 0) 
    Matrix[indMin] = 0
  return(Matrix)
}
