# tiff plot code for publication

tiff_plot <- function(filename){
  tiff(paste(filename,".tiff",sep=""),width=2000,height=1400,compression="lzw",res=200)
}
