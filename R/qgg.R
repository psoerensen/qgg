
.onLoad <- function (libname=find.package("qgg"), pkgname="qgg") {
     OS <- .Platform$OS.type
     if(OS=="windows") dll <- paste(find.package("qgg"),"/libs/x64/qgg.dll",sep="")    
     #if(OS=="windows") dll <- paste(find.package("qgg"),"/libs/qgg.dll",sep="")    
     if(!OS=="windows") dll <- paste(find.package("qgg"),"/libs/qgg.so",sep="")
     if(file.exists(dll)) dyn.load(dll)
}

.onUnload <- function (libname=find.package("qgg"), pkgname="qgg") {
     OS <- .Platform$OS.type
     if(OS=="windows") dll <- paste(find.package("qgg"),"/libs/x64/qgg.dll",sep="")    
     #if(OS=="windows") dll <- paste(find.package("qgg"),"/libs/qgg.dll",sep="")    
     if(!OS=="windows") dll <- paste(find.package("qgg"),"/libs/qgg.so",sep="")
     if(file.exists(dll)) dyn.unload(dll)
}


srun --pty -c 1 --mem=6g --partition=short --job-name=Rint bash
source /com/extra/R/3.4/load.sh

install.packages("devtools")
library(devtools)
install_github("psoerensen/qgg")

