
.onLoad <- function (libname=find.package("qgg"), pkgname="qgg") {
     OS <- .Platform$OS.type
     if(OS=="windows") dll <- paste(find.package("qgg"),"/libs/qgg.dll",sep="")    
     if(!OS=="windows") dll <- paste(find.package("qgg"),"/libs/qgg.so",sep="")
     if(file.exists(dll)) dyn.load(dll)
}

.onUnload <- function (libname=find.package("qgg"), pkgname="qgg") {
     OS <- .Platform$OS.type
     if(OS=="windows") dll <- paste(find.package("qgg"),"/libs/qgg.dll",sep="")    
     if(!OS=="windows") dll <- paste(find.package("qgg"),"/libs/qgg.so",sep="")
     if(file.exists(dll)) dyn.unload(dll)
}
