auc <- function(yObs=NULL, yPred=NULL, plot=FALSE, fOut=NULL){
     n0 <- length(yObs[yObs==0])
     n1 <- length(yObs[yObs==1])
     y <- cbind(yObs, yPred)
     y <- y[order(y[,2], decreasing=TRUE),]
     y <- cbind(y, seq(from=nrow(y), to=1))
     rd <- mean(y[y[,1]==1,][,3])
     auc <- (1/n0)*(rd-(n1/2)-(1/2))
     auc 
}

