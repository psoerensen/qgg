auc <- function(yobs=NULL, ypred=NULL) {
     n0 <- length(yobs[yobs==0])
     n1 <- length(yobs[yobs==1])
     y <- cbind(yobs, ypred)
     y <- y[order(y[,2], decreasing=TRUE),]
     y <- cbind(y, seq(from=nrow(y), to=1))
     rd <- mean(y[y[,1]==1,][,3])
     auc <- (1/n0)*(rd-(n1/2)-(1/2))
     auc 
}


rnag <- function(yobs=NULL,ypred=NULL) {
     fit0 <- glm(yobs~1,family=binomial(link='logit'))
     fit1 <- glm(yobs~1+ypred,family=binomial(link='logit'))
     n <- length(yobs)
     LR <- anova(fit1)$Deviance[2]
     L0 <-  as.numeric(logLik(fit0))
     r2nag <- (1-exp(-LR/n))/(1-exp(-(-2*L0)/n))
     return(r2nag)
}

qcpred <- function(yobs=NULL,ypred=NULL) {
     r2 <- summary(lm(yobs ~ ypred))$r.squared
     pa <- cor(ypred, yobs)
     mspe <- sum((ypred - yobs)^2)/length(yobs)
     intercept <- lm(yobs ~ ypred )$coef[1]
     slope <- lm(yobs ~ ypred)$coef[2]
     aurc <- r2nag <- NA
     if(nlevels(factor(yobs))==2) aurc <- auc(yobs=yobs,ypred=ypred)
     if(nlevels(factor(yobs))==2) r2nag <- rnag(yobs=yobs,ypred=ypred)
     res <- round(c(pa,r2,r2nag,aurc,intercept,slope,mspe),3)
     names(res) <- c("Corr","R2","Nagel R2", "AUC", "intercept", "slope", "MSPE")
     return(res)
}

