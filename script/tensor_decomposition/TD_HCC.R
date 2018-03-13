#HCC
require(rTensor)
x <- read.csv("../../data/cancer_data/cancer_expression_data/HCC/expression.csv",sep="\t")
Z <- array(NA,c(dim(x)[1],20,3))
Z[,,1] <- data.matrix(x[,-1][,1:20])
Z[,,2] <- data.matrix(x[,-1][,21:40])
Z[,,3] <- data.matrix(x[,-1][,41:60])
Z <- apply(Z,c(2,3),scale)
HOSVD <- hosvd(as.tensor(Z))
plot(HOSVD$U[[2]][,1],ylim=c(-0.5,0),type="h");abline(0,0,col=2,lty=2)
plot(HOSVD$U[[3]][,2],type="h");abline(0,0,col=2,lty=2)
order(-abs(HOSVD$Z@data[,1,2]))[1:10]
HOSVD$Z@data[,1,2][order(-abs(HOSVD$Z@data[,1,2]))[1:10]]
P <- pchisq(rowSums(scale(HOSVD$U[[1]][,1:2])^2),2,lower.tail=F)

table(p.adjust(P,"BH")<0.01)

table(p.adjust(P,"BH")<0.05)

table(p.adjust(P,"BH")<0.1)
