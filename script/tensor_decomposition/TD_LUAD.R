#LUAD
require(rTensor)
x <- read.csv("../../data/cancer_data/cancer_expression_data/LUAD/expression.csv",sep="\t")
Z <- array(NA,c(dim(x)[1],27,2))
Z[,,1] <- data.matrix(x[,-1][,seq(1,54,by=2)])
Z[,,2] <- data.matrix(x[,-1][,seq(2,54,by=2)])
Z <- apply(Z,c(2,3),scale)
HOSVD <- hosvd(as.tensor(Z))
plot(HOSVD$U[[2]][,1],ylim=c(-0.5,0),type="h");abline(0,0,col=2,lty=2)
plot(HOSVD$U[[3]][,2],type="h");abline(0,0,col=2,lty=2)
order(-abs(HOSVD$Z@data[,1,2]))[1:10]
HOSVD$Z@data[,1,2][order(-abs(HOSVD$Z@data[,1,2]))[1:10]]
P <- pchisq(scale(HOSVD$U[[1]][,2])^2,1,lower.tail=F)
table(p.adjust(P,"BH")<0.01)

table(p.adjust(P,"BH")<0.05)

table(p.adjust(P,"BH")<0.1)
