#New_data
require(rTensor)
#BDC
x <- read.csv("../../data/cancer_data/cancer_expression_data/BDC/bdc_expression.csv",sep=",")
Z <- array(NA,c(dim(x)[1],7,2))
Z[,,1] <- data.matrix(x[,-1][,1:7])
Z[,,2] <- data.matrix(x[,-1][,8:14])
Z <- apply(Z,c(2,3),scale)
HOSVD <- hosvd(as.tensor(Z))
plot(HOSVD$U[[2]][,1],ylim=c(-0.5,0),type="h");abline(0,0,col=2,lty=2)
plot(HOSVD$U[[3]][,2],type="h");abline(0,0,col=2,lty=2)
order(-abs(HOSVD$Z@data[,1,2]))[1:10]
# [1]  2  1  3  4  8  6  9 10  7 14
HOSVD$Z@data[,1,2][order(-abs(HOSVD$Z@data[,1,2]))[1:10]]
# [1] 165.5347435  10.7837727  -4.2657159   3.2511170  -1.4915595  -0.8796976
# [7]  -0.7844165  -0.7237482   0.7067796   0.6982120
P <- pchisq(scale(HOSVD$U[[1]][,2])^2,1,lower.tail=F)
table(p.adjust(P,"BH")<0.01)

#FALSE  TRUE 
#17417   199 

table(p.adjust(P,"BH")<0.05)

#FALSE  TRUE 
#17306   310 

table(p.adjust(P,"BH")<0.1)

#FALSE  TRUE 
#17238   378  
