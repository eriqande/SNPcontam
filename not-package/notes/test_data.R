zeros <- t(matrix(c(40,160,50,80,100,20,30,15,140),nrow=3))
ones <- 200 - zeros
genos <- t(matrix(c(1,2,0,2,1,0,1,0,0,1,2,2,1,1,1),nrow=5))
data <- matrix(c(1,2,1,0,2,0,2,0,0,0,1,2,0,0,1),nrow=3)
clean_data <- P_likelihood(zeros,ones,genos,.5)
contam_data <- Pcontam(zeros,ones,genos,.5) 