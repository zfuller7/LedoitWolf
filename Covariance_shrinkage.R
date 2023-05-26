#Ledoit Wolf Covariance Matrix Shrinkage Example
#Zach Fuller
#STAT 616
#Run on R version 4.2.2

#First import excel file with pct change data
library("readxl")
library("psych")
data <- read_excel("filepath") #File should be stock return data. I used a file with returns from 5 stocks.
#If you want to use a different number of stocks, then some edits will be necessary.

#Create Covariance Matrix and find standard deviations
cm <- cov(data)
c <- c(cm[1,1], cm[2,2], cm[3,3], cm[4,4], cm[5,5])
s <- sqrt(c)

#Correlation Matrix F 
Fi <- cor(data)
F2 <- as.vector(Fi)
a <- mean(F2)

# Create a matrix of the same size as Fi filled with the mean of F2
mean_matrix <- matrix(data = a, nrow = nrow(Fi), ncol = ncol(Fi))

# Create a matrix where each row is s multiplied by s[i]
s_matrix <- s %*% t(s)

# Element-wise multiplication of mean_matrix and s_matrix
constantcorr_matrix <- mean_matrix * s_matrix

# Set the diagonal elements to the diagonal elements of Fi
diag(constantcorr_matrix) <- diag(Fi)

# Convert the matrix to a vector
constantcorr <- as.vector(constantcorr_matrix)

TargetMatrix <- matrix(data=constantcorr, nrow=5, byrow=TRUE)


#The next steps lead up to the shrinkage constant, delta.
#First find pihat, rhohat, and gammahat
pihat <- sum(var(sqrt(length(data$ChangeSPY))*cm))

rhohat1 <- tr(var(sqrt(length(data$ChangeSPY))*cm))
rhohat2 <- rhohat1 #I'm still working on this part, so I put a placeholder.
rhohat <- rhohat1+rhohat2

gammahat <- sum((F2-as.vector(cm))^2)

#Now, delta
delta <- (pihat-rhohat)/gammahat

#Finally, the revised covariance matrix, covshrink
covshrink <- delta*Fi+(1-delta)*cm


#Compare results to another shrinkage estimator
install.packages("cvCovEst")
install.packages("rbibutils")
library(cvCovEst)

l <- linearShrinkLWEst(data)

matlabmatrix <- read_excel("C:\\Users\\zachary.fuller\\Desktop\\Matlabmatrix.xlsx")
matlabmatrix <- as.matrix(matlabmatrix)
#Compute Frobenius Norm of difference between methods
norm((l-matlabmatrix), type="F")
