# load the EPIC whole blood cell type count library
library(FlowSorted.Blood.EPIC)

data(IDOLOptimizedCpGs)
data(IDOLOptimizedCpGs.compTable)

# the beta matrix is Noob normalized
cellFrac <- projectCellType_CP(getBeta(preprocessNoob(rgset))[IDOLOptimizedCpGs,],
                               IDOLOptimizedCpGs.compTable,lessThanOne = TRUE)

# correct the beta matrix based on the cell count ( from ChAMP )

# The idea is to apply those cell fractions (except the smallest ones) as covariates into regression model, and then get the residual value.
# In theory, these residual contains information of other factors ( like cancer/normal, etc) in original matrix apart from these cell fractions.
# Basically, we get the corrected beta matrix by adding the mean beta values of each probe to the regression residual values.
# The residual values are from linear model that was built on five largest cell fractions.
lm.o1 <- lm(t(b1) ~ cellFrac[,-1*which.min(colMeans(cellFrac))])
tmp.m1 <- t(lm.o1$res)+rowMeans(b1); 
tmp.m1[tmp.m1 <= 0] <- min(tmp.m1[which(tmp.m1 > 0)])
tmp.m1[tmp.m1 >= 1] <- max(tmp.m1[which(tmp.m1 < 1)])

# name the corrected beta matrix as myCell1
myCell1 <- tmp.m1