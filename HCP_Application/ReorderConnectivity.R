rm(list = ls())
require(R.matlab)
RV=read.csv("/Users/michellemiranda/Documents/CompositeHybrid_Toolbox/HCP_Application/Results/RV2.txt",header=FALSE)
library(gplots)
x=as.matrix(RV)
dis_con=dist(1-x)
aa <- hclust(dis_con,method = "complete", members = NULL)
x2=x[aa$order,aa$order]
write.table(x2, file = "/Users/michellemiranda/Documents/CompositeHybrid_Toolbox/HCP_Application/Results/ConnectivityOrdered.txt", append = FALSE, sep = " ", dec = ".",
            row.names = FALSE, col.names = FALSE)
ConnecOrder=aa$order
write.table(ConnecOrder, file = "/Users/michellemiranda/Documents/CompositeHybrid_Toolbox/HCP_Application/Results/ConnectivityROIorder.txt", append = FALSE, sep = " ", dec = ".",
            row.names = FALSE, col.names = FALSE)


