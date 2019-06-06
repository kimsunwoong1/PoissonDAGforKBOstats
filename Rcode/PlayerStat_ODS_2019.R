# set your work directory for source code
# setwd("C:/Users/uos/Dropbox/Group study/2019OngoingProject/ParkHyewon/KBO_Stat_Relationship/rcode")
source("ODS_Algorithm_2019.R")

#####library#########################
library(igraph)
library(graph)
library(RBGL)
library(Rgraphviz)
library(corrplot)
library(xtable)
library(jmuOutlier)

#####Data ##############################

# set your work directory for data
# setwd("C:/Users/uos/Dropbox/Group study/2019OngoingProject/ParkHyewon/KBO_Stat_Relationship/KBOdata")
namelist = c("KimJaeHwan", "Romak", "LeeJeongHu", "ParkHaeMin", "Bernadina", "ChoiJeong")

##### 2016-2018 Data ##
data_all = read.csv(file = "3yearsData.csv")
data_all2 = data_all[,-(2:3)]
aggregate(data_all2, by = list(data_all2$name), mean)

data_pa = read.csv(file = "minPAData.csv", na.strings = "-")
apply(data_pa[,4:17], 2, mean)
apply(data_pa[,4:17], 2, sd)
##### Corr plot: 600 / 500 ##############################
for(name in namelist){
  data = read_data(name)
  corr = cor(data, method = "spearman")
  corrplot.mixed(corr, upper= "color",number.cex = 0.8 ,lower.col = "black", na.label = "0", na.label.col = "white")
}

##### Permutation spearman correlation test - It will take some time ####

corr_perm = function(name){
  data = read_data(name)
  pmat = matrix(nrow = ncol(data), ncol = ncol(data))
  for(i in 1:ncol(data)){
    j = i+1
    while(j <= ncol(data)){
      pmat[i,j] = pmat[j,i] = perm.cor.test(data[,c(i,j)], method = "spearman")$p.value
      j = j+1
    }
  }
  return(pmat)
}

cort_all = lapply(namelist,corr_perm)
names(cort_all) = namelist

for(name in namelist){
  data = read_data(name)
  corr = cor(data, method = "spearman")
  corp = cort_all[[which(name == namelist)]]
  corrplot.mixed(corr, upper= "color",number.cex = 0.8 ,lower.col = "black", 
                 na.label = "0", na.label.col = "white", p.mat = corp)
}

# Get positions & plot formatted p-values

for(name in namelist){
  data = read_data(name)
  corr = cor(data, method = "spearman")
  corp = cort_all[[which(name == namelist)]]
  
  corrplot(corr, type = "upper", p.mat = corp, method = "color",
           insig = "p-value", sig.level = -0.1, number.cex = 0.8, tl.pos = "d", na.label.col = "white")
  
  corrplot(corr, type = "lower", add = T, method = "number", tl.pos = "d", col = "black", 
           cl.pos = "n",  number.cex = 0.9, na.label.col = "white")
}

#### DAG for 6 players *ODS_Algorithm(20181005).R: nfolds = n - It will take some time####
par(mai = c(0,0,0,0))
for(name in namelist[1]){
  set.seed(7) # Some seeds may cause an error.
  ODSresult = Alg_fun(name)
  plot_fun(ODSresult$DAG)
}


####figure function ##########################
plot_fun = function(estimated_DAG, edgemode="directed"){
  estimated_DAG = t(estimated_DAG)
  am.graph<-new("graphAM", adjMat= estimated_DAG, edgemode=edgemode)
  plot(am.graph, attrs = list(node = list(fillcolor = "lightblue", fontsize = "10", height = "3", width = "3"), edge = list(arrowsize=1)))
}
?plot
#### function for reading data ####
read_data = function(name){
  data = read.csv(file = paste0(name, "Data.csv"), na.strings = "-")
  data = na.omit(data)
  data = data[, 5:17]
  colnames(data)[4:5] = c("2B", "3B")
  data = sapply(data, as.numeric)
  return(data)
}

####Algorithm function###############################
Alg_fun = function(name){
  data = read_data(name)
  id = which(apply(data, 2, function(x) max(table(x)/length(x))) > 0.99 )
  if( length(id)>0 )data = data[,-id]
  
  data_res = ODS_Algorithm(data, sparsity_level = 2 )
  
  colnames(data_res$DAG) = rownames(data_res$DAG)  = colnames(data)
  return(data_res)
}

