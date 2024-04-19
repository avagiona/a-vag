library("dplyr")
library("igraph")
library("NetHypGeom")
library("PSICQUIC")
library("stringr")
#set pathway
setwd("~/R-4.2.0")
#import hippie database
#hippie_v2.3<-read.delim("http://cbdm-01.zdv.uni-mainz.de/~mschaefer/hippie/hippie_current.txt", 
#                     sep="\t",header = FALSE, stringsAsFactors = FALSE)
#it takes a bit of time to import database online, I saved and the file, simplicity reasons. 
#saveRDS(hippie_v2.3,file="hippie_v2.3.RData")
hippie_v2.3<-readRDS("hippie_v2.3.RData") 
#calculate percentage of edges supported by more that one experiment
hippie_v2.3$V6 <- sub("experiments:", "", hippie_v2.3$V6)
hippie_v2.3$V6 <- sub("\\pmids.*", "", hippie_v2.3$V6)
hippie_v2.3$V6 <- sub("in vivo,", "", hippie_v2.3$V6)
hippie_v2.3$V6 <- sub("in vivo;", "", hippie_v2.3$V6)
hippie_v2.3$V6 <- sub("in vitro,", "", hippie_v2.3$V6)
hippie_v2.3$V6 <- sub("in vitro;", "", hippie_v2.3$V6)
hippie_v2.3$V6 <- sub(";", ",", hippie_v2.3$V6)
#0.69#
thres_0.69<-hippie_v2.3 %>% filter(V5 >= 0.69)
thres_0.69$number <- lengths(strsplit(thres_0.69$V6, ","))
exp_0.69<-sum(thres_0.69$number>1)
percentage_0.69<-exp_0.69/nrow(thres_0.69)
#0.70#
thres_0.70<-hippie_v2.3 %>% filter(V5 >= 0.70)
thres_0.70$number <- lengths(strsplit(thres_0.70$V6, ","))
exp_0.70<-sum(thres_0.70$number>1)
percentage_0.70<-exp_0.70/nrow(thres_0.70)
#0.71#
thres_0.71<-hippie_v2.3 %>% filter(V5 >= 0.71)
thres_0.71$number <- lengths(strsplit(thres_0.71$V6, ","))
exp_0.71<-sum(thres_0.71$number>1)
percentage_0.71<-exp_0.71/nrow(thres_0.71)
#0.72#
thres_0.72<-hippie_v2.3 %>% filter(V5 >= 0.72)
thres_0.72$number <- lengths(strsplit(thres_0.72$V6, ","))
exp_0.72<-sum(thres_0.72$number>1)
percentage_0.72<-exp_0.72/nrow(thres_0.72)
#0.73#
thres_0.73<-hippie_v2.3 %>% filter(V5 >= 0.73) 
thres_0.73$number <- lengths(strsplit(thres_0.73$V6, ","))
exp_0.73<-sum(thres_0.73$number>1)
percentage_0.73<-exp_0.73/nrow(thres_0.73)
#0.74#
thres_0.74<-hippie_v2.3 %>% filter(V5 >= 0.74) 
thres_0.74$number <- lengths(strsplit(thres_0.74$V6, ","))
exp_0.74<-sum(thres_0.74$number>1)
percentage_0.74<-exp_0.74/nrow(thres_0.74)

#edgelists for different confidence scores
edgelist_0.69<-thres_0.69 %>% dplyr::select (V1,V3)
edgelist_0.70<-thres_0.70 %>% dplyr::select (V1,V3)
edgelist_0.71<-thres_0.71 %>% dplyr::select (V1,V3)
edgelist_0.72<-thres_0.72 %>% dplyr::select (V1,V3)
edgelist_0.73<-thres_0.73 %>% dplyr::select (V1,V3)
edgelist_0.74<-thres_0.74 %>% dplyr::select (V1,V3)

#graphs from edgelists
net_0.69 <- graph_from_data_frame(edgelist_0.69, directed = FALSE)
net_0.70 <- graph_from_data_frame(edgelist_0.70, directed = FALSE)
net_0.71 <- graph_from_data_frame(edgelist_0.71, directed = FALSE)
net_0.72 <- graph_from_data_frame(edgelist_0.72, directed = FALSE)
net_0.73 <- graph_from_data_frame(edgelist_0.73, directed = FALSE)
net_0.74 <- graph_from_data_frame(edgelist_0.74, directed = FALSE)

#simplify to remove self and parallel interactions 
net_0.69<-simplify(net_0.69)
net_0.70<-simplify(net_0.70)
net_0.71<-simplify(net_0.71)
net_0.72<-simplify(net_0.72)
net_0.73<-simplify(net_0.73)
net_0.74<-simplify(net_0.74)

#parameter estimation netork 0.69
nodes_net_0.69<-data.frame(name = igraph::V(net_0.69)$name) #Nodes
edges_net_0.69<-igraph::E(net_0.69) #Edges
clusters_net_0.69<-decompose.graph(net_0.69)#get the clusters
giant_net_0.69<-clusters_net_0.69[[ which.max(sapply(clusters_net_0.69, vcount)) ]] #get the LCC
nodes_giant_net_0.69<-data.frame(name = igraph::V(giant_net_0.69)$name) #Nodes LCC
edges_giant_net_0.69<-igraph::E(giant_net_0.69) #Edges LCC

gma_0.69<-fit_power_law(degree(giant_net_0.69))$alpha #estimate_gma
c_0.69<-transitivity(giant_net_0.69,type = "localaverage") #estimate clustering coefficient
mean_degree_0.69<-mean(degree(giant_net_0.69))#estimate average node degree for artificial networks

#estimate Temperature
psmodels_0.69=list()
clustering_psmodels_0.69=list()
mean_clustering_models_0.69=NULL
for (i in 1:10){
  psmodels_0.69[[i]]<-ps_model(N = 15639, avg.k = 24.46, gma = 2.98, Temp = 0)
  clustering_psmodels_0.69[[i]]<-transitivity(psmodels_0.69[[i]][["network"]],type = "localaverage")
  mean_clustering_models_0.69<-(mean(clustering_psmodels_0.69[[i]]))
}

# Some maths by hand here, based on SF3 (LaBNE+HM paper) network 0.69
# c=a*T+b
# when T=0, c=0.77
# when T=1, c=0 so the equation is: c=-0.77*T + 0.77
# when c=0.13 (clustering coefficient of the LCC), T=0.83

#apply LaBNE+HM 
#coordinates_giant_net_0.69<-labne_hm(giant_net_0.69,gma = 2.98, Temp = 0.83, w = 2*pi)
#saveRDS(coordinates_giant_net_0.69,file="coordinates_0.69.RData")
#Coord_G_net_0.69<-readRDS("coordinates_0.69.RData")

#parameter estimation netork 0.70
nodes_net_0.70<-data.frame(name = igraph::V(net_0.70)$name) #Nodes
edges_net_0.70<-igraph::E(net_0.70) #Edges
clusters_net_0.70<-decompose.graph(net_0.70)#get the clusters
giant_net_0.70<-clusters_net_0.70[[ which.max(sapply(clusters_net_0.70, vcount)) ]] #get the LCC
nodes_giant_net_0.70<-data.frame(name = igraph::V(giant_net_0.70)$name) #Nodes LCC
edges_giant_net_0.70<-igraph::E(giant_net_0.70) #Edges LCC

gma_0.70<-fit_power_law(degree(giant_net_0.70))$alpha #estimate_gma
c_0.70<-transitivity(giant_net_0.70,type = "localaverage") #estimate clustering coefficient
mean_degree_0.70<-mean(degree(giant_net_0.70))#estimate average node degree for artificial networks

#estimate Temperature
psmodels_0.70=list()
clustering_psmodels_0.70=list()
mean_clustering_models_0.70=NULL
for (i in 1:10){
  psmodels_0.70[[i]]<-ps_model(N = 15633, avg.k = 24.43, gma = 2.98, Temp = 0)
  clustering_psmodels_0.70[[i]]<-transitivity(psmodels_0.70[[i]][["network"]],type = "localaverage")
  mean_clustering_models_0.70<-(mean(clustering_psmodels_0.70[[i]]))
}

# Some maths by hand here, based on SF3 (LaBNE+HM paper) network 0.70
# c=a*T+b
# when T=0, c=0.77
# when T=1, c=0 so the equation is: c=-0.77*T + 0.77
# when c=0.13 (clustering coefficient of the LCC), T=0.83

#apply LaBNE+HM 
#coordinates_giant_net_0.70<-labne_hm(giant_net_0.70,gma = 2.98, Temp = 0.83, w = 2*pi)
#saveRDS(coordinates_giant_net_0.70,file="coordinates_0.70.RData")
#Coord_G_net_0.70<-readRDS("coordinates_0.70.RData")

#parameter estimation netork 0.71
nodes_net_0.71<-data.frame(name = igraph::V(net_0.71)$name) #Nodes
edges_net_0.71<-igraph::E(net_0.71) #Edges
clusters_net_0.71<-decompose.graph(net_0.71)#get the clusters
giant_net_0.71<-clusters_net_0.71[[ which.max(sapply(clusters_net_0.71, vcount)) ]] #get the LCC
nodes_giant_net_0.71<-data.frame(name = igraph::V(giant_net_0.71)$name) #Nodes LCC
edges_giant_net_0.71<-igraph::E(giant_net_0.71) #Edges LCC
egde<-as.data.frame(get.edgelist(giant_net_0.71))# get the edgelist as dataframe as input for mercator

gma_0.71<-fit_power_law(degree(giant_net_0.71))$alpha #estimate_gma
c_0.71<-transitivity(giant_net_0.71,type = "localaverage") #estimate clustering coefficient
mean_degree_0.71<-mean(degree(giant_net_0.71))#estimate average node degree for artificial networks

#estimate Temperature
psmodels_0.71=list()
clustering_psmodels_0.71=list()
mean_clustering_models_0.71=NULL
for (i in 1:10){
  psmodels_0.71[[i]]<-ps_model(N = 15588, avg.k = 23.94, gma = 2.97, Temp = 0)
  clustering_psmodels_0.71[[i]]<-transitivity(psmodels_0.71[[i]][["network"]],type = "localaverage")
  mean_clustering_models_0.71<-(mean(clustering_psmodels_0.71[[i]]))
}

# Some maths by hand here, based on SF3 (LaBNE+HM paper) network 0.71
# c=a*T+b
# when T=0, c=0.77
# when T=1, c=0 so the equation is: c=-0.77*T + 0.77
# when c=0.13 (clustering coefficient of the LCC), T=0.83

#apply LaBNE+HM 
#coordinates_giant_net_0.71<-labne_hm(giant_net_0.71,gma = 2.97, Temp = 0.83, w = 2*pi)
#saveRDS(coordinates_giant_net_0.71,file="coordinates_0.71.RData")
Coord_G_net_0.71<-readRDS("coordinates_0.71.RData")

#parameter estimation netork 0.72
nodes_net_0.72<-data.frame(name = igraph::V(net_0.72)$name) #Nodes
edges_net_0.72<-igraph::E(net_0.72) #Edges
clusters_net_0.72<-decompose.graph(net_0.72)#get the clusters
giant_net_0.72<-clusters_net_0.72[[ which.max(sapply(clusters_net_0.72, vcount)) ]] #get the LCC
nodes_giant_net_0.72<-data.frame(name = igraph::V(giant_net_0.72)$name) #Nodes LCC
edges_giant_net_0.72<-igraph::E(giant_net_0.72) #Edges LCC

gma_0.72<-fit_power_law(degree(giant_net_0.72))$alpha #estimate_gma
c_0.72<-transitivity(giant_net_0.72,type = "localaverage") #estimate clustering coefficient
mean_degree_0.72<-mean(degree(giant_net_0.72))#estimate average node degree for artificial networks

#estimate Temperature
psmodels_0.72=list()
clustering_psmodels_0.72=list()
mean_clustering_models_0.72=NULL
for (i in 1:10){
  psmodels_0.72[[i]]<-ps_model(N = 15585, avg.k = 23.91, gma = 2.97, Temp = 0)
  clustering_psmodels_0.72[[i]]<-transitivity(psmodels_0.72[[i]][["network"]],type = "localaverage")
  mean_clustering_models_0.72<-(mean(clustering_psmodels_0.72[[i]]))
}

# Some maths by hand here, based on SF3 (LaBNE+HM paper) network 0.72
# c=a*T+b
# when T=0, c=0.77
# when T=1, c=0 so the equation is: c=-0.77*T + 0.77
# when c=0.127 (clustering coefficient of the LCC), T=0.83


#apply LaBNE+HM 
#coordinates_giant_net_0.72<-labne_hm(giant_net_0.72,gma = 2.97, Temp = 0.83, w = 2*pi)
#saveRDS(coordinates_giant_net_0.72,file="coordinates_0.72.RData")
#Coord_G_net_0.72<-readRDS("coordinates_0.72.RData")

#parameter estimation netςork 0.73
nodes_net_0.73<-data.frame(name = igraph::V(net_0.73)$name) #Nodes
edges_net_0.73<-igraph::E(net_0.73) #Edges
clusters_net_0.73<-decompose.graph(net_0.73)#get the clusters
giant_net_0.73<-clusters_net_0.73[[ which.max(sapply(clusters_net_0.73, vcount)) ]] #get the LCC
nodes_giant_net_0.73<-data.frame(name = igraph::V(giant_net_0.73)$name) #Nodes LCC
edges_giant_net_0.73<-igraph::E(giant_net_0.73) #Edges LCC

gma_0.73<-fit_power_law(degree(giant_net_0.73))$alpha #estimate_gma
c_0.73<-transitivity(giant_net_0.73,type = "localaverage") #estimate clustering coefficient
mean_degree_0.73<-mean(degree(giant_net_0.73))#estimate average node degree for artificial networks

#estimate Temperature
psmodels_0.73=list()
clustering_psmodels_0.73=list()
mean_clustering_models_0.73=NULL
for (i in 1:10){
  psmodels_0.73[[i]]<-ps_model(N = 14367, avg.k = 17.36, gma = 3, Temp = 0)
  clustering_psmodels_0.73[[i]]<-transitivity(psmodels_0.73[[i]][["network"]],type = "localaverage")
  mean_clustering_models_0.73<-(mean(clustering_psmodels_0.73[[i]]))
}

# Some maths by hand here, based on SF3 (LaBNE+HM paper) network 0.73
# c=a*T+b
# when T=0, c=0.76
# when T=1, c=0 so the equation is: c=-0.76*T + 0.76
# when c=0.15 (clustering coefficient of the LCC), T=0.80


#apply LaBNE+HM 
#coordinates_giant_net_0.73<-labne_hm(giant_net_0.73,gma = 3, Temp = 0.80, w = 2*pi)
#saveRDS(coordinates_giant_net_0.73,file="coordinates_0.73.RData")
#Coord_G_net_0.73<-readRDS("coordinates_0.73.RData")

#parameter estimation netork 0.74
nodes_net_0.74<-data.frame(name = igraph::V(net_0.74)$name) #Nodes
edges_net_0.74<-igraph::E(net_0.74) #Edges
clusters_net_0.74<-decompose.graph(net_0.74)#get the clusters
giant_net_0.74<-clusters_net_0.74[[ which.max(sapply(clusters_net_0.74, vcount)) ]] #get the LCC
nodes_giant_net_0.74<-data.frame(name = igraph::V(giant_net_0.74)$name) #Nodes LCC
edges_giant_net_0.74<-igraph::E(giant_net_0.74) #Edges LCC

gma_0.74<-fit_power_law(degree(giant_net_0.74))$alpha #estimate_gma
c_0.74<-transitivity(giant_net_0.74,type = "localaverage") #estimate clustering coefficient
mean_degree_0.74<-mean(degree(giant_net_0.74))#estimate average node degree for artificial networks

#estimate Temperature
psmodels_0.74=list()
clustering_psmodels_0.74=list()
mean_clustering_models_0.74=NULL
for (i in 1:10){
  psmodels_0.74[[i]]<-ps_model(N = 13729, avg.k = 15.01, gma = 3, Temp = 0)
  clustering_psmodels_0.74[[i]]<-transitivity(psmodels_0.74[[i]][["network"]],type = "localaverage")
  mean_clustering_models_0.74<-(mean(clustering_psmodels_0.74[[i]]))
}

# Some maths by hand here, based on SF3 (LaBNE+HM paper) network 0.74
# c=a*T+b
# when T=0, c=0.73
# when T=1, c=0 so the equation is: c=-0.76*T + 0.76
# when c=0.15 (clustering coefficient of the LCC), T=0.80

#apply LaBNE+HM 
#coordinates_giant_net_0.74<-labne_hm(giant_net_0.74,gma = 3, Temp = 0.80, w = 2*pi)
#saveRDS(coordinates_giant_net_0.74,file="coordinates_0.74.RData")
#Coord_G_net_0.74<-readRDS("coordinates_0.74.RData")




















