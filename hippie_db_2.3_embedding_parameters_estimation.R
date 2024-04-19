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

#0.71#
thres_0.71<-hippie_v2.3 %>% filter(V5 >= 0.71)
thres_0.71$number <- lengths(strsplit(thres_0.71$V6, ","))
exp_0.71<-sum(thres_0.71$number>1)
percentage_0.71<-exp_0.71/nrow(thres_0.71)

#edgelists for different confidence scores
edgelist_0.71<-thres_0.71 %>% dplyr::select (V1,V3)

#graphs from edgelists
net_0.71 <- graph_from_data_frame(edgelist_0.71, directed = FALSE)

#simplify to remove self and parallel interactions 
net_0.71<-simplify(net_0.71)

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

#apply LaBNE+HM and embed
coordinates_giant_net_0.71<-labne_hm(giant_net_0.71,gma = 2.97, Temp = 0.83, w = 2*pi)
saveRDS(coordinates_giant_net_0.71,file="coordinates_0.71.RData")

#import the file 
Coord_G_net_0.71<-readRDS("coordinates_0.71.RData")
















