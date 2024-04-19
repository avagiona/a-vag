#REQUIREMENTS
library("NetHypGeom") 
library("igraph")
library("dplyr")
library("ggplot2")
library("ggpubr")
library("enrichR")
library("lattice") 
library("RColorBrewer")
library("tidyr")
library("writexl")
library("enrichR")
library("plyr")
library("tibble")
library("ggrepel")
library("org.Hs.eg.db")
#set pathway
setwd("~/R/R-4.1.0")

#import hippie database version 2.2
edgelist<-read.delim("hippie_v2.2.txt",sep="\t",header = TRUE,stringsAsFactors = FALSE)

#subset 0.71 network
thres_0.71<-edgelist %>% 
            dplyr::filter(V5 >= 0.71)%>%
            dplyr::select (V1,V3) #0.71 edgelist
net_0.71 <- simplify(graph_from_data_frame(thres_0.71, directed = FALSE)) #graph from edgelist
clusters_net_0.71<-decompose.graph(net_0.71) #get the clusters
giant_net_0.71<-clusters_net_0.71[[ which.max(sapply(clusters_net_0.71, vcount)) ]] #get the LCC
nodes_giant_net_0.71<-data.frame(name = igraph::V(giant_net_0.71)$name) #Nodes LCC
edge<-as.data.frame(get.edgelist(giant_net_0.71)) #get the edgelist of LCC
edge<-read.csv("egdelist_net_0.71.csv",fileEncoding = "UTF-8-BOM", header=TRUE)

#estimate parameters
gma<-fit_power_law(degree(giant_net_0.71))$alpha #estimate_gma
c<-transitivity(giant_net_0.71,type = "localaverage") #estimate clustering coefficient
avg_k<-mean(degree(giant_net_0.71))#estimate average node degree for artificial networks

#construct 10 artificial networks
#estimate Temperature
psmodels=list()
clustering_psmodels=list()
mean_clustering_models=NULL
for (i in 1:10){
  psmodels[[i]]<-ps_model(N = 13077, avg.k = 14.24, gma = 2.74, Temp = 0)
  clustering_psmodels[[i]]<-transitivity(psmodels[[i]][["network"]],type = "localaverage")
  mean_clustering_models<-(mean(clustering_psmodels[[i]]))
  return(mean_clustering_models)
}
# Some maths by hand here, based on SF3 (LaBNE+HM paper)
# when T=0, c=0.77
# when T=1, c=0 so the equation is: c=-0.77*T + 0.77
# when c=0.15 (clustering coefficient of the LCC), T=0.80

#apply LaBNE+HM and embed
coordinates_giant_net_0.71<-labne_hm(giant_net_0.71,gma = 2.74, Temp = 0.80, w = 2*pi)
saveRDS(coordinates_giant_net_0.71,file="coordinates_0.71.RData")
coor_0.71<-readRDS("coordinates_0.71.RData")

#evaluate the embedding to the hyperbolic disc
N <- vcount(giant_net_0.71)
beta <- 1 / (gma - 1) # Parameter controlling popularity fading
m <- round(avg_k/2) # Parameter controlling average node degree
Temp=0.80

# Connection probability
conn <- get_conn_probs(giant_net_0.71, coor_0.71$polar, bins = 20)
theo <- get_theoretical_conn_probs(conn$dist, N, avg_k, gma, Temp)
conn <- rbind(conn, theo)
conn$label <- rep(c("LaBNE+HM", "Theory"), each = 20)

p_conn <- ggplot(conn, aes(dist, prob, colour = label, shape = label)) + 
  geom_point(size = 2) + geom_line() +
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", 
                                              scales::math_format(10^.x))) + 
  annotation_logticks(sides = "l") +
  scale_colour_manual(values = c("#339e2b", "#e3196a")) + 
  labs(x = "Hyperbolic distance", y = "Connection probability") + 
  theme_bw() + theme(legend.title = element_blank(), 
                     legend.background = element_blank(), 
                     legend.justification = c(0, 0), legend.position = c(0, 0))
p_conn

# Real degrees vs expected degrees
degs <- tibble(k = degree(giant_net_0.71), exp_k = numeric(N))

R <- 2*log(N) - 
  2*log((2*Temp*(1 - exp(-0.5*(1 - beta)*2*log(N))))/(sin(Temp*pi)*m*(1 - beta)))

for(i in 1:N){
  # Compute the hyperbolic distance between a node and all others
  d <- hyperbolic_dist(coor_0.71$polar[i, ], coor_0.71$polar)
  # Compute the probability that the node will connect with all others
  prob <- 1 / (1 + exp((d - R)/(2*Temp)))
  # Compute the expected node degree
  degs$exp_k[i] <- sum(prob)
}

p_deg <- ggplot(degs, aes(k, round(exp_k))) + geom_point(size = 0.2) +
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x), 
                labels = scales::trans_format("log10", scales::math_format())) +
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x), 
                labels = scales::trans_format("log10", scales::math_format())) +
  annotation_logticks() + labs(x = "Node degree", y = "Expected node degree") + 
  theme_bw()
p_deg

# Clustering
epochs <- 10
cc <- transitivity(giant_net_0.71, "average")
theo_cc <- numeric(epochs)

for(i in 1:epochs){
  ps_net <- ps_model(N = N, avg.k = avg_k, gma = gma, Temp = Temp)
  theo_cc[i] <- transitivity(ps_net$network, "average")
}

saveRDS(ps_net,file="ps_net_clustering.RData")
ps_net<-readRDS("ps_net_clustering.RData")
theo_cc<-readRDS("theo_cc.RData")

clust <- tibble(label = c("Real", "Theory"), cc = c(cc, mean(theo_cc)),
                err = c(0, sd(theo_cc)))
dodge <- position_dodge(width = 0.9)
p_cc <- ggplot(clust, aes(label, cc)) + geom_col(position = dodge, width = 0.5) +
  geom_errorbar(aes(ymin = cc - err, ymax = cc + err), 
                position = dodge, width = 0.25) + 
  labs(x = "", y = "Clustering coefficient") + theme_bw()
p_cc


SF_1<- ggarrange(p_conn, p_deg, p_cc,nrow = 2, ncol = 2,
                        labels = c("a", "b", "c"))
SF_1

#import nodes of hPIN with entrez ID
nodes_hPIN<-read.csv("nodes_mapped.csv", fileEncoding = "UTF-8-BOM")

#function to calculate differences
clusters <- function(df){
  df<-df[order(df$theta),]
  differences<-diff(df$theta, lag=1)
  df$differences<-append(differences,0)
  df<-df[order((df$differences),decreasing = TRUE),]
  return(df)
}
hPIN_diff<-clusters(nodes_hPIN)

#plot_differences_g_gap_htt_net
SF_2 <- ggplot(hPIN_diff, aes(x=theta,y=differences)) +
  geom_bar(stat = "identity",width = 0.01)+
  geom_hline(yintercept=0.011344,color="red")+
  annotate("text",x=0.1,y=0.0117,color="red",size=3,label=c('g=0.011344'))+
  theme_bw()+ylab("Difference between consecutive angles")+
  xlab("Angular dimension θ")
SF_2

#create list with the clusters to enrich
split_tibble <- function(tibble, col = 'col') tibble %>% split(., .[, col])
final_nodes_hPIN<-split_tibble(nodes_hPIN, "cluster")
#load the "GO_Biological_Process_2021" 
dbs_GO_BP <- c("GO_Biological_Process_2021" )

#enrichment analysis for the different clusters 
enrichment_GO_BP_hPIN<-list()
for (i in 1:length(final_nodes_hPIN)) {
  enrichment_analysis_GO_BP_hPIN<-enrichr(c(final_nodes_hPIN[[i]]$Symbol),dbs_GO_BP)
  enrichment_GO_BP_hPIN[[paste0("cluster", i)]] <- enrichment_analysis_GO_BP_hPIN
}

#plot the network hPIN with the interactions
rainbowcol<-rainbow(12)

F_1<-ggplot(nodes_hPIN,aes(x=x,y=y))+
  scale_colour_manual(values = rainbowcol,labels=c(labels_hPIN$term))+
  theme_classic()+
  theme(axis.title.y = element_blank(), axis.title.x = element_blank(), 
        axis.text.x = element_blank(), axis.text.y = element_blank(), 
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank()) + 
  guides(colour = guide_legend(ncol=2))+
  geom_segment(data=edge_in,aes(x=x_a,y=y_a,xend=x,yend=y),color='grey90',size=0.2,alpha=0.7)+
  geom_point(aes(colour = as.factor(cluster)),size=1)+
  theme(legend.position="bottom",legend.text=element_text(size=12),
        legend.title = element_blank(),panel.background = element_blank(),
        strip.background = element_rect(colour=NA, fill=NA),
        panel.border = element_rect(fill = NA, color = "black"))
F_1


#Htt interactors 
htt_interactors<-read.csv("htt_dataset.csv",fileEncoding = "UTF-8-BOM")
#DF with interactors, coordinates, paralogs
htt_dataset<-merge(nodes_net,htt_interactors,by.x='id',by.y='values_V3')
htt_diff<-clusters(htt_dataset)

#plot_differences_g_gap_htt_net
SF_3 <- ggplot(htt_diff, aes(x=theta,y=differences)) +
  geom_bar(stat = "identity",width = 0.05)+
  geom_hline(yintercept=0.059198411,color="red")+
  annotate("text",x=0.025,y=0.062,color="red",size=3,label=c('g=0.059198'))+
  theme_bw()+ylab("Difference between consecutive angles")+
  xlab("Angular dimension θ")
SF_3

colors<-c("1" = "#66C2A5", "2" = "#CF948C","3" = "#D58EC4","4" = "#B7D84C",
          "5" = "#EFCC6B","6" = "#FF62BC","7" = "coral1","8" = "steelblue1",
          "9" = "yellow2","10" = "hotpink2","11" = "palegreen1","12" = "maroon",
          "13" = "azure4","14"="red", "15"="blue","16"="yellow", "17"="pink")

#plot network of 17 clusters after enrichment analysis 
F_2<-ggplot(htt_dataset,aes(x=x,y=y))+
  geom_point(aes(colour = as.factor(cluster)),size=1.2)+
  scale_colour_manual(values = colors,labels=c(labels$`X[[i]]`))+
  theme_classic()+
  theme(axis.title.y = element_blank(), axis.title.x = element_blank(), 
        axis.text.x = element_blank(), axis.text.y = element_blank(), 
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank()) + 
  geom_text_repel(aes(label=SYMBOL),box.padding   = 0.2,
                  point.padding = 1,segment.color = 'black',hjust=0.6,vjust=0.5,
                  segment.size = 0.1,size=2, color='black',force=1,max.overlaps = 2000)+
  guides(colour = guide_legend(ncol=2))+
  theme(legend.position="bottom",panel.background = element_blank(),
        strip.background = element_rect(colour=NA, fill=NA),
        panel.border = element_rect(fill = NA, color = "black"),
        legend.text=element_text(size=7.5),legend.title = element_blank())
F_2

