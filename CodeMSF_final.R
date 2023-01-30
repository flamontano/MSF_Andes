#######################
# MSF South America
########################

## 1.  LOADING PACKAGES ## -----

rm(list=ls())

# required packages -------
library(asnipe)
library(igraph)
library(survival)
library(tnet)
library(ggplot2)
library(dplyr)
library(tidyr)
library(raster)
library(car)
library(tidyverse)
library(furrr)
library(visreg)
library(betareg)
library(assortnet)


### LOADING DATA ------------

flocks.df<-read.csv("data_MS/Andeanflocks.csv", header = T)
row.names(flocks.df)<-flocks.df[,1]

#predictors per network
environ<-aggregate(cbind(latitude, elev, Forest50, hfp_venter) ~ network_id, flocks.df, mean )

### creating matrices for networks -----------------
mats<-list()
out <- split(flocks.df, f = flocks.df$network)
names(out) 

for (i in 1:length(out)){
  fx<-as.data.frame(out[[i]][,2:649])# 
  fx[] <- lapply(fx, function(x) as.numeric(as.character(x)))
  fx<-fx[,which(colSums(fx)!=0)]#eliminating species that do not occur 
  rownames(fx)<-out[[i]][,1]
  fx<-fx[which(rowSums(fx)>1),]
  mats[[i]]<-t(fx)
}
names(mats)<-names(out)


### OBSERVED NETWORKS, NETWORK PERMUTATIONS AND NETWORK COVARIANCE ------------------------

permutations = 3500

networks<-list()# observed networks
network.cv<-data.frame(obs_cov= rep(NA, length(mats)), 
                       RandomQt2.5 = rep(NA, length(mats)), 
                       RandomQ97.5= rep(NA, length(mats)), 
                       P_cov=rep(NA, length(mats))) # network CV
networks.random<-list()#randomized networks per site, a list of arrays
random.cv<-list()# network CV for randomized networks

for (i in 1:length(mats)){
  bgs<-t(mats[[i]])
  network<-asnipe::get_network(bgs)
  observed.cv<-raster::cv(network)
  network.random <- network_permutation(bgs, permutations=permutations, association_matrix=network)
  
  random.cor <- rep(NA, permutations)
 
  for (j in 1:permutations) {
    random.cor[j] <- cv(network.random[j,,])
    }
  network.cv[i,1]<-observed.cv
  network.cv[i,2:3]<-quantile(random.cor,c(0.025,0.975))
  random.cv[[i]]<-random.cor
  networks[[i]]<-network
  networks.random[[i]]<-network.random
  network.cv[i,4]<-sum(observed.cv < random.cor)/permutations
}

rownames(network.cv)<-names(mats)

network.cv # If observed CV is larger than 95% of CVrandom, the observed network contains more preferred/avoided relationships than expected 
            # P val = number of times the value from the random networks was larger than the value from the observed network, divided by the number of permutations (Farine, 2017). 


### NETWORK LEVEL METRICS ------------------------- 

#function for connectance
connec.fun=function(x)(sum(x)/((length(x)^2)))

#function for all calculations

get_network_property = function(network){
  gs = igraph::graph_from_adjacency_matrix (network, "undirected" , weighted = T )
  cluster = tnet::clustering_w(network)
  com = igraph::cluster_fast_greedy (gs) # detecting communities with Clauset et al. (2004) algorithm
  modularity = igraph::modularity(com)# Q index, based on community assignments 
  degrees= igraph::degree(gs, normalized = T) # number of connections per node divided by the number of nodes n minus 1
  av.degree=mean(degrees)
  connectance = connec.fun(degrees) #
  tibble::tibble(cluster = cluster, modularity = modularity, av.degree= av.degree,
                 connectance = connectance)
}

n_perm = permutations 

all_metrics = vector("list", length = length(networks)) # 12:03 pm
for (i in 1:length(networks)){
  cat("i = ", i, "\n")
  
  #observed values
  obs = get_network_property(network = networks[[i]])
  obs$type = "obs"
  
  #random values
  network.random<-networks.random[[i]]
  
  future::plan("multicore", workers = 40)
  
  rand = furrr::future_map_dfr(1:n_perm, .f = function(j){
    get_network_property(network.random[j, ,])
  }, .progress = T )
  
  head(rand)
  rand = mutate(rand, type = "rand")
  rand = as_tibble(rand)
  
  all_i = bind_rows(obs, rand)
  all_metrics[[i]] = all_i
  # if(i %% 10 == 0) saveRDS(all_metrics, file = "file.rds")
}

names(all_metrics)<-names(mats)

# comparing obs and randomized values to see significance of each observed metric

get_ci = function(network_met){
  rand_ci = filter(network_met, type == "rand") %>% 
    dplyr::select(-type) %>% 
    pivot_longer(cols = everything()) %>% 
    group_by(name) %>% 
    summarise(ci_low = quantile(value, 0.025),
              ci_high = quantile(value, 0.975))
  obs_ci = filter(network_met, type == "obs") %>% 
    dplyr::select(-type) %>% 
    pivot_longer(cols = everything(), values_to = "obs")
  left_join(obs_ci, rand_ci)
  
}

all_CI = lapply(all_metrics, get_ci)
network_level = bind_rows(all_CI, .id = "network_id")

#significance (column sig) is binary: 1 = NOT significant (value is between the 95%CI), 0  = significant
network_level$sig<-NA
for (i in 1: nrow(network_level)){
  network_level$sig[i] = between(network_level[i,3], network_level[i,4], network_level[i,5])*1
}

#dataframe for regressions --

metrics = as.data.frame(network_level %>% 
  dplyr::select(-ci_low, -ci_high)  %>% 
  pivot_wider(names_from = name, values_from = c(obs, sig)))
rownames(metrics)<-metrics[,1]

metrics<- merge(metrics, network.cv[,c(1,4)], by=0); metrics = metrics[,-1]
  richness=lapply(mats,colSums)#observed richness
metrics$av.richness=sapply(richness,mean)# mean per flock richness
metrics$richperflock.min=t(sapply(richness,range))[,1]# max per flock richness
metrics$richperflock.max=t(sapply(richness,range))[,2]# min per flock richness
metrics$network_rich = sapply(mats, dim)[1,]
metrics$Nflocks = sapply(mats, dim)[2,]

# correlation between SR per network and number of flocks
cor(metrics$Nflocks, metrics$network_rich)#0.346
    
## DATA FRAMES FROM RANDOMIZED NETWORKS----------------------

cluster_data = function(network_met){
  rand_data = filter(network_met, type == "rand") %>% 
    dplyr::select(-type)%>% dplyr::select(cluster)
}

mod_data = function(network_met){
  rand_data = filter(network_met, type == "rand") %>% 
    dplyr::select(-type)%>% dplyr::select(modularity)
}
avdeg_data = function(network_met){
  rand_data = filter(network_met, type == "rand") %>% 
    dplyr::select(-type)%>% dplyr::select(av.degree)
}

conn_data = function(network_met){
  rand_data = filter(network_met, type == "rand") %>% 
    dplyr::select(-type)%>% dplyr::select(connectance)
}

cluster_r = lapply(all_metrics, cluster_data)
mod_r = lapply(all_metrics, mod_data)
avdeg_r = lapply(all_metrics, avdeg_data)
conn_r = lapply(all_metrics, conn_data)

cluster.df= do.call("cbind", cluster_r); names(cluster.df) = names(mats)
mod.df= do.call("cbind", mod_r); names(mod.df) = names(mats)
avdeg.df= do.call("cbind", avdeg_r); names(avdeg.df) = names(mats)
conn.df = do.call("cbind", conn_r); names(conn.df) = names(mats)
conn.df<-log(conn.df) #log of connectance for regressions

random.networks.df<-list(cluster.df, mod.df, avdeg.df, conn.df)
names(random.networks.df)<-c("cluster", "modul", "avdeg", "connec")


#### calculating r_com for modularity

#Function to calculate r_c, with default number of bootstraps = 100, needs sites in rows, spp in columns

t_Mats <- lapply(mats,function(x){
  t(x)})

#function modified from Shizuka & Farine
calc_rc=function(data, n.bootstraps=100, plot.result=F){
  network.community <- matrix(0,ncol(data),ncol(data))
  network.present <- matrix(0,ncol(data),ncol(data))
  network <- get_network(data,data_format="GBI", association_index="SRI")
  community.observed <- fastgreedy.community(graph.adjacency(network,mode="undirected",weighted=TRUE))
  
  for (i in 1:n.bootstraps) {
    gbi.boot <- data[sample(1:nrow(data),nrow(data),replace=TRUE),]
    network.boot <- get_network(gbi.boot,data_format="GBI", association_index="SRI")
    community.boot <- fastgreedy.community(graph.adjacency(network.boot,mode="undirected",weighted=TRUE))
    network.community <- network.community + outer(community.boot$membership, community.boot$membership,"==")
    network.present <- network.present + outer((rowSums(network.boot)>0),(rowSums(network.boot)>0),"*")
  }
  P <- network.community/network.present
  P[!is.finite(P)] <- 0
  rc <- assortment.discrete(P,community.observed$membership)$r
  return(rc)
}
#end function from Shizuka & Farine

r_com_all = lapply(t_Mats, calc_rc)
r_com_df = as.data.frame(t(bind_rows(r_com_all, .id = "network_id"))); names(r_com_df)<-"r_com"


### REGRESSIONS ###

network.df<-left_join(metrics, environ, by="network_id");network.df<-network.df[,-1]
network.df$latitude<-abs(network.df$latitude)

network.df$obs_connectance<-log(network.df$obs_connectance)# log for observed connectance

# input data for analyses - modularity and connectace

inp <- network.df[,c(11, 9, 1:4)]# 
inp$r_com<-r_com_df$r_com
inp$Forest <- as.numeric(scale(network.df$Forest50)) #
inp$HFP <-as.numeric(scale(network.df$hfp_venter)) #
inp$Elev<-as.numeric(scale(network.df$elev))
inp$Lat<-as.numeric(scale(network.df$latitude))

### regression models with and without interaction terms

# mean species richness

mod.wint<-glm(inp[,1]~Elev*Lat+Forest+HFP, data = inp, family = Gamma (link = "log"))
mod<-glm(inp[,1]~Elev+Lat+Forest+HFP, data = inp, family = Gamma (link = "log"))

#covariance

mod.wint<-glm(inp[,2]~Elev*Lat+Forest+HFP, data = inp, family = Gamma (link = "log"))
mod<-glm(inp[,2]~Elev+Lat+Forest+HFP, data = inp, family = Gamma (link = "log"))

# average degree

mod.wint<-betareg(inp[,5]~Elev*Lat+Forest+HFP, data = inp)
# without interaction term
mod<-betareg(inp[,5]~Elev+Lat+Forest+HFP, data = inp)

#connectance 
mod.wint<-glm(inp[,6]~Elev*Lat+Forest+HFP, data = inp)
mod<-glm(inp[,6]~Elev+Lat+Forest+HFP, data = inp)

#modularity

mod.wint<-glm(inp[,4]~Elev*Lat+Forest+HFP, data = inp)
mod<-glm(inp[,4]~Elev+Lat+Forest+HFP, data = inp)

# clustering 

mod.wint<-betareg(inp[,3]~Elev*Lat+Forest+HFP, data = inp)
mod<-betareg(inp[,3]~Elev+Lat+Forest+HFP, data = inp)


## subsets of data ##

# 1. taking only networks with modularity r_com values >0.4

hist(inp$r_com)

inp_rcom<-inp[(inp$r_com>0.4),] #50 networks, might change because of bootstrapping

mod.rcom.wint<-glm(inp_rcom$obs_modularity ~Elev*Lat+Forest+HFP, data = inp_rcom)
mod.rcom<-  glm(inp_rcom$obs_modularity~Elev+Lat+Forest+HFP, data = inp_rcom)


#2.  Datasets with 20+ flocks (N = 72)

inp$Nflocks<-network.df$Nflocks
inp20 <- inp[inp$Nflocks>19,]

# mean species richness

mod.wint<-glm(inp20[,1]~Elev*Lat+Forest+HFP, data = inp20, family = Gamma (link = "log"))
mod<-glm(inp20[,1]~Elev+Lat+Forest+HFP, data = inp20, family = Gamma (link = "log"))

#covariance

mod.wint<-glm(inp20[,2]~Elev*Lat+Forest+HFP, data = inp20, family = Gamma (link = "log"))
mod<-glm(inp20[,2]~Elev+Lat+Forest+HFP, data = inp20, family = Gamma (link = "log"))

# average degree

mod.wint<-betareg(inp20[,5]~Elev*Lat+Forest+HFP, data = inp20)
# without interaction term
mod<-betareg(inp20[,5]~Elev+Lat+Forest+HFP, data = inp20)

#connectance 
mod.wint<-glm(inp20[,6]~Elev*Lat+Forest+HFP, data = inp20)
mod<-glm(inp20[,6]~Elev+Lat+Forest+HFP, data = inp20)

#modularity

mod.wint<-glm(inp20[,4]~Elev*Lat+Forest+HFP, data = inp20)
mod<-glm(inp20[,4]~Elev+Lat+Forest+HFP, data = inp20)

# clustering 

mod.wint<-betareg(inp20[,3]~Elev*Lat+Forest+HFP, data = inp20)
mod<-betareg(inp20[,3]~Elev+Lat+Forest+HFP, data = inp20)

#3.  Datasets with 10+ flocks, below 3500 m asl (eliminating flocks in Polylepis dominated forests)

inp$altura<-network.df$elev
inp3<- inp[inp$altura<3500,] #81 networks

# mean species richness

mod.wint<-glm(inp3[,1]~Elev*Lat+Forest+HFP, data = inp3, family = Gamma (link = "log"))
mod<-glm(inp3[,1]~Elev+Lat+Forest+HFP, data = inp3, family = Gamma (link = "log"))

#covariance

mod.wint<-glm(inp3[,2]~Elev*Lat+Forest+HFP, data = inp3, family = Gamma (link = "log"))
mod<-glm(inp3[,2]~Elev+Lat+Forest+HFP, data = inp3, family = Gamma (link = "log"))

# average degree

mod.wint<-betareg(inp3[,5]~Elev*Lat+Forest+HFP, data = inp3)
# without interaction term
mod<-betareg(inp3[,5]~Elev+Lat+Forest+HFP, data = inp3)

#connectance 
mod.wint<-glm(inp3[,6]~Elev*Lat+Forest+HFP, data = inp3)
mod<-glm(inp3[,6]~Elev+Lat+Forest+HFP, data = inp3)

#modularity

mod.wint<-glm(inp3[,4]~Elev*Lat+Forest+HFP, data = inp3)
mod<-glm(inp3[,4]~Elev+Lat+Forest+HFP, data = inp3)

# clustering 

mod.wint<-betareg(inp3[,3]~Elev*Lat+Forest+HFP, data = inp3)
mod<-betareg(inp3[,3]~Elev+Lat+Forest+HFP, data = inp3)




