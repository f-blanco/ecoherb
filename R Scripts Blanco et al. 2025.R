######R Scripts Blanco et al. 2025#######




####libraries#######
library(divDyn)
library(epm)
library(ggalluvial)

library(ggpubr)
library(AICcmodavg)
library(segmented)
library(bipartite)
library(gdata)
library(funModeling)
library(betapart)
library(doSNOW)
library(scales)
library(picante)
library(plotly)
library(binovisualfields)
library(googlesheets4)
library(googledrive)
library(tidyverse)
library(ggplot2)
library(ggdendro)
library(openxlsx)
library(nlme)
library(paleoTS)
library(cluster)
library(plotrix)
library(wesanderson)
library(munsell)
library(RColorBrewer)
library(stringr)
library(dplyr)
library(tidyr)
library(vegan)
library(ape)
library(geiger)
library(phytools)
library(ade4)
library(doParallel)
library(picante)
library(fastmatch)
library(foreach)
library(factoextra)
library(countrycode)
library(HDInterval)
library(devtools)


###### Functions Input Network Analysis #########


#### Functions to create the Infomap input for the network analysis

mat.link.input.bipartite<-function(mat){
  edg<-list()
  
  for(i in 1:ncol(mat)){
    
    edg[[i]]<-cbind(i,(which(mat[,i]>0)+ncol(mat)),mat[which(mat[,i]>0),i])
    
  }
  do.call(rbind,edg)
}
write.input.info.bipartite<-function(mat,path,file){
  link.node<-mat.link.input.bipartite(mat)
  vert<-cbind(c(1:(ncol(mat)+nrow(mat))),c(colnames(mat),rownames(mat))) 
  input<-list(cbind("*Vertices",nrow(vert)),vert,cbind("*Edges",nrow(link.node)),link.node)####old input
  #input<-list(cbind("*Vertices",nrow(vert)),vert,cbind("*Bipartite",ncol(occ_analysis_functional_types)+1),link.node)###bipartite
  write.table(input[[1]], paste(path,file,sep=""),row.names=FALSE,col.names=FALSE,quote=FALSE)
  write.table(input[[2]], paste(path,file,sep=""), append= TRUE,row.names=FALSE,col.names=FALSE,quote=c(2))
  write.table(input[[3]], paste(path,file,sep=""), append= TRUE,row.names=FALSE,col.names=FALSE,quote=FALSE )
  write.table(input[[4]], paste(path,file,sep=""), append= TRUE,row.names=FALSE,col.names=FALSE,quote=FALSE)
  return(input)
}


run.infomap <- function(path.info, path.in, name, extension, path.out, parameters){
  system(paste("bash -c", shQuote(paste("cd ", path.info ," && ./Infomap ", path.in, name , extension  ," ", path.out, " ", parameters, sep=""),type="sh")),wait=T)
  out.win<-gsub("/mnt/c","C:",path.out)
  par<-read.table(paste(out.win,name,".tree",sep="")) 
  scan_csv <- scan(paste(out.win,name,".tree",sep=""), what="character", sep=NULL)
  ####Activate to have one output for every possible combination of relaxation and percolation value
  #par<-read.table(paste(out.win,name,"_","states",".tree",sep="")) 
  #scan_csv <- scan(paste(out.win,name,"_","states",".tree",sep=""), what="character", sep=NULL)
  mod_value <- scan_csv[grep("savings",scan_csv)+1]
  quality <- as.numeric(gsub("%","",mod_value))
  list(par,quality)
}



sam.val.fun <- function(x,path.unix,grep,treshold,Ntrain,Ntest,Nsam){
  path.win<-gsub("/mnt/c","C:",path.unix)
  setwd(path.win)
  write.table(x,"input_partitions.txt",sep=" ",row.names=F,col.names=F,quote=F)
  k<- paste(Ntrain,Ntest,Nsam,sep=" ")
  system(paste("bash -c", shQuote(paste("cd ",path.unix," && ./partition-validation -s ", runif(1,0,1000)," -t ", treshold, " --validation-sampling ", k, " input_partitions.txt ", grep, sep=""),type="sh")),wait=T)
  read.table(paste(grep,"_validation",sep=""))
}

jac.fun<-function(x,y){sum(x%in%y)/length(unique(c(x,y)))}


##### Load data  ####
#Load raw data
large_herb_matrix<- read.xlsx()
#Load fdata
no_na_fdata_art <- read.xlsx()
#Load the occ species
occ_analysis_fdata <- read.xlsx()
#Load ftypes
functional_types <- read.xlsx()
#Load occ ftypes
occ_analysis_functional_types<- read.xlsx()

###########  Functional distances ###### 

###Calculate the distance matrix

####Order the traits we can

hyp_ordered <- c("bra","mes","hyp")
no_na_fdata_art$hypsodonty <- ordered(no_na_fdata_art$hypsodonty, levels=hyp_ordered)
buccal_ordered <- c("2","3","M")
no_na_fdata_art$buccal_cusp <- ordered(no_na_fdata_art$buccal_cusp, levels=buccal_ordered)
lingual_ordered <- c("1","2","3","M")
no_na_fdata_art$lingual_cusp <- ordered(no_na_fdata_art$lingual_cusp, levels=lingual_ordered)
longitudinal_ordered <- c("0","1","2")
no_na_fdata_art$longitudinal_lophs <- ordered(no_na_fdata_art$longitudinal_lophs, levels=longitudinal_ordered)
transverse_ordered <- c("0","1","2","3","M")
no_na_fdata_art$transverse_lophs <- ordered(no_na_fdata_art$transverse_lophs, levels=transverse_ordered)
horizodonty_ordered <- c("0","1","2","3")
no_na_fdata_art$horizodonty <- ordered(no_na_fdata_art$horizodonty, levels=horizodonty_ordered)
bs_ordered <- c("B","C","D","E","F","G","H","I","J")
no_na_fdata_art$body_size_categories <- ordered(no_na_fdata_art$body_size_categories, levels= bs_ordered)

# dental traits
weights <- rep(1/6, ncol(no_na_fdata_art))

# Body size
weights[colnames(no_na_fdata_art) %in% c("body_size_categories")] <- bs_weight <- 1

# tooth shape
weights[colnames(no_na_fdata_art) %in% c("tooth_shape_multicusp")] <- 1


#### Calculate the dissimilarity matrix
distances <- NULL
#which variables are ordered
ordratio <- c(2,4,5,6,7,8,14)
distances <- daisy(no_na_fdata_art, metric="gower", type = list(ordratio = ordratio),weights = weights)

distance_matrix <- as.matrix(distances)

pcoa_sp_distances <- NULL
pcoa_sp_distances <- pcoa(distance_matrix, correction="none")
pcoa_vectors <-as.data.frame (pcoa_sp_distances$vectors) 
pcoa_disim_analy <- pcoa_vectors[,1:3]


#########################Betapart to calculate the functional distances between sites


pcoa_disim_analy <- range01(pcoa_disim_analy)


###To use a matrix as occ_analysis

occ_analysis <- occ_analysis_fdata


t_occ_analysis<- t(occ_analysis) 
t_occ_analysis <- t_occ_analysis[rowSums(t_occ_analysis)>5,]
t_occ_analysis <- t_occ_analysis[,colSums(t_occ_analysis)>0]



##### First, a loop to pairwise calculate the sites functional distance and create the Infomap input 

combinations_sites <- combn(rownames(t_occ_analysis),2,simplify = FALSE)
nestedness_matrix <- matrix(data=0,ncol = length(rownames(t_occ_analysis)),nrow =length(rownames(t_occ_analysis)),dimnames = list(c(rownames(t_occ_analysis)),c(rownames(t_occ_analysis))))
turnover_matrix<- matrix(data=0,ncol = length(rownames(t_occ_analysis)),nrow =length(rownames(t_occ_analysis)),dimnames = list(c(rownames(t_occ_analysis)),c(rownames(t_occ_analysis))))
betadiv_matrix<- matrix(data=0,ncol = length(rownames(t_occ_analysis)),nrow =length(rownames(t_occ_analysis)),dimnames = list(c(rownames(t_occ_analysis)),c(rownames(t_occ_analysis))))


for (i in 1:length(combinations_sites)) {
  
  
  
  focus_sites<- combinations_sites[[i]] 
  

  
  occ_functional_tryal <- as.data.frame(t_occ_analysis[focus_sites,])
  occ_beta_pair <- occ_functional_tryal[,colSums(occ_functional_tryal)>0]
  pcoa_traits <- pcoa_disim_analy[colnames(occ_beta_pair),]

  functiona_dissim_values <- functional.beta.pair(occ_beta_pair, pcoa_traits, index.family="sorensen")
  

  nestedness_matrix[focus_sites[1],focus_sites[2]] <- functiona_dissim_values$funct.beta.sne
  turnover_matrix[focus_sites[1],focus_sites[2]] <- functiona_dissim_values$funct.beta.sim
  betadiv_matrix[focus_sites[1],focus_sites[2]] <- functiona_dissim_values$funct.beta.sor



  print(i/length(combinations_sites))
  
}



class(nestedness_matrix)


##Second, copy the upper triangle into the lower to have a symmetric matrix
lowerTriangle(nestedness_matrix) = upperTriangle(nestedness_matrix, byrow=TRUE)
lowerTriangle(turnover_matrix) = upperTriangle(turnover_matrix, byrow=TRUE)
lowerTriangle(betadiv_matrix) = upperTriangle(betadiv_matrix, byrow=TRUE)


turnover_matrix[turnover_matrix<0] <- 0
 
##We calculate the nmds to plot it

nmds_nested <-  metaMDS(nestedness_matrix, wascores = T, try=100, autotransform = FALSE, k=2, distance = "euclidean", trace=2)
nmds_turnover <-  metaMDS(turnover_matrix, wascores = T, try=100, autotransform = FALSE, k=2, distance = "euclidean", trace=2)
nmds_betadiv <-  metaMDS(betadiv_matrix, wascores = T, try=100, autotransform = FALSE, k=2, distance = "euclidean", trace=2)


##Plot evolution of betadiv parameters through time (previous bin with next bin)

selected_continent <- unique(sapply(strsplit(rownames(turnover_matrix),"_"),function(x)x[2]))


for (z in selected_continent) {
  
  continent_names<- rownames(nestedness_matrix)[sapply(strsplit(rownames(nestedness_matrix),"_"),function(x)x[2]) %in% z]
  
  nested_continent <- nestedness_matrix[continent_names,continent_names]
  turnover_continent <- turnover_matrix[continent_names,continent_names]
  betadiv_continent <- betadiv_matrix[continent_names,continent_names]
  
  name <- NULL
  nestedness_value <- NULL
  turnover <- NULL
  betadiv <- NULL
  bin_value_nest <- NULL
  bin_value_tur <- NULL
  bin_value_beta <- NULL
  bin_name <- NULL
  
  for (i in 1:length(rownames(nested_continent))-1) {
    
    
    bin_value_nest <- nested_continent[i,i+1]
    bin_value_tur <- turnover_continent[i,i+1]
    bin_value_beta <- betadiv_continent[i,i+1]
    bin_name <- rownames(nested_continent)[i]
    
    
    name<- c(name,bin_name)
    nestedness_value<- c(nestedness_value,bin_value_nest)
    turnover<- c(turnover,bin_value_tur)
    betadiv <- c(betadiv,bin_value_beta)
    
    
  }
  
  
  
  betapart_plot <- cbind(name,nestedness_value,turnover,betadiv)
  betapart_plot <- as.data.frame(betapart_plot)
  
  betapart_plot$mean_age <- bins_loc$mean_age[rownames(bins_loc) %in% betapart_plot$name]
  
  
  quartz(width = 16, height = 12)
  plot(NULL, xlim= rev(range(betapart_plot$mean_age)),ylim = c(0,1),  ylab = z, xlab="time")
  points(betapart_plot$mean_age,betapart_plot$nestedness_value,col="grey",cex=3,pch=19)
  points(betapart_plot$mean_age,betapart_plot$turnover,col="#F8485E",cex=3,pch=19)
  points(betapart_plot$mean_age,betapart_plot$betadiv,col="#00C1D4",cex=3,pch=19)
  lines(betapart_plot$mean_age,betapart_plot$nestedness_value,col="grey",lwd=5)
  lines(betapart_plot$mean_age,betapart_plot$turnover,col="#F8485E",lwd=5)
  lines(betapart_plot$mean_age,betapart_plot$betadiv,col="#00C1D4",lwd=5)
  legend("topright",legend = c("nestedness","turnover","betadiv"),pch=17,col=c("grey","#F8485E","#00C1D4"))
  
  
}






########### Run Network Analysis #######

# input #
# ::::: #
path <-  ""

write.input.info.bipartite(occ_analysis_functional_types,path,file="Infomap.net")



# Run infomap #
# :::::::::::::: #
path.info<-"/Applications/infomap-1.1.4" 
path.in<-"/Users/fernando/Desktop/" 
name<-"Infomap"
extension<-".net"
path.out<-"/Users/fernando/Desktop/"


par <- list()
q <- c()
for( i in 1:100){
  seed <- sample(1:1000,1) 
  parameters<-paste("-N 100","-s", seed) 
  res.i <- run.infomap (path.info, path.in, name, extension, path.out, parameters)
  par[[i]] <- res.i[[1]] 
  q[[i]] <- res.i[[2]]
}



  
# Search completeness  #
# ::::::::::::::::::::: #
# Download "partition-validation": https://github.com/mapequation/partition-validation

par.c <- lapply(par,function(x)as.character(x[order(x[,3]),1]))
par.c <- do.call(cbind,lapply(par.c,function(x)sapply(strsplit(x,":"),function(x)x[1])))####change the number to get the modules or submodules
par.c <- par.c[,order(1-as.numeric(q))]


# Module significance  #
# :::::::::::::::::::::: #

id.nodes <- as.character(par[[1]][order(par[[1]][,3]),3])
mod.ref <- tapply(id.nodes,par.c[,1],list)
mod.alt <- apply(par.c[,-1],2,function(x)tapply(id.nodes,x,list))
id.mod <- names(mod.ref)

#similarity to the most similar
#::::::::::::::::::::::::::::::

jac.max<-list()
for(i in id.mod){
  sim.avg<-c()
  for(j in 1:length(mod.alt)){
    sim <- unlist(lapply(mod.alt[[j]],function(x){jac.fun(mod.ref[[i]],x)}))
    sim.avg[j] <- max(sim)
    names(sim.avg)[j]	<- names(sim)[sim==max(sim)]
    print(c(i,j))
  }
  jac.max[[i]]<- sim.avg
}


#proportion more similar than a threshold (0.75)
#::::::::::::::::::::::::::::::::::::::::
# Module robustness
p.sim <- sapply(jac.max,function(x)sum(x>=0.75)/length(x)) 


#Node probability module assignation
#::::::::::::::::::::::::::
res.prob <- list()
for(i in 1:length(jac.max)){
  id.mod.i <- names(jac.max[[i]])
  nodes.mod.i <-c()
  for(j in 1:length(id.mod.i)){
    nodes.mod.i<-c(nodes.mod.i,mod.alt[[j]][[id.mod.i[j]]])
  }
  prob <- data.frame(table(nodes.mod.i)/length(id.mod.i))			
  res.prob[[i]] <- prob[prob[,1]%in% mod.ref[[names(jac.max)[i]]],]
  print(i)
}
names(res.prob) <- names(jac.max)

head(res.prob)

#module similarity matrix
#:::::::::::::::::::::::::
sim.best.mod <- data.frame(matrix(ncol=length(jac.max),nrow=length(jac.max)))
colnames (sim.best.mod) <- rownames(sim.best.mod) <- names(jac.max)
for(i in 1:length(jac.max)){
  for(j in 1:length(jac.max)){
    sim.best.mod[i,j] <- sum(names(jac.max[[i]])==names(jac.max[[j]]))/length(mod.alt)
  }
}
module_similarity<- round(sim.best.mod,3)


#module similarity matrix only based on the best partition. 
#:::::::::::::::::::::::::

par <- cbind(id.nodes,par.c[,1])
id.gf <- rownames(occ_analysis_functional_types)
id.lc <- colnames(occ_analysis_functional_types)
par.gf <- par[par[,1]%in%id.gf,]
par.lc <- par[par[,1]%in%id.lc,]
id.mod <- unique(par.lc[,2])

sim.best.mod.B <- data.frame(matrix(ncol=length(jac.max),nrow=length(jac.max)))
colnames (sim.best.mod.B) <- rownames(sim.best.mod.B) <- id.mod
for(i in 1:length(id.mod)){
  for(j in 1:length(id.mod)){
    lc.mod.i <- as.character(par.lc[par.lc[,2]==id.mod[i],1])
    lc.mod.j <- as.character(par.lc[par.lc[,2]==id.mod[j],1])
    mat.jac <- cbind(rowSums(occ_analysis_functional_types[lc.mod.i]),rowSums(occ_analysis_functional_types[lc.mod.j]))
    sim.best.mod.B[i,j] <- sum(apply(mat.jac,1,min))/sum(apply(mat.jac,1,max))
  }
}
module_similarity.B <- round(sim.best.mod.B,3)



# affinity, fidelity and indval #
# :::::::::::::::::::::::::::::::: #

par <- cbind(id.nodes,par.c[,1])
id.gf <- rownames(occ_analysis_functional_types)
id.lc <- colnames(occ_analysis_functional_types)
par.gf <- par[par[,1]%in%id.gf,]
par.lc <- par[par[,1]%in%id.lc,]


specificity.lc <- c()
afinity.lc <- c()
for(i in 1:nrow(par.lc)){
  mod.i <- par.lc[i,2]
  gf.mod.i <- as.character(par.gf[par.gf[,2]==mod.i,1])
  lc.mod.i <- as.character(par.lc[par.lc[,2]==mod.i,1])
  val.mod <- sum(occ_analysis_functional_types[rownames(occ_analysis_functional_types)%in% gf.mod.i,colnames(occ_analysis_functional_types)==par.lc[i,1]])
  val.tot.lc.i <- sum(occ_analysis_functional_types[,colnames(occ_analysis_functional_types)==par.lc[i,1]])
  val.mod.gf <- sum(occ_analysis_functional_types[rownames(occ_analysis_functional_types)%in% gf.mod.i,colnames(occ_analysis_functional_types)==par.lc[i,1]]>0)
  specificity.lc[i]<-val.mod/val.tot.lc.i
  afinity.lc [i] <- val.mod.gf/length(gf.mod.i)
}


specificity.gf <- c()
afinity.gf <- c()
for(i in 1:nrow(par.gf)){
  mod.i <- par.gf[i,2]
  gf.mod.i <- as.character(par.gf[par.gf[,2]==mod.i,1])
  lc.mod.i <- as.character(par.lc[par.lc[,2]==mod.i,1])
  val.mod <- sum(occ_analysis_functional_types[rownames(occ_analysis_functional_types)==par.gf[i,1] ,colnames(occ_analysis_functional_types)%in%lc.mod.i])
  val.tot.gf.i <- sum(occ_analysis_functional_types[rownames(occ_analysis_functional_types)==par.gf[i,1],])
  val.mod.lc <- sum(occ_analysis_functional_types[rownames(occ_analysis_functional_types)==par.gf[i,1],colnames(occ_analysis_functional_types)%in%lc.mod.i]>0)
  specificity.gf[i]<-val.mod/val.tot.gf.i
  afinity.gf [i] <- val.mod.lc/length(lc.mod.i)
}

res_fin <- data.frame(ID = c(as.character(par.lc[,1]),as.character(par.gf[,1])),specificidad = c(specificity.lc,specificity.gf), afinidad= c(afinity.lc,afinity.gf))
res_fin$indval <- res_fin[,2]* res_fin[,3]
res_fin[,2:4]<-apply(res_fin[,2:4],2,as.numeric)



# Final table		     #
# :::::::::::::::::::::::::::::::: #

res.prob<-do.call(rbind,res.prob)
colnames(res.prob) <- c("ID","P_class")
colnames(par)<-c("ID","module")
par <- as.data.frame(par)
RES <- merge(par,res_fin,by="ID")
RES <- merge(RES,res.prob,by="ID")
head(RES)
unique(RES$module)

network_stats_table <- RES

network_stats_table$new_module <- bins_loc$module[match(network_stats_table$module,bins_loc$old_module)]


network_stats_table$new_module <- network_stats_table$new_module-1


########### Sensitivity analysis #######

#Sensitivity to the selection of kmeans groups
kmeans_groups <- seq(100,500,100)

detectCores()
cl <- makeCluster(12)
registerDoSNOW(cl)
sensi_loop <- NULL


sensi_loop<- foreach(i=kmeans_groups, .options.snow = opts, .packages=c("stats")) %dopar% {
  
  randomizations <- 100
  
  f_groups_list <- NULL
  module_number_list <- NULL
  links_number_list <- NULL
  node_number_list <- NULL
  porportion_links_list <- NULL
  quality_list <- NULL
  ran_list <- NULL
  


for (ran in 1:randomizations) {
    
   
    km <- kmeans(x = pcoa_vectors, centers=i,nstart = 100, iter.max = 100000)
    functional_types <- km$cluster
      
      
    occ_analysis_fdata <- occ_analysis_fdata[names(functional_types),]
    
    
    occ_sensitivity <- aggregate(occ_analysis_fdata, by=list(functional_types),sum)
    rownames(occ_sensitivity) <- paste("FT",occ_sensitivity$Group.1,sep="_")
    occ_sensitivity<- occ_sensitivity[,-1]
    
    
    file_name <-  paste(i,"_",ran,"_","Infomap.net",sep = "")
    path_folder <-  "/Users/fernando/Desktop/large_herb_sesnsitivity_analysis/"
    write.input.info.bipartite(occ_sensitivity,file=file_name,path=path_folder)
    
    
    path.info<-"/Applications/infomap-1.1.4" ##if we run it in the iMacPro
    #path.info<-"/Applications/Infomap"
    path.in<-"/Users/fernando/Desktop/large_herb_sesnsitivity_analysis/"
    name<-paste(i,"_",ran,"_","Infomap",sep = "")
    extension<-".net"
    path.out<-"/Users/fernando/Desktop/large_herb_sesnsitivity_analysis/tree/"
    
    
    seed <- sample(1:1000,1)
    parameters<-paste("-N 100", "-s",seed)
    res.i <- run.infomap (path.info, path.in, name, extension, path.out, parameters)
    
    
    f_groups <- i
    
    new_network <- res.i[[1]]
    module_number <- length(unique(sapply(strsplit(new_network$V1,":"),function(x)x[1])))
    
    link.node<-mat.link.input.bipartite(occ_sensitivity)
    vert<-cbind(c(1:(ncol(occ_sensitivity)+nrow(occ_sensitivity))),c(colnames(occ_sensitivity),rownames(occ_sensitivity))) 
    links_number <- nrow(link.node)
    node_number <- nrow(vert)
    
    proportion_links <- links_number/node_number
    quality_number <- res.i[[2]]
    
    f_groups_list <- c(f_groups_list,f_groups)
    module_number_list <- c(module_number_list,module_number)
    links_number_list <- c(links_number_list,links_number)
    node_number_list <- c(node_number_list,node_number)
    porportion_links_list <- c(porportion_links_list,proportion_links)
    quality_list <- c(quality_list,quality_number)
    ran_list <- c(ran_list,ran)
  
  
  
  
}
  
  
  return(list(f_groups_list,module_number_list,links_number_list,node_number_list,porportion_links_list,quality_list,ran_list))
  
}




fg<- unlist(lapply(sensi_loop, function(x) x[[1]]))
mod<- unlist(lapply(sensi_loop, function(x) x[[2]]))
link<- unlist(lapply(sensi_loop, function(x) x[[3]]))
node<- unlist(lapply(sensi_loop, function(x) x[[4]]))
prop<- unlist(lapply(sensi_loop, function(x) x[[5]]))
qua<- unlist(lapply(sensi_loop, function(x) x[[6]]))
r<- unlist(lapply(sensi_loop, function(x) x[[7]]))


sensitivity_network_matrix <- do.call(rbind, Map(data.frame, functional_groups=fg,modules=mod,Links=link,nodes=node,proportion=prop,quality_code_length=qua,randomization=r))
mean_matrix <- aggregate(sensitivity_network_matrix,by=list(sensitivity_network_matrix$functional_groups),mean)


dim(sensitivity_network_matrix)


head(sensitivity_network_matrix)
sensitivity_network_matrix$proportion

p <- ggplot(sensitivity_network_matrix, aes(functional_groups, quality_code_length)) +
  geom_point() +
  geom_smooth(method = "loess")

quartz()
p


####Sensitivity Analysis to sampling



detectCores()
cl <- makeCluster(12)
registerDoSNOW(cl)

iteraciones <- 100

pb <- txtProgressBar(max = length(1:iteraciones), style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

sampling_loop <- NULL



sampling_loop<- foreach(i_sampling_loop=1:iteraciones, .options.snow = opts, .packages=c("stats","betapart")) %dopar% {
  
  
  
  functional_types_sampling <-functional_types[names(functional_types) %in% rownames(occ_large_herb)]
  sampling_occ_fdata <- occ_large_herb[names(functional_types_sampling),]
  sampling_occ_def <- sampling_occ_fdata[,colSums(sampling_occ_fdata)>0]
  sampling_occ_def <- as.data.frame(sampling_occ_def)
  sampling_occ_def <- sampling_occ_def[,colSums(sampling_occ_def)>=20]
  
  
  for (z_sampling in 1:length(colnames(sampling_occ_def))) {
    
    bin_selected <- colnames(sampling_occ_def)[z_sampling]
    sub_matrix_bin_selected <- sampling_occ_def[bin_selected]
    sp_list_bin_selected <- names(which(rowSums(sub_matrix_bin_selected)==1))
    sampled_spp<- sample(sp_list_bin_selected,20)
    sampled_spp_selected <- rownames(sampling_occ_def)[rownames(sampling_occ_def) %in% sampled_spp]
    sampling_occ_def[sampled_spp_selected,bin_selected]<- 0
    
  }
  
  
  
  FT_sampling_sensitivity <- aggregate(sampling_occ_def, by=list(functional_types_sampling),sum)
  rownames(FT_sampling_sensitivity) <- paste("FT",FT_sampling_sensitivity$Group.1,sep="_")
  FT_sampling_sensitivity<- FT_sampling_sensitivity[,-1]
  
  
  
  t_sampling_occ_def<- t(sampling_occ_def) 
  t_sampling_occ_def <- t_sampling_occ_def[,colSums(t_sampling_occ_def)>0]
  combinations_sites_sampling <- combn(rownames(t_sampling_occ_def),2,simplify = F)
  sampling_nestedness_matrix <- matrix(data=0,ncol = length(rownames(t_sampling_occ_def)),nrow =length(rownames(t_sampling_occ_def)),dimnames = list(c(rownames(t_sampling_occ_def)),c(rownames(t_sampling_occ_def))))
  sampling_turnover_matrix<- matrix(data=0,ncol = length(rownames(t_sampling_occ_def)),nrow =length(rownames(t_sampling_occ_def)),dimnames = list(c(rownames(t_sampling_occ_def)),c(rownames(t_sampling_occ_def))))
  sampling_betadiv_matrix<- matrix(data=0,ncol = length(rownames(t_sampling_occ_def)),nrow =length(rownames(t_sampling_occ_def)),dimnames = list(c(rownames(t_sampling_occ_def)),c(rownames(t_sampling_occ_def))))
  
  
  for (x_sampling_combinations in 1:length(combinations_sites_sampling)) {
    
    focus_sites_sampling<- combinations_sites_sampling[[x_sampling_combinations]] 
    
    occ_sampling_tryal <- as.data.frame(t_sampling_occ_def[focus_sites_sampling,])
    occ_beta_pair_sampling <- occ_sampling_tryal[,colSums(occ_sampling_tryal)>0]
    pcoa_traits_sampling <- pcoa_disim_analy[colnames(occ_beta_pair_sampling),]
    
    sampling_functiona_dissim_values <- functional.beta.pair(occ_beta_pair_sampling, pcoa_traits_sampling, index.family="sorensen")
    
    sampling_nestedness_matrix[focus_sites_sampling[1],focus_sites_sampling[2]] <- sampling_functiona_dissim_values$funct.beta.sne
    sampling_turnover_matrix[focus_sites_sampling[1],focus_sites_sampling[2]] <- sampling_functiona_dissim_values$funct.beta.sim
    sampling_betadiv_matrix[focus_sites_sampling[1],focus_sites_sampling[2]] <- sampling_functiona_dissim_values$funct.beta.sor
    
  }
  
  
  
  
  return(list(sampling_nestedness_matrix,sampling_turnover_matrix,sampling_betadiv_matrix,FT_sampling_sensitivity,sampling_occ_def))
  
}


###We extract the 100 matrices for each matrix
sampling_nested_100_it<- lapply(sampling_loop, function(x) x[[1]])
sampling_turnover_100_it<- lapply(sampling_loop, function(x) x[[2]])
sampling_betadiv_100_it<- lapply(sampling_loop, function(x) x[[3]])
FT_sampling_100_it<- lapply(sampling_loop, function(x) x[[4]])
occ_def_sampling_100_it<- lapply(sampling_loop, function(x) x[[5]])




###We merge the 100 matrices into 1 and make the mean matrix, which is the one we are going to use
nested_mean_matrix_sampling <- Reduce("+",sampling_nested_100_it)/length(sampling_nested_100_it)
turnover_mean_matrix_sampling <- Reduce("+",sampling_turnover_100_it)/length(sampling_turnover_100_it)
betadiv_mean_matrix_sampling <- Reduce("+",sampling_betadiv_100_it)/length(sampling_betadiv_100_it)


lowerTriangle(turnover_mean_matrix_sampling) = upperTriangle(turnover_mean_matrix_sampling, byrow=TRUE)
sampling_nmds_turnover <-  metaMDS(turnover_mean_matrix_sampling, wascores = T, try=100, autotransform = FALSE, k=2, distance = "euclidean", trace=2)




####### Null model network 


random_q <- list()


for (n in 1:100) {
  
  occ_null <- aggregate(occ_analysis_fdata, by=list(sample(functional_types)),sum)
  
  rownames(occ_null) <- paste("FT",occ_null$Group.1,sep="_")
  occ_null<- occ_null[,-1]
  
 
  file <-  "Infomap_random.net"
  path <-  "/Users/fernando/Desktop/"
  write.input.info.bipartite(occ_null,path,file)
  
  
  # Run infomap #
  # :::::::::::::: #
 
 
  path.info<-"/Applications/infomap-1.1.4" 
  path.in<-"/Users/fernando/Desktop/" 
  name<-"Infomap_random"
  extension<-".net"
  path.out<-"/Users/fernando/Desktop/"
  
  
  
  
  par <- list()
  q <- c()
  for( i in 1:100){
    seed <- sample(1:1000,1)
    parameters<-paste("-N 100", "-s",seed)
    res.i <- run.infomap (path.info, path.in, name, extension, path.out, parameters)
    par[[i]] <- res.i[[1]]
    q[i] <- res.i[[2]]
  }
  
  random_q[[n]] <- q
  
}



null_orders_q <- random_q; save(null_orders_q,file= "~/Desktop/null_order_fe_q.RData")


max_null <- unlist(lapply(null_orders_q, max))



max_null[max_null<0] <- 0


## OBSERVED Q

occ_fe <- aggregate(occ_analysis_fdata, by=list((functional_types)),sum)
rownames(occ_fe) <- paste("FT",occ_fe$Group.1,sep="_")
occ_fe<- occ_fe[,-1]


file <-  "Infomap_observed.net"
path <-  "/Users/fernando/Desktop/"
write.input.info.bipartite(occ_fe,path,file)


# Run infomap #
# :::::::::::::: #

path.info<-"/Applications/infomap-1.1.4" 
path.in<-"/Users/fernando/Desktop/" 
name<-"Infomap_observed" 
extension<-".net"
path.out<-"/Users/fernando/Desktop/"

par <- list()
q <- c()
for( i in 1:100){
  seed <- sample(1:1000,1)
  parameters<-paste("-N 100", "-s",seed)
  res.i <- run.infomap (path.info, path.in, name, extension, path.out, parameters)
  par[[i]] <- res.i[[1]]
  q[i] <- res.i[[2]]
}





observed_orders <- 11.8665 ## max(q)





sum(max_null>observed_orders)/100

quartz(width = 8, height = 4)
cex.axis <- 0.9
col="grey50"
plot(density(max_null),xlim=c(0,observed_orders/100),type="n", main="", axes=F, ylab="", xlab="")
polygon(density(max_null), col=alpha("slategray", 0.2),  border=alpha("slategray", 0.9))
abline(v=observed_orders/100, col="orangered", lwd=2)
axis(1,tcl=-0.5,cex.axis=cex.axis,col=col,col.axis=col)
axis(2,tcl=-0.5,cex.axis=cex.axis,col=col,col.axis=col, las=1)
mtext("null modularity",side=1,col=col,line=2,cex=cex.axis)
mtext("density",side=2,col=col,line=2.3,cex=cex.axis)



########### Plot modules ######
  
##### Normal PLot
  
par_original  <- read.table("/Users/fernando/Desktop/Infomap.tree")
par_original$p2 <- sapply(strsplit(par_original$V1,":"),function(x)x[1])
par.1 <-par_original[par_original$V3%in%colnames(occ_analysis_functional_types),c("V3","p2")]
colnames(par.1) <- c("ID","module")

par.1 <-par.1[match(colnames(occ_analysis_functional_types),par.1$ID),]




  RES_original <- par.1



  mod_loc <- RES_original$ID
  
  bins_loc <- matrix(0,ncol =0 ,nrow =length(RES_original$ID))
  rownames(bins_loc) <- RES_original$ID
  bins_loc <- as.data.frame(bins_loc)

 dim(bins_loc)
 
  #### Create a new columns for the plot with the bins names and different values
  bins_loc$module<- RES_original$module
  bins_loc$old_module<- RES_original$module
  bins_loc$sp_richness <- colSums(occ_analysis_fdata)
  fe_ocurrence <- occ_analysis_functional_types
  fe_ocurrence[fe_ocurrence!=0] <- 1
  bins_loc$fe_richness <- colSums(fe_ocurrence)[match(rownames(bins_loc), colnames(fe_ocurrence))]
  bins_fe_diversity <- diversity(t(occ_analysis_functional_types),"shannon")
  bins_loc$fe_diversity <- bins_fe_diversity[names(bins_fe_diversity) %in% rownames(bins_loc)]
  bins_loc$bin<- sapply(strsplit(rownames(bins_loc),"_"),function(x)x[1])
  bins_loc$continent<- sapply(strsplit(rownames(bins_loc),"_"),function(x)x[2])

  ##to calculate the mean age of the bins
  bin_age_list <- list()
for (i in 1:length(breaks_cenozoic)-1) {
  age_1 <- rev(breaks_cenozoic)[i]
  age_2 <- rev(breaks_cenozoic)[1+i]
  bin_age <- (age_1+age_2)/2
  bin_age_list[i] <- bin_age
} 
  
  ##correct to put age 0 to the last bin
  bin_age_list[1] <- 0
  names(bin_age_list) <- c(1:23)
  

  bins_loc$mean_age <-as.numeric(bin_age_list)[match(bins_loc$bin,names(bin_age_list))]
  order_age <- aggregate(bins_loc$mean_age,by=list(bins_loc$module), median)
  
  #### Order the modules by age for the plot
  RES$module <- match(RES$module,order(order_age$x, decreasing = TRUE))
  
  
  bins_loc$module <- match(as.numeric(bins_loc$module),order(order_age$x, decreasing = TRUE))

  bins_loc <- bins_loc[order(bins_loc$mean_age,decreasing = TRUE),]
  
 ##PLot adjustments
  y_position_indval <-  as.numeric(bins_loc$module)*10-5+0.5*8+1
  line_africa <- as.numeric(bins_loc$module[bins_loc$continent=="Africa"])*10-5+0.5*8+1
  line_europe <- as.numeric(bins_loc$module[bins_loc$continent=="Europe"])*10-5+0.5*8+1
  line_america <- as.numeric(bins_loc$module[bins_loc$continent=="America"])*10-5+0.5*8+1
  line_asia <- as.numeric(bins_loc$module[bins_loc$continent=="Asia"])*10-5+0.5*8+1


  ####To calculate the mean age module replacement
  mean_age_module_replacement <- diff(order_age$x)
  
  mean(mean_age_module_replacement)
  std.error(mean_age_module_replacement)


  lwd=2
  line <- 2
  tcl=0.15
  line.ylabs <- -0.25
  
  bins_loc$number_continent[bins_loc$continent=="America"] <- 1
  bins_loc$number_continent[bins_loc$continent=="Europe"] <- 2
  bins_loc$number_continent[bins_loc$continent=="Africa"] <- 3
  bins_loc$number_continent[bins_loc$continent=="Asia"] <- 4

  length(unique(bins_loc$module))
  
  line_colours <- c("#e76f51","#e9c46a","#2a9d8f","#875B52")
  
 
  point_colours <- line_colours[bins_loc$number_continent]
  

  unique(bins_loc$module)
  

  points <- c(15,16,17,18)
  points <- points[bins_loc$number_continent]

  
  n_module <- max(as.numeric(bins_loc$module))
  y.lab_line <- 1.3
  cex.axis <- 1
  col="#313131"
  col_var <- "3414141"
  point.transp <- 0.1
  

  
  quartz(width = 16, height = 12)

  ########corregido para RES (Indval)
  plot(NULL, xlim=rev(range(bins_loc$mean_age)),ylim = c(5,(n_module+1)*10), axes = FALSE, ylab = "", xlab="")
  abline(h=1:n_module*10-5, col="#525252", lty=1, lwd=0.5)
  points(bins_loc$mean_age,y_position_indval, pch=points, col=alpha(point_colours,0.6), cex=4)
  lines(bins_loc$mean_age[bins_loc$continent=="America"],line_america,col=alpha(line_colours[1],0.7),lwd=20)
  lines(bins_loc$mean_age[bins_loc$continent=="Europe"],line_europe,col=alpha(line_colours[2],0.7),lwd=20)
  lines(bins_loc$mean_age[bins_loc$continent=="Africa"],line_africa,col=alpha(line_colours[3],0.7),lwd=20)
  lines(bins_loc$mean_age[bins_loc$continent=="Asia"],line_asia,col=alpha(line_colours[4],0.7),lwd=20)
  axis(1,labels=FALSE,cex.axis=cex.axis,col=col,col.axis=col, las=1, line=0, tcl=tcl)
  axis(1,cex.axis=cex.axis,col=col,col.axis=col, las=1, line= line.ylabs-0.25, tick=FALSE,tcl=tcl)
  axis(2,at=1:n_module*10,labels = c(paste("M",c(1:n_module), sep="")),cex.axis=cex.axis,col=module_colours,col.axis=col, las=1, line= line.ylabs, tick=FALSE,tcl=tcl)
  mtext("time (Ma)",side=1,col= col,line= 3,cex=cex.axis, font=2)
  legend("topleft",legend = c("America","Europe","Africa","Asia"),pch=unique(points),col=line_colours)

  
########### Plot NMDS-Modules and continent distance time #########
  
  continent_colour <- c("#e76f51","#e9c46a","#2a9d8f","#875B52")
  module_cols <- c("#9e0142","#d53e4f","#f46d43","#fdae61","#fee08b","#e6f598","#abdda4","#66c2a5","#3288bd","#5e4fa2","#875B52")

  ge_cols <- c("#FFCB45","#ffba08","#e85d04","#d00000","#6a040f", "#1d3557")
  new_palette <- colorRampPalette(ge_cols,bias= 0.9)
  ge_cols <- new_palette(20)
  
  new_age_cols <- c("#FFCB45","#ffba08","#e85d04","#d00000","#6a040f", "#1d3557")
  new_palette <- colorRampPalette(new_age_cols,bias=1.5)
  
  nmds_df_plot <- data.frame(nmds_turnover$points)
  
  ###### For the species sampling
  #nmds_df_plot <- data.frame(sampling_nmds_turnover$points)
  

  nmds_df_plot$module <- bins_loc$module[match(rownames(nmds_df_plot),rownames(bins_loc))]
  nmds_df_plot$continent <- bins_loc$continent[match(rownames(nmds_df_plot),rownames(bins_loc))]
  nmds_df_plot$number_continent <- bins_loc$number_continent[match(rownames(nmds_df_plot),rownames(bins_loc))]
  nmds_df_plot$continent_colour <-  continent_colour[nmds_df_plot$number_continent]
  nmds_df_plot$age <- bins_loc$mean_age[match(rownames(nmds_df_plot),rownames(bins_loc))]
  nmds_df_plot$bin_plot <- as.numeric(bins_loc$bin[match(rownames(nmds_df_plot),rownames(bins_loc))])-1
  nmds_df_plot$bin_colour <- ge_cols[as.numeric(nmds_df_plot$bin_plot)]
  nmds_df_plot$new_age <- round(nmds_df_plot$age, digits = 2)*100
  new_age_cols <- new_palette(max(nmds_df_plot$new_age))
  nmds_df_plot$new_age_colour <- new_age_cols[nmds_df_plot$new_age]

  nmds_df_plot <- nmds_df_plot[order(nmds_df_plot$module, decreasing = T),]
  nmds_df_plot$module_colour <-  module_cols[nmds_df_plot$module]
 
  
  
  
  font_family <- "sans"
  quartz(width = 20, height = 10,family = font_family)
  par(mfcol=c(2,4))
  #quartz(width = cm(20),height = cm(13),family = font_family)
  
 
for (plot_continent in unique(nmds_df_plot$continent)) {
  points_to_plot_2 <- nmds_df_plot$MDS2[nmds_df_plot$continent==plot_continent]
  points_to_plot_1 <- nmds_df_plot$MDS1[nmds_df_plot$continent==plot_continent]
  
  plot(NULL,xlim=range(nmds_df_plot$MDS2),ylim=range(nmds_df_plot$MDS1),xlab = "NMDS2", ylab= "NMDS1", main=plot_continent)
  for (i in 1:length(points_to_plot_2)) {
    colfunc <- colorRampPalette(c(nmds_df_plot$new_age_colour[nmds_df_plot$MDS2==points_to_plot_2[i]],nmds_df_plot$new_age_colour[nmds_df_plot$MDS2==points_to_plot_2[i+1]]))
    segments(points_to_plot_2[i],points_to_plot_1[i],points_to_plot_2[i+1],points_to_plot_1[i+1],
             col = colfunc(1000),lwd=9)
  }
  points(nmds_df_plot$MDS2,nmds_df_plot$MDS1,col=ifelse(nmds_df_plot$continent==plot_continent,alpha(nmds_df_plot$new_age_colour[nmds_df_plot$continent==plot_continent],0.1),alpha(nmds_df_plot$new_age_colour,0.1)),pch=19,cex=4)

  
  plot(NULL,xlim=range(nmds_df_plot$MDS2),ylim=range(nmds_df_plot$MDS1),xlab = "NMDS2", ylab= "NMDS1", main=plot_continent)
  lines(points_to_plot_2,points_to_plot_1,col=alpha("#979797",0.3),lwd=7)
  points(nmds_df_plot$MDS2,nmds_df_plot$MDS1,col=ifelse(nmds_df_plot$continent==plot_continent,nmds_df_plot$module_colour,alpha(nmds_df_plot$module_colour,0.1)),pch=19,cex=4)

}
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  plot(0, 0, type = 'n', bty = 'n', xaxt = 'n', yaxt = 'n')
  legend(0.97,0.96,legend = round(unique(nmds_df_plot$age[order(nmds_df_plot$age)]),digits = 2),pch=19,col=unique(nmds_df_plot$new_age_colour),xpd = TRUE, horiz = F, cex = 1.3, seg.len=1, bty = 'n',title ="Age")
  legend(0.97,-0.34,legend = unique(nmds_df_plot$module),pch=19,col=unique(nmds_df_plot$module_colour),xpd = TRUE, horiz = F, cex = 1.3, seg.len=1, bty = 'n',title="Module")
  
 
  

  ###continent distances
  
  distance_matrix_for_continents <- turnover_matrix
  

  
  continent_bins_total <- unique(nmds_df_plot$continent)

  
  continent_loop <- combn(continent_bins_total,2,simplify = FALSE)
  
  
  quartz(width = 20, height = 10)
  par(mfcol=c(2,3), oma = c(0, 0, 0, 6))
  
  
  for (i in 1:length(continent_loop)) {
    
    continent_selected_funcdist <- continent_loop[[i]]
    
  
    continent_1 <-rownames(nmds_df_plot) [nmds_df_plot$continent %in%continent_selected_funcdist[1]]
    continent_2 <- rownames(nmds_df_plot) [nmds_df_plot$continent %in%continent_selected_funcdist[2]]
    

    continent_distnace_matrix <- distance_matrix_for_continents[continent_1,continent_2]
    if(ncol(continent_distnace_matrix)!=nrow(continent_distnace_matrix)){
        dimension_matrix <- min(ncol(continent_distnace_matrix),nrow(continent_distnace_matrix))
        continent_distnace_matrix <- continent_distnace_matrix[1:dimension_matrix,1:dimension_matrix]
    }
   
    plot(NULL,xlim=rev(range(nmds_df_plot$age[rownames(nmds_df_plot)%in%rownames(continent_distnace_matrix)])),ylim=range(continent_distnace_matrix),xlab = "Time" , ylab= "Functional distance", main=paste(continent_selected_funcdist[1],continent_selected_funcdist[2],sep= " vs "))
    lines(nmds_df_plot$age[rownames(nmds_df_plot)%in%rownames(continent_distnace_matrix)],diag(continent_distnace_matrix),col=alpha("#99A799",0.7),lwd=10)
  
  }
  
  
  tamaño_linea <- 2
  ajuste_loes <- 0.3
  intervalo_confianza <- FALSE
 
  continent_selected_funcdist <- c("America","Europe")
  continent_1 <-rownames(nmds_df_plot) [nmds_df_plot$continent %in%continent_selected_funcdist[1]]
  continent_2 <- rownames(nmds_df_plot) [nmds_df_plot$continent %in%continent_selected_funcdist[2]]
  continent_distnace_matrix <- distance_matrix_for_continents[continent_1,continent_2]
  if(ncol(continent_distnace_matrix)!=nrow(continent_distnace_matrix)){
    dimension_matrix <- min(ncol(continent_distnace_matrix),nrow(continent_distnace_matrix))
    continent_distnace_matrix <- continent_distnace_matrix[1:dimension_matrix,1:dimension_matrix]
  }
  ggplot_matrix_distance_continents <- data.frame(age=nmds_df_plot$age[rownames(nmds_df_plot)%in%rownames(continent_distnace_matrix)],continent_distance=diag(continent_distnace_matrix))
  
  residuals_continent_distance <- lm(continent_distance ~ age,data=ggplot_matrix_distance_continents)
  seg_breakpoints <- segmented(residuals_continent_distance, seg.Z = ~ age,npsi = 3) 
  fitted_breakpoints <- fitted(seg_breakpoints)
  model_breakpoints <- data.frame(mean_age_breakpoints=ggplot_matrix_distance_continents$age,fit_break=fitted_breakpoints)
  lines_break_points <- seg_breakpoints$psi[,2]
  
  a <- ggplot(ggplot_matrix_distance_continents, aes(age, continent_distance)) +
    geom_point(color=alpha("#525252",0.1),size=4) +
    geom_smooth(method="loess", color=alpha("#e89a5e",0.8), fill=alpha("#e89a5e",0.3), se=intervalo_confianza, span=ajuste_loes, size=tamaño_linea) +
    theme_niwot() +
    labs(y = "Functional distance\n", x = "Time") +
    ggtitle(paste(continent_selected_funcdist[1],continent_selected_funcdist[2],sep= " vs "))+
    xlim(rev(range(nmds_df_plot$age[rownames(nmds_df_plot)%in%rownames(continent_distnace_matrix)]))) +
    ylim(0,1)+ 
    geom_vline(xintercept = lines_break_points, linetype = "dashed",color=alpha("#525252",0.7))
  
  
  
  continent_selected_funcdist <- c("America","Asia")
  continent_1 <-rownames(nmds_df_plot) [nmds_df_plot$continent %in%continent_selected_funcdist[1]]
  continent_2 <- rownames(nmds_df_plot) [nmds_df_plot$continent %in%continent_selected_funcdist[2]]
  continent_distnace_matrix <- distance_matrix_for_continents[continent_1,continent_2]
  if(ncol(continent_distnace_matrix)!=nrow(continent_distnace_matrix)){
    dimension_matrix <- min(ncol(continent_distnace_matrix),nrow(continent_distnace_matrix))
    continent_distnace_matrix <- continent_distnace_matrix[1:dimension_matrix,1:dimension_matrix]
  }
  ggplot_matrix_distance_continents <- data.frame(age=nmds_df_plot$age[rownames(nmds_df_plot)%in%rownames(continent_distnace_matrix)],continent_distance=diag(continent_distnace_matrix))
  
  residuals_continent_distance <- lm(continent_distance ~ age,data=ggplot_matrix_distance_continents)
  seg_breakpoints <- segmented(residuals_continent_distance, seg.Z = ~ age,npsi = 3) 
  fitted_breakpoints <- fitted(seg_breakpoints)
  model_breakpoints <- data.frame(mean_age_breakpoints=ggplot_matrix_distance_continents$age,fit_break=fitted_breakpoints)
  lines_break_points <- seg_breakpoints$psi[,2]
  
 b <- ggplot(ggplot_matrix_distance_continents, aes(age, continent_distance)) +
    geom_point(color=alpha("#525252",0.1),size=4) +
    geom_smooth(method="loess", color=alpha("#b76552",0.8), fill=alpha("#b76552",0.3), se=intervalo_confianza, span=ajuste_loes, size=tamaño_linea) +
    theme_niwot() +
    labs(y = "Functional distance\n", x = "Time") +
    ggtitle(paste(continent_selected_funcdist[1],continent_selected_funcdist[2],sep= " vs "))+
    xlim(rev(range(nmds_df_plot$age[rownames(nmds_df_plot)%in%rownames(continent_distnace_matrix)]))) +
    ylim(0,1)+ 
    geom_vline(xintercept = lines_break_points, linetype = "dashed",color=alpha("#525252",0.7))
 
 
  continent_selected_funcdist <- c("America","Africa")
  continent_1 <-rownames(nmds_df_plot) [nmds_df_plot$continent %in%continent_selected_funcdist[1]]
  continent_2 <- rownames(nmds_df_plot) [nmds_df_plot$continent %in%continent_selected_funcdist[2]]
  continent_distnace_matrix <- distance_matrix_for_continents[continent_1,continent_2]
  if(ncol(continent_distnace_matrix)!=nrow(continent_distnace_matrix)){
    dimension_matrix <- min(ncol(continent_distnace_matrix),nrow(continent_distnace_matrix))
    continent_distnace_matrix <- continent_distnace_matrix[1:dimension_matrix,1:dimension_matrix]
  }
  ggplot_matrix_distance_continents <- data.frame(age=nmds_df_plot$age[rownames(nmds_df_plot)%in%rownames(continent_distnace_matrix)],continent_distance=diag(continent_distnace_matrix))
  
  residuals_continent_distance <- lm(continent_distance ~ age,data=ggplot_matrix_distance_continents)
  seg_breakpoints <- segmented(residuals_continent_distance, seg.Z = ~ age,npsi = 1) 
  fitted_breakpoints <- fitted(seg_breakpoints)
  model_breakpoints <- data.frame(mean_age_breakpoints=ggplot_matrix_distance_continents$age,fit_break=fitted_breakpoints)
  lines_break_points <- seg_breakpoints$psi[,2]
  
  c <- ggplot(ggplot_matrix_distance_continents, aes(age, continent_distance)) +
    geom_point(color=alpha("#525252",0.1),size=4) +
    geom_smooth(method="loess", color=alpha("#898670",0.8), fill=alpha("#898670",0.3), se=intervalo_confianza, span=ajuste_loes, size=tamaño_linea) +
    theme_niwot() +
    labs(y = "Functional distance\n", x = "Time") +
    ggtitle(paste(continent_selected_funcdist[1],continent_selected_funcdist[2],sep= " vs "))+
    xlim(rev(range(nmds_df_plot$age[rownames(nmds_df_plot)%in%rownames(continent_distnace_matrix)]))) +
    ylim(0,1)+ 
    geom_vline(xintercept = lines_break_points, linetype = "dashed",color=alpha("#525252",0.7))
  
  
  continent_selected_funcdist <- c("Europe","Asia")
  continent_1 <-rownames(nmds_df_plot) [nmds_df_plot$continent %in%continent_selected_funcdist[1]]
  continent_2 <- rownames(nmds_df_plot) [nmds_df_plot$continent %in%continent_selected_funcdist[2]]
  continent_distnace_matrix <- distance_matrix_for_continents[continent_1,continent_2]
  if(ncol(continent_distnace_matrix)!=nrow(continent_distnace_matrix)){
    dimension_matrix <- min(ncol(continent_distnace_matrix),nrow(continent_distnace_matrix))
    continent_distnace_matrix <- continent_distnace_matrix[1:dimension_matrix,1:dimension_matrix]
  }
  ggplot_matrix_distance_continents <- data.frame(age=nmds_df_plot$age[rownames(nmds_df_plot)%in%rownames(continent_distnace_matrix)],continent_distance=diag(continent_distnace_matrix))
  
  residuals_continent_distance <- lm(continent_distance ~ age,data=ggplot_matrix_distance_continents)
  seg_breakpoints <- segmented(residuals_continent_distance, seg.Z = ~ age,npsi = 3) 
  fitted_breakpoints <- fitted(seg_breakpoints)
  model_breakpoints <- data.frame(mean_age_breakpoints=ggplot_matrix_distance_continents$age,fit_break=fitted_breakpoints)
  lines_break_points <- seg_breakpoints$psi[,2]
  
  d <- ggplot(ggplot_matrix_distance_continents, aes(age, continent_distance)) +
    geom_point(color=alpha("#525252",0.1),size=4) +
    geom_smooth(method="loess", color=alpha("#b8905e",0.8), fill=alpha("#b8905e",0.3), se=intervalo_confianza, span=ajuste_loes, size=tamaño_linea) +
    theme_niwot() +
    labs(y = "Functional distance\n", x = "Time") +
    ggtitle(paste(continent_selected_funcdist[1],continent_selected_funcdist[2],sep= " vs "))+
    xlim(rev(range(nmds_df_plot$age[rownames(nmds_df_plot)%in%rownames(continent_distnace_matrix)]))) +
    ylim(0,1)+ 
    geom_vline(xintercept = lines_break_points, linetype = "dashed",color=alpha("#525252",0.7))
  
  
  continent_selected_funcdist <- c("Europe","Africa")
  continent_1 <-rownames(nmds_df_plot) [nmds_df_plot$continent %in%continent_selected_funcdist[1]]
  continent_2 <- rownames(nmds_df_plot) [nmds_df_plot$continent %in%continent_selected_funcdist[2]]
  continent_distnace_matrix <- distance_matrix_for_continents[continent_1,continent_2]
  if(ncol(continent_distnace_matrix)!=nrow(continent_distnace_matrix)){
    dimension_matrix <- min(ncol(continent_distnace_matrix),nrow(continent_distnace_matrix))
    continent_distnace_matrix <- continent_distnace_matrix[1:dimension_matrix,1:dimension_matrix]
  }
  ggplot_matrix_distance_continents <- data.frame(age=nmds_df_plot$age[rownames(nmds_df_plot)%in%rownames(continent_distnace_matrix)],continent_distance=diag(continent_distnace_matrix))
  
  residuals_continent_distance <- lm(continent_distance ~ age,data=ggplot_matrix_distance_continents)
  seg_breakpoints <- segmented(residuals_continent_distance, seg.Z = ~ age,npsi = 2) 
  fitted_breakpoints <- fitted(seg_breakpoints)
  model_breakpoints <- data.frame(mean_age_breakpoints=ggplot_matrix_distance_continents$age,fit_break=fitted_breakpoints)
  lines_break_points <- seg_breakpoints$psi[,2]
  
 e <- ggplot(ggplot_matrix_distance_continents, aes(age, continent_distance)) +
    geom_point(color=alpha("#525252",0.1),size=4) +
    geom_smooth(method="loess", color=alpha("#8ab17d",0.8), fill=alpha("#8ab17d",0.3), se=intervalo_confianza, span=ajuste_loes, size=tamaño_linea) +
    theme_niwot() +
    labs(y = "Functional distance\n", x = "Time") +
    ggtitle(paste(continent_selected_funcdist[1],continent_selected_funcdist[2],sep= " vs "))+
    xlim(rev(range(nmds_df_plot$age[rownames(nmds_df_plot)%in%rownames(continent_distnace_matrix)]))) +
    ylim(0,1)+ 
    geom_vline(xintercept = lines_break_points, linetype = "dashed",color=alpha("#525252",0.7))
 
 
  continent_selected_funcdist <- c("Asia","Africa")
  continent_1 <-rownames(nmds_df_plot) [nmds_df_plot$continent %in%continent_selected_funcdist[1]]
  continent_2 <- rownames(nmds_df_plot) [nmds_df_plot$continent %in%continent_selected_funcdist[2]]
  continent_distnace_matrix <- distance_matrix_for_continents[continent_1,continent_2]
  if(ncol(continent_distnace_matrix)!=nrow(continent_distnace_matrix)){
    dimension_matrix <- min(ncol(continent_distnace_matrix),nrow(continent_distnace_matrix))
    continent_distnace_matrix <- continent_distnace_matrix[1:dimension_matrix,1:dimension_matrix]
  }
  ggplot_matrix_distance_continents <- data.frame(age=nmds_df_plot$age[rownames(nmds_df_plot)%in%rownames(continent_distnace_matrix)],continent_distance=diag(continent_distnace_matrix))
  
  residuals_continent_distance <- lm(continent_distance ~ age,data=ggplot_matrix_distance_continents)
  seg_breakpoints <- segmented(residuals_continent_distance, seg.Z = ~ age,npsi = 3) 
  fitted_breakpoints <- fitted(seg_breakpoints)
  model_breakpoints <- data.frame(mean_age_breakpoints=ggplot_matrix_distance_continents$age,fit_break=fitted_breakpoints)
  lines_break_points <- seg_breakpoints$psi[,2]
  
  f <- ggplot(ggplot_matrix_distance_continents, aes(age, continent_distance)) +
    #geom_line( color=alpha("#90be6d",0.8),size=3)+
    geom_point(color=alpha("#525252",0.1),size=4) +
    geom_smooth(method="loess", color=alpha("#597c71",0.8), fill=alpha("#597c71",0.3), se=intervalo_confianza, span=ajuste_loes, size=tamaño_linea) +
    theme_niwot() +
    labs(y = "Functional distance\n", x = "Time") +
    ggtitle(paste(continent_selected_funcdist[1],continent_selected_funcdist[2],sep= " vs "))+
    xlim(rev(range(nmds_df_plot$age[rownames(nmds_df_plot)%in%rownames(continent_distnace_matrix)]))) +
    ylim(0,1)+ 
    geom_vline(xintercept = lines_break_points, linetype = "dashed",color=alpha("#525252",0.7))
  
 


  font_family <- "sans"
  require(gridExtra)
  quartz(width = cm(8),height = cm(5),family = font_family)
  grid.arrange(arrangeGrob(a,b,c,ncol = 3),arrangeGrob(d,e,f, ncol = 3), ncol=1)
  
########### Module Characterization #########
  
  network_stats_table_FT <- network_stats_table[network_stats_table$ID %in% rownames(occ_analysis_functional_types),]

  
  
  #### FT richness of modules
  
  network_stats_table_FT$fe_richness <- rowSums(occ_analysis_functional_types)[match(network_stats_table_FT$ID,rownames(occ_analysis_functional_types))]

  network_stats_table_FT$FT_number <- as.numeric(sapply(strsplit(network_stats_table_FT$ID,"_"),function(x)x[2]))
  network_stats_table_FT<- network_stats_table_FT[order(network_stats_table_FT$FT_number),]
  network_stats_table_FT$fe_diversity <- diversity(occ_analysis_functional_types,"shannon")
  
  network_stats_table_FT$old_module <- as.numeric(network_stats_table_FT$module)
  network_stats_table_FT<- network_stats_table_FT[order(network_stats_table_FT$module),]
  
  length(functional_types)
  

  dim(network_table_no_singletone)
  
  ####Plot first to have the correct module number ordered by age RUN MODULE PLOT 
 
  bins_loc$mean_age <-as.numeric(bin_age_list)[match(bins_loc$bin,names(bin_age_list))]
  order_age <- aggregate(bins_loc$mean_age,by=list(bins_loc$module), median)
  
  module_number_changes <- aggregate(bins_loc$module,by=list(bins_loc$old_module), median)
  
  network_stats_table_FT$module <- module_number_changes$x[match(network_stats_table_FT$old_module,as.numeric(module_number_changes$Group.1))]
  
  
  #### Most abundant traits in the highest indval FT per module

  ##we withdraw module one because only has 1 FT
 
  network_table_no_singletone <- network_stats_table_FT[network_stats_table_FT$module!=1,]
  
  data_traits_modules <- NULL


  
  for (i in 1:length(unique(network_table_no_singletone$module))+1) {
    

    module_1 <- network_table_no_singletone[network_table_no_singletone$module==i,]
    module_1_indval <- module_1[order(module_1$indval, decreasing = T),]
    
    most_representative_FT_module<- module_1[module_1$indval>=quantile(module_1_indval$indval,probs=0.95),]
    sp_FT_representatives_modules <- names(functional_types)[functional_types %in% most_representative_FT_module$FT_number]
    traits_sp_FT_representatives_modules<- no_na_fdata_art[rownames(no_na_fdata_art) %in% sp_FT_representatives_modules,]
    
    table_stats_traits <- lapply(traits_sp_FT_representatives_modules, function(x)table(x))
    
    data_traits_modules[[i]]<- c(table_stats_traits)
  
    
  }
  

  
  ##we plot the results of the most abundant traits
  
  
trait_palette <- hcl.colors(14,"viridis")
  
  for (z in colnames(no_na_fdata_art) ) {
    quartz(width = 20, height = 10)
    par(mfcol=c(3,3), oma = c(0, 0, 0, 4),family = font_family)
      for (trait_names in z) {
      for (i in 1:length(unique(network_table_no_singletone$module))+1){
        selected_trait_to_plot <-  data_traits_modules[[i]][names(data_traits_modules[[i]])==trait_names]
        plot(selected_trait_to_plot[[1]],lwd=20,col=trait_palette[colnames(no_na_fdata_art)==z],main=paste(c("Module",i-1)),bty="n",
             ylab = names(selected_trait_to_plot),col.axis="#525252",col.lab="#525252",col.main="#525252")
      }
    }
  }
  
  
  



  ##### PLot diversity throuhg time 
  
  
  modules_plot_lines <- unique(bins_loc$module[order(bins_loc$module,decreasing = F)])
  module_cols <-  hcl.colors(10,"viridis")
  bins_loc$module_color_lines <- module_cols[bins_loc$module]
  diversity_plot <- bins_loc[order(bins_loc$mean_age,decreasing = T),]
  diversity_plot$fe_diversity <- round(diversity_plot$fe_diversity,2)

  ggplot_diversity_plot <- diversity_plot
  ggplot_diversity_plot <- ggplot_diversity_plot[order(ggplot_diversity_plot$module,decreasing = F),]
  ggplot_diversity_plot$module <- as.character(ggplot_diversity_plot$module)
  level_order_diversity_plot <- unique(ggplot_diversity_plot$module)
  ggplot_diversity_plot$evenness <- round(ggplot_diversity_plot$fe_diversity/log(ggplot_diversity_plot$fe_richness),2)
  ggplot_diversity_plot$evenness[is.na(ggplot_diversity_plot$evenness)] <- 0
  

 corrected_ggplot_diversity <- ggplot_diversity_plot
 corrected_ggplot_diversity$occ <- log(colSums(occ_occurences)[match(rownames(corrected_ggplot_diversity),colnames(occ_occurences))])
 corrected_ggplot_diversity$normal_occ <- colSums(occ_occurences)[match(rownames(corrected_ggplot_diversity),colnames(occ_occurences))]
 corrected_ggplot_diversity$log_richness <- log(corrected_ggplot_diversity$fe_richness)
 corrected_ggplot_diversity$quadratic_occ <- corrected_ggplot_diversity$occ^2
 corrected_ggplot_diversity <- corrected_ggplot_diversity[corrected_ggplot_diversity$sp_richness>20,]
 
 
 corrected_ggplot_diversity$number_continent[corrected_ggplot_diversity$continent=="America"] <- 1
 corrected_ggplot_diversity$number_continent[corrected_ggplot_diversity$continent=="Europe"] <- 2
 corrected_ggplot_diversity$number_continent[corrected_ggplot_diversity$continent=="Africa"] <- 3
 corrected_ggplot_diversity$number_continent[corrected_ggplot_diversity$continent=="Asia"] <- 4
 continent_cols <- c("#e76f51","#e9c46a","#2a9d8f","#875B52") 
 corrected_ggplot_diversity$continent_color_lines <- continent_cols[corrected_ggplot_diversity$number_continent]
 level_order_diversity_continent <- unique(corrected_ggplot_diversity$continent)
 
 rownames(corrected_ggplot_diversity)
 

 
 ###RAW DIVERSITY
   a <- ggplot(corrected_ggplot_diversity, aes(mean_age, fe_diversity)) +
    geom_point(color=alpha("#525252",0.1),size=4) +
    geom_smooth(method="loess", color=alpha("#90be6d",0.8), fill=alpha("#90be6d",0.3), se=TRUE) +
    theme_niwot() +
    labs(y = "Raw Shannon diversity", x = "Time") +
    xlim(rev(range(corrected_ggplot_diversity$mean_age))) +
    ylim(0,max(corrected_ggplot_diversity$fe_diversity))
  
   b <- ggplot(corrected_ggplot_diversity, aes(mean_age, evenness)) +
     geom_point(color=alpha("#525252",0.1),size=4) +
     geom_smooth(method="loess", color=alpha("#90AACB",0.8), fill=alpha("#90AACB",0.3), se=TRUE) +
     theme_niwot() +
     labs(y = "Raw Evenness diversity", x = "Time") +
     xlim(rev(range(corrected_ggplot_diversity$mean_age))) +
     ylim(0,max(corrected_ggplot_diversity$evenness))
  
  c <- ggplot(corrected_ggplot_diversity, aes(mean_age, log_richness)) +
    geom_point(color=alpha("#525252",0.1),size=4) +
    geom_smooth(method="loess", color=alpha("#f9c74f",0.8), fill=alpha("#f9c74f",0.3), se=TRUE) +
    theme_niwot() +
    labs(y = "Raw Functional richness", x = "Time") +
    xlim(rev(range(corrected_ggplot_diversity$mean_age))) +
    ylim(0,max(corrected_ggplot_diversity$log_richness))
  
  
  ###OCC per continent
  

  shannon_lm_continent <- lm(fe_diversity ~ occ*continent,data=corrected_ggplot_diversity)
  corrected_ggplot_diversity$residuals_shannon <- shannon_lm_continent$residuals
  
  d <- ggplot(corrected_ggplot_diversity, aes(occ, fe_diversity,colour=factor(continent))) +
    geom_point(size=4, alpha=0.2) +
    stat_smooth(method = "lm",se=FALSE,alpha=0.05)+
    theme_niwot() +
    theme(legend.position = "none",legend.key.size = unit(0.0001,"cm"))+
    labs(y = "Shannon diversity", x = "Occurences") +
    xlim(min(corrected_ggplot_diversity$occ),max(corrected_ggplot_diversity$occ)) +
    ylim(min(corrected_ggplot_diversity$fe_diversity),max(corrected_ggplot_diversity$fe_diversity))+
    scale_color_manual(values=continent_cols)
  
  
  evenness_lm_continent <- lm(evenness ~ occ*continent,data=corrected_ggplot_diversity)
  corrected_ggplot_diversity$residuals_evenness <- evenness_lm_continent$residuals
  
  summary(evenness_lm_continent)
  
  e <- ggplot(corrected_ggplot_diversity, aes(occ, evenness,colour=factor(continent))) +
    geom_point(size=4, alpha=0.2) +
    stat_smooth(method = "lm",se=FALSE,alpha=0.05)+
    theme_niwot() +
    theme(legend.position = "none",legend.key.size = unit(0.0001,"cm"))+
    labs(y = "Evenness diversity", x = "Occurences") +
    xlim(min(corrected_ggplot_diversity$occ),max(corrected_ggplot_diversity$occ)) +
    ylim(min(corrected_ggplot_diversity$evenness),max(corrected_ggplot_diversity$evenness))+
    scale_color_manual(values=continent_cols)
  
  
  richness_lm_continent <- lm(fe_richness ~ occ*continent,data=corrected_ggplot_diversity)
  corrected_ggplot_diversity$residuals_richness <- richness_lm_continent$residuals
  
  f <- ggplot(corrected_ggplot_diversity, aes(occ, fe_richness,colour=factor(continent))) +
    geom_point(size=4, alpha=0.2) +
    stat_smooth(method = "lm",se=F,alpha=0.05)+
    theme_niwot() +
    theme(legend.position = "right",legend.key.size = unit(0.0001,"cm"))+
    labs(y = "Functional richnness", x = "Occurences") +
    xlim(min(corrected_ggplot_diversity$occ),max(corrected_ggplot_diversity$occ)) +
    ylim(min(corrected_ggplot_diversity$fe_richness),max(corrected_ggplot_diversity$fe_richness))+
    scale_color_manual(values=continent_cols)
  
  
  
  ###RESIDUALS diversity total
  
  residuals_lm_shannon <- lm(residuals_shannon ~ mean_age,data=corrected_ggplot_diversity)
  seg_breakpoints <- segmented(residuals_lm_shannon, seg.Z = ~ mean_age, npsi = 3) 
  fitted_breakpoints <- fitted(seg_breakpoints)
  model_breakpoints <- data.frame(mean_age_breakpoints=corrected_ggplot_diversity$mean_age,fit_break=fitted_breakpoints)
  lines_break_points <- seg_breakpoints$psi[,2]

  g <- ggplot(corrected_ggplot_diversity, aes(mean_age, residuals_shannon)) +
    geom_point(color=alpha("#525252",0.1),size=4) +
    geom_smooth(method="loess", color=alpha("#90be6d",0.8), fill=alpha("#90be6d",0.3), se=TRUE) +
    theme_niwot() +
    labs(y = "Functional Shannon", x = "Time") +
    xlim(rev(range(corrected_ggplot_diversity$mean_age))) +
    ylim(min(corrected_ggplot_diversity$residuals_shannon),max(corrected_ggplot_diversity$residuals_shannon))+ 
    geom_line(data = model_breakpoints, aes(x = mean_age_breakpoints, y = fit_break), colour = "tomato")+
    geom_vline(xintercept = lines_break_points, linetype = "dashed")
  
  
  residuals_lm_evenness <- lm(residuals_evenness ~ mean_age,data=corrected_ggplot_diversity)
  seg_breakpoints <- segmented(residuals_lm_evenness, seg.Z = ~ mean_age, npsi = 1) 
  fitted_breakpoints <- fitted(seg_breakpoints)
  model_breakpoints <- data.frame(mean_age_breakpoints=corrected_ggplot_diversity$mean_age,fit_break=fitted_breakpoints)
  lines_break_points <- seg_breakpoints$psi[, 2]
  
  h <- ggplot(corrected_ggplot_diversity, aes(mean_age, residuals_evenness)) +
    geom_point(color=alpha("#525252",0.1),size=4) +
    geom_smooth(method="loess", color=alpha("#90AACB",0.8), fill=alpha("#90AACB",0.3), se=TRUE) +
    theme_niwot() +
    labs(y = "Functional Evenness", x = "Time") +
    xlim(rev(range(corrected_ggplot_diversity$mean_age))) +
    ylim(min(corrected_ggplot_diversity$residuals_evenness),max(corrected_ggplot_diversity$residuals_evenness))+ 
    geom_vline(xintercept = lines_break_points, linetype = "dashed")
  

  richnness_lm_residualss <- lm(residuals_richness ~ mean_age,data=corrected_ggplot_diversity)
  seg_breakpoints <- segmented(richnness_lm_residualss, seg.Z = ~ mean_age, npsi = 1) 
  fitted_breakpoints <- fitted(seg_breakpoints)
  model_breakpoints <- data.frame(mean_age_breakpoints=corrected_ggplot_diversity$mean_age,fit_break=fitted_breakpoints)
  lines_break_points <- seg_breakpoints$psi[, 2]
  
  i <- ggplot(corrected_ggplot_diversity, aes(mean_age, residuals_richness)) +
    geom_point(color=alpha("#525252",0.1),size=4) +
    geom_smooth(method="loess", color=alpha("#f9c74f",0.8), fill=alpha("#f9c74f",0.3), se=TRUE) +
    theme_niwot() +
    labs(y = "Functional richnness", x = "Time") +
    xlim(rev(range(corrected_ggplot_diversity$mean_age))) +
    ylim(min(corrected_ggplot_diversity$residuals_richness),max(corrected_ggplot_diversity$residuals_richness))+
    #geom_line(data = model_breakpoints, aes(x = mean_age_breakpoints, y = fit_break), colour = "tomato")+
    geom_vline(xintercept = lines_break_points, linetype = "dashed")
  

  
  ##Diversity per module

  module_cols_diversity <- c("#9e0142","#d53e4f","#f46d43","#fdae61","#fee08b","#e6f598","#abdda4"#,"#66c2a5"
                             ,"#3288bd","#5e4fa2")

  level_order_diversity_plot <- rev(unique(corrected_ggplot_diversity$module))
  
 
  
  j <- ggplot(corrected_ggplot_diversity, aes(mean_age, residuals_shannon,colour=factor(module,levels=level_order_diversity_plot))) +
    geom_point(size=4,alpha=0.2) +
    geom_smooth(method = "lm", se=F,alpha=0.7) +
    theme_niwot() +
    theme(legend.position = "none",legend.key.size = unit(0.05,"cm"))+
    labs(y = "Modules Shannon", x = "Time") +
    xlim(rev(range(corrected_ggplot_diversity$mean_age))) +
    ylim(min(corrected_ggplot_diversity$residuals_shannon),max(corrected_ggplot_diversity$residuals_shannon))+
    scale_color_manual(values=rev(module_cols_diversity))
  
  k <- ggplot(corrected_ggplot_diversity, aes(mean_age, residuals_evenness,colour=factor(module,levels=level_order_diversity_plot))) +
    geom_point(size=4,alpha=0.2) +
    geom_smooth(method = "lm", se=F,alpha=0.7) +
    theme_niwot() +
    theme(legend.position = "none",legend.key.size = unit(0.05,"cm"))+
    labs(y = "Modules Evenness", x = "Time") +
    xlim(rev(range(corrected_ggplot_diversity$mean_age))) +
    ylim(min(corrected_ggplot_diversity$residuals_evenness),max(corrected_ggplot_diversity$residuals_evenness))+
    scale_color_manual(values=rev(module_cols_diversity))
  
  
  
  l <- ggplot(corrected_ggplot_diversity, aes(mean_age, residuals_richness,colour=factor(module,levels=level_order_diversity_plot))) +
    geom_point(size=4,alpha=0.2) +
    geom_smooth(method = "lm", se=F,alpha=0.7) +
    theme_niwot() +
    theme(legend.position = "right",legend.key.size = unit(0.05,"cm"))+
    labs(y = "Modules richness", x = "Time") +
    xlim(rev(range(corrected_ggplot_diversity$mean_age))) +
    ylim(min(corrected_ggplot_diversity$residuals_richness),max(corrected_ggplot_diversity$residuals_richness))+
    scale_color_manual(values=rev(module_cols_diversity))
  

  ###Diversity by continents
  
  

  m <- ggplot(corrected_ggplot_diversity, aes(mean_age, residuals_shannon,colour=factor(continent,levels=level_order_diversity_continent),fill=factor(continent,levels=level_order_diversity_continent))) +
    geom_point(size=4,alpha=0.2) +
    geom_smooth(method = "loess", se=T,alpha=0.05) +
    theme_niwot() +
    theme(legend.position = "none",legend.key.size = unit(0.0001,"cm"))+
    labs(y = "Continent Shannon", x = "Time") +
    xlim(rev(range(corrected_ggplot_diversity$mean_age))) +
    ylim(min(corrected_ggplot_diversity$residuals_shannon),max(corrected_ggplot_diversity$residuals_shannon))+
    scale_color_manual(values=continent_cols)
  
  n <- ggplot(corrected_ggplot_diversity, aes(mean_age, residuals_evenness,colour=factor(continent,levels=level_order_diversity_continent),fill=factor(continent,levels=level_order_diversity_continent))) +
    geom_point(size=4,alpha=0.2) +
    geom_smooth(method = "loess", se=T,alpha=0.05) +
    theme_niwot() +
    theme(legend.position = "none",legend.key.size = unit(0.0001,"cm"))+
    labs(y = "Continent Evenness", x = "Time") +
    xlim(rev(range(corrected_ggplot_diversity$mean_age))) +
    ylim(min(corrected_ggplot_diversity$residuals_evenness),max(corrected_ggplot_diversity$residuals_evenness))+
    scale_color_manual(values=continent_cols)
  
  
  
  o <- ggplot(corrected_ggplot_diversity, aes(mean_age, residuals_richness,colour=factor(continent,levels=level_order_diversity_continent),fill=factor(continent,levels=level_order_diversity_continent))) +
    geom_point(size=4,alpha=0.2) +
    geom_smooth(method = "loess", se=T,alpha=0.05) +
    theme_niwot() +
    theme(legend.position = "none",legend.key.size = unit(0.0001,"cm"))+
    labs(y = "Continent richness", x = "Time") +
    xlim(rev(range(corrected_ggplot_diversity$mean_age))) +
    ylim(min(corrected_ggplot_diversity$residuals_richness),max(corrected_ggplot_diversity$residuals_richness))+
    scale_color_manual(values=continent_cols)
  
  
  
  
  font_family <- "sans"
  require(gridExtra)
  quartz(width = cm(8),height = cm(5),family = font_family)
  grid.arrange(arrangeGrob(a,b,c,ncol = 3),arrangeGrob(d,e,f, ncol = 3),arrangeGrob(g,h,i,ncol=3),arrangeGrob(j,k,l,ncol=3),arrangeGrob(m,n,o,ncol=3), ncol=1)
  
  

  font_family <- "sans"
  require(gridExtra)
  quartz(width = cm(8),height = cm(2),family = font_family)
  grid.arrange(arrangeGrob(g,h,i,ncol=3), ncol=1)
  
  

  font_family <- "sans"
  require(gridExtra)
  quartz(width = cm(8),height = cm(3),family = font_family)
  grid.arrange(arrangeGrob(g,h,i,ncol=3),arrangeGrob(m,n,o,ncol=3), ncol=1)
  
  

  font_family <- "sans"
  require(gridExtra)
  quartz(width = cm(8),height = cm(3),family = font_family)
  grid.arrange(arrangeGrob(a,b,c,ncol = 3),arrangeGrob(d,e,f, ncol = 3), ncol=1)
  
  
  
  
  
  ###AICc to compare which model fits better to explain occ
  shannon_lm_quadratic <- lm(fe_diversity ~ occ*quadratic_occ,data=corrected_ggplot_diversity)
  shannon_lm_continent <- lm(fe_diversity ~ occ*continent,data=corrected_ggplot_diversity)
  models <- list(shannon_lm_continent, shannon_lm_quadratic)
  
  #specify model names
  mod.names <- c("shannon_lm_continent", "shannon_lm_quadratic")
  
  #calculate AIC of each model
  aictab(cand.set = models, modnames = mod.names)
  
  
  
 ####PERMANOVA to test the accuracy of network partiontion compared to the functional distances
  
  ##We take the same bins as in the turnover, which has a "subsample"
  permanova_bins_loc <- bins_loc[rownames(bins_loc) %in% rownames(turnover_matrix),]
  permanova_bins_loc <- permanova_bins_loc[order(permanova_bins_loc$module, decreasing = T),]
  
  permanova_modules_distances <- adonis(turnover_matrix ~ module,data = permanova_bins_loc)
 

  permanova_modules_distances$aov.tab
  

  
  
  
  



  