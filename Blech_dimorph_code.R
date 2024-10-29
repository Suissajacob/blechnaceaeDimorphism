#This script will reproduce figures from Suissa and Smith 2024 Evol. 

#Load libraries
library(tidyverse)
library(phytools)
library(ape)
library(geiger)
library(phylotools)
library(factoextra)
library(FactoMineR)
library(ggfortify)
library(ggsci)

##### 
#Part 1 read in the quantitative data and analyze dimorphic leaves in and build PCA
#####

#read in dimorphic_4.df

dimorphic_4.df<- read.csv(path/to/merged_dimorphic_v4.csv)


#subset the data
df <- dimorphic_4.df%>%
  select(-state)%>%
  select(-stat2)%>%
  column_to_rownames(var= "species")

# for the sake of the data group species by moomorphic an non-monomorphic
dimorphic_3.5.df<- dimorphic_3.5.df%>%
  mutate(., stat2 = ifelse(state == "Mono" , "Mono", "Non-mono"))


df_pca<- prcomp(df, scale. = TRUE, center = TRUE)

summary(df_pca)

##This will give us the loadings matrix
df_pca$rotation
df_pca$x
df_pca$scale

fviz_pca_biplot(df_pca, axes = c(1, 2), label ="all",  palette = "Dark2",  habillage=dimorphic_4.df$stat2, pointsize = 3, invisible="quali" , addEllipses=FALSE, ellipse.level=0.95, repel= TRUE)+theme_bw()

fviz_pca_biplot(df_pca, axes = c(1, 2),  label = "off",  palette = "Dark2",  habillage=dimorphic_4.df$stat2, pointsize = 3, invisible="quali" , addEllipses=FALSE, ellipse.level=0.95, repel= TRUE)+theme_bw()


#################
# visualize eigenvalues (scree plot). Show the percentage of variances explained by each principal component.
#################


fviz_eig(df_pca)


fviz_pca_var(df_pca,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)


#Test clustering

set.seed(123)

# function to compute total within-cluster sum of square 
wss <- function(k) {
  kmeans(scores, k, nstart = 10 )$tot.withinss
}

# Compute and plot wss for k = 1 to k = 15
k.values <- 1:15

# extract wss for 2-15 clusters
wss_values <- map_dbl(k.values, wss)

plot(k.values, wss_values,
     type="b", pch = 19, frame = FALSE, 
     xlab="Number of clusters K",
     ylab="Total within-clusters sum of squares")



############################################################ 
#Phylo pca
############################################################ 

#bring in modified data

df<- read.csv(path/to/merged_dimorphic_v4.csv, row.names = 1)

# bring in tree
tree<- read.tree(path/to/blechnaceae_Tree_4pnt_lam0.tre)


# check names in tree and in dataframe
chk<-name.check(tree,df)
summary(chk)

#check which are in data nand not tree
chk$data_not_tree

# now prune the tre based on the data available
tree.pruned<-drop.tip(tree,chk$tree_not_data)
fern.Data.pruned<-df[!(rownames(df)%in%
                         chk$data_not_tree),,drop=FALSE]

fern.Data.pruned1<- fern.Data.pruned%>%
  select(!c(stat2, state))


#Run phylo pca
dimPhyPCA<- phyl.pca(tree.pruned,fern.Data.pruned1, method = "lambda", mode = "corr")

dd<- dimPhyPCA$S[,1:2]

dd<- cbind(dd, fern.Data.pruned)

dd <- dd%>%
  select(stat2, PC1, PC2)

dd2<- dimPhyPCA$S[,2:3]

dd2<- cbind(dd2, fern.Data.pruned)

dd2 <- dd2%>%
  select(stat2, PC2, PC3)

#generate the variance explained by each PC
diag(dimPhyPCA$Eval)/sum(dimPhyPCA$Eval)*100

cols<-setNames(c("#57768F","#F9C750"),c("Mono", "Non-mono"))


# run a phylomorphospace without colors
phylomorphospace(X = dd[,2:3], colors=cols, tree = dim_simmap,node.by.map=TRUE, label = "off", control=list(lwd=2), node.size=c(0,1.2))


#DD2

phylomorphospace(X = dd2[,2:3], colors=cols, tree = dim_simmap,node.by.map=TRUE, label = "off", control=list(lwd=2), node.size=c(0,1.2))



##############################
#Model fitting
##############################
# Read in the tree
mytimetree_0.rel <- read.tree("/Users/jacobsuissa/Documents/Students/Makaleh_smith/Dimorphism_project/trees/Time_Trees/blechnaceae_Tree_4pnt_lam0.tre")

# Bring in the data
leaf.dat <- read.csv("/Users/jacobsuissa/Documents/Students/Makaleh_smith/Dimorphism_project/CSVs/Dimorphism_data/dimorphism_4k_species_Blechnaceae.csv")

# Drop NAs and select just the data of interest
leaf.dat_NAomit <- leaf.dat %>%
  filter(!is.na(dimorphism)) %>%
  dplyr::select(species_name, dimorphism)

# Change names
leaf.dat_NAomit$dimorphism <- str_replace_all(leaf.dat_NAomit$dimorphism, "monomorphic", "A")
leaf.dat_NAomit$dimorphism <- str_replace_all(leaf.dat_NAomit$dimorphism, "NELD", "C")
leaf.dat_NAomit$dimorphism <- str_replace_all(leaf.dat_NAomit$dimorphism, "ELD", "B")

leaf.dat_NAomit <- leaf.dat_NAomit %>%
  arrange(dimorphism)

# Filter tips in the tree to just match tips in the dataframe
rownames(leaf.dat_NAomit) <- leaf.dat_NAomit[, 1]

# Check species in tree
chk <- name.check(mytimetree_0.rel, leaf.dat_NAomit)

# Now prune the tree based on the data available
tree.pruned <- drop.tip(mytimetree_0.rel, chk$tree_not_data)
fern.Data.pruned <- leaf.dat_NAomit[!(rownames(leaf.dat_NAomit) %in% chk$data_not_tree), , drop = FALSE]

# Run ancestral reconstruction stochastic character mapping for discrete characters
# Format the data for stochastic character mapping
leaf.mode <- setNames(fern.Data.pruned[, 2], rownames(fern.Data.pruned))


#ARD mimic
ARD_mimic <- matrix(c(0, 1, 4,
                      2, 0, 0,
                      3, 0, 0),
                    3,3,byrow = TRUE,
                    dimnames = list(0:2, 0:2))


#set data and  tree
rownames(fern.Data.pruned) <- NULL


corHMM_output_R1_ARD_mimic <- corHMM(phy = tree.pruned, data = fern.Data.pruned, rate.cat = 1, rate.mat =ARD_mimic, node.states = "marginal", fixed.nodes=FALSE, get.tip.states =TRUE, root.p=c(1,0,0), nstarts=10, upper.bound=500)


#Precurser ARD

ordered <- matrix(c(0, 1, 0,
                    2, 0, 3,
                    0, 4, 0),
                  3,3, byrow = TRUE,
                  dimnames = list(0:2, 0:2))

corHMM_output_R1_ordered <- corHMM(phy = tree.pruned, data = fern.Data.pruned, rate.cat = 1, rate.mat =ordered, node.states = "marginal", fixed.nodes=FALSE, get.tip.states =TRUE, root.p=c(1,0,0), nstarts=10, upper.bound=500)



#Precurser ER
ordered_ER <- matrix(c(0, 1, 0,
                       1, 0, 1,
                       0, 1, 0),
                     3,3,byrow = TRUE,
                     dimnames = list(0:2, 0:2))


corHMM_output_R1_ordered_ER <- corHMM(phy = tree.pruned, data = fern.Data.pruned, rate.cat = 1, rate.mat =ordered_ER, node.states = "marginal", fixed.nodes=FALSE, get.tip.states =TRUE, root.p=c(1,0,0), nstarts=10, upper.bound=500)




## Ordered SYM

ordered_SYM <- matrix(c(0, 1, 0,
                        1, 0, 2,
                        0, 2, 0),
                      3,3,byrow = TRUE,
                      dimnames = list(0:2, 0:2))

corHMM_output_R1_ordered_SYM <- corHMM(phy = tree.pruned, data = fern.Data.pruned, rate.cat = 1, rate.mat =ordered_SYM, node.states = "marginal", fixed.nodes=FALSE, get.tip.states =TRUE, root.p=c(1,0,0), nstarts=10, upper.bound=500)


#Single rate ER

corHMM_output_R1_full_mat_ER<- corHMM(phy = tree.pruned, data = fern.Data.pruned, rate.cat = 1, model = "ER", node.states = "marginal", fixed.nodes=FALSE, get.tip.states =TRUE, root.p=c(1,0,0), nstarts=10, upper.bound=500)



#Single_rate SYM

corHMM_output_R1_full_mat_SYM<- corHMM(phy = tree.pruned, data = fern.Data.pruned, rate.cat = 1, model = "SYM", node.states = "marginal", fixed.nodes=FALSE, get.tip.states =TRUE, root.p=c(1,0,0), nstarts=10, upper.bound=500)

#Single_rate ARD

corHMM_output_R1_full_mat_ARD<- corHMM(phy = tree.pruned, data = fern.Data.pruned, rate.cat = 1, model = "ARD", node.states = "marginal", fixed.nodes=FALSE, get.tip.states =TRUE, root.p=c(1,0,0), nstarts=10, upper.bound=500)



#Hidden rates ER

corHMM_output_R2_full_mat_ER<- corHMM(phy = tree.pruned, data = fern.Data.pruned, rate.cat = 2, model = "ER", node.states = "marginal", fixed.nodes=FALSE, get.tip.states =TRUE, root.p=c(1,0,0,0,0,0), nstarts=10, upper.bound=500)

#plot
par( mar=c(1,1,1,1))
plotMKmodel(corHMM_output_R2_full_mat_ER, display = "row", rate.cat = 2)




#Hidden rates SYM

corHMM_output_R2_full_mat_SYM<- corHMM(phy = tree.pruned, data = fern.Data.pruned, rate.cat = 2, model = "SYM", node.states = "marginal", fixed.nodes=FALSE, get.tip.states =TRUE, root.p=c(1,0,0,0,0,0), nstarts=10, upper.bound=500)

#Plot
par( mar=c(1,1,1,1))
plotMKmodel(corHMM_output_R2_full_mat_SYM, display = "row", rate.cat = 2)


#Hidden rates true ARD

corHMM_output_R2_full_mat_ARD_true<- corHMM(phy = tree.pruned, data = fern.Data.pruned, rate.cat = 2, model = "ARD", node.states = "marginal", fixed.nodes=FALSE, get.tip.states =TRUE, root.p=c(1,0,0,0,0,0), nstarts=10, upper.bound=500)

#Plot
par( mar=c(1,1,1,1))
plotMKmodel(corHMM_output_R2_full_mat_ARD_true, display = "row", rate.cat = 2)



#examine hidden rates ARD-ARD
ARD_mimic <- matrix(c(0, 1, 4,
                      2, 0, 0,
                      3, 0, 0),
                    3,3,byrow = TRUE,
                    dimnames = list(0:2, 0:2))

ARD_mimic_2 <- matrix(c(0, 1, 4,
                        2, 0, 0,
                        3, 0, 0),
                      3,3,byrow = TRUE,
                      dimnames = list(0:2, 0:2))



StateMats <- list(ARD_mimic, ARD_mimic_2) 

StateMats

RateClassMat <- getRateCatMat(2) 


full_mat_ARD_ARD<-getFullMat(StateMats, RateClassMat)


#set data and  tree
rownames(fern.Data.pruned) <- NULL


corHMM_output_R2_full_mat_ARD_ARD <- corHMM(phy = tree.pruned, data = fern.Data.pruned, rate.cat = 2, rate.mat =full_mat_ARD_ARD, node.states = "marginal", fixed.nodes=FALSE, get.tip.states =TRUE, root.p=c(1,0,0,0,0,0), nstarts=100, upper.bound=500, n.cores = 3)

#plot
par( mar=c(1,1,1,1))
plotMKmodel(corHMM_output_R2_full_mat_ARD_ARD, display = "row", rate.cat = 2, )



plotMKmodel(corHMM_output_R1_ARD_mimic, display = "row", rate.cat = 1)


#Ordered ARD hidden
# Make my own matrix
ordered <- matrix(c(0, 1, 0,
                    2, 0, 3,
                    0, 4, 0),
                  3,3, byrow = TRUE,
                  dimnames = list(0:2, 0:2))

ordered_2 <- matrix(c(0, 1, 0,
                      2, 0, 3,
                      0, 4, 0),
                    3,3, byrow = TRUE,
                    dimnames = list(0:2, 0:2))


StateMats_ARD_ordered <- list(ordered, ordered_2) 


RateClassMat_ordered <- getRateCatMat(2) 


full_mat_ordered<-getFullMat(StateMats_ARD_ordered, RateClassMat)


#set data and  tree
# rownames(fern.Data.pruned) <- NULL


corHMM_output_R2_ordered <- corHMM(phy = tree.pruned, data = fern.Data.pruned, rate.cat = 2, rate.mat =full_mat_ordered, node.states = "marginal", fixed.nodes=FALSE, get.tip.states =TRUE, root.p=c(1,0,0,0,0,0), nstarts=10, upper.bound=1000)

par( mar=c(1,1,1,1))
plotMKmodel(corHMM_output_R2_ordered, display = "row", rate.cat = 2)


#Ordered ER hidden
# Make my own matrix
ordered <- matrix(c(0, 1, 0,
                    1, 0, 1,
                    0, 1, 0),
                  3,3, byrow = TRUE,
                  dimnames = list(0:2, 0:2))

ordered_2 <- matrix(c(0, 1, 0,
                      1, 0, 1,
                      0, 1, 0),
                    3,3, byrow = TRUE,
                    dimnames = list(0:2, 0:2))


StateMats_ARD_ordered <- list(ordered, ordered_2) 


RateClassMat_ordered <- getRateCatMat(2) 


full_mat_ordered_ER<-getFullMat(StateMats_ARD_ordered, RateClassMat)


#set data and  tree
# rownames(fern.Data.pruned) <- NULL


corHMM_output_R2_ordered_ER <- corHMM(phy = tree.pruned, data = fern.Data.pruned, rate.cat = 2, rate.mat =full_mat_ordered_ER, node.states = "marginal", fixed.nodes=FALSE, get.tip.states =TRUE, root.p=c(1,0,0,0,0,0), nstarts=10, upper.bound=500)

par( mar=c(1,1,1,1))
plotMKmodel(corHMM_output_R2_ordered_ER, display = "row", rate.cat = 2)

#Ordered SYM hidden
# Make my own matrix
ordered <- matrix(c(0, 1, 0,
                    1, 0, 2,
                    0, 2, 0),
                  3,3, byrow = TRUE,
                  dimnames = list(0:2, 0:2))

ordered_2 <- matrix(c(0, 1, 0,
                      1, 0, 2,
                      0, 2, 0),
                    3,3, byrow = TRUE,
                    dimnames = list(0:2, 0:2))


StateMats_ARD_ordered <- list(ordered, ordered_2) 


RateClassMat_ordered <- getRateCatMat(2) 


full_mat_ordered_sym<-getFullMat(StateMats_ARD_ordered, RateClassMat)


#set data and  tree
# rownames(fern.Data.pruned) <- NULL


corHMM_output_R2_ordered_sym <- corHMM(phy = tree.pruned, data = fern.Data.pruned, rate.cat = 2, rate.mat =full_mat_ordered_sym, node.states = "marginal", fixed.nodes=FALSE, get.tip.states =TRUE, root.p=c(1,0,0,0,0,0), nstarts=10, upper.bound=500)
corHMM_output_R2_ordered

########
aic <- c(
  corHMM_output_R1_full_mat_ER$AIC,
  corHMM_output_R1_full_mat_SYM$AIC,
  corHMM_output_R1_full_mat_ARD$AIC,
  corHMM_output_R1_ordered$AIC,
  corHMM_output_R1_ordered_SYM$AIC,
  corHMM_output_R1_ordered_ER$AIC,
  corHMM_output_R1_ARD_mimic$AIC,
  corHMM_output_R2_full_mat_ARD_ARD$AIC,
  corHMM_output_R2_full_mat_ARD_true$AIC,
  corHMM_output_R2_ordered$AIC,
  corHMM_output_R2_ordered_sym$AIC,
  corHMM_output_R2_ordered_ER$AIC,
  corHMM_output_R2_full_mat_ER$AIC,
  corHMM_output_R2_full_mat_SYM$AIC
)

model.fit <- data.frame(
  row.names =  c(
    "corHMM_output_R1_full_mat_ER",
    "corHMM_output_R1_full_mat_SYM",
    "corHMM_output_R1_full_mat_ARD",
    "corHMM_output_R1_ordered",
    "corHMM_output_R1_ordered_SYM",
    "corHMM_output_R1_ordered_ER",
    "corHMM_output_R1_ARD_mimic",
    "corHMM_output_R2_full_mat_ARD_ARD",
    "corHMM_output_R2_full_mat_ARD_true",
    "corHMM_output_R2_ordered",
    "corHMM_output_R2_ordered_sym",
    "corHMM_output_R2_ordered_ER",
    "corHMM_output_R2_full_mat_ER",
    "corHMM_output_R2_full_mat_SYM"
  ),
  log.l = c(
    corHMM_output_R1_full_mat_ER$loglik,
    corHMM_output_R1_full_mat_SYM$loglik,
    corHMM_output_R1_full_mat_ARD$loglik,
    corHMM_output_R1_ordered$loglik,
    corHMM_output_R1_ordered_SYM$loglik,
    corHMM_output_R1_ordered_ER$loglik,
    corHMM_output_R1_ARD_mimic$loglik,
    corHMM_output_R2_full_mat_ARD_ARD$loglik,
    corHMM_output_R2_full_mat_ARD_true$loglik,
    corHMM_output_R2_ordered$loglik,
    corHMM_output_R2_ordered_sym$loglik,
    corHMM_output_R2_ordered_ER$loglik,
    corHMM_output_R2_full_mat_ER$loglik,
    corHMM_output_R2_full_mat_SYM$loglik
  ),
  d.f. = c(
    sum(!is.na(corHMM_output_R1_full_mat_ER$solution)),
    sum(!is.na(corHMM_output_R1_full_mat_SYM$solution)),
    sum(!is.na(corHMM_output_R1_full_mat_ARD$solution)),
    sum(!is.na(corHMM_output_R1_ordered$solution)),
    sum(!is.na(corHMM_output_R1_ordered_SYM$solution)),
    sum(!is.na(corHMM_output_R1_ordered_ER$solution)),
    sum(!is.na(corHMM_output_R1_ARD_mimic$solution)),
    sum(!is.na(corHMM_output_R2_full_mat_ARD_ARD$solution)),
    sum(!is.na(corHMM_output_R2_full_mat_ARD_true$solution)),
    sum(!is.na(corHMM_output_R2_ordered$solution)),
    sum(!is.na(corHMM_output_R2_ordered_sym)),
    sum(!is.na(corHMM_output_R2_ordered_ER)),
    sum(!is.na(corHMM_output_R2_full_mat_ER$solution)),
    sum(!is.na(corHMM_output_R2_full_mat_SYM$solution))
  ),
  
  AIC = aic,
  delta.AIC = aic - min(aic),
  weight = exp(-0.5 * aic - min(aic)) / sum(exp(-0.5 * aic - min(aic)))
)



#get model averages

obj <- list(corHMM_output_R1_full_mat_ER,
            corHMM_output_R1_full_mat_SYM,
            corHMM_output_R1_full_mat_ARD,
            corHMM_output_R1_ordered,
            corHMM_output_R1_ordered_SYM,
            corHMM_output_R1_ordered_ER,
            corHMM_output_R1_ARD_mimic,
            corHMM_output_R2_full_mat_ARD_ARD,
            corHMM_output_R2_full_mat_ARD_true,
            corHMM_output_R2_full_mat_ER,
            corHMM_output_R2_full_mat_SYM)

AICcs <- unlist(lapply(obj, function(x) x$AICc))
AICwt <- exp(-0.5 * AICcs - min(AICcs))/sum(exp(-0.5 * AICcs - min(AICcs)))
res <- matrix(0, dim(obj[[1]]$states)[1], dim(obj[[1]]$states)[2])
for(i in 1:length(obj)){
  States <- colnames(obj[[i]]$solution)
  if(length(grep("R2", States)) == 0){
    ASR_i <- obj[[i]]$states[,grep("R1", States)]
  }else{
    ASR_i <- obj[[i]]$states[,grep("R1", States)] + obj[[i]]$states[,grep("R2", States)]
  }
  res <- res + (ASR_i * AICwt[i])
}


# Run simmap
model = corHMM_output_R2_full_mat_ARD_ARD$solution
model[is.na(model)] <- 0 
diag(model) <- -rowSums(model) 
phy= corHMM_output_R2_full_mat_ARD_ARD$phy
data = corHMM_output_R2_full_mat_ARD_ARD$data

tips = corHMM_output_R2_full_mat_ARD_ARD$tip.states

# Check if tip labels of the tree match the row names of the data

corhmm_simmap <-corHMM::makeSimmap(tree = phy, data = data, model = model, rate.cat = 2)


cols<-setNames(c("#57768F", "#41AA8B", "#F9C750", "#57768F", "#41AA8B", "#F9C750"), c("1", "2", "3","4","5","6"))

#to plot the rate tree use these colors
cols<-setNames(c("red", "red", "red", "blue", "blue", "blue"), c("1", "2", "3","4","5","6"))

sampled_simmaps <-sample(dimorph_hidden_ARD_simmap, 100, replace=FALSE)


# next plot an outline tree:
plotTree(dimorph_hidden_ARD_simmap[[1]],ftype="i",fsize=0.45,offset=0.5,
         lwd=2)
par(fg="transparent",lend=1)

plotTree(dimorph_hidden_ARD_simmap[[1]],ftype="i",fsize=0.45,offset=0.5,
         lwd=1,color="white",add=TRUE)

## now plot our 100 stochastic map trees with transparency


for(i in 1:length(sampled_simmaps)) plot(sampled_simmaps[[i]],
                                         colors=sapply(cols,make.transparent,alpha=0.01),
                                         add=TRUE,lwd=2,ftype="i",fsize=0.45,offset=0.5)
par(fg="black")
nodelabels(pie=res,piecol=cols,
           cex=0.25, outlines= FALSE, frame = "none")




# Run Hisse

# Read in the tree
mytimetree_0.rel <- read.tree("/Users/jacobsuissa/Documents/Students/Makaleh_smith/Dimorphism_project/trees/Time_Trees/blechnaceae_Tree_4pnt_lam0.tre")

# Bring in the data
leaf.dat <- read.csv("/Users/jacobsuissa/Documents/Students/Makaleh_smith/Dimorphism_project/CSVs/Dimorphism_data/dimorphism_4k_species_Blechnaceae.csv")

# Drop NAs and select just the data of interest
leaf.dat_NAomit <- leaf.dat %>%
  filter(!is.na(dimorphism)) %>%
  dplyr::select(species_name, dimorphism)

# Change names
leaf.dat_NAomit$dimorphism <- str_replace_all(leaf.dat_NAomit$dimorphism, "monomorphic", "0")
leaf.dat_NAomit$dimorphism <- str_replace_all(leaf.dat_NAomit$dimorphism, "NELD", "1")
leaf.dat_NAomit$dimorphism <- str_replace_all(leaf.dat_NAomit$dimorphism, "ELD", "1")

leaf.dat_NAomit <- leaf.dat_NAomit %>%
  arrange(dimorphism)

# Filter tips in the tree to just match tips in the dataframe
rownames(leaf.dat_NAomit) <- leaf.dat_NAomit[, 1]

# Check species in tree
chk <- name.check(mytimetree_0.rel, leaf.dat_NAomit)

# Now prune the tree based on the data available
tree.pruned <- drop.tip(mytimetree_0.rel, chk$tree_not_data)
fern.Data.pruned <- leaf.dat_NAomit[!(rownames(leaf.dat_NAomit) %in% chk$data_not_tree), , drop = FALSE]

#change character to numeric

fern.Data.pruned$dimorphism <- as.numeric(fern.Data.pruned$dimorphism)

#reorder the species
identical(fern.Data.pruned$species_name, tree.pruned$tip.label)

fern.Data.pruned <- fern.Data.pruned[match(tree.pruned$tip.label, fern.Data.pruned$species_name), ]
identical(fern.Data.pruned$species_name, tree.pruned$tip.label)


#Get states in right format
states <- data.frame(fern.Data.pruned$dimorphism, fern.Data.pruned$dimorphism, 
                     row.names=fern.Data.pruned$species_name)

states <- states[fern.Data.pruned$species_name,]

states.trans <- states


for(i in 1:Ntip(tree.pruned)){
  if(states[i,1] == 0){
    states.trans[i,1] = 1
    states.trans[i,2] = 0
  }
  if(states[i,1] == 1){
    states.trans[i,1] = 0
    states.trans[i,2] = 1
  }
}

#Make the rownames a column name
states.trans1 <- states.trans%>%
  rownames_to_column(var = "Species")


#reorder the species
identical(states.trans1$Species, tree.pruned$tip.label)

#turn to numeric
states.trans2$fern.Data.pruned.dimorphism <- as.numeric(states.trans2$fern.Data.pruned.dimorphism)

states.trans2$fern.Data.pruned.dimorphism.1 <- as.numeric(states.trans2$fern.Data.pruned.dimorphism.1)


#Run dull null model
turnover <- c(1,1)
extinction.fraction <- c(1,1)
f = c(1,1)

#musse no hidden rate
trans.rates.bisse <-  TransMatMakerHiSSE(hidden.traits=0)

#Now, we can call MuHiSSE and estimate the parameters under this model using the default settings:
datmat<- as.matrix(fern.Data.pruned)

tree1<-multi2di(tree.pruned)

dull.null <-  hisse(phy=tree1, data=datmat, f=f, turnover=turnover, 
                    eps=extinction.fraction, hidden.states=FALSE, 
                    trans.rate=trans.rates.bisse)

turnover <- c(1,2)
extinction.fraction <- c(1,1)
true_BiSSE <- hisse(phy=tree1, data=datmat, f=f, turnover=turnover, 
                    eps=extinction.fraction, hidden.states=FALSE, 
                    trans.rate=trans.rates.bisse)



#standard hisse model

turnover <- c(1,2,3,4)
extinction.fraction <- rep(1, 4) 
f = c(1,1)

trans.rate.hisse <- TransMatMakerHiSSE(hidden.traits=1)

#hisse
HiSSE <- hisse(phy=tree1, data=datmat, f=f, turnover=turnover, 
               eps=extinction.fraction, hidden.states=TRUE, 
               trans.rate=trans.rate.hisse)


#CID 2 model

turnover <- c(1, 1, 2, 2)
extinction.fraction <- rep(1, 4) 
f = c(1,1)
trans.rate_cid2 <- TransMatMakerHiSSE(hidden.traits=1, make.null=TRUE)

cid2_HiSSE <- hisse(phy=tree1, data=datmat, f=f, turnover=turnover, 
                    eps=extinction.fraction, hidden.states=TRUE, 
                    trans.rate=trans.rate_cid2)


#CID 4 model
turnover <- c(1, 1, 2, 2, 3, 3, 4, 4)
extinction.fraction <- rep(1, 8) 
trans.rate_cid4 <- TransMatMakerHiSSE(hidden.traits=3, make.null=TRUE)

cid4_HiSSE <- hisse(phy=tree1, data=datmat, f=f, turnover=turnover, 
                    eps=extinction.fraction, hidden.states=TRUE, 
                    trans.rate=trans.rate_cid4)


#Get model averaged rates across each model
hisse_marg_recon<-MarginReconHiSSE(tree1, states.trans1, f, HiSSE$solution, hidden.states=2, condition.on.survival=TRUE, root.type="madfitz", root.p=NULL, includes.fossils = FALSE,
                                   k.samples = NULL, AIC=NULL, get.tips.only=FALSE, verbose=TRUE, n.cores=NULL, dt.threads=2)

#bisse
bisse_marg_recon<-MarginReconHiSSE(tree1, states.trans1, f, true_BiSSE$solution, hidden.states=2, condition.on.survival=TRUE, root.type="madfitz", root.p=NULL, includes.fossils = FALSE,
                                   k.samples = NULL, AIC=NULL, get.tips.only=FALSE, verbose=TRUE, n.cores=NULL, dt.threads=2)

#null dull bisse
null_marg_recon<-MarginReconHiSSE(tree1, states.trans1, f, dull.null$solution, hidden.states=2, condition.on.survival=TRUE, root.type="madfitz", root.p=NULL, includes.fossils = FALSE,
                                  k.samples = NULL, AIC=NULL, get.tips.only=FALSE, verbose=TRUE, n.cores=NULL, dt.threads=2)

#cid2
cid2_marg_recon<-MarginReconHiSSE(tree1, states.trans1, f, cid2_HiSSE$solution, hidden.states=2, condition.on.survival=TRUE, root.type="madfitz", root.p=NULL, includes.fossils = FALSE,
                                  k.samples = NULL, AIC=NULL, get.tips.only=FALSE, verbose=TRUE, n.cores=NULL, dt.threads=2)

#cid4
cid4_marg_recon<-MarginReconHiSSE(tree1, states.trans1, f, cid4_HiSSE$solution, hidden.states=2, condition.on.survival=TRUE, root.type="madfitz", root.p=NULL, includes.fossils = FALSE,
                                  k.samples = NULL, AIC=NULL, get.tips.only=FALSE, verbose=TRUE, n.cores=NULL, dt.threads=2)

#Compute AIC weights
aic_values<- c(cid4_HiSSE$AIC, cid2_HiSSE$AIC, dull.null$AIC, true_BiSSE$AIC, HiSSE$AIC)

aicweights <-exp(-0.5 * (aic_values - min(aic_values))) / sum(exp(-0.5 * (aic_values - min(aic_values))))

AIDcweights <- setNames(as.numeric(aicweights), 
                        c("cid4_marg_recon", "cid2_marg_recon", "null_marg_recon",
                          "bisse_marg_recon", "hisse_marg_recon"))


hissMods<- list(cid4_marg_recon,cid2_marg_recon, null_marg_recon, bisse_marg_recon, hisse_marg_recon )

hiss_mod_Avg <- hisse::GetModelAveRates(x= hissMods, AIC.weights = aicweights, type= "both" )



#estimate node states
hisse_marg_recon<-MarginReconHiSSE(tree1, states.trans1, f, HiSSE$solution, hidden.states=2, condition.on.survival=TRUE, root.type="madfitz", root.p=NULL, includes.fossils = FALSE,k.samples = NULL, AIC=NULL, get.tips.only=FALSE, verbose=TRUE, n.cores=NULL, dt.threads=2)

#plot
hisse::plot.hisse.states(hisse_marg_recon, rate.param = "net.div", type = "phylogram", show.tip.label = TRUE ,fsize = 0.3, legend.kernel= "gaussian"  )


# plot model averages
for (i in seq_along(hissMods)) {
  hissMods[[i]]$AIC <- aic_values[i]
}

hisse::plot.hisse.states(hissMods, rate.param = "net.div", type = "phylogram", show.tip.label = TRUE ,fsize = 0.3, legend.kernel= "gaussian"  )



#make plots

# Summary statistics of turnover, net.div, speciation, and extinction rates
summary_stats_tips <- hiss_mod_Avg$tips %>%
  summarise(mean_turnover = mean(turnover, na.rm = TRUE),
            mean_net_div = mean(net.div, na.rm = TRUE),
            mean_speciation = mean(speciation, na.rm = TRUE),
            mean_extinction = mean(extinction, na.rm = TRUE),
            extinction_fraction = mean(extinct.frac, na.rm = TRUE))


hiss_mod_Avg$tips%>%
  group_by(state)%>%
  mutate(mean_net_div = mean(net.div, na.rm = TRUE))%>%
  distinct(state, .keep_all = TRUE)

ggplot(hiss_mod_Avg$tips, aes(x=as.character(state), y=net.div, color=as.factor(state))) +
  geom_point(position = "jitter") + geom_violin(alpha=0.5, draw_quantiles = c(0.25, 0.5, 0.75))+theme_classic()+ guides(color=FALSE)+ xlab("Leaf dimorphism") +ylab ("Net diversification") +scale_color_manual(values= c("darkgray", "black"))




#Hygromorphic data

#Read in phys data
dat <- read.csv(/path/to/Final_trial_measurements_total.csv)

dat%>%
  group_by(species_and_herb, pinnule,  humidity_perc)%>%
  mutate(mean_dat = mean(length_mm)) %>%
  distinct(species_and_herb, pinnule, humidity_perc, .keep_all = TRUE)%>%
  mutate(humidity_perc= as.character(humidity_perc))%>%
  arrange(Species)%>%
  ggplot(., aes(x = humidity_perc , y = log10(mean_dat),  group=humidity_perc)) + 
  geom_jitter(color="black", alpha=0.6)+
  geom_boxplot(fill="gray", alpha=0.5) +
  labs(x = "Relative humidity (%)", y = "log(pinnule length (mm))")+
  theme_classic() + guides(color=FALSE) + guides(fill=FALSE) +facet_wrap(~Species, scales = "free")+stat_compare_means(comparisons = list(c("0", "100")), label = "p.signif", method='t.test', p.adjust.method = "bonferroni")



#run t-test with bonferoni correction

aa<-dat%>%
  group_by(species_and_herb, pinnule,  humidity_perc)%>%
  mutate(mean_dat = mean(length_mm)) %>%
  distinct(species_and_herb, pinnule, humidity_perc, .keep_all = TRUE)%>%
  mutate(humidity_perc= as.character(humidity_perc))%>%
  arrange(Species)

t.test(x = aa$humidity_perc, y = aa$mean_dat, p.adjust.method="bonferroni")

summary(lm(aa$mean_dat~aa$humidity_perc + aa$Species))


#get mean and SD for each species:
dat%>%
  group_by(species_and_herb, pinnule,  humidity_perc)%>%
  mutate(mean_dat = mean(length_mm)) %>%
  distinct(species_and_herb, pinnule, humidity_perc, .keep_all = TRUE)%>%
  mutate(humidity_perc= as.character(humidity_perc))%>%
  arrange(Species)%>%
  group_by(Species,humidity_perc)%>%
  mutate(mean_length= mean(mean_dat))%>%
  mutate(SD_length= sd(mean_dat))%>%
  select(Species,mean_length,SD_length )%>%
  distinct(Species,.keep_all = TRUE)

