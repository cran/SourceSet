## ---- eval=FALSE---------------------------------------------------------
#  install.packages("SourceSet")
#  help(package = "SourceSet")

## ---- message=FALSE, warning=FALSE, eval = TRUE--------------------------
library(SourceSet)
#?sourceSet

## ---- message=FALSE, warning=FALSE---------------------------------------
data("simulation")
names(simulation)

## ---- message=FALSE, warning=FALSE---------------------------------------
library(mvtnorm)
set.seed(111)
# sample size
n<-50

# parameters of control condition
param.cond1<-simulation$condition1

# parameters of perturbed condition (`5`: true source; `2`: dysregulation intensity, strong)
param.cond2<-simulation$condition2$`5`$`2`

# condition 1 
data.cond1<-rmvnorm(n = n,mean =param.cond1$mu ,sigma =param.cond1$S )
# condition 2
data.cond2<-rmvnorm(n = n,mean =param.cond2$mu ,sigma=param.cond2$S)

# Input arguments for the sourceSet function
data<-rbind(data.cond1,data.cond2)
classes<-c(rep(1,nrow(data.cond1)),rep(2,nrow(data.cond2)))
graphs<-list("source.node5"=simulation$graph)

## ---- message=FALSE, warning=FALSE---------------------------------------
simulation$graph
ripped<-gRbase::rip(simulation$graph)

# number of cliques
length(ripped$cliques)
# size of larger clique
max(sapply(ripped$cliques,length))

## ---- echo=TRUE, message=FALSE, warning=FALSE----------------------------
result<-sourceSet(graphs ,data ,classes ,seed = 123 ,permute =FALSE ,shrink =FALSE, alpha=0.05  )
class(result)

## ------------------------------------------------------------------------
names(result$source.node5)

## ------------------------------------------------------------------------
# source set: primary disregulation
result$source.node5$primarySet

# secondary disregulation
result$source.node5$secondarySet

# all affected variables
unique(unlist(result$source.node5$orderingSet))

## ------------------------------------------------------------------------
# number of orderings
length(result$source.node5$Decompositions)
# number of unique components
nrow(result$source.node5$Components)
# ordering with root clique 'C3'
result$source.node5$Decompositions$C5

## ------------------------------------------------------------------------
# alpha and corrected threshold
result$source.node5$Threshold[c("alpha","value")]

# ordering source set when 'C5' is used as root clique
result$source.node5$orderingSet$C5

# manual indentification of the ordering source set
union(result$source.node5$Elements$C2,result$source.node5$Elements$C5)

## ------------------------------------------------------------------------
# source set of each ordering
result$source.node5$orderingSet

# primary disregulation
result$source.node5$primarySet

# secondary disregulation
result$source.node5$secondarySet

# all affaceted variables
unique(unlist(result$source.node5$orderingSet))

## ------------------------------------------------------------------------
info<-infoSource(result)
names(info)

## ---- echo=FALSE, message=FALSE, warning=FALSE---------------------------
knitr::kable(info$graph,caption = "> info$graph")

## ---- echo=FALSE, message=FALSE, warning=FALSE---------------------------
knitr::kable(info$variable[,-1],caption=" > info$variable")

## ---- fig.align="center", fig.height=2, message=FALSE, warning=FALSE-----
easyLookSource(result)

## ---- message=FALSE, warning=FALSE---------------------------------------
set.seed(222)
data2.cond1<-rmvnorm(n = n,mean =simulation$condition1$mu ,sigma =simulation$condition1$S )
data2.cond2<-rmvnorm(n = n,mean =simulation$condition2$`10`$`2`$mu ,sigma =simulation$condition2$`10`$`2`$S)

# Input arguments for the sourceSet function
data2<-rbind(data2.cond1,data2.cond2)
classes<-c(rep(1,nrow(data2.cond1)),rep(2,nrow(data2.cond2)))
graphs<-list("source.node10"=simulation$graph)

result2<-sourceSet(graphs,data2,classes,seed=222,permute = FALSE,shrink = FALSE)

## ---- eval=FALSE---------------------------------------------------------
#  sourceSankeyDiagram(result2,height = 100,width = 800)

## ---- echo=FALSE, message=FALSE, warning=FALSE---------------------------
sourceSankeyDiagram(result2,height = 150,width = 800)

## ---- eval=FALSE, message=FALSE, warning=FALSE---------------------------
#  # Lunch cytoscape and run the following commands
#  
#  # simulation 1: sourceset composed by variable 5
#  cytoID.5<-sourceCytoscape(result,collection.name = "Simulation")
#  # simulation 2: sourceset composed by variable 10,9,8,5
#  cytoID.10<-sourceCytoscape(result2,collection.name = "Simulation")

## ---- message=FALSE, warning=FALSE---------------------------------------
library(ALL); data("ALL")
ALL

## ---- echo=TRUE, message=FALSE, warning=FALSE----------------------------
# Retrieve expression matrix and phenotypes
library(Biobase)
data.all<-Biobase::exprs(ALL)
pheno.all<-Biobase::pData(ALL)

# Select individuals:
# BT: type and stage of the disease ('B' indicates B-cell ALL)
# mol.biol: molecular biology of cancer ('BCR/ABL' and 'NEG')
ind<- intersect(grep("B",pheno.all$BT), which(pheno.all$mol.biol %in% c("BCR/ABL","NEG")))
code<- rownames(pheno.all)[ind]
group<- paste(pheno.all$mol.biol[ind])
data.all<-data.all[,code]

## ---- echo=TRUE, message=FALSE, warning=FALSE----------------------------
# Convert identifiers:
# map between manufacturer identifiers and Entrez Gene identifiers
library(hgu95av2.db)
mapped_probes <- mappedkeys(hgu95av2ENTREZID)
# remove not mapped probes
data.all<-data.all[rownames(data.all) %in% mapped_probes,]
# convert identifiers
entrez.id<-paste(hgu95av2ENTREZID[rownames(data.all)])
# merge not unique code on mean value
data.all<-apply(data.all,2,function(x,f) { tapply(x,f,mean) }, f=entrez.id)

## ---- echo=TRUE, message=FALSE, warning=FALSE----------------------------
# sourceSet function arguments:
data<-t(data.all)
classes<-sapply(group,function(x) switch(x,"BCR/ABL"=2,"NEG"=1))
table(classes)
ncol(data)

## ---- eval=TRUE, message=FALSE, warning=FALSE----------------------------
library(graphite); library(graph)
# pathways selection
names<-c("Axon guidance","Cell cycle","Chronic myeloid leukemia","ErbB signaling pathway",
"Neurotrophin signaling pathway","Pathways in cancer","Ras signaling pathway","Viral myocarditis")
# retrieve a list of pathways from a database for a given species
pathways  <- pathways("hsapiens", "kegg")[names]
# convert the node identifiers of pathways
pathways<-convertIdentifiers(pathways,"entrez")

# For each pathway, build a graphNEL object representing its topology
graphs<-lapply(pathways,function(p) pathwayGraph(p))

# Match node IDs with the names of the data matrix columns (delete the prefix 'ENTREZID:')
# (graphite version 1.24.1)
graph::nodes(graphs[[1]])[1:3]
colnames(data)[1:3]
for(i in 1:length(graphs)) graph::nodes(graphs[[i]])<-gsub("ENTREZID:","",graph::nodes(graphs[[i]]))

## ------------------------------------------------------------------------
graphs$`Chronic myeloid leukemia`

## ---- message=FALSE, warning=FALSE, include=FALSE------------------------
path<-system.file("extdata","ALLsourceresult.RData",package = "SourceSet")
load(file = path)

## ---- eval=FALSE, echo=TRUE, warning=FALSE-------------------------------
#  # It requires about 16 minutes:
#  # run instead: load(file=system.file("extdata","ALLsourceresult.RData",package = "SourceSet"))
#  results.all<-sourceSet(graphs,data,classes,seed =111 ,permute =TRUE ,shrink =TRUE )

## ---- fig.height=3, fig.width=10-----------------------------------------
# Convert identifiers:
# map between Entrez Gene identifiers and gene symbols
library(org.Hs.eg.db)
mapped.genes.symbol <- as.list(org.Hs.egSYMBOL[rownames(data.all)])

# Lists of primary genes for the analyzed graphs
primary<-lapply(results.all,function(x) x$primarySet)
# Number of primary genes 
n.primary<-length(unique(unlist(primary)))

easyLookSource(sourceObj = results.all, map.name.variable = mapped.genes.symbol,
               maxnum.variable = n.primary,
               label.variable = "Genes",label.graph = "Pathways")

## ---- message=FALSE, warning=FALSE---------------------------------------
info.all<-infoSource(results.all,map.name.variable = mapped.genes.symbol)

## ---- echo=FALSE, message=FALSE, warning=FALSE---------------------------
knitr::kable(info.all$variable[1:5,],caption = "> info.all$variable")

## ---- echo=FALSE, message=FALSE, warning=FALSE---------------------------
knitr::kable(info.all$graph,caption = "> info.all$graph")

## ---- eval=FALSE, message=FALSE, warning=FALSE---------------------------
#  sourceSankeyDiagram(results.all,height = 600,width = 800,map.name.variable = mapped.genes.symbol)

## ---- echo=FALSE---------------------------------------------------------
sourceSankeyDiagram(results.all,height = 600,width = 800,map.name.variable = mapped.genes.symbol)

## ---- eval=FALSE---------------------------------------------------------
#  # NB: Remember to launch cytoscape before running the following commands
#  
#  # Create two collections of pathways to visualize the results
#  graph.signaling<-names(results.all)[grep("signaling",names(results.all))]
#  graph.other<-setdiff(names(results.all),graph.signaling)
#  
#  # Signaling collection + union source set
#  cytoID.signaling<-sourceCytoscape(results.all, name.graphs = graph.signaling,
#        collection.name ="SignalingPathway", map.name.variable = mapped.genes.symbol)
#  cytoID.signaling.union<-sourceUnionCytoscape(results.all, name.graphs =graph.signaling ,
#        collection.name ="SignalingPathway" ,network.name ="SignalingUnion",
#        map.name.variable =mapped.genes.symbol)
#  
#  # Other collection + union source set
#  cytoID.other<-sourceCytoscape(results.all, name.graphs = graph.other,
#        collection.name ="OtherPathway", map.name.variable = mapped.genes.symbol)
#  cytoID.other.union<-sourceUnionCytoscape(results.all ,name.graphs =graph.other,
#        collection.name ="OtherPathway" ,network.name ="OtherUnion",
#        map.name.variable =mapped.genes.symbol)

## ---- message=FALSE, warning=FALSE, eval=FALSE---------------------------
#  # Install from Bioconductor
#  source("https://bioconductor.org/biocLite.R")
#  biocLite("r2cytoscape")
#  
#  # or.. Install from GitHub
#  library(devtools)
#  install_github("cytoscape/r2cytoscape")

