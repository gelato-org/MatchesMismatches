
#######################################
# SUBSET ANALYSIS COMPARING GELATO WITH LANGUAGE DIVERGENCE TIMES
#######################################
##******************************************
#******************************************
### Analysis Paper "A global analysis of matches and mismatches between human genetic and linguistic histories"
# GeLaTo database
# second part: phylogenetic comparisons
# Chiara BArbieri June 2022

#******************************************
#*
#*
#*

### read the two main files
# list of 404 populations :
perpopRED<-read.table("PerpopRED_MaMi2022.txt", header = T, sep = "\t", as.is=T)

# list of pairwise comparisons :
FstListinfo<-read.csv("DatasetS2_pairwise_FstListREDinfo_MaMi2022.csv", header=T, sep=",")


# color palette for major language families
#
Lfamil<-table(perpopRED$glottolog.NAME)
MainFamilies<-unlist(labels(Lfamil[which(Lfamil>4)])) # minimum 5 populations per Lang Family
sizes<-as.numeric(unlist((Lfamil[which(Lfamil>4)])))
perpopREDfamily<-perpopRED[which(perpopRED$glottolog.NAME %in% MainFamilies),]
colorchoice<-c( "darkorange4" ,"#93032E"  ,   "#33A02C"    , "#A6D854"  ,   "#377EB8"  ,   "#E31A1C"   ,  "#FFD92F" ,    "#FF7F00"  ,   "#666666" ,   
                "cyan4"  ,     "#BC80BD"   ,  "#FED9A6" ,    "tan3" ,       "#6A3D9A" ,    "deeppink"   )
MainFamilies2<-as.data.frame(MainFamilies)
MainFamilies2$COLOR<-colorchoice
MainFamilies2$names<-paste0(MainFamilies, " (",sizes, ")" )

library(ggplot2)
#******************************************
#******************************************
#******************************************
#******************************************

meltFstREDinfo<-FstListinfo

library(phylobase)
library(pegas)
library(ggplot2)

### Prepare Lang Distances corresponding to the pairs i have in GeLaTo
#*******
#*
## AUSTRONESIAN

TARGETFamily<-"Austronesian"
listAUS<-unique(perpopRED$Gray2009)
listAUS<-listAUS[-which(listAUS=="0")]

AUS<-read.nexus("withinFamilyTMRCA/Gray2009.tree")
# from https://github.com/D-PLACE/dplace-data/blob/master/phylogenies/gray_et_al2009/summary.trees
taxaAUS<-read.csv("withinFamilyTMRCA/taxaGray2009.csv")
# from https://github.com/D-PLACE/dplace-data/blob/master/phylogenies/gray_et_al2009/taxa.csv

#*******

listAUSreverse<-as.character(taxaAUS$taxon[match(listAUS,taxaAUS$glottocode)])
## remove duplicated languages

g1 <- as(AUS, "phylo4")
bi_subset <- subset(g1, tips.include=listAUSreverse) # select only the tips i have. note, i have duplicated tips.
alberoLang<-as(bi_subset,"phylo")

alberoLang<-makeNodeLabel(alberoLang, method = "number", prefix = "")

alberoLang$edge.length<-alberoLang$edge.length*1000  # Austronesian tree is in fractions of year
alberoLangPhylo4 <- as(alberoLang, "phylo4")
alberoLang$tip.label<-taxaAUS$glottocode[match(alberoLang$tip.label,taxaAUS$taxon)] # names in glottocode
write.tree(alberoLang,paste0(TARGETFamily ,"alberoLangPhylo4.phy"))

#****************************************
# INDO EUROPEAN
#****************************************


TARGETFamily<-"Indo-European"
listAUS<-unique(perpopRED$Bouckaert2012)
listAUS<-listAUS[-which(listAUS=="0")]
AUS<-read.nexus("withinFamilyTMRCA/Bouckaert2012.trees")
# from "https://github.com/D-PLACE/dplace-data/blob/master/phylogenies/bouckaert_et_al2012/summary.trees"

taxaAUS<-read.table("withinFamilyTMRCA/Bouckaert2012_taxa.csv", sep=";", header=T, quote = "")
# from "https://github.com/D-PLACE/dplace-data/blob/master/phylogenies/bouckaert_et_al2012/taxa.csv"

listAUSreverse<-as.character(taxaAUS$taxon[match(listAUS,taxaAUS$glottocode)])
## remove duplicated languages

g1 <- as(AUS, "phylo4")
bi_subset <- subset(g1, tips.include=listAUSreverse) # select only the tips i have. note, i have duplicated tips.
alberoLang<-as(bi_subset,"phylo")

alberoLang<-makeNodeLabel(alberoLang, method = "number", prefix = "")
alberoLangPhylo4 <- as(alberoLang, "phylo4")
alberoLang$tip.label<-taxaAUS$glottocode[match(alberoLang$tip.label,taxaAUS$taxon)] # names in glottocode
write.tree(alberoLang,paste0(TARGETFamily ,"alberoLangPhylo4.phy")) # manually check and adapt

#******************************************
# Turkic
#******************************************
#*
#* Hrushka 2015
#*  

TARGETFamily<-"Turkic"

listAUS<-unique(perpopRED$Hruschka2015)
listAUS<-listAUS[-which(listAUS=="0")]
AUS<-read.nexus("withinFamilyTMRCA/Hrushka_summary2.trees")
#from https://github.com/D-PLACE/dplace-data/blob/master/phylogenies/hruschka_et_al2015/summary.trees
taxaAUS<-read.csv("withinFamilyTMRCA/Hrushka_taxa2.csv", header=T, quote = "")
# from https://github.com/D-PLACE/dplace-data/blob/master/phylogenies/hruschka_et_al2015/taxa.csv
colnames(taxaAUS)[1]<-"code"
colnames(taxaAUS)[2]<-"taxon"

listAUSreverse<-as.character(taxaAUS$taxon[match(listAUS,taxaAUS$glottocode)])
## remove duplicated languages

g1 <- as(AUS, "phylo4")
bi_subset <- subset(g1, tips.include=listAUSreverse) # select only the tips i have. note, i have duplicated tips.
alberoLang<-as(bi_subset,"phylo")

alberoLang<-makeNodeLabel(alberoLang, method = "number", prefix = "")
alberoLangPhylo4 <- as(alberoLang, "phylo4")
alberoLang$tip.label<-taxaAUS$glottocode[match(alberoLang$tip.label,taxaAUS$taxon)] # names in glottocode
write.tree(alberoLang,paste0(TARGETFamily ,"alberoLangPhylo4.phy"))


#########



library("reshape")
library("nodiv")

library('Quartet')
library("phytools")
library("phylogram")

my.max <- function(x) ifelse( !all(is.na(x)), max(x, na.rm=T), NA)
meltFstREDinfo$threefamiliesNEW<-NA

### LOOP for each of the 3 lang families


for (k in 1:3){
  TARGETFamily<-ThreeFamilies[k]
  
  #### Indo-European
  
  if(TARGETFamily=="Indo-European"){
    listAUS<-unique(perpopRED$Bouckaert2012)
    listAUS<-listAUS[-which(listAUS=="0")]
    AUS<-read.nexus("withinFamilyTMRCA/Bouckaert2012.trees")
    taxaAUS<-read.table("withinFamilyTMRCA/Bouckaert2012_taxa.csv", sep=";", header=T, quote = "")
    meltFstREDinfoTEMP<-meltFstREDinfo
    
    meltFstREDinfoTEMP$glottocodeBase11<-perpopRED$Bouckaert2012[match(meltFstREDinfoTEMP$Pop1,perpopRED$PopName)]
    meltFstREDinfoTEMP$glottocodeBase22<-perpopRED$Bouckaert2012[match(meltFstREDinfoTEMP$Pop2,perpopRED$PopName)]
    ## root with Armenian as outgroup, according to the linguistic tree IE
    rootalo<-c("Armenian_Hemsheni","Armenian")
    
    
  }
  
  #***********************
  #*  #### Austronesian
  
  if(TARGETFamily=="Austronesian"){
    listAUS<-unique(perpopRED$Gray2009)
    listAUS<-listAUS[-which(listAUS=="0")]
    AUS<-read.nexus("withinFamilyTMRCA/Gray2009.tree")
    taxaAUS<-read.csv("withinFamilyTMRCA/taxaGray2009.csv")
    meltFstREDinfoTEMP<-meltFstREDinfo
    
    meltFstREDinfoTEMP$glottocodeBase11<-perpopRED$Gray2009[match(meltFstREDinfoTEMP$Pop1,perpopRED$PopName)]
    meltFstREDinfoTEMP$glottocodeBase22<-perpopRED$Gray2009[match(meltFstREDinfoTEMP$Pop2,perpopRED$PopName)]
    ## root with Ami as outgroup, according to the linguistic tree AUSTRONESIAN
    rootalo<-c("Ami")
  }
  
  #***********************
  #*  #*  #### Turkic
  
  if(TARGETFamily=="Turkic"){
    listAUS<-unique(perpopRED$Hruschka2015)
    listAUS<-listAUS[-which(listAUS=="0")]
    AUS<-read.nexus("withinFamilyTMRCA/Hrushka_summary2.trees")
    taxaAUS<-read.csv("withinFamilyTMRCA/Hrushka_taxa2.csv", header=T, quote = "")
    colnames(taxaAUS)[1]<-"code"
    colnames(taxaAUS)[2]<-"taxon"
    meltFstREDinfoTEMP<-meltFstREDinfo
    
    meltFstREDinfoTEMP$glottocodeBase11<-perpopRED$Hruschka2015 [match(meltFstREDinfoTEMP$Pop1,perpopRED$PopName)]
    meltFstREDinfoTEMP$glottocodeBase22<-perpopRED$Hruschka2015[match(meltFstREDinfoTEMP$Pop2,perpopRED$PopName)]
    
    ## root with two Chuvash speaking as outgroup, according to the linguistic tree Turkic
    rootalo<-c("Chuvash","Chuvash_Tatarstan")
    
  }
  
  
  #***********************
  #*
  alberoLang<-read.tree(paste0(TARGETFamily ,"alberoLangPhylo4.phy"))  # the tree selection i saved in the first part of the script above
  alberoLangPhylo4 <- as(alberoLang, "phylo4")
  alberoLang$tip.label<-taxaAUS$glottocode[match(alberoLang$tip.label,taxaAUS$glottocode)] # names in glottocode
  distancealberoLang<-cophenetic.phylo(alberoLang) # names in glottocode
  distancealberoLang<-distancealberoLang/2 # remember the distances are double!!
  MELTdistancealberoLang<-melt(distancealberoLang,varnames=c('glottocodeBase11', 'glottocodeBase22'))
  colnames(MELTdistancealberoLang)[3]<-"threefamilies"
  
  ## assign the time divergence for the three main families to the master info FST pairwise file
  for (q in 1:nrow(meltFstREDinfo)){
    
    coso<-as.numeric(MELTdistancealberoLang$threefamilies)[which(MELTdistancealberoLang$glottocodeBase11==meltFstREDinfoTEMP$glottocodeBase11[q]&
                                                                   MELTdistancealberoLang$glottocodeBase22==meltFstREDinfoTEMP$glottocodeBase22[q])]
    if(length(coso)!=0){
      meltFstREDinfo$threefamiliesNEW[q]<-coso
    }
    
  }
  
  distanceTMRCAred<-distancealberoLang
  
  for (i in 1:nrow(distanceTMRCAred)){
    for (j in 1:ncol(distanceTMRCAred)){
      temprows<-meltFstREDinfoTEMP[which(meltFstREDinfoTEMP$glottocodeBase11==rownames(distanceTMRCAred)[i]
                                         &meltFstREDinfoTEMP$glottocodeBase22==colnames(distanceTMRCAred)[j]),]
      distanceTMRCAred[i,j]<-mean(na.omit(temprows$TMRCA_doubleNe))  # note i use the mean in case of more gen pop for the same language
    }
  }
  diag(distanceTMRCAred)<-NA
  
  meltFstREDinfoTEMPdRED<-meltFstREDinfoTEMP[!is.na(meltFstREDinfoTEMP$TMRCA_doubleNe),] # exclude the pairs for which i do not have genetic TMRCA reconstructed
  glottredlist<-unique(meltFstREDinfoTEMPdRED$glottocodeBase11)
  
  listAUSreverseInGenTMRCA<-taxaAUS$taxon[match(glottredlist,taxaAUS$glottocode)]
  
  
  # proportion time Lang/Gen matrix distance
  distancePROPORTIONred<- distancealberoLang/distanceTMRCAred
  
  
  nodevalues<-matrix(NA, Nnode(alberoLang),6)
  colnames(nodevalues)<-c("nodename","maxDivergenceTimeLang",
                          "maxDivergenceTimeGen","meanDivergenceTimeGen", 
                          "MaxCoupleproportionDivergenceTime","meanproportionDivergenceTime")
  nodevalues[,1]<-c(1:Nnode(alberoLang))
  
  for (j in 1:Nnode(alberoLang)){
    settemp<- labels(descendants(alberoLangPhylo4,which(attributes(alberoLangPhylo4)$label==j),type = "all"))
     
    timetemp<-melt(distancealberoLang[which(colnames(distancealberoLang)%in%settemp),which(rownames(distancealberoLang)%in%settemp)])
    timetemp$value<-round(timetemp$value)
    nodevalues[j,2]<-(max(timetemp$value))
    rowsmaxdivergence<-which(timetemp$value==max(timetemp$value))
    maxdivergence<-timetemp[rowsmaxdivergence,] # i cannot use only the maximum divergence time otherwise i do not have useful matches with genetic data, so i pick up the maximum from genetic divergence of all the derived nodes
    
    gentimetemp<-melt(distanceTMRCAred[which(colnames(distanceTMRCAred)%in%settemp),which(rownames(distanceTMRCAred)%in%settemp)])[rowsmaxdivergence,]
    gentimetempProportion<- melt( distancePROPORTIONred[which(colnames(distancePROPORTIONred)%in%settemp),which(rownames(distancePROPORTIONred)%in%settemp)])[rowsmaxdivergence,]
    nodevalues[j,3]<-my.max(gentimetemp$value)
    nodevalues[j,4]<-mean(gentimetemp$value, na.rm = T)
    nodevalues[j,6]<-mean(gentimetempProportion$value, na.rm = T)
  }
  
  nodevalues<-as.data.frame(nodevalues)
  nodevalues$mainFamily<-TARGETFamily
 
  listAUSreverseInGenTMRCA<-taxaAUS$taxon[match(glottredlist,taxaAUS$glottocode)]
  
  alberoLangNamesLang<-alberoLang
  alberoLangNamesLang$tip.label<-taxaAUS$taxon [match(alberoLangNamesLang$tip.label,taxaAUS$glottocode)]
  
  alberoLangNamesLang$tip.label[which(alberoLangNamesLang$tip.label%in%listAUSreverseInGenTMRCA)]<-paste0(alberoLangNamesLang$tip.label[which(alberoLangNamesLang$tip.label%in%listAUSreverseInGenTMRCA)], "_GEN_TIME") # mark the names of languages for which i do have the genetic divergence time
  
  proportionGenLang <- nodevalues$meanproportionDivergenceTime
  
  cexplay=1.5
  
  pdf(paste0(TARGETFamily,"_2022.pdf"),width=12, height=5,useDingbats=FALSE)
  par(mfrow=c(1,4))
  plot_nodes_phylo(round(nodevalues$maxDivergenceTimeLang), alberoLangNamesLang, cex = cexplay, main= "maximum language divergence time")
  plot_nodes_phylo(round(nodevalues$maxDivergenceTimeGen), alberoLangNamesLang, cex = cexplay, main= "maximum genetic divergence time")
  plot_nodes_phylo(round(nodevalues$meanDivergenceTimeGen), alberoLangNamesLang, cex = cexplay, main= "mean genetic divergence time")
  plot_nodes_phylo(proportionGenLang, alberoLangNamesLang, cex = cexplay, main= "Mean Lang/Gen proportion")
  dev.off()
  
  
  #**********************************************************
  ### FST tree of the selected languages
  # compare phylogeny fst and phylogeny language time tree
  #**********************************************************
  if(TARGETFamily=="Indo-European"){
    listpopLang<-perpopRED$PopName[which(perpopRED$Bouckaert2012%in% alberoLang$tip.label)]
  }
  if(TARGETFamily=="Austronesian"){
    listpopLang<-perpopRED$PopName[which(perpopRED$Gray2009%in% alberoLang$tip.label)]
  }
  if(TARGETFamily=="Turkic"){
    listpopLang<-perpopRED$PopName[which(perpopRED$Hruschka2015%in% alberoLang$tip.label)]
  }
  
  FstREDlangMatrix<-matrix(NA, length(listpopLang),length(listpopLang))
  rownames(FstREDlangMatrix)<-listpopLang
  colnames(FstREDlangMatrix)<-listpopLang
  
  
  for (i in 1:nrow(FstREDlangMatrix)){
    for (j in 1:ncol(FstREDlangMatrix)){
      temp<-which(meltFstREDinfo$Pop1==rownames(FstREDlangMatrix)[i]&meltFstREDinfo$Pop2==colnames(FstREDlangMatrix)[j])
      if(length(temp)>0){
        FstREDlangMatrix[i,j]<-meltFstREDinfo$FstLinear[temp]
      }
    }
  }
  
  diag(FstREDlangMatrix)<-0
  
  # make the language time matrix with the populations present in gelato
  
  timetreelangMatrix<-matrix(NA, length(listpopLang),length(listpopLang))
  rownames(timetreelangMatrix)<-listpopLang
  colnames(timetreelangMatrix)<-listpopLang
  
  for (i in 1:nrow(timetreelangMatrix)){
    for (j in 1:ncol(timetreelangMatrix)){
      temp<-which(meltFstREDinfo$Pop1==rownames(FstREDlangMatrix)[i]&meltFstREDinfo$Pop2==colnames(FstREDlangMatrix)[j])
      if(length(temp)>0){
        timetreelangMatrix[i,j]<-as.numeric(meltFstREDinfo$threefamiliesNEW[temp])
      }
    }
  }
  
  diag(timetreelangMatrix)<-0
  timetreelangMatrix<-timetreelangMatrix/2
  
  phy1 <- nj(FstREDlangMatrix)
  phy2 <- nj(timetreelangMatrix)
  
  phy1 <- root(phy1, rootalo)
  phy2 <- root(phy2, rootalo)
  phy1$edge.length[phy1$edge.length < 0] = 0.002
  
  
  write.tree(phy1,paste0(TARGETFamily, "_FST_2022.tree"))
  write.tree(phy2,paste0(TARGETFamily, "_lang_2022.tree"))
  
  ### QUARTET measurements
  statuses <- QuartetStatus(phy1, phy2)
  QuartetDivergence(statuses, similarity = FALSE) 
  print (QuartetDivergence(statuses, similarity = FALSE) )
  SimilarityMetrics(statuses, similarity = TRUE)
  print(SimilarityMetrics(statuses, similarity = TRUE))
  
  pdf(paste0(TARGETFamily,"Quartet_2022.pdf"),useDingbats=FALSE, height = 10, width = 10)
  VisualizeQuartets(phy2, phy1, scale=0.6)
  dev.off() # i cannot group the quartet plots in a single figure!
  
  ### compare phylogenies with Phytools
  
  #  TARGETFamily<-ThreeFamilies[k]
  #  phy1<-read.tree(paste0(TARGETFamily, "_FST.tree"))
  #  phy2<-read.tree(paste0(TARGETFamily, "_lang.tree"))
  
  # plot compared phylogenies
  pdf(paste0(TARGETFamily,"Cophilo2022.pdf"),width=15, height=7,useDingbats=FALSE)
  
  plot(cophylo(phy2,phy1,rotate=T), fsize=0.6)
  # nodelabels.cophylo(phy2$node.label)
  #  nodelabels.cophylo(which="right",phy1$node.label )
  
  dev.off()
  
}

#****************************************#****************************************
#*#****************************************
#*#****************************************
#*#****************************************
#*#****************************************
## prepare correlation between nodes, with 95% credible intervals from the original Bayesian linguistic phylogeny
#****************************************
#### with the support of Simon Greenhill
#****************************************

library(tidyverse)


# function written by Balthasar Bickel

get_smallest_clade_age <- function(language1, language2, age_table) {
  d <- filter(age_table, grepl(language1, nodesGlotto) & grepl(language2, nodesGlotto)) %>%
    slice_min(order_by = clade_size)
  return(data.frame(Language1 = language1,
                    Language2 = language2,
                    LinguisticDivergenceTime_median = d$median,
                    LinguisticDivergenceTime_lower = d$hpdlower,
                    LinguisticDivergenceTime_upper = d$hpdupper
  ))
}

# --------------------------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------
#### PLOT CORRELATIONS SIMPLE
# -----------------------------------------------------------------
## AUSTRONESIAN
# 
TARGETFamily<-"Austronesian"

listAUS<-unique(perpopRED$Gray2009)
listAUS<-listAUS[-which(listAUS=="0")]

taxaAUS<-read.csv("withinFamilyTMRCA/taxaGray2009.csv")
### get the minimum age between pairs of languages, internal clades that are parents to several children languages
ie_AUSTR.df <- read.csv("AustronesianalberoLangPhylo4.csv") ## file provided by Simon G. Nodes Time divergence and 95% probability. 

# get the clade size by counting the languages in the node names:
ie_AUSTR.df <-  mutate(ie_AUSTR.df, clade_size = str_count(node, fixed(",")))

listsplit<-strsplit(as.character(ie_AUSTR.df$node),",")

#returnglotto<-function(x){(match(x,taxaAUS$taxon))}
returnglotto<-function(x){taxaAUS$glottocode[match(x,taxaAUS$taxon)]}
# names in glottocode

ie_AUSTR.df$nodesGlotto<-lapply(listsplit, returnglotto)
# ie_AUSTR.df$nodesGlotto<-lapply(ie_AUSTR.df$nodesGlotto, as.character)

get_smallest_clade_age("Rennellese", "Ogan",ie_AUSTR.df )

# listaglot<-unique(levels(ie_AUSTR.df$node[[1]]))

meltFstREDinfo$Gray2009_DivTime<-as.numeric(meltFstREDinfo$Gray2009_DivTime)
meltFstREDGray2009<-meltFstREDinfo[!is.na(meltFstREDinfo$Gray2009_DivTime),]
meltFstREDGray2009<-meltFstREDGray2009[!is.na(meltFstREDGray2009$TMRCA_doubleNe),] # take only the comparisons where i have a linguistic div time and a genetic div time


addendumLang<-c()
for (i in 1:nrow(meltFstREDGray2009)){
  pop1<-meltFstREDGray2009$Pop1[i]
  pop2<-meltFstREDGray2009$Pop2[i]
  glottmatch1<-perpopRED$Gray2009[match(pop1,perpopRED$PopName)]
  glottmatch2<-perpopRED$Gray2009[match(pop2,perpopRED$PopName)]
  addendumLang1<-  get_smallest_clade_age(glottmatch1, glottmatch2,ie_AUSTR.df )
  addendumLang<-rbind(addendumLang,addendumLang1)
}

for (i in 1:nrow(meltFstREDGray2009)){
  pop1<-meltFstREDGray2009$Pop1[i]
  pop2<-meltFstREDGray2009$Pop2[i]
  glottmatch1<-perpopRED$Gray2009[match(pop1,perpopRED$PopName)]
  glottmatch2<-perpopRED$Gray2009[match(pop2,perpopRED$PopName)]
  temptime<- addendumLang[intersect( grep(glottmatch1, addendumLang$Language1) , grep(glottmatch2, addendumLang$Language2)),]
  meltFstREDGray2009$LinguisticDivergenceTime_median[i]<-mean(temptime$LinguisticDivergenceTime_median)*1000
  meltFstREDGray2009$LinguisticDivergenceTime_lower[i]<-mean(temptime$LinguisticDivergenceTime_lower)*1000
  meltFstREDGray2009$LinguisticDivergenceTime_upper[i]<-mean(temptime$LinguisticDivergenceTime_upper)*1000
}


## now plot
outliers<-c( "Mamanwa", "Rennell_and_Bellona", "Mamanwa1")
meltFstREDGray2009<-meltFstREDGray2009[-which(meltFstREDGray2009$Pop1%in%outliers),]
meltFstREDGray2009<-meltFstREDGray2009[-which(meltFstREDGray2009$Pop2%in%outliers),]

#adjust for CI which expand out of the limit of the y axis
maxTMRCA<-20000
meltFstREDGray2009$TMRCA_doubleNe_95[which(meltFstREDGray2009$TMRCA_doubleNe_95>maxTMRCA)]<-maxTMRCA
meltFstREDGray2009<-meltFstREDGray2009[which(meltFstREDGray2009$case=="single"),] # don't need double values to plot

colorino<-MainFamilies2$COLOR[which(MainFamilies2$MainFamilies==meltFstREDGray2009$FAMILY[1])]

gg<-ggplot(meltFstREDGray2009,aes(LinguisticDivergenceTime_median,TMRCA_doubleNe))
AUSTR<-gg+
  ylim(0,20000)+
  xlim(0,6000)+
  geom_errorbar(aes(ymin=TMRCA_doubleNe_5, ymax=TMRCA_doubleNe_95,),size=3,width=3,
                alpha=0.1)+
  geom_errorbarh(aes(xmin=LinguisticDivergenceTime_lower, xmax=LinguisticDivergenceTime_upper),size=3,height=3,
                 alpha=0.1)+
  
  geom_point(size=3,alpha=0.7, fill=colorino, shape=21, color="black")+
  # geom_text(aes(label=popslistemp), size=1)+
  xlab("Time distance from language tree - years ago")+
  ylab("Time distance from genetic data - years ago")+
  geom_abline(slope=1, intercept = 0, alpha=0.5)+theme_light()
# ggtitle(meltFstREDGray2009$FAMILY[1])+theme(plot.title = element_text(color = colorino))

# ggsave("correlationTimeGray2009_Austronesian_noOutlierMamanwa_RennellBelloneMINI_doubleBAR.pdf", useDingbats=FALSE, height = 5, width = 5)

# #****************************************
# # INDO EUROPEAN
# #****************************************
# 

TARGETFamily<-"Indo-European"

listAUS<-unique(perpopRED$Bouckaert2012)
listAUS<-listAUS[-which(listAUS=="0")]
taxaAUS<-read.table("withinFamilyTMRCA/Bouckaert2012_taxa.csv", sep=";", header=T, quote = "")
# from "https://github.com/D-PLACE/dplace-data/blob/master/phylogenies/bouckaert_et_al2012/taxa.csv"
# 


ie_IE.df <- read.csv("Indo-EuropeanalberoLangPhylo4BOUCKAERT.csv") 
# get the clade size by counting the languages in the node names:
ie_IE.df <-ie_IE.df %>% mutate(ie_IE.df, clade_size = str_count(node, fixed(",")))

listsplit<-strsplit(as.character(ie_IE.df$node),",")

returnglotto<-function(x){taxaAUS$glottocode[match(x,taxaAUS$taxon)]}
# names in glottocode


ie_IE.df$nodesGlotto<-lapply(listsplit, returnglotto)
ie_IE.df$nodesGlotto<-lapply(ie_IE.df$nodesGlotto, as.character)


listaglot<-unique(levels(ie_IE.df$node[[1]]))


meltFstREDinfo$bouckaert2012_DivTime<-as.numeric(meltFstREDinfo$bouckaert2012_DivTime)
meltFstREDBouckaert<-meltFstREDinfo[!is.na(meltFstREDinfo$bouckaert2012_DivTime),]
meltFstREDBouckaert<-meltFstREDBouckaert[!is.na(meltFstREDBouckaert$TMRCA_doubleNe),] # take only the comparisons where i have a linguistic div time and a genetic div time


addendumLang<-c()
for (i in 1:nrow(meltFstREDBouckaert)){
  pop1<-meltFstREDBouckaert$Pop1[i]
  pop2<-meltFstREDBouckaert$Pop2[i]
  glottmatch1<-perpopRED$Bouckaert2012 [match(pop1,perpopRED$PopName)]
  glottmatch2<-perpopRED$Bouckaert2012[match(pop2,perpopRED$PopName)]
  addendumLang1<-  get_smallest_clade_age(glottmatch1, glottmatch2,ie_IE.df )
  addendumLang<-rbind(addendumLang,addendumLang1)
}

for (i in 1:nrow(meltFstREDBouckaert)){
  pop1<-meltFstREDBouckaert$Pop1[i]
  pop2<-meltFstREDBouckaert$Pop2[i]
  glottmatch1<-perpopRED$Bouckaert2012[match(pop1,perpopRED$PopName)]
  glottmatch2<-perpopRED$Bouckaert2012[match(pop2,perpopRED$PopName)]
  temptime<- addendumLang[intersect( grep(glottmatch1, addendumLang$Language1) , grep(glottmatch2, addendumLang$Language2)),]
  meltFstREDBouckaert$LinguisticDivergenceTime_median[i]<-mean(temptime$LinguisticDivergenceTime_median)
  meltFstREDBouckaert$LinguisticDivergenceTime_lower[i]<-mean(temptime$LinguisticDivergenceTime_lower)
  meltFstREDBouckaert$LinguisticDivergenceTime_upper[i]<-mean(temptime$LinguisticDivergenceTime_upper)
}



## now plot

outliers<-c( "Sardinian") # Sardinians are too genetically divergent
meltFstREDBouckaert<-meltFstREDBouckaert[-which(meltFstREDBouckaert$Pop1%in%outliers),]
meltFstREDBouckaert<-meltFstREDBouckaert[-which(meltFstREDBouckaert$Pop2%in%outliers),]

maxTMRCA<-11500 # for the plot area
meltFstREDBouckaert$TMRCA_doubleNe_95[which(meltFstREDBouckaert$TMRCA_doubleNe_95>maxTMRCA)]<-maxTMRCA
meltFstREDBouckaert<-meltFstREDBouckaert[which(meltFstREDBouckaert$case=="single"),] # don't need double values to plot

colorino<-MainFamilies2$COLOR[which(MainFamilies2$MainFamilies==meltFstREDBouckaert$FAMILY[1])]

gg<-ggplot(meltFstREDBouckaert,aes(LinguisticDivergenceTime_median,TMRCA_doubleNe))
BOUCKAERT<-gg+
  xlim(0,7500)+
  ylim(0,maxTMRCA)+
  geom_errorbar(aes(ymin=TMRCA_doubleNe_5, ymax=TMRCA_doubleNe_95,),size=3,width=3,
                alpha=0.1)+
  geom_errorbarh(aes(xmin=LinguisticDivergenceTime_lower, xmax=LinguisticDivergenceTime_upper),size=3,height=3,
                 alpha=0.1)+
  
  geom_point(size=3,alpha=0.7, fill=colorino, shape=21, color="black")+
  # geom_text(aes(label=popslistemp), size=1)+
  xlab("Time distance from language tree - years ago")+
  ylab("Time distance from genetic data - years ago")+
  geom_abline(slope=1, intercept = 0, alpha=0.5)+theme_light()
# ggtitle(TARGETFamily)+theme(plot.title = element_text(color = colorino))

# ggsave("correlationTimeBouckaertIE_noSardiniaMINI_doubleBAR.pdf", useDingbats=FALSE, height = 5, width = 5)


#******************************************
# Chang Data
#*******************************************************************************************************
#*# CHANG ET AL 2015
#*************************

listAUS<-unique(perpopRED$Chang2015)
listAUS<-listAUS[-which(listAUS=="0")]

AUS<-read.nexus("withinFamilyTMRCA/Chang2015.trees")
# from https://github.com/D-PLACE/dplace-data/blob/master/phylogenies/chang_et_al2015/summary.trees

taxaAUS<-read.table("withinFamilyTMRCA/Chang2015_taxa.csv", sep=";", header=T, quote = "")
# from https://github.com/D-PLACE/dplace-data/blob/master/phylogenies/chang_et_al2015/taxa.csv


#*************************
# manually elaborate
meltFstREDChang<-read.table("meltFstREDChang.txt",sep="\t", as.is=T, header=T)
## patch CI lang

gg<-ggplot(meltFstREDChang,aes(chang2015_DivTime,TMRCA_doubleNe))

CHANG<-gg+
  ylim(0,maxTMRCA)+
  xlim(0,7500)+
  geom_errorbar(aes(ymin=TMRCA_doubleNe_5, ymax=TMRCA_doubleNe_95),size=3,width=3,
                alpha=0.1)+
  geom_errorbarh(aes(xmin=LinguisticDivergenceTime_lower, xmax=LinguisticDivergenceTime_upper),size=3,height=3,
                 alpha=0.1)+
  
  geom_point(size=3,alpha=0.7, fill=colorino, shape=21, color="black")+
  geom_text(aes(label=popslistemp), size=1)+
  xlab("Time distance from language tree - years ago")+
  ylab("Time distance from genetic data - years ago")+
  geom_abline(slope=1, intercept = 0, alpha=0.5)+theme_light()+
  ggtitle(meltFstREDChang$FAMILY[1])+theme(plot.title = element_text(color = colorino))



# #******************************************
# # Turkic
# #******************************************
# #*
# #* Hrushka 2015
# #*  
# 
TARGETFamily<-"Turkic"

listAUS<-unique(perpopRED$Hruschka2015)
listAUS<-listAUS[-which(listAUS=="0")]
taxaAUS<-read.csv("withinFamilyTMRCA/Hrushka_taxa2.csv", header=T, quote = "")
# from https://github.com/D-PLACE/dplace-data/blob/master/phylogenies/hruschka_et_al2015/taxa.csv
colnames(taxaAUS)[1]<-"code"
colnames(taxaAUS)[2]<-"taxon"
#
# 

ie_TURK.df <- read.csv("TurkicalberoLangPhylo4.csv") 
# get the clade size by counting the languages in the node names:
ie_TURK.df <-ie_TURK.df %>% mutate(ie_TURK.df, clade_size = str_count(node, fixed(",")))

listsplit<-strsplit(as.character(ie_TURK.df$node),",")
returnglotto<-function(x){taxaAUS$glottocode[match(x,taxaAUS$code)]}

ie_TURK.df$nodesGlotto<-lapply(listsplit, returnglotto)
ie_TURK.df$nodesGlotto<-lapply(ie_TURK.df$nodesGlotto, as.character)

meltFstREDinfo$Hruschka2015_DivTime<-as.numeric(meltFstREDinfo$Hruschka2015_DivTime)
meltFstREDHrushka<-meltFstREDinfo[!is.na(meltFstREDinfo$Hruschka2015_DivTime),]
meltFstREDHrushka<-meltFstREDHrushka[!is.na(meltFstREDHrushka$TMRCA_doubleNe),]


addendumLang<-c()
for (i in 1:nrow(meltFstREDHrushka)){
  pop1<-meltFstREDHrushka$Pop1[i]
  pop2<-meltFstREDHrushka$Pop2[i]
  glottmatch1<-perpopRED$Hruschka2015 [match(pop1,perpopRED$PopName)]
  glottmatch2<-perpopRED$Hruschka2015[match(pop2,perpopRED$PopName)]
  addendumLang1<-  get_smallest_clade_age(glottmatch1, glottmatch2,ie_TURK.df )
  addendumLang<-rbind(addendumLang,addendumLang1)
}

for (i in 1:nrow(meltFstREDHrushka)){
  pop1<-meltFstREDHrushka$Pop1[i]
  pop2<-meltFstREDHrushka$Pop2[i]
  glottmatch1<-perpopRED$Hruschka2015 [match(pop1,perpopRED$PopName)]
  glottmatch2<-perpopRED$Hruschka2015[match(pop2,perpopRED$PopName)]
  temptime<- addendumLang[intersect( grep(glottmatch1, addendumLang$Language1) , grep(glottmatch2, addendumLang$Language2)),]
  meltFstREDHrushka$LinguisticDivergenceTime_median[i]<-mean(temptime$LinguisticDivergenceTime_median)
  meltFstREDHrushka$LinguisticDivergenceTime_lower[i]<-mean(temptime$LinguisticDivergenceTime_lower)
  meltFstREDHrushka$LinguisticDivergenceTime_upper[i]<-mean(temptime$LinguisticDivergenceTime_upper)
}



#adjust for CI which expand out of the limit of the y axis
maxTMRCA<-22000
meltFstREDHrushka$TMRCA_doubleNe_95[which(meltFstREDHrushka$TMRCA_doubleNe_95>maxTMRCA)]<-maxTMRCA
meltFstREDHrushka<-meltFstREDHrushka[which(meltFstREDHrushka$case=="single"),] # don't need double values to plot
colorino<-MainFamilies2$COLOR[which(MainFamilies2$MainFamilies==meltFstREDHrushka$FAMILY[1])]

gg<-ggplot(meltFstREDHrushka,aes(LinguisticDivergenceTime_median,TMRCA_doubleNe))

HRUS<-gg+
  xlim(0,3000)+
  
  geom_errorbar(aes(ymin=TMRCA_doubleNe_5, ymax=TMRCA_doubleNe_95,),size=3,width=3,
                alpha=0.1)+
  geom_errorbarh(aes(xmin=LinguisticDivergenceTime_lower, xmax=LinguisticDivergenceTime_upper),size=3,height=3,
                 alpha=0.1)+
  
  geom_point(size=3,alpha=0.7, fill=colorino, shape=21, color="black")+
  # geom_text(aes(label=popslistemp), size=1)+
  xlab("Time distance from language tree - years ago")+
  ylab("Time distance from genetic data - years ago")+
  geom_abline(slope=1, intercept = 0, alpha=0.5)+theme_light()
#  ggtitle(TARGETFamily)+theme(plot.title = element_text(color = colorino))

# ggsave("correlationTimeturkicMINI_doubleBAR.pdf", useDingbats=FALSE, height = 5, width = 5)


#******************************************
#* Savalyev 2020
#* 
#* TARGETFamily<-"Turkic"
#******************************************
#* Savalyev and Robbeets 2020, manual adding from Simon's file
#* 

meltFstREDSavalyev<-read.table("meltFstREDSavalyev.txt",sep="\t", as.is=T, header=T)
## patch CI lang

#adjust for CI which expand out of the limit of the y axis
maxTMRCA<-22000
meltFstREDSavalyev$TMRCA_doubleNe_95[which(meltFstREDSavalyev$TMRCA_doubleNe_95>maxTMRCA)]<-maxTMRCA
colorino<-MainFamilies2$COLOR[which(MainFamilies2$MainFamilies==TARGETFamily)]

gg<-ggplot(meltFstREDSavalyev,aes(Savelyev2020_DivTime,TMRCA_doubleNe))


SAVAL<-gg+
  xlim(0,3000)+
  geom_errorbar(aes(ymin=TMRCA_doubleNe_5, ymax=TMRCA_doubleNe_95),size=3,width=3,
                alpha=0.1)+
  geom_errorbarh(aes(xmin=LinguisticDivergenceTime_lower, xmax=LinguisticDivergenceTime_upper),size=3,height=3,
                 alpha=0.1)+
  
  geom_point(size=3,alpha=0.7, fill=colorino, shape=21, color="black")+
  geom_text(aes(label=popslistemp), size=1)+
  xlab("Time distance from language tree - years ago")+
  ylab("Time distance from genetic data - years ago")+
  geom_abline(slope=1, intercept = 0, alpha=0.5)+
  theme_light()+
  ggtitle(TARGETFamily)+theme(plot.title = element_text(color = colorino))


#*****************************************************
## MAIN FIGURE  COMPARISON IE, AUSTR, TURKIC 
# - FIGURE 4
#*****************************************************
#*
library(ggpubr)
ggarrange(BOUCKAERT, AUSTR, HRUS  + rremove("x.text"), 
          labels = c("D", "E", "F"),
          ncol = 1, nrow = 3)
ggsave("combined3LangFamiliesCorrelation_Fig4_2022.pdf", useDingbats=FALSE, height = 12, width = 4)

#*****************************************************
## SUPPLEMENTARY COMPARISON CHANG AND SAVALYEV
#*****************************************************
#*
library(ggpubr)
ggarrange(CHANG, SAVAL + rremove("x.text"), 
          labels = c("A.", "B."),
          ncol = 2, nrow = 1)
ggsave("Fig_S12_combinedChangAndSavelyev_2022_errorbar.pdf", useDingbats=FALSE, height = 5, width = 10)



