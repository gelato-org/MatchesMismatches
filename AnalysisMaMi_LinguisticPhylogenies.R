
#######################################
# SUBSET ANALYSIS COMPARING GELATO WITH LANGUAGE DIVERGENCE TIMES
#######################################
#
#
#
# Chiara Barbieri 2021

#list of populations
perpopRED<-read.table("/Users/chiarabarbieri/Library/Mobile Documents/com~apple~CloudDocs/GeLaTo/perpopREDMaMi2021.txt",   sep = "\t", header=T, as.is=T)
#list of pairwise distances
FstListREDinfo<-read.table("/Users/chiarabarbieri/Library/Mobile Documents/com~apple~CloudDocs/GeLaTo/FstListREDinfo2021.txt",header = T, as.is=T, sep="\t")
# color for major language families

Lfamil<-table(perpopRED$glottolog.NAME)
MainFamilies<-unlist(labels(Lfamil[which(Lfamil>4)])) # minimum 5 populations per Lang Family
#MainFamilies<-MainFamilies[-which(MainFamilies=="LI")]   # exclude the Language Isolates LI
perpopREDfamily<-perpopRED[which(perpopRED$glottolog.NAME %in% MainFamilies),]
colorchoice<-c( "darkorange4" ,"#A6CEE3"  ,   "#33A02C"    , "#A6D854"  ,   "#377EB8"  ,   "#E31A1C"   ,  "#FFD92F" ,    "#FF7F00"  ,   "#666666" ,   
                "cyan4"  ,     "#BC80BD"   ,  "#FED9A6" ,    "tan3" ,       "#6A3D9A" ,    "deeppink"   )
MainFamilies2<-as.data.frame(MainFamilies)
MainFamilies2$COLOR<-colorchoice

library(phylobase)
library(pegas)
library(ggplot2)


meltFstREDinfo<-FstListREDinfo

#******************************************
#*


# -----------------------------------------------------------------
### put all the dates from all the trees in the FST file

#******************************************************************************
# INDO EUROPEAN
#******************************************************************************
#* Bouckaert 2012
#* 
#* 
TARGETFamily<-"Indo-European"

listAUS<-unique(perpopRED$Bouckaert2012)
listAUS<-listAUS[-which(listAUS=="")]
AUS<-read.nexus("withinFamilyTMRCA/Bouckaert2012.trees")
 # from "https://github.com/D-PLACE/dplace-data/blob/master/phylogenies/bouckaert_et_al2012/summary.trees"

taxaAUS<-read.table("withinFamilyTMRCA/Bouckaert2012_taxa.csv", sep=";", header=T, quote = "")
# from "https://github.com/D-PLACE/dplace-data/blob/master/phylogenies/bouckaert_et_al2012/taxa.csv"

#*************************

listAUSreverse<-as.character(taxaAUS$taxon[match(listAUS,taxaAUS$glottocode)])
## remove duplicated languages

g1 <- as(AUS, "phylo4")
bi_subset <- subset(g1, tips.include=listAUSreverse) # select only the tips i have. note, i have duplicated tips.
alberoLang<-as(bi_subset,"phylo")

alberoLang<-makeNodeLabel(alberoLang, method = "number", prefix = "")
#write.tree(alberoLang,paste0(TARGETFamily ,"alberoLangPhylo4.phy"))
alberoLangPhylo4 <- as(alberoLang, "phylo4")
alberoLang$tip.label<-taxaAUS$glottocode[match(alberoLang$tip.label,taxaAUS$taxon)] # names in glottocode
write.tree(alberoLang,paste0(TARGETFamily ,"alberoLangPhylo4.phy"))

distancealberoLang<-cophenetic.phylo(alberoLang) # names in glottocode
distancealberoLang<-distancealberoLang/2 # remember the distances are double!!

meltdistanceAUSunique<-melt(distancealberoLang)
colnames(meltdistanceAUSunique)<-c("proxy1", "proxy2", "LangDivergenceTime")

meltdistanceAUSunique$pop1<-perpopRED$PopName[match(meltdistanceAUSunique$proxy1,perpopRED$Bouckaert2012)]
meltdistanceAUSunique$pop2<-perpopRED$PopName[match(meltdistanceAUSunique$proxy2,perpopRED$Bouckaert2012)]

meltFstREDinfo$langproxy1<-perpopRED$Bouckaert2012[match(meltFstREDinfo$Pop1,perpopRED$PopName)]
meltFstREDinfo$langproxy2<-perpopRED$Bouckaert2012[match(meltFstREDinfo$Pop2,perpopRED$PopName)]

keeperNOLANG<-union(which(meltFstREDinfo$langproxy1==""),which(meltFstREDinfo$langproxy2==""))
meltFstREDinfoTESt<-meltFstREDinfo[-keeperNOLANG,]

# loop to assign language divergence time corresponding to genetic glottocodes pairs
for (i in 1:nrow(meltFstREDinfoTESt)){
  meltFstREDinfoTESt$lang_DivTime[i]<-meltdistanceAUSunique[which(meltdistanceAUSunique$proxy1==meltFstREDinfoTESt$langproxy1[i]&
                                                                    meltdistanceAUSunique$proxy2==meltFstREDinfoTESt$langproxy2[i]),3]
}
meltFstREDinfo$bouckaert2012_DivTime<-"NA"
meltFstREDinfo$bouckaert2012_DivTime[-keeperNOLANG]<-as.numeric(meltFstREDinfoTESt$lang_DivTime)  # create the column lang divergence time in the list of fst file
meltFstREDinfo<-meltFstREDinfo[,-which(colnames(meltFstREDinfo)%in%c("langproxy1","langproxy2"))]


#*******************************************************************************************************
#*# CHANG ET AL 2015
#*************************

listAUS<-unique(perpopRED$Chang2015)
listAUS<-listAUS[-which(listAUS=="")]

AUS<-read.nexus("withinFamilyTMRCA/Chang2015.trees")
# from https://github.com/D-PLACE/dplace-data/blob/master/phylogenies/chang_et_al2015/summary.trees

taxaAUS<-read.table("withinFamilyTMRCA/Chang2015_taxa.csv", sep=";", header=T, quote = "")
 # from https://github.com/D-PLACE/dplace-data/blob/master/phylogenies/chang_et_al2015/taxa.csv

#*************************

listAUSreverse<-as.character(taxaAUS$taxon[match(listAUS,taxaAUS$glottocode)])
## remove duplicated languages

g1 <- as(AUS, "phylo4")
bi_subset <- subset(g1, tips.include=listAUSreverse) # select only the tips i have. note, i have duplicated tips.
alberoLang<-as(bi_subset,"phylo")

alberoLang<-makeNodeLabel(alberoLang, method = "number", prefix = "")

alberoLangPhylo4 <- as(alberoLang, "phylo4")

alberoLang$tip.label<-taxaAUS$glottocode[match(alberoLang$tip.label,taxaAUS$taxon)] # names in glottocode
write.tree(alberoLang,paste0(TARGETFamily ,"CHANGalberoLangPhylo4.phy"))

distancealberoLang<-cophenetic.phylo(alberoLang) # names in glottocode
distancealberoLang<-distancealberoLang/2 # remember the distances are double!!
#distancealberoLang<-distancealberoLang*1000 # remember the distances are double!! Austronesian case

meltdistanceAUSunique<-melt(distancealberoLang)
colnames(meltdistanceAUSunique)<-c("proxy1", "proxy2", "LangDivergenceTime")

meltdistanceAUSunique$pop1<-perpopRED$PopName[match(meltdistanceAUSunique$proxy1,perpopRED$Chang2015)]
meltdistanceAUSunique$pop2<-perpopRED$PopName[match(meltdistanceAUSunique$proxy2,perpopRED$Chang2015)]
meltFstREDinfo$langproxy1<-perpopRED$Chang2015[match(meltFstREDinfo$Pop1,perpopRED$PopName)]
meltFstREDinfo$langproxy2<-perpopRED$Chang2015[match(meltFstREDinfo$Pop2,perpopRED$PopName)]

keeperNOLANG<-union(which(meltFstREDinfo$langproxy1==""),which(meltFstREDinfo$langproxy2==""))
meltFstREDinfoTESt<-meltFstREDinfo[-keeperNOLANG,]

# loop to assign language divergence time corresponding to genetic glottocodes pairs
for (i in 1:nrow(meltFstREDinfoTESt)){
  meltFstREDinfoTESt$lang_DivTime[i]<-meltdistanceAUSunique[which(meltdistanceAUSunique$proxy1==meltFstREDinfoTESt$langproxy1[i]&
                                                                    meltdistanceAUSunique$proxy2==meltFstREDinfoTESt$langproxy2[i]),3]
}
meltFstREDinfo$chang2015_DivTime<-"NA"
meltFstREDinfo$chang2015_DivTime[-keeperNOLANG]<-as.numeric(meltFstREDinfoTESt$lang_DivTime)  # create the column lang divergence time in the list of fst file
meltFstREDinfo<-meltFstREDinfo[,-which(colnames(meltFstREDinfo)%in%c("langproxy1","langproxy2"))]


#******************************************
# Austronesian
#******************************************

TARGETFamily<-"Austronesian"

listAUS<-unique(perpopRED$Gray2009)
listAUS<-listAUS[-which(listAUS=="")]

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

alberoLang$edge.length<-alberoLang$edge.length*1000
#write.tree(alberoLang,paste0(TARGETFamily ,"alberoLangPhylo4.phy"))

alberoLangPhylo4 <- as(alberoLang, "phylo4")
alberoLang$tip.label<-taxaAUS$glottocode[match(alberoLang$tip.label,taxaAUS$taxon)] # names in glottocode
write.tree(alberoLang,paste0(TARGETFamily ,"alberoLangPhylo4.phy"))

distancealberoLang<-cophenetic.phylo(alberoLang) # names in glottocode
distancealberoLang<-distancealberoLang/2 # remember the distances are double!!
#distancealberoLang<-distancealberoLang*1000 # multiply the distances to get years -  Austronesian case

meltdistanceAUSunique<-melt(distancealberoLang)
colnames(meltdistanceAUSunique)<-c("proxy1", "proxy2", "LangDivergenceTime")

meltdistanceAUSunique$pop1<-perpopRED$PopName[match(meltdistanceAUSunique$proxy1,perpopRED$Gray2009)]
meltdistanceAUSunique$pop2<-perpopRED$PopName[match(meltdistanceAUSunique$proxy2,perpopRED$Gray2009)]

meltFstREDinfo$langproxy1<-perpopRED$Gray2009[match(meltFstREDinfo$Pop1,perpopRED$PopName)]
meltFstREDinfo$langproxy2<-perpopRED$Gray2009[match(meltFstREDinfo$Pop2,perpopRED$PopName)]

keeperNOLANG<-union(which(meltFstREDinfo$langproxy1==""),which(meltFstREDinfo$langproxy2==""))
meltFstREDinfoTESt<-meltFstREDinfo[-keeperNOLANG,]

# loop to assign language divergence time corresponding to genetic glottocodes pairs
for (i in 1:nrow(meltFstREDinfoTESt)){
  meltFstREDinfoTESt$lang_DivTime[i]<-meltdistanceAUSunique[which(meltdistanceAUSunique$proxy1==meltFstREDinfoTESt$langproxy1[i]&
                                                                    meltdistanceAUSunique$proxy2==meltFstREDinfoTESt$langproxy2[i]),3]
}
meltFstREDinfo$Gray2009_DivTime<-"NA"
meltFstREDinfo$Gray2009_DivTime[-keeperNOLANG]<-as.numeric(meltFstREDinfoTESt$lang_DivTime)  # create the column lang divergence time in the list of fst file
meltFstREDinfo<-meltFstREDinfo[,-which(colnames(meltFstREDinfo)%in%c("langproxy1","langproxy2"))]


#******************************************
# Turkic
#******************************************
#*
#* Hrushka 2015
#*  

TARGETFamily<-"Turkic"

listAUS<-unique(perpopRED$Hruschka2015)
listAUS<-listAUS[-1]
AUS<-read.nexus("withinFamilyTMRCA/Hrushka_summary2.trees")
#from https://github.com/D-PLACE/dplace-data/blob/master/phylogenies/hruschka_et_al2015/summary.trees
taxaAUS<-read.csv("withinFamilyTMRCA/Hrushka_taxa2.csv", header=T, quote = "")
# from https://github.com/D-PLACE/dplace-data/blob/master/phylogenies/hruschka_et_al2015/taxa.csv
colnames(taxaAUS)[1]<-"code"
colnames(taxaAUS)[2]<-"taxon"
AUS$tip.label<-as.character(taxaAUS$taxon[match(AUS$tip.label, taxaAUS$code)])

##

listAUSreverse<-as.character(taxaAUS$taxon[match(listAUS,taxaAUS$glottocode)])
## remove duplicated languages

g1 <- as(AUS, "phylo4")
bi_subset <- subset(g1, tips.include=listAUSreverse) # select only the tips i have. note, i have duplicated tips.
alberoLang<-as(bi_subset,"phylo")

alberoLang<-makeNodeLabel(alberoLang, method = "number", prefix = "")
#write.tree(alberoLang,paste0(TARGETFamily ,"alberoLangPhylo4.phy"))

alberoLangPhylo4 <- as(alberoLang, "phylo4")

alberoLang$tip.label<-taxaAUS$glottocode[match(alberoLang$tip.label,taxaAUS$taxon)] # names in glottocode
write.tree(alberoLang,paste0(TARGETFamily ,"alberoLangPhylo4.phy"))

distancealberoLang<-cophenetic.phylo(alberoLang) # names in glottocode
distancealberoLang<-distancealberoLang/2 # remember the distances are double!!

meltdistanceAUSunique<-melt(distancealberoLang)
colnames(meltdistanceAUSunique)<-c("proxy1", "proxy2", "LangDivergenceTime")

meltdistanceAUSunique$pop1<-perpopRED$PopName[match(meltdistanceAUSunique$proxy1,perpopRED$Hruschka2015)]
meltdistanceAUSunique$pop2<-perpopRED$PopName[match(meltdistanceAUSunique$proxy2,perpopRED$Hruschka2015)]

meltFstREDinfo$langproxy1<-perpopRED$Hruschka2015[match(meltFstREDinfo$Pop1,perpopRED$PopName)]
meltFstREDinfo$langproxy2<-perpopRED$Hruschka2015[match(meltFstREDinfo$Pop2,perpopRED$PopName)]

keeperNOLANG<-union(which(meltFstREDinfo$langproxy1==""),which(meltFstREDinfo$langproxy2==""))
meltFstREDinfoTESt<-meltFstREDinfo[-keeperNOLANG,]

# loop to assign language divergence time corresponding to genetic glottocodes pairs
for (i in 1:nrow(meltFstREDinfoTESt)){
  meltFstREDinfoTESt$lang_DivTime[i]<-meltdistanceAUSunique[which(meltdistanceAUSunique$proxy1==meltFstREDinfoTESt$langproxy1[i]&
                                                                    meltdistanceAUSunique$proxy2==meltFstREDinfoTESt$langproxy2[i]),3]
}
meltFstREDinfo$Hruschka2015_DivTime<-"NA"
meltFstREDinfo$Hruschka2015_DivTime[-keeperNOLANG]<-as.numeric(meltFstREDinfoTESt$lang_DivTime) # create the column lang divergence time in the list of fst file
meltFstREDinfo<-meltFstREDinfo[,-which(colnames(meltFstREDinfo)%in%c("langproxy1","langproxy2"))]

#******************************************
#* Savalyev and Robbeets 2020
#* 

TARGETFamily<-"Turkic"

listAUS<-unique(perpopRED$savelyev2020)
listAUS<-listAUS[-1]
AUS<-read.nexus("withinFamilyTMRCA/Savalyev_summary.trees")
taxaAUS<-read.csv("withinFamilyTMRCA/savalyev_taxa.csv", header=T, quote = "")

##

listAUSreverse<-as.character(taxaAUS$taxon[match(listAUS,taxaAUS$glottocode)])
## remove duplicated languages

g1 <- as(AUS, "phylo4")
bi_subset <- subset(g1, tips.include=listAUSreverse) # select only the tips i have. note, i have duplicated tips.
alberoLang<-as(bi_subset,"phylo")

alberoLang<-makeNodeLabel(alberoLang, method = "number", prefix = "")
#write.tree(alberoLang,paste0(TARGETFamily ,"SAVALYEValberoLangPhylo4.phy"))

alberoLangPhylo4 <- as(alberoLang, "phylo4")

alberoLang$tip.label<-taxaAUS$glottocode[match(alberoLang$tip.label,taxaAUS$taxon)] # names in glottocode
write.tree(alberoLang,paste0(TARGETFamily ,"SAVALYEValberoLangPhylo4.phy"))

distancealberoLang<-cophenetic.phylo(alberoLang) # names in glottocode
distancealberoLang<-distancealberoLang/2 # remember the distances are double!!
distancealberoLang<-distancealberoLang*1000 # adjust to make the years

meltdistanceAUSunique<-melt(distancealberoLang)
colnames(meltdistanceAUSunique)<-c("proxy1", "proxy2", "LangDivergenceTime")

meltdistanceAUSunique$pop1<-perpopRED$PopName[match(meltdistanceAUSunique$proxy1,perpopRED$savelyev2020)]
meltdistanceAUSunique$pop2<-perpopRED$PopName[match(meltdistanceAUSunique$proxy2,perpopRED$savelyev2020)]

meltFstREDinfo$langproxy1<-perpopRED$savelyev2020[match(meltFstREDinfo$Pop1,perpopRED$PopName)]
meltFstREDinfo$langproxy2<-perpopRED$savelyev2020[match(meltFstREDinfo$Pop2,perpopRED$PopName)]

keeperNOLANG<-union(which(meltFstREDinfo$langproxy1==""),which(meltFstREDinfo$langproxy2==""))
meltFstREDinfoTESt<-meltFstREDinfo[-keeperNOLANG,]

# loop to assign language divergence time corresponding to genetic glottocodes pairs
for (i in 1:nrow(meltFstREDinfoTESt)){
  meltFstREDinfoTESt$lang_DivTime[i]<-meltdistanceAUSunique[which(meltdistanceAUSunique$proxy1==meltFstREDinfoTESt$langproxy1[i]&
                                                                    meltdistanceAUSunique$proxy2==meltFstREDinfoTESt$langproxy2[i]),3]
}
meltFstREDinfo$Savelyev2020_DivTime<-"NA"
meltFstREDinfo$Savelyev2020_DivTime[-keeperNOLANG]<-as.numeric(meltFstREDinfoTESt$lang_DivTime)  # create the column lang divergence time in the list of fst file
meltFstREDinfo<-meltFstREDinfo[,-which(colnames(meltFstREDinfo)%in%c("langproxy1","langproxy2"))]


# --------------------------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------
#### PLOT CORRELATIONS SIMPLE
# -----------------------------------------------------------------


meltFstREDinfoSINGLE<-meltFstREDinfo[which(meltFstREDinfo$case=="single"),]  # i need only one point for each pair
meltFstREDinfoSINGLE<-as.data.frame(meltFstREDinfoSINGLE)
meltFstREDinfoSINGLE$TMRCA_doubleNe<-as.numeric(meltFstREDinfoSINGLE$TMRCA_doubleNe)
meltFstREDinfoSINGLE$TMRCA_doubleNe_5<-as.numeric(meltFstREDinfoSINGLE$TMRCA_doubleNe_5)
meltFstREDinfoSINGLE$TMRCA_doubleNe_95<-as.numeric(meltFstREDinfoSINGLE$TMRCA_doubleNe_95)


#* # ***************************************************
### Indo-European comparisons divergence time gene language
# ***************************************************
#* Bouckaert 2012

meltFstREDinfoSINGLE$bouckaert2012_DivTime<-as.numeric(meltFstREDinfoSINGLE$bouckaert2012_DivTime)
meltFstREDBouckaert<-meltFstREDinfoSINGLE[!is.na(meltFstREDinfoSINGLE$bouckaert2012_DivTime),]
meltFstREDBouckaert<-meltFstREDBouckaert[!is.na(meltFstREDBouckaert$TMRCA_doubleNe),] # take only the comparisons where i have a linguistic div time and a genetic div time

outliers<-c( "Sardinian") # Sardinians are too genetically divergent
meltFstREDBouckaert<-meltFstREDBouckaert[-which(meltFstREDBouckaert$Pop1%in%outliers),]
meltFstREDBouckaert<-meltFstREDBouckaert[-which(meltFstREDBouckaert$Pop2%in%outliers),]

maxTMRCA<-11500 # for the plot area
meltFstREDBouckaert$TMRCA_doubleNe_95[which(meltFstREDBouckaert$TMRCA_doubleNe_95>maxTMRCA)]<-maxTMRCA

colorino<-MainFamilies2$COLOR[which(MainFamilies2$MainFamilies==meltFstREDBouckaert$FAMILY[1])]

gg<-ggplot(meltFstREDBouckaert,aes(bouckaert2012_DivTime,TMRCA_doubleNe))
BOUCKAERT<-gg+
  xlim(0,7500)+
  ylim(0,maxTMRCA)+
  geom_errorbar(aes(ymin=TMRCA_doubleNe_5, ymax=TMRCA_doubleNe_95),size=3,width=3,
                alpha=0.1)+
  geom_point(size=3,alpha=0.7, fill=colorino, shape=21, color="black")+
#  geom_text(aes(label=popslistemp), size=1)+
  xlab("Time distance from language tree - years ago")+
  ylab("Time distance from genetic data - years ago")+
  geom_abline(slope=1, intercept = 0, alpha=0.5)+theme_light()+
  ggtitle(meltFstREDBouckaert$FAMILY[1])+theme(plot.title = element_text(color = colorino))

# ggsave("correlationTimeBouckaertIE_noSardiniaMINI.pdf", useDingbats=FALSE, height = 5, width = 5)

#******************************************
# Chang Data
meltFstREDinfoSINGLE$chang2015_DivTime<-as.numeric(meltFstREDinfoSINGLE$chang2015_DivTime)

meltFstREDChang<-meltFstREDinfoSINGLE[!is.na(meltFstREDinfoSINGLE$chang2015_DivTime),]
meltFstREDChang<-meltFstREDChang[!is.na(meltFstREDChang$TMRCA_doubleNe),]

gg<-ggplot(meltFstREDChang,aes(chang2015_DivTime,TMRCA_doubleNe))

CHANG<-gg+
  ylim(0,maxTMRCA)+
  xlim(0,7500)+
  geom_errorbar(aes(ymin=TMRCA_doubleNe_5, ymax=TMRCA_doubleNe_95),size=3,width=3,
                alpha=0.1)+
  geom_point(size=3,alpha=0.7, fill=colorino, shape=21, color="black")+
  geom_text(aes(label=popslistemp), size=1)+
  xlab("Time distance from language tree - years ago")+
  ylab("Time distance from genetic data - years ago")+
  geom_abline(slope=1, intercept = 0, alpha=0.5)+theme_light()+
  ggtitle(meltFstREDChang$FAMILY[1])+theme(plot.title = element_text(color = colorino))

# ggsave("correlationTimeChangIE_MINI.pdf", useDingbats=FALSE, height = 5, width = 5)

#* # ***************************************************
### Turkic comparisons divergence time gene language
# ***************************************************
#* Hrushka 2012
meltFstREDinfoSINGLE$Hruschka2015_DivTime<-as.numeric(meltFstREDinfoSINGLE$Hruschka2015_DivTime)

meltFstREDHrushka<-meltFstREDinfoSINGLE[!is.na(meltFstREDinfoSINGLE$Hruschka2015_DivTime),]
meltFstREDHrushka<-meltFstREDHrushka[!is.na(meltFstREDHrushka$TMRCA_doubleNe),]

#adjust for CI which expand out of the limit of the y axis
maxTMRCA<-22000
meltFstREDHrushka$TMRCA_doubleNe_95[which(meltFstREDHrushka$TMRCA_doubleNe_95>maxTMRCA)]<-maxTMRCA
#  perpopREDIE<-perpopRED[which(perpopRED$Bouckaert2012!=""),]
colorino<-MainFamilies2$COLOR[which(MainFamilies2$MainFamilies==meltFstREDHrushka$FAMILY[1])]

gg<-ggplot(meltFstREDHrushka,aes(Hruschka2015_DivTime,TMRCA_doubleNe))

HRUS<-gg+
  xlim(0,3000)+
  geom_errorbar(aes(ymin=TMRCA_doubleNe_5, ymax=TMRCA_doubleNe_95),size=3,width=3,
                alpha=0.1)+
  geom_point(size=3,alpha=0.7, fill=colorino, shape=21, color="black")+
 # geom_text(aes(label=popslistemp), size=1)+
  xlab("Time distance from language tree - years ago")+
  ylab("Time distance from genetic data - years ago")+
  geom_abline(slope=1, intercept = 0, alpha=0.5)+ theme_light()+
  ggtitle(meltFstREDHrushka$FAMILY[1])+theme(plot.title = element_text(color = colorino))

#******************************************
#* Savalyev 2020
meltFstREDinfoSINGLE$Savelyev2020_DivTime<-as.numeric(meltFstREDinfoSINGLE$Savelyev2020_DivTime)

meltFstREDSavalyev<-meltFstREDinfoSINGLE[!is.na(meltFstREDinfoSINGLE$Savelyev2020_DivTime),]
meltFstREDSavalyev<-meltFstREDSavalyev[!is.na(meltFstREDSavalyev$TMRCA_doubleNe),]
#adjust for CI which expand out of the limit of the y axis
maxTMRCA<-22000
meltFstREDSavalyev$TMRCA_doubleNe_95[which(meltFstREDSavalyev$TMRCA_doubleNe_95>maxTMRCA)]<-maxTMRCA
colorino<-MainFamilies2$COLOR[which(MainFamilies2$MainFamilies==TARGETFamily)]

gg<-ggplot(meltFstREDSavalyev,aes(Savelyev2020_DivTime,TMRCA_doubleNe))

SAVAL<-gg+
  xlim(0,3000)+
  geom_errorbar(aes(ymin=TMRCA_doubleNe_5, ymax=TMRCA_doubleNe_95),size=3,width=3,
                alpha=0.1)+
  geom_point(size=3,alpha=0.7, fill=colorino, shape=21, color="black")+
  geom_text(aes(label=popslistemp), size=1)+
  xlab("Time distance from language tree - years ago")+
  ylab("Time distance from genetic data - years ago")+
  geom_abline(slope=1, intercept = 0, alpha=0.5)+
  theme_light()+
  ggtitle(TARGETFamily)+theme(plot.title = element_text(color = colorino))

####


#* # ***************************************************
### Austronesian comparisons divergence time gene language
# ***************************************************

meltFstREDinfoSINGLE$Gray2009_DivTime<-as.numeric(meltFstREDinfoSINGLE$Gray2009_DivTime)
meltFstREDGray2009<-meltFstREDinfoSINGLE[!is.na(meltFstREDinfoSINGLE$Gray2009_DivTime),]
meltFstREDBouckaert<-meltFstREDBouckaert[!is.na(meltFstREDBouckaert$TMRCA_doubleNe),] # take only the comparisons where i have a linguistic div time and a genetic div time

outliers<-c( "Mamanwa", "Rennell_and_Bellona", "Mamanwa1")
meltFstREDGray2009<-meltFstREDGray2009[-which(meltFstREDGray2009$Pop1%in%outliers),]
meltFstREDGray2009<-meltFstREDGray2009[-which(meltFstREDGray2009$Pop2%in%outliers),]

#adjust for CI which expand out of the limit of the y axis
maxTMRCA<-20000
meltFstREDGray2009$TMRCA_doubleNe_95[which(meltFstREDGray2009$TMRCA_doubleNe_95>maxTMRCA)]<-maxTMRCA

colorino<-MainFamilies2$COLOR[which(MainFamilies2$MainFamilies==meltFstREDGray2009$FAMILY[1])]

gg<-ggplot(meltFstREDGray2009,aes(Gray2009_DivTime,TMRCA_doubleNe))
AUSTR<-gg+
  ylim(0,20000)+
  xlim(0,5000)+
  geom_errorbar(aes(ymin=TMRCA_doubleNe_5, ymax=TMRCA_doubleNe_95),size=3,width=3,
                alpha=0.1)+
  geom_point(size=3,alpha=0.7, fill=colorino, shape=21, color="black")+
  # geom_text(aes(label=popslistemp), size=1)+
  xlab("Time distance from language tree - years ago")+
  ylab("Time distance from genetic data - years ago")+
  geom_abline(slope=1, intercept = 0, alpha=0.5)+theme_light()+
  ggtitle(meltFstREDGray2009$FAMILY[1])+theme(plot.title = element_text(color = colorino))

# ggsave("correlationTimeGray2009_Austronesian_noOutlierMamanwa_RennellBelloneMINI.pdf", useDingbats=FALSE, height = 5, width = 5)




#*****************************************************
## MAIN FIGURE  COMPARISON IE, AUSTR, TURKIC 
# - FIGURE 4
#*****************************************************
#*
library(ggpubr)
ggarrange(BOUCKAERT, AUSTR, HRUS  + rremove("x.text"), 
          labels = c("A", "B", "C"),
          ncol = 3, nrow = 1)
ggsave("combined3LangFamiliesCorrelation_Fig4_2021.pdf", useDingbats=FALSE, height = 4, width = 12)

library(patchwork)
BOUCKAERT+HRUS+AUSTR
ggsave("combined3LangFamiliesCorrelation_Fig4_2021.pdf", useDingbats=FALSE, height = 4, width = 12)

#*****************************************************
## SUPPLEMENTARY COMPARISON CHANG AND SAVALYEV
#*****************************************************
#*
library(ggpubr)
ggarrange(CHANG, SAVAL + rremove("x.text"), 
          labels = c("A", "B"),
          ncol = 2, nrow = 1)
ggsave("combinedChangAndSavelyev_2021.pdf", useDingbats=FALSE, height = 5, width = 10)





#**********************************************************************************
#**********************************************************************************
#**********************************************************************************
#**********************************************************************************
## supplementary Figure
## phylogenetic comparison
## ADD COLOR NODE PROPORTIONAL TO MEDIAN PROPORTION TIMES LANG AND GENES
# need phylo4 object to make a vector of the descendants tips median divergences
#**********************************************************************************
#* based on the linguistic tree, i check correspondence genetic divergence time on the nodes of the linguistic tree
#* 
library(reshape)
library(nodiv)

ThreeFamilies<-c("Indo-European", "Austronesian","Turkic")
my.max <- function(x) ifelse( !all(is.na(x)), max(x, na.rm=T), NA)

for (k in 1:3){
  TARGETFamily<-ThreeFamilies[k]
if(TARGETFamily=="Indo-European"){
listAUS<-unique(perpopRED$Bouckaert2012)
listAUS<-listAUS[-which(listAUS=="")]
AUS<-read.nexus("withinFamilyTMRCA/Bouckaert2012.trees")
taxaAUS<-read.table("withinFamilyTMRCA/Bouckaert2012_taxa.csv", sep=";", header=T, quote = "")
LISTTIME<-meltFstREDinfo$bouckaert2012_DivTime
meltFstREDinfoTEMP<-meltFstREDinfo

meltFstREDinfoTEMP$glottocodeBase11<-perpopRED$Bouckaert2012[match(meltFstREDinfoTEMP$glottocodeBase1,perpopRED$glottocodeBase)]
meltFstREDinfoTEMP$glottocodeBase22<-perpopRED$Bouckaert2012[match(meltFstREDinfoTEMP$glottocodeBase2,perpopRED$glottocodeBase)]
}
#***********************
  if(TARGETFamily=="Austronesian"){
listAUS<-unique(perpopRED$Gray2009)
listAUS<-listAUS[-which(listAUS=="")]
AUS<-read.nexus("withinFamilyTMRCA/Gray2009.tree")
taxaAUS<-read.csv("withinFamilyTMRCA/taxaGray2009.csv")
meltFstREDinfoTEMP<-meltFstREDinfo

meltFstREDinfoTEMP$glottocodeBase11<-perpopRED$Gray2009[match(meltFstREDinfoTEMP$glottocodeBase1,perpopRED$glottocodeBase)]
meltFstREDinfoTEMP$glottocodeBase22<-perpopRED$Gray2009[match(meltFstREDinfoTEMP$glottocodeBase2,perpopRED$glottocodeBase)]
}
#***********************
 if(TARGETFamily=="Turkic"){
listAUS<-unique(perpopRED$Hruschka2015)
listAUS<-listAUS[-1]
AUS<-read.nexus("withinFamilyTMRCA/Hrushka_summary2.trees")
taxaAUS<-read.csv("withinFamilyTMRCA/Hrushka_taxa2.csv", header=T, quote = "")
colnames(taxaAUS)[1]<-"code"
colnames(taxaAUS)[2]<-"taxon"
meltFstREDinfoTEMP<-meltFstREDinfo

meltFstREDinfoTEMP$glottocodeBase11<-perpopRED$Hruschka2015 [match(meltFstREDinfoTEMP$glottocodeBase1,perpopRED$glottocodeBase)]
meltFstREDinfoTEMP$glottocodeBase22<-perpopRED$Hruschka2015[match(meltFstREDinfoTEMP$glottocodeBase2,perpopRED$glottocodeBase)]
}
#***********************
#*
alberoLang<-read.tree(paste0(TARGETFamily ,"alberoLangPhylo4.phy"))
alberoLangPhylo4 <- as(alberoLang, "phylo4")
alberoLang$tip.label<-taxaAUS$glottocode[match(alberoLang$tip.label,taxaAUS$taxon)] # names in glottocode
distancealberoLang<-cophenetic.phylo(alberoLang) # names in glottocode
distancealberoLang<-distancealberoLang/2 # remember the distances are double!!


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


nodevalues<-matrix(NA, Nnode(alberoLang),7)
colnames(nodevalues)<-c("nodename","maxDivergenceTimeLang","meanDivergenceTimeLang",
                        "maxDivergenceTimeGen","meanDivergenceTimeGen", 
                        "MaxCoupleproportionDivergenceTime","meanproportionDivergenceTime")
nodevalues[,1]<-c(1:Nnode(alberoLang))

for (j in 1:Nnode(alberoLang)){
  settemp<- labels(descendants(alberoLangPhylo4,which(attributes(alberoLangPhylo4)$label==j),type = "all"))
  settemp<- taxaAUS$glottocode[which (taxaAUS$taxon%in%settemp)] # glottocodes
  timetemp<-(distancealberoLang[which(colnames(distancealberoLang)%in%settemp),which(rownames(distancealberoLang)%in%settemp)])
  
  timetemp<-melt(distancealberoLang[which(colnames(distancealberoLang)%in%settemp),which(rownames(distancealberoLang)%in%settemp)])
  maxdivergence<-timetemp[which(timetemp$value==max(timetemp$value)),] # i cannot use only the maximum divergence time otherwise i do not have useful matches with genetic data, so i pick up the maximum from genetic divergence of all the derived nodes
  nodevalues[j,2]<-max(timetemp$value)
  nodevalues[j,3]<-mean(timetemp$value[timetemp$value>0]) # exclude zero which is the same population match
  
  maxcouplestimenode<-c(as.character(maxdivergence$X1),as.character(maxdivergence$X2))
  gentimetemp<-distanceTMRCAred[which(colnames(distanceTMRCAred)%in%settemp),which(rownames(distanceTMRCAred)%in%settemp)]
  gentimetempProportion<- distancePROPORTIONred[which(colnames(distancePROPORTIONred)%in%settemp),which(rownames(distancePROPORTIONred)%in%settemp)]
  nodevalues[j,4]<-my.max(gentimetemp)
  nodevalues[j,5]<-mean(gentimetemp, na.rm = T)
  nodevalues[j,6]<-mean(melt(gentimetempProportion)[which(timetemp$value==max(timetemp$value)),]$value, na.rm = T)
  nodevalues[j,7]<-mean(gentimetempProportion, na.rm = T)
}

nodevalues<-as.data.frame(nodevalues)
nodevalues$mainFamily<-TARGETFamily
#nodevalues11<-rbind(nodevalues11,nodevalues)
# from here
#nodevalues<-nodevalues11[which(nodevalues11$mainFamily==TARGETFamily),]
listAUSreverseInGenTMRCA<-taxaAUS$taxon[match(glottredlist,taxaAUS$glottocode)]

alberoLangNamesLang<-alberoLang
alberoLangNamesLang$tip.label<-taxaAUS$taxon [match(alberoLangNamesLang$tip.label,taxaAUS$glottocode)]

alberoLangNamesLang$tip.label[which(alberoLangNamesLang$tip.label%in%listAUSreverseInGenTMRCA)]<-paste0(alberoLangNamesLang$tip.label[which(alberoLangNamesLang$tip.label%in%listAUSreverseInGenTMRCA)], "_GEN_TIME")

proportionGenLang <- nodevalues$meanproportionDivergenceTime
proportionGenLang[!is.na(nodevalues$MaxCoupleproportionDivergenceTime)] <- nodevalues$MaxCoupleproportionDivergenceTime[!is.na(nodevalues$MaxCoupleproportionDivergenceTime)] #
# replace with max divergence time when available 

cexplay=1.5

pdf(paste0(TARGETFamily,"_2021.pdf"),width=12, height=5,useDingbats=FALSE)
par(mfrow=c(1,4))
plot_nodes_phylo(round(nodevalues$maxDivergenceTimeLang), alberoLangNamesLang, cex = cexplay, main= "maximum language divergence time")
plot_nodes_phylo(round(nodevalues$maxDivergenceTimeGen), alberoLangNamesLang, cex = cexplay, main= "maximum genetic divergence time")
plot_nodes_phylo(round(nodevalues$meanDivergenceTimeGen), alberoLangNamesLang, cex = cexplay, main= "mean genetic divergence time")
plot_nodes_phylo(proportionGenLang, alberoLangNamesLang, cex = cexplay, main= "Mean Lang/Gen proportion")
dev.off()
}

#write.table(nodevalues11, "breakdowninternalNodeTreesFamilies.txt", sep="\t", row.names = F, quote=F) 
#nodevalues11<-read.table("breakdowninternalNodeTreesFamilies.txt", sep="\t", header=T)


#**********************************************************************************
#**********************************************************************************
#**********************************************************
### FST tree of the selected languages
# compare phylogeny fst and phylogeny language time tree
#**********************************************************

library('Quartet')
library("phytools")
library("phylogram")


meltFstREDinfo$threefamilies<-NA
for (i in 1:nrow(meltFstREDinfo)){
  lista<-c(meltFstREDinfo$Gray2009_DivTime[i],meltFstREDinfo$bouckaert2012_DivTime[i],meltFstREDinfo$Hruschka2015_DivTime[i])
pick<-which(lista!="NA")
if(length(pick)>0){
meltFstREDinfo$threefamilies[i]<-lista[pick]
}
}

ThreeFamilies<-c("Indo-European", "Austronesian","Turkic")

for (k in 1:3){
  TARGETFamily<-ThreeFamilies[k]
  if(TARGETFamily=="Indo-European"){
    listAUS<-unique(perpopRED$Bouckaert2012)
    listAUS<-listAUS[-which(listAUS=="")]
    AUS<-read.nexus("withinFamilyTMRCA/Bouckaert2012.trees")
    taxaAUS<-read.table("withinFamilyTMRCA/Bouckaert2012_taxa.csv", sep=";", header=T, quote = "")
    alberoLang<-read.tree(paste0(TARGETFamily ,"alberoLangPhylo4.phy"))
    alberoLangPhylo4 <- as(alberoLang, "phylo4")
    alberoLang$tip.label<-taxaAUS$glottocode[match(alberoLang$tip.label,taxaAUS$taxon)] # names in glottocode
    listpopLang<-perpopRED$PopName[which(perpopRED$Bouckaert2012%in% alberoLang$tip.label)]
    ## root with Armenian as outgroup, according to the linguistic tree IE
    rootalo<-c("Armenian_Hemsheni","Armenian")
    
     }
  #***********************
  if(TARGETFamily=="Austronesian"){
    listAUS<-unique(perpopRED$Gray2009)
    listAUS<-listAUS[-which(listAUS=="")]
    AUS<-read.nexus("withinFamilyTMRCA/Gray2009.tree")
    taxaAUS<-read.csv("withinFamilyTMRCA/taxaGray2009.csv")
    alberoLang<-read.tree(paste0(TARGETFamily ,"alberoLangPhylo4.phy"))
    alberoLangPhylo4 <- as(alberoLang, "phylo4")
    alberoLang$tip.label<-taxaAUS$glottocode[match(alberoLang$tip.label,taxaAUS$taxon)] # names in glottocode
    listpopLang<-perpopRED$PopName[which(perpopRED$Gray2009%in% alberoLang$tip.label)]
    ## root with Ami as outgroup, according to the linguistic tree AUSTRONESIAN
    rootalo<-c("Ami")
    
  }
  #***********************
  if(TARGETFamily=="Turkic"){
    listAUS<-unique(perpopRED$Hruschka2015)
    listAUS<-listAUS[-1]
    AUS<-read.nexus("withinFamilyTMRCA/Hrushka_summary2.trees")
    taxaAUS<-read.csv("withinFamilyTMRCA/Hrushka_taxa2.csv", header=T, quote = "")
    colnames(taxaAUS)[1]<-"code"
    colnames(taxaAUS)[2]<-"taxon"
    alberoLang<-read.tree(paste0(TARGETFamily ,"alberoLangPhylo4.phy"))
    alberoLangPhylo4 <- as(alberoLang, "phylo4")
    alberoLang$tip.label<-taxaAUS$glottocode[match(alberoLang$tip.label,taxaAUS$taxon)] # names in glottocode
    listpopLang<-perpopRED$PopName[which(perpopRED$Hruschka2015%in% alberoLang$tip.label)]
    ## root with two Chuvash speaking as outgroup, according to the linguistic tree Turkic
    rootalo<-c("Chuvash","Chuvash_Tatarstan")
    
  }
  #***********************
  #*
# make the fst matrix with the populations present in the language tree

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
       timetreelangMatrix[i,j]<-as.numeric(meltFstREDinfo$threefamilies[temp])
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

write.tree(phy1,paste0(TARGETFamily, "_FST.tree"))
write.tree(phy2,paste0(TARGETFamily, "_lang.tree"))

pdf(paste0(TARGETFamily,"treesComparison_2021.pdf"),width=12, height=12,useDingbats=FALSE)

#par(mfrow=c(1,3))
layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))
# cophilo comparison
plot(cophylo(phy2,phy1,rotate=T), fsize=0.6)

# lang tree
#pdf(paste0(TARGETFamily,"Langtree_2021.pdf"),width=5, height=7,useDingbats=FALSE)
plot.phylo(phy2,align.tip.label = TRUE, cex=0.5)
axisPhylo(side = 3)
#dev.off()

# FST tree
#pdf(paste0(TARGETFamily,"FSTtree_2021.pdf"),width=5, height=7,useDingbats=FALSE)
plot.phylo(phy1,align.tip.label = TRUE,cex=0.4)
axisPhylo(side = 3)
#dev.off()

dev.off()

### QUARTET measurements
statuses <- QuartetStatus(phy1, phy2)
QuartetDivergence(statuses, similarity = FALSE)
SimilarityMetrics(statuses, similarity = TRUE)
ncol(AllQuartets(Ntip(phy2)))

pdf(paste0(TARGETFamily,"Quartet.pdf"),useDingbats=FALSE, height = 10, width = 10)
VisualizeQuartets(phy2, phy1, scale=0.6)
dev.off() # i cannot group the quartet plots in a single figure!

### compare phylogenies with Phytools

#  TARGETFamily<-ThreeFamilies[k]
#  phy1<-read.tree(paste0(TARGETFamily, "_FST.tree"))
#  phy2<-read.tree(paste0(TARGETFamily, "_lang.tree"))
  
  # plot compared phylogenies
  pdf(paste0(TARGETFamily,"Cophilo.pdf"),width=15, height=7,useDingbats=FALSE)
  
  plot(cophylo(phy2,phy1,rotate=T), fsize=0.6)
  # nodelabels.cophylo(phy2$node.label)
  #  nodelabels.cophylo(which="right",phy1$node.label )
  
  dev.off()

}



