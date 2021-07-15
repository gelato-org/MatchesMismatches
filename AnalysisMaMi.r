#******************************************
#******************************************
### analysis Paper MaMi
# Chiara BArbieri 30 March 2021

#******************************************
#*
#*
#*

#write.table(FstListREDinfo, "/Users/chiarabarbieri/Library/Mobile Documents/com~apple~CloudDocs/GeLaTo/FstListREDinfo2021.txt", sep="\t", row.names = F, quote=F) 

#write.table(perpopRED, "/Users/chiarabarbieri/Library/Mobile Documents/com~apple~CloudDocs/GeLaTo/perpopREDMaMi2021.txt", row.names = F,  sep = "\t", quote = F)


#source("/Users/chiarabarbieri/Library/Mobile Documents/com~apple~CloudDocs/ggworld.R")

# setwd("/Users/chiarabarbieri/switchdrive/GeLaTo")

### read the two main files
# list of 404 populations :
perpopRED<-read.table("/Users/chiarabarbieri/Library/Mobile Documents/com~apple~CloudDocs/GeLaTo/perpopREDMaMi2021.txt",   sep = "\t", header=T, as.is=T)
# list of pairwise comparisons :
FstListREDinfo<-read.table("/Users/chiarabarbieri/Library/Mobile Documents/com~apple~CloudDocs/GeLaTo/FstListREDinfo2021.txt",header = T, as.is=T, sep="\t")

# color palette for major language families

Lfamil<-table(perpopRED$glottolog.NAME)
MainFamilies<-unlist(labels(Lfamil[which(Lfamil>4)])) # minimum 5 populations per Lang Family
#MainFamilies<-MainFamilies[-which(MainFamilies=="LI")]   # exclude the Language Isolates LI
perpopREDfamily<-perpopRED[which(perpopRED$glottolog.NAME %in% MainFamilies),]
colorchoice<-c( "darkorange4" ,"#A6CEE3"  ,   "#33A02C"    , "#A6D854"  ,   "#377EB8"  ,   "#E31A1C"   ,  "#FFD92F" ,    "#FF7F00"  ,   "#666666" ,   
                "cyan4"  ,     "#BC80BD"   ,  "#FED9A6" ,    "tan3" ,       "#6A3D9A" ,    "deeppink"   )
MainFamilies2<-as.data.frame(MainFamilies)
MainFamilies2$COLOR<-colorchoice


library(ggplot2)


#******************************************
### MAP
#******************************************

library(maps)
library('maps')
library('geosphere')
library(rworldmap)


source("/Users/chiarabarbieri/Library/Mobile Documents/com~apple~CloudDocs/ggworld.R")
perpopRED$lat<-as.numeric(perpopRED$lat)
perpopRED$lon<-as.numeric(perpopRED$lon)


perpopRED2<-perpopRED[!is.na(perpopRED$lat),] # exclude populations for which i do not have geographic coordinates
perpopREDSHIFTMAP<-MoveDataCoord(perpopRED2)   # perform coordinate shift to plot Pacific centered maps

## map population location and language families
#### MAP WITH MAJOR LANGUAGE FAMILIES ASSIGNED TO COLOR CODE,  
# FIGURE 1a
#******************************************
#*

base_plot + geom_point(data = perpopREDSHIFTMAP,
                       aes(x = lon.1, y = lat.1), 
                       color="black",shape=3,size=0.5)+
  theme(legend.position="bottom")+
  geom_point(data = perpopREDSHIFTMAP[which(perpopREDSHIFTMAP$MainFamilies%in%MainFamilies),],
             aes(x = lon.1, y = lat.1, color=MainFamilies), 
             size = 2, alpha=0.5)+
  scale_color_manual(values=colorchoice)

ggsave("WholeGelatoMap_Pacific_LittleCrossesMajorFamilies.pdf", useDingbats=FALSE, width = 15, height = 8, dpi = 300, units = "in")







#******************************************
  #*#******************************************
  #*   genetic and linguistic enclaves
  #*#******************************************
  
  # ***************************************************
  # evaluate the incidence of neighbors from the SAME AND THE different L famiy for each population
  # closer neighbor from same and different L family
  
  perpopRED$hasOtherMatches<-"YES"  
perpopRED$closeFstSameFamily<-NA  
perpopRED$geodistSameFamily<-NA  
perpopRED$closeFstDIFFFamily<-NA  
perpopRED$geodistDIFFFamily<-NA  
perpopRED$DIFFFamilyClosestpop<-NA
perpopRED$SameFamilyMostDistantClosestFst<-NA


# how far the same family is gen closer than other families?


for (i in 1:nrow(perpopRED)){
  targetpop<-perpopRED$PopName[i]
  temp<-FstListREDinfo[which(FstListREDinfo$Pop1==targetpop),]
  
  samefamily<-temp[which(temp$FAMILY!="DIVERSE"),]
  escludiniVicini<-which(samefamily$GEOdist<10&samefamily$glottocodeBase1==samefamily$glottocodeBase2)
  if(length(escludiniVicini)!=0){
    samefamily<-samefamily[-escludiniVicini,] # exclude when there is a neighbor too close (LESS THAN 10 KM) from same exact language as DUPLICATED SAMPLE
  } 
  if(nrow(samefamily)==0){
    perpopRED$hasOtherMatches[i]<-"NO" 
  }
  else{
    perpopRED$closeFstSameFamily[i]<-sort(samefamily$Fst)[1]  # the closest Fst
    perpopRED$geodistSameFamily[i]<-samefamily$GEOdist[order(samefamily$Fst)][1]    # the geographic distance from the closest Fst
  }
  DIFFfamily<-temp[which(temp$FAMILY=="DIVERSE"),]
  perpopRED$closeFstDIFFFamily[i]<-sort(DIFFfamily$Fst)[1]  # the closest Fst from a different language family
  perpopRED$geodistDIFFFamily[i]<-DIFFfamily$GEOdist[order(DIFFfamily$Fst)][1]    # the geographic distance from the closest Fst of a different language family
  perpopRED$DIFFFamilyClosestpop[i]<-DIFFfamily[order(DIFFfamily$Fst),][1,]$Pop2  # the pop which closest fst in mismatch
  
  if(length(which(sort(samefamily$Fst)<perpopRED$closeFstDIFFFamily[i]))==0){
    perpopRED$SameFamilyMostDistantClosestFst[i]<-"NONE"    }
  else{
    perpopRED$SameFamilyMostDistantClosestFst[i]<-tail(samefamily$GEOdist[which(sort(samefamily$Fst)<perpopRED$closeFstDIFFFamily[i])],1) # the geographic distance from the closest Fst of the same language family that is less than the closest fst from diff L family
  }
}


perpopRED$proportionFST_diff_sameFamily<- perpopRED$closeFstDIFFFamily/perpopRED$closeFstSameFamily 
perpopRED$proportionGeoDistSameDiffFamily<- perpopRED$geodistDIFFFamily/perpopRED$geodistSameFamily 
perpopRED$Mismatch3a<-NA
#perpopRED$Mismatchsingular2[which(perpopRED$closeFstSameFamil==0&perpopRED$geodistSameFamily>100&perpopRED$closeFstDIFFFamily!=0)] <- "secondaryMATCH_FSTzero_geodist100"
perpopRED$Mismatch3a[which(perpopRED$closeFstDIFFFamily==0&perpopRED$closeFstSameFamil!=0)] <- "secondaryMISMATCH_FSTzeroDiffFamily"

perpopRED$Mismatch3a[which(perpopRED$proportionFST_diff_sameFamily<1&perpopRED$proportionGeoDistSameDiffFamily>1)] <- "MISMATCH"
perpopRED$Mismatch3a[which(perpopRED$proportionFST_diff_sameFamily>1&perpopRED$proportionGeoDistSameDiffFamily<1)] <- "MATCH"

perpopRED$Mismatch3a[which(perpopRED$hasOtherMatches!="YES"  )]<-"ZeroSameFamilyNeighbors"

table(perpopRED$Mismatch3a)/nrow(perpopRED)

### table to visualize how many cases

#*** 2021 

MATCH                            MISMATCH 
0.128712871                         0.069306931 
secondaryMISMATCH_FSTzeroDiffFamily             ZeroSameFamilyNeighbors 
0.004950495                         0.049504950 
#***


## genetic and linguistic enclaves to exclude from further analysis

ListEnclaves<-perpopRED$PopName [which(perpopRED$Mismatch3a%in%c("MISMATCH","secondaryMISMATCH_FSTzeroDiffFamily"))]

> ListEnclaves
[1] "Yoruba"            "Mengen"            "Bengali"           "Hazara"           
[5] "Kharia"            "Gui"               "Khwe"              "Nama"             
[9] "Mongola"           "Avar_Gunibsky"     "Aleut"             "Dai"              
[13] "Jew_Georgian"      "Khomani"           "Spanish_PaisVasco" "Yaquis"           
[17] "Yukagir_Forest"    "Yukagir_Tundra"    "Zapotec"           "Wayku"            
[21] "Han-NChina"        "Evenk_FarEast"     "Cocama"            "Guarani"          
[25] "Guarani_GN"        "Karitiana"         "Surui"             "Azeri_Azerbajan"  
[29] "Hungarian1"        "Hungarian2"       

enclavesByMistake<-c("Yoruba"  ,"Nama", "Han-NChina" , "Bengali", "Evenk_FarEast") # these pops are not real enclaves, their linguistically unrelated pair is the one driving this genetic proximity effect
ListEnclaves<-ListEnclaves[-which(ListEnclaves%in%(enclavesByMistake))]

perpopRED$ListEnclaves<-NA
perpopRED$ListEnclaves[which(perpopRED$PopName%in%ListEnclaves)]<-"Enclave"


#*#******************************************
#*
#* potential FIGURE 1 different distribution FST closest same family and different family
#* 

aa<-perpopRED[,c(1,39,41,16)] 
meltperpop<-melt(aa)

aa$difference<-aa$closeFstSameFamily-aa$closeFstDIFFFamily
length(which(aa$difference>0))/(404-length(which(is.na(aa$difference))))
[1] 0.1822917

# 18 % of the pops who have another genetic population of the same language family do have closer FST with speakers of another language family


#**************************************************
### table of cases mismatch enclaves

mismatch3a<-perpopRED[grep("MISMATCH",perpopRED$Mismatch3a),]
perpopRED$DIFFFamilyClosestFAM<-perpopRED$glottolog.NAME[match(perpopRED$DIFFFamilyClosestpop,perpopRED$PopName)]
colnameschoice<-c("PopName","glottolog.NAME", "geodistSameFamily","geodistDIFFFamily","DIFFFamilyClosestpop", "DIFFFamilyClosestFAM")
mismatch3a<-mismatch3a[,colnameschoice]
#write.table(mismatch3a[order(mismatch3a$geodistDIFFFamily),],"tableSupplMismatches3a.txt", sep="\t", row.names = F, quote=F) 




### 
#*#******************************************
#*  #*#******************************************
#*    #*#******************************************
#*    
#*    


## exclude neighbors who speak the same language as DUPLICATE POPS
escludiniVicini<-FstListREDinfo[which(FstListREDinfo$GEOdist<10&FstListREDinfo$glottocodeBase1==FstListREDinfo$glottocodeBase2),]

Pop1             Pop2                       popslistemp         Fst       family1       family2        region2
7803       Ju_hoan_North              San                  Ju_hoan_NorthSan 0.003176030           Kxa           Kxa         AFRICA
9447                Kove    Kove_Tamuniai                 KoveKove_Tamuniai 0.000000000  Austronesian  Austronesian        OCEANIA
10395            Mamanwa         Mamanwa1                   MamanwaMamanwa1 0.050776000  Austronesian  Austronesian SOUTHEAST_ASIA
13579        Sulka_Ganai     Sulka_Watwat           Sulka_GanaiSulka_Watwat 0.000000000            LI            LI        OCEANIA
25464           BedouinA         BedouinB                  BedouinABedouinB 0.015446000  Afro-Asiatic  Afro-Asiatic        EURASIA
34834             Cusco2     Peru_Quechua                Cusco2Peru_Quechua 0.000427919      Quechuan      Quechuan       AMERICAS
36912              Druze      Palestinian                  DruzePalestinian 0.008987530  Afro-Asiatic  Afro-Asiatic        EURASIA
45009       Greek_Athens        Greek_WGA             Greek_AthensGreek_WGA 0.003147220 Indo-European Indo-European        EURASIA
46714          GujaratiA        GujaratiB                GujaratiAGujaratiB 0.001633070 Indo-European Indo-European        EURASIA
46731          GujaratiA        GujaratiC                GujaratiAGujaratiC 0.002983290 Indo-European Indo-European        EURASIA
46899          GujaratiB        GujaratiC                GujaratiBGujaratiC 0.003168390 Indo-European Indo-European        EURASIA
48305         Han-NChina         Han_HGDP                Han-NChinaHan_HGDP 0.002480980  Sino-Tibetan  Sino-Tibetan        EURASIA
49975         Hungarian1       Hungarian2              Hungarian1Hungarian2 0.000000000        Uralic        Uralic        EURASIA
60866               Kinh Vietnamese_North              KinhVietnamese_North 0.001434990 Austroasiatic Austroasiatic SOUTHEAST_ASIA
63575 Lebanese_Christian  Lebanese_Muslim Lebanese_ChristianLebanese_Muslim 0.001025050  Afro-Asiatic  Afro-Asiatic        EURASIA
81071              Uzbek   Uzbek_Tashkent               UzbekUzbek_Tashkent 0.000917217        Turkic        Turkic        EURASIA


FstListREDinfo_noDuplicateNeighbors<-FstListREDinfo[-which(FstListREDinfo$GEOdist<10&FstListREDinfo$glottocodeBase1==FstListREDinfo$glottocodeBase2),]



#******************************************
#*#******************************************
#*#******************************************
#*#******************************************
#*#******************************************
#*#******************************************
#*    COMPARING MEDIAN FST BETWEEN AND WITHIN
# POPULATION HEURISTIC 2
#* FST distribution between and within language families
#*#******************************************
#*#******************************************

library(pairwiseCI)
FstListREDinfo_noDuplicateNeighbors<-FstListREDinfo[-which(FstListREDinfo$GEOdist<10&FstListREDinfo$glottocodeBase1==FstListREDinfo$glottocodeBase2),]

perpopRED$MedianWithinSMALLER<-NA
  
for (i in 1:nrow(perpopRED)){
  TARGET<-perpopRED$PopName[i]
  meltini<-FstListREDinfo_noDuplicateNeighbors[which(FstListREDinfo_noDuplicateNeighbors$Pop1==TARGET),]
  meltiniSAME<-meltini[which(meltini$SameFamily=="YES"),]
  maxgeofam<-max(meltiniSAME$GEOdist, na.rm = T)
  if(maxgeofam<500){
    maxgeofam<-500
  }
  meltini<-meltini[which(meltini$GEOdist<=maxgeofam),]
  meltiniNO<-meltini[which(meltini$SameFamily!="YES"),]
  
  if(length(which(meltini$SameFamily=="YES"))>2&
     length(which(meltini$SameFamily=="NO"))>2){   # need to have at least 3 comparisons between and within family
    
    perpopRED$MedDiff[i]=Median.diff(meltiniNO$FstLinear, meltiniSAME$FstLinear, conf.level=0.95, alternative="lesser",R=10000)$estimate
    perpopRED$MedDiffCIlower[i]=Median.diff(meltiniNO$FstLinear, meltiniSAME$FstLinear, conf.level=0.95, alternative="lesser",R=10000)$conf.int[1]
    perpopRED$MedDiffCIupper[i]=Median.diff(meltiniNO$FstLinear, meltiniSAME$FstLinear, conf.level=0.95, alternative="lesser",R=10000)$conf.int[2]
     perpopRED$MedianWithinSMALLER[i]<-median(meltiniSAME$FstLinear)<median(meltiniNO$FstLinear) # if TRUE, the median FST within is smaller than the FST between
  }
}

perpopRED$MedDiffCI<-perpopRED$MedDiffCIupper-perpopRED$MedDiffCIlower

########################################################
# Figure S8 supplementary - comparison of difference of medians and CI, colored per language family
########################################################

library(ggrepel)
 colorchoice2<-colorchoice
 colorchoice2[10]<-"gray20"
 colorchoice2[11:16]<-colorchoice[10:15]

 ggplot(perpopRED,aes(MedDiff,MedDiffCI, color=MainFamilies))+
  geom_point()+
  geom_text_repel(aes(label=PopName, color=MainFamilies), size=0.5)+
scale_color_manual(values = colorchoice2)+
  theme_light()+ xlab("Difference Median FST between - within") +
  ylab("Confidence Interval median difference") 
ggsave("DifferenceMedianwithText.pdf", height = 12, width = 14, useDingbats=FALSE)


fst1<-ggplot(perpopRED,aes(MedDiff,MedDiffCI, color=MainFamilies))+
  geom_point(alpha=0.5, size=3,show.legend = FALSE)+
  #geom_text_repel(aes(label=PopName, color=MainFamilies), size=0.5)+
  scale_color_manual(values = colorchoice2)+
  theme_light() +
xlab("Difference Median FST between - within") +
  ylab("Confidence Interval median difference") +
  geom_vline(xintercept=0)+
  scale_y_log10()

ggsave("DifferenceMedianLogScale.pdf", height = 12, width = 14, useDingbats=FALSE)


fst2<-ggplot(perpopREDfamily ,aes(MedDiff,MedDiffCI, color=MainFamilies))+
  geom_point(alpha=0.5, size=3,show.legend = FALSE)+
  geom_text_repel(aes(label=PopName, color=MainFamilies), size=0.5, show.legend = FALSE)+
  scale_color_manual(values = colorchoice)+
  theme_light()+ xlab("Difference Median FST between - within") +
  ylab("Confidence Interval median difference") +
  geom_vline(xintercept=0)+
  scale_y_log10()+
  facet_wrap(~MainFamilies,ncol=3)
ggsave("DifferenceMedianLogScaleMajorLangFam.pdf", height = 12, width = 14, useDingbats=FALSE)

library(patchwork)
(fst1 + fst2)
ggsave("DifferenceMedianLogScaleGLOBALandMajorLangFam.pdf", height = 6, width = 12, useDingbats=FALSE)



### 
# mark the mismatches in fst distribution
ListEnclaves<-perpopRED$PopName[which(perpopRED$ListEnclaves=="Enclave")] 

Misalligned<-perpopRED$PopName[which(perpopRED$MedDiff<0&perpopRED$MedDiffCI<0.01)]
Misalligned<-c(Misalligned,"Hoan", "Naro")
Misalligned<-Misalligned[-which(Misalligned%in%ListEnclaves)]


## all the problematic cases together
library(reshape) 

ListEnclaves<-perpopRED$PopName[which(perpopRED$ListEnclaves=="Enclave")] 
intersecti<-intersect(DRIFTONI, ListEnclaves)

perpopRED$singlepops<-"ND"
perpopRED$singlepops[which(perpopRED$Mismatch3a =="MATCH")]<-"Match"
perpopRED$singlepops[which(perpopRED$PopName %in% DRIFTONI)]<-"Drifted"
perpopRED$singlepops[which(perpopRED$ListEnclaves=="Enclave")]<-"Enclave"
perpopRED$singlepops[which(perpopRED$PopName %in% Misalligned)]<-"Misalligned"


table(perpopRED$singlepops)

#**************************************************
### stacked bar proportion cases to exclude
# supplementary S12
#**************************************************
library(reshape) 
perpopREDfamily<-perpopRED[which(perpopRED$glottolog.NAME %in% MainFamilies),]

size<-as.numeric(table(perpopREDfamily$MainFamilies))
proportionSingleCases<-table(perpopREDfamily$MainFamilies,perpopREDfamily$singlepops)/size

proportionSingleMELT<-melt(proportionSingleCases)
colnames(proportionSingleMELT)<-c("Family", "variable", "proportion")

gi<-ggplot(proportionSingleMELT, aes(Family, proportion))

ggplot(proportionSingleMELT, aes(fill=variable, y=proportion, x=Family)) + 
  geom_bar(position="fill", stat="identity", alpha=0.8)+
  scale_fill_manual(values = c("cyan", "red", "orange", "yellow","gray20"), name="")+
theme(axis.text.x = element_text(angle = 45, vjust=1, hjust=1)) +xlab("") +ylab("")

  ggsave("proportionSingleCasesEachFamily.pdf", useDingbats=FALSE, width = 7, height = 7)

########
## the map with singular mismatches of linguistic and genetic migrants, misaligned and drifted
# Figure 1 C
########
  
library(ggrepel)

#enclavesByMistake<-c("Yoruba"  ,"Nama", "Han-NChina" , "Bengali", "Evenk_FarEast") # these pops are not real enclaves, their linguistically unrelated pair is the one driving this genetic proximity effect

perpopRED2<-perpopRED[!is.na(perpopRED$lat),]
perpopREDSHIFTMAP<-MoveDataCoord(perpopRED2)
perpopREDSHIFTMAPinterest<-perpopREDSHIFTMAP[-which(perpopREDSHIFTMAP$singlepops=="ND"),]
perpopREDSHIFTMAPmismatch<-perpopREDSHIFTMAPinterest[which(perpopREDSHIFTMAPinterest$singlepops%in%c("Misalligned","Enclave")),]
perpopREDSHIFTMAPnodata<-perpopREDSHIFTMAP[which(perpopREDSHIFTMAP$Mismatch3a=="ZeroSameFamilyNeighbors"),]
#perpopREDSHIFTMAPinterest<-perpopREDSHIFTMAPinterest[-which(perpopREDSHIFTMAPinterest$PopName%in% enclavesByMistake),]
#perpopREDSHIFTMAPmismatch<-perpopREDSHIFTMAPinterest[which(perpopREDSHIFTMAPinterest$Mismatch3a%in%c("MISMATCH","secondaryMISMATCH_FSTzeroDiffFamily")),]

#perpopREDSHIFTMAPNoMismatch<-perpopREDSHIFTMAP[which(perpopREDSHIFTMAP$hasOtherMatches=="NO"),]

base_plot + geom_point(data = perpopREDSHIFTMAP,
                       aes(x = lon.1, y = lat.1),
                       color = "black",
                       size = 0.5, shape=3)+
  geom_point(data=perpopREDSHIFTMAPinterest, aes(x = lon.1, y = lat.1, 
                                                 color=singlepops),alpha=0.8)+
  geom_point(data=perpopREDSHIFTMAPnodata, aes(x = lon.1, y = lat.1),
                                                 color="white", shape= 13,alpha=0.8)+
  geom_text_repel(data=perpopREDSHIFTMAPmismatch, aes(x = lon.1, y = lat.1, 
                                                      label=PopName), size=2) +
  xlab("") + ylab("") +
  scale_colour_manual(values = c("pink", "darkred", "blue", "orange"), name="")

                   #   labels=c("Match", "Mismatch Genetic Enclave","Mismatch Linguistic Enclave", "no neighbors from the same Language Family"))

ggsave("MapMismatchpopulationsAndDrifted_.pdf", useDingbats=FALSE, width = 15, height = 8, dpi = 300, units = "in")






# ***************************************************
#### SUPPLEMENTARY FIGURES 
# ***************************************************
# ***************************************************
# ***************************************************




# ***************************************************
# sensitivity test with different geo thresholds


### all the close FST distances on a map
base_plot + geom_point(data = perpopREDSHIFTMAP,
                       aes(x = lon.1, y = lat.1),
                       color = "black",
                       size = 0.5, shape=3)+
  geom_path(data=betweenfamilSMALLgeoplottabileSHIFTMAP, aes(x = lon.1, y = lat.1, 
                                                             group=index, alpha=percentileFST, color=percentileFST))+
  scale_alpha(range = c(0.01,0.2))+
  scale_colour_gradient2(midpoint=0.955, low="white", mid="yellow",
                         high="purple")+
  xlab("") + ylab("") 

ggsave("MapPairFSTPercentileDistribution.pdf", useDingbats=FALSE, width = 15, height = 8, dpi = 300, units = "in")


threshold<-c(500,1000,1500,2000,2500,3000)
blockneighborthresholdGeo<-c()
for (i in 1:length(threshold)){
  testino2<-ddply(FstListREDinfo, "Pop1", function(x) length(which(x$GEOdist<threshold[i]&x$FAMILY=="DIVERSE"))) # no continental filter

  colnames(testino2)<-c("PopName","GeoCloseUnrelated1")
  testinomerged<-testino2
  testinomerged$GeoCloseUnrelated<-testinomerged$GeoCloseUnrelated1/2
  testinomerged$threshold<-threshold[i]
  
  blockneighborthresholdGeo<-rbind(blockneighborthresholdGeo,testinomerged)
}


blockneighborthresholdGeoINFO<-merge(perpopRED,blockneighborthresholdGeo)


# violin plots with number of populations from a different language family
# fig S3B


blockneighborthresholdGeoINFO$threshold<-as.character(blockneighborthresholdGeoINFO$threshold)
blockneighborthresholdGeoINFO$threshold[which(blockneighborthresholdGeoINFO$threshold=="500")]<-" 500"
ga<-ggplot(blockneighborthresholdGeoINFO,aes(x=threshold,y=GeoCloseUnrelated1))
ga+ geom_jitter(shape=16, position=position_jitter(0.2))+
  geom_violin(trim=FALSE, alpha=0.4)+
  stat_summary(fun.data="mean_sdl", mult=1, 
               geom="pointrange", color="orange")+
  theme_light()+
  ylim(0,150)+xlab("Geographic radius in km")+ylab("number of populations from a different language family")

ggsave("Suppl_DensityDistributionMismatches_vsGeoDistanceThresholds_2021.pdf", useDingbats=FALSE, height = 4, width = 6, units = "in")


# frequencies at least one neighbor from diff family
freqOneneighbor<-c()
for (i in 1:length(threshold)){
  thresholdtemp<-threshold[i]
popsNumberofNeighbors<-ddply(FstListREDinfo, "Pop1", function(x) length(which(x$GEOdist<thresholdtemp&x$FAMILY=="DIVERSE"))) # no continental filter
freqOneneighbor[i]<-length(which(popsNumberofNeighbors$V1>0))/404
}
round(freqOneneighbor*100)
[1]  57  85  93  98  99 100  # add manually on the violin plot figure



# ***************************************************
# evaluate the incidence of neighbors from a different L famiy for each population

# count how many populations from different L families in a radius of 1000 km in the same continent

popsNumberofNeighbors<-ddply(FstListREDinfo, "Pop1", function(x) length(which(x$GEOdist<1000&x$FAMILY=="DIVERSE"))) # no continental filter
perpopRED$nNeighborsDiffFamily<-popsNumberofNeighbors
colnames(popsNumberofNeighbors)<-c("PopName","GeoCloseUnrelated")
perpopRED3<-merge(perpopRED,popsNumberofNeighbors)
perpopRED4<-perpopRED3[-which(perpopRED3$GeoCloseUnrelated==0),]
perpopRED4_SHIFTMAP<-MoveDataCoord(perpopRED4)

# plot on a map the density of unrelated pop geographically close
# fig S3A


base_plot + geom_point(data = perpopREDSHIFTMAP,
                       aes(x = lon.1, y = lat.1),
                       color = "black",
                       size = 0.5, shape=3)+
  geom_point(data=perpopRED4_SHIFTMAP,  aes(x = lon.1, y = lat.1,color=as.numeric(GeoCloseUnrelated)), alpha=0.5, size=2)+
  # ggtitle("number of neighbors from a different language family")+
  scale_color_gradient(low = "blue", high = "red") +
  xlab("") + ylab("") +
  theme(legend.title = element_blank())
ggsave("DensityPopulationsWithDifferentLFamilyNeighbors_noContinentFilter_map_PacificCenter_1000km_2021_.pdf", useDingbats=FALSE, width = 15, height = 8, dpi = 300, units = "in")

# count how many LANGUAGE FAMILIES  in a radius of 1000 km in the same continent

popsNumberofNeighborFamilies<-ddply(FstListREDinfo, "Pop1", function(x) length(unique(x[which(x$GEOdist<1000&x$FAMILY=="DIVERSE"),]$family2))) # no continental filter
perpopRED$nDiffFamilies1000km<-popsNumberofNeighbors

colnames(popsNumberofNeighborFamilies)<-c("PopName","GeoCloseUnrelatedFamilies")
perpopRED3<-merge(perpopRED,popsNumberofNeighborFamilies)
perpopRED4<-perpopRED3[-which(perpopRED3$GeoCloseUnrelatedFamilies==0),]
perpopRED4_SHIFTMAP<-MoveDataCoord(perpopRED4)

# plot on a map the density of unrelated Language Families geographically close
# fig S3C


base_plot + geom_point(data = perpopREDSHIFTMAP,
                       aes(x = lon.1, y = lat.1),
                       color = "black",
                       size = 0.5, shape=3)+
  geom_point(data=perpopRED4_SHIFTMAP,  aes(x = lon.1, y = lat.1,color=as.numeric(GeoCloseUnrelatedFamilies)), alpha=0.5, size=2)+
  # ggtitle("number of neighbors from a different language family")+
  scale_color_gradient(low = "blue", high = "red") +
  xlab("") + ylab("") +
  theme(legend.title = element_blank())
ggsave("DensityDifferentLFamilyNeighbors_noContinentFilter_map_PacificCenter_1000km_2021.pdf", useDingbats=FALSE, width = 15, height = 8, dpi = 300, units = "in")



#******************************************
#******************************************
#******************************************
#******************************************
# ANALYSIS PERCENTILES
#******************************************
# OVERVIEW OF MISMATCHES WITH CLOSE FST DISTANCES
#******************************************
#******************************************
#* Figure S6

# ------------------------------------------------------
# CONTINENT MEDIAN FST and a range of percentile threshold values



## adjust dataset without very Drifted and within and between regions
FstListGlotto_infowithinREgion<-FstListREDinfo[which(FstListREDinfo$region1==FstListREDinfo$region2),]

## exclude drifted pops or the Fst averages will be higher than normal
DRIFTONI<-perpopRED[which(perpopRED$medianFSTregion>0.1&perpopRED$medianFST>0.1),]$PopName

FstListGlotto_infowithinREgionNoDrif<-FstListGlotto_infowithinREgion[-c(which(FstListGlotto_infowithinREgion$Pop%in%DRIFTONI), which(FstListGlotto_infowithinREgion$Pop1%in%DRIFTONI)),]
FstListGlottoIBD_infoNoDrift<-FstListREDinfo[-c(which(FstListREDinfo$Pop2%in%DRIFTONI), which(FstListREDinfo$Pop1%in%DRIFTONI)),]

# set a series of percentiles up to 0.5
percentiles<-seq(from = .002, to = .5, by = .002)
#percentilesCoarse<-seq(from = .005, to = .5, by = .005)

# calculate thresholds per continents

library(plyr)
medianGEO<-ddply(FstListGlotto_infowithinREgionNoDrif, "region1", function(x) median(x$Fst))
valuespercentileRegions<-data.frame(row.names =medianGEO$region1)
for (i in 1:length(percentiles)){
  percTarget<-percentiles[i]
  temp<-ddply(FstListGlotto_infowithinREgionNoDrif, "region1", function(x) quantile(x$Fst, percTarget))
  valuespercentileRegions[,i]<-temp[,2]
  #colnames(valuespercentileRegions[i])<-colnames(temp)[2]
}
colnames(valuespercentileRegions)<-percentiles
medianGEO<-cbind(medianGEO,valuespercentileRegions)
colnames(medianGEO)[2]<-"medianRegion"

aggiuntaWorld<-c()  # the list of percentile threshold on global distribution

for (i in 1:length(percentiles)){
  percTarget<-percentiles[i]
  aggiuntaWorld[i]<-quantile(FstListGlottoIBD_infoNoDrift$Fst, percTarget)
}
aggiuntaWorldline<-c("WORLD", median(FstListGlottoIBD_infoNoDrift$Fst),aggiuntaWorld) 
medianGEO<-rbind(medianGEO,aggiuntaWorldline)


# assign the lowest percentile to each pair fst

FstListREDinfobeforetest<-FstListREDinfo
FstListREDinfobeforetesttoday<-FstListREDinfo
FstListREDinfo$percentileFST<-NA

for (i in 1:nrow(FstListREDinfo)){
  if(FstListREDinfo$region1[i]!=FstListREDinfo$region2[i] ){
    FstListREDinfo$percentileFST[i]<- percentiles[which(FstListREDinfo$Fst[i]<aggiuntaWorld)][1]
  }
  else {
    reftemp<-  valuespercentileRegions[which( medianGEO$region1== FstListREDinfo$region2[i]),] # list of percentiles corresponding to the continent of the two populations
    FstListREDinfo$percentileFST[i]<- percentiles[which(FstListREDinfo$Fst[i]<reftemp)][1]
    
  }
}

# prepare plottable file for pairwise connections

betweenfamilSMALLgeoplottabile1<-FstListREDinfo
betweenfamilSMALLgeoplottabile1$index<-c(1:nrow(FstListREDinfo))
betweenfamilSMALLgeoplottabile1<-betweenfamilSMALLgeoplottabile1[,-which(colnames(betweenfamilSMALLgeoplottabile1)%in%c("lat2","lon2"))]

betweenfamilSMALLgeoplottabile2<-FstListREDinfo
betweenfamilSMALLgeoplottabile2$index<-c(1:nrow(FstListREDinfo))
betweenfamilSMALLgeoplottabile2<-betweenfamilSMALLgeoplottabile2[,-which(colnames(betweenfamilSMALLgeoplottabile2)%in%c("lat1","lon1"))]

betweenfamilSMALLgeoplottabile<-rbind(betweenfamilSMALLgeoplottabile1,setNames(betweenfamilSMALLgeoplottabile2,names(betweenfamilSMALLgeoplottabile1)))
betweenfamilSMALLgeoplottabile$lat=as.numeric(as.character(betweenfamilSMALLgeoplottabile$lat1))
betweenfamilSMALLgeoplottabile$lon=as.numeric(as.character(betweenfamilSMALLgeoplottabile$lon1))
betweenfamilSMALLgeoplottabile$index=as.numeric(as.character(betweenfamilSMALLgeoplottabile$index))

betweenfamilSMALLgeoplottabile<-betweenfamilSMALLgeoplottabile[!is.na(betweenfamilSMALLgeoplottabile$lat),]
betweenfamilSMALLgeoplottabile<-betweenfamilSMALLgeoplottabile[!is.na(betweenfamilSMALLgeoplottabile$percentileFST),]
#betweenfamilSMALLgeoplottabile$percentileFST<-1-betweenfamilSMALLgeoplottabile$percentileFST
betweenfamilSMALLgeoplottabileSHIFTMAP<-MoveDataCoord(betweenfamilSMALLgeoplottabile) 

# one map with increasing color connections for increasing percentile
# from percentile 0.2 to 0.01

betweenfamilSMALLgeoplottabileSHIFTMAP2<-betweenfamilSMALLgeoplottabileSHIFTMAP[which(betweenfamilSMALLgeoplottabileSHIFTMAP$percentileFST<0.1),]

base_plot + geom_point(data = perpopREDSHIFTMAP,
                       aes(x = lon.1, y = lat.1),
                       color = "black",
                       size = 0.5, shape=3)+
  geom_path(data=betweenfamilSMALLgeoplottabileSHIFTMAP2, aes(x = lon.1, y = lat.1, 
                                                              group=index, alpha=percentileFST, size=percentileFST, color=percentileFST))+
  scale_alpha(range = c(0.05,0.1))+
  scale_size(range = c(0.5,.01))+
  scale_colour_gradient2(midpoint=0.1, low="darkred", mid="darkorange",
                         high="white", name="Percentile FST distribution")+
  #ggtitle("weight of mismatches according to FST percentile distribution")+
  xlab("") + ylab("") +
  guides(alpha=FALSE, size=FALSE)


### now in Supplementary!!!
ggsave("MapPairMismatchesPercentileDistribution_until10percentile_2021.pdf", useDingbats=FALSE, width = 15, height = 8, dpi = 300, units = "in")



##-------------------------------------
# figure  area proportion of linguistically unrelated in close FSTs
# now supplementary Fig S6B
##-------------------------------------



#FstListREDinfo_diverse<-FstListREDinfo[which(FstListREDinfo$FAMILY=="DIVERSE"),]
FstListGlotto_infowithinREgion<-FstListREDinfo[which(FstListREDinfo$region1==FstListREDinfo$region2),]
FstListREDinfo_BetweenRegion<-FstListREDinfo[which(FstListREDinfo$region1!=FstListREDinfo$region2),]



numberTotalSmallFST<-c()
numberTotalSmallFSTDIVERSE<-c()
percentageMismatchesinTotalSmallFST<-c()

#write.table(medianGEO, "medianGEO.txt", sep = " ")

meltFstGlottoIBD_infoWithPercentageMismatch<-data.frame(row.names = FstListREDinfo$popslistemp)

for (j in 1:length(percentiles)){
  mismatchlist<-c()
  percTarget<-percentiles[j]
  for (i in 1:5){  # the five macro continents
    mismatchlist1<-FstListGlotto_infowithinREgion[which(FstListGlotto_infowithinREgion$region1==medianGEO[i,1]),] # select continent
    thresholdtemp<-as.numeric(medianGEO[i,which(colnames(medianGEO)==percTarget)]) # select percentile value fst
    mismatchlist1<-mismatchlist1[which(mismatchlist1$Fst<= thresholdtemp),]
    mismatchlist<-rbind(mismatchlist,mismatchlist1)
  }
  thresholdtempworld<-as.numeric(medianGEO[6,which(colnames(medianGEO)==percTarget)])
  mismatchlistBEtweenregion<-FstListREDinfo_BetweenRegion[which(FstListREDinfo_BetweenRegion$Fst<= thresholdtempworld),] # use the target quantile all over the world Fst WORLD
  TheMismatches<-rbind(mismatchlist, mismatchlistBEtweenregion)
  numberTotalSmallFST[j]<-nrow(TheMismatches) # number of pairs within FST percentile
  numberTotalSmallFSTDIVERSE[j]<-length(which(TheMismatches$FAMILY=="DIVERSE"))  # number of pairs in mismatch
  
  percentageMismatchesinTotalSmallFST[j]<-length(which(TheMismatches$FAMILY=="DIVERSE"))/nrow(TheMismatches)
  meltFstGlottoIBD_infoWithPercentageMismatch[,j]<-NA
  meltFstGlottoIBD_infoWithPercentageMismatch[TheMismatches$popslistemp,j]<-percTarget
}
colnames(meltFstGlottoIBD_infoWithPercentageMismatch)<-percentiles
#meltFstGlottoIBD_infoWithPercentageMismatchSort<-
# meltFstGlottoIBD_infoWithPercentageMismatch[FstListREDinfo$popslistemp,]

listperpercentileNumberCases<-rbind(numberTotalSmallFST,numberTotalSmallFSTDIVERSE,percentageMismatchesinTotalSmallFST)
colnames(listperpercentileNumberCases)<-percentiles

# count the cases below each percentile
f<-function(x, output){sort(x)[1]}
FstListREDinfo$percentileFSTdoublecheck<-apply(meltFstGlottoIBD_infoWithPercentageMismatch,
                                               1, f)

plot(FstListREDinfo$percentileFSTdoublecheck, FstListREDinfo$percentileFST) # check if the two slightly different method give the same percentiles

#plotproportion<-data.frame(percentiles,percentageMismatchesinTotalSmallFST)
plotproportion<-data.frame(percentiles,nrow(FstListREDinfo)-numberTotalSmallFST,numberTotalSmallFST-numberTotalSmallFSTDIVERSE,numberTotalSmallFSTDIVERSE)
colnames(plotproportion)<-c("percentiles","pairsOutside","pairsInsideSameFamily","pairsInsideDiffFamily")

ggplot(plotproportion,aes(x=as.numeric(percentiles), y=percentageMismatchesinTotalSmallFST, color=as.numeric(percentiles)))+
  geom_segment(aes(xend=as.numeric(percentiles), yend=0, color=as.numeric(percentiles)), alpha=0.9)+
  geom_line(size=1)+
  # geom_area(aes(x=as.numeric(percentiles), y=percentageMismatchesinTotalSmallFST), alpha=0.5, fill="darkorange", color="darkorange")+
  #ggtitle("proportion of couples from Diff L Families over number of pairs genetically close - FST percentile distrib")+
  labs(x="Percentile global Fst distribution", y="proportion of pairs from Different Language Families")+
  #geom_vline(xintercept = 0.2, color="darkorange", size=1)+
  scale_colour_gradient2(midpoint=0.2, low="darkred", mid="darkorange",
                         high="white", name="")+
  theme_light() +
  theme(legend.position = "none", text=element_text(size=5))

ggsave("ProportionaMismatchesSensitivityThresholdFST_2021_better.pdf", useDingbats=FALSE, width = 4, height = 2)



### for each population, i annotate the percentile FST threshold in which they appear in a mismatch pair
perpopRED$percentileMismatch<-NA
singlepopinmismatch<-c()

for (i in length(percentiles):1){
  percTarget<-percentiles[i]
  temp<-FstListREDinfo[which(FstListREDinfo$percentileFST<=percTarget),]
  listoni<-  unique(temp$Pop1)
  perpopRED$percentileMismatch[which(perpopRED$PopName%in%listoni)]<-percTarget
  singlepopinmismatch[i]<-length(listoni)
}



# CONNECT MISMATCHES BETWEEN MAJOR FAMILIES
# with a big circle

##### FIGURE S6C  ****************************************

perpopREDfamily<-perpopRED[which(perpopRED$glottolog.NAME %in% MainFamilies),]

meltFstGlottoIBD_infoLessPercentile<-FstListREDinfo[which(FstListREDinfo$percentileFST<0.1),] # the closest 0.1 percentile


meltFstREDinfoMAJORFAM<-meltFstGlottoIBD_infoLessPercentile[which(meltFstGlottoIBD_infoLessPercentile$FAMILY=="DIVERSE"),]
meltFstREDinfoMAJORFAM<-meltFstREDinfoMAJORFAM[which(meltFstREDinfoMAJORFAM$family1%in%MainFamilies&meltFstREDinfoMAJORFAM$family2%in%MainFamilies),]

tablonMismatchMajorFamilies<-table(meltFstREDinfoMAJORFAM$family1,meltFstREDinfoMAJORFAM$family2)


mat<-as.matrix(tablonMismatchMajorFamilies)

orderfamily<-c("Kxa","Khoe-Kwadi" ,"Atlantic-Congo","Afro-Asiatic","Turkic"  ,"Indo-European" ,"Uralic",
               "Abkhaz-Adyge" ,"Nakh-Daghestanian", "Tungusic", "Mongolic-Khitan","Sino-Tibetan",
               "Austronesian" ,"Quechuan","Tupian")
orderspecial<-orderfamily[which(orderfamily%in%rownames(mat))]

matorder<-mat[orderspecial,]
matorder2<-matorder[,orderspecial]

rownames(MainFamilies2)<-MainFamilies2$MainFamilies
MainFamiliesRED<-MainFamilies2[orderspecial,]

grid.col<-(MainFamiliesRED$COLOR)
names(grid.col)<-MainFamiliesRED$MainFamilies # create a color combination with the color palette

library(circlize)

# adjust for sample size TOTAL DATASET between major families
meltFstREDinfo222<-FstListREDinfo[which(FstListREDinfo$family1%in%MainFamilies),]
meltFstREDinfo222<-meltFstREDinfo222[which(meltFstREDinfo222$family2%in%MainFamilies),]

#samplesize1<-table(c(meltFstREDinfo222$family1,meltFstREDinfo222$family2)) # all the possible pairs for each language family but only between other main language families
samplesize1<-table(meltFstREDinfo222$family1)

orderspecialsamplesize<-samplesize1[orderspecial]

matAdjustSize<-matorder2
for(j in 1:nrow(matorder2)){
  for (k in 1:ncol(matorder2)){
    matAdjustSize[j,k]<-((matorder2[j,k]/(orderspecialsamplesize[j]*orderspecialsamplesize[k])))
  }
}
matAdjustSizeROUND<-round(matAdjustSize/min(matAdjustSize[-which(matAdjustSize==0)]))
matAdjustSizeROUND<-as.matrix(matAdjustSizeROUND)

# set color proportional to average lower percentile FST


matFSTpercentile<-matAdjustSizeROUND
for(j in 1:nrow(matAdjustSizeROUND)){
  for (k in 1:ncol(matAdjustSizeROUND)){
    coppia<-c(rownames(matAdjustSizeROUND)[j],colnames(matAdjustSizeROUND)[k])
    templist<-meltFstREDinfoMAJORFAM[which(meltFstREDinfoMAJORFAM$family1%in%coppia),]
    #    templist2<-meltFstREDinfoMAJORFAM[which(meltFstREDinfoMAJORFAM$family2==rownames(matorder2)[j]&meltFstREDinfoMAJORFAM$family1==colnames(matorder2)[k]),]
    #   templist<-rbind(templist, templist2)
    matFSTpercentile[j,k]<-median(templist$percentileFST)
  }
}

matFSTpercentilemelt<-melt(matFSTpercentile)
listvalues<-sort(unique(matFSTpercentile))

sort(unique(matFSTpercentile))


[1] 0.046 0.048 0.049 0.050 0.052 0.056 0.057 0.058 0.059 0.060 0.062 0.063 0.064 0.065 0.066 0.068 0.070 0.071 0.074 0.076 0.078
[22] 0.086


colfunc <- colorRampPalette(c("darkred", "darkorange"))
listcolors<-colfunc(length((listvalues)))
listvalues<-cbind(listvalues,listcolors)

matFSTpercentilemeltcolor<- listvalues[,2][match(matFSTpercentilemelt[,3],listvalues[,1])]

matpercentilecolor<-matrix(matFSTpercentilemeltcolor,
                           nrow(matAdjustSizeROUND),
                           ncol(matAdjustSizeROUND))


pdf("circlize_Mismatches_families_adjustSampleSize_below10percentile.pdf", width = 8, height = 8)
ff<-chordDiagramFromMatrix(matAdjustSizeROUND, directional = 0, 
                           transparency = 0.1, symmetric=T, order = orderspecial, 
                           col=matpercentilecolor,grid.col = grid.col)
dev.off()




#******************************************
## figure distribution FST SAME OR DIFFERENT FAMILY for each population
##supplementary Figure S9
#******************************************


pdf("EACHPOP_violin_YESorNO_GEOradius_2021_V2.pdf")  

for (i in 1:nrow(perpopRED)){
  TARGET<-perpopRED$PopName[i]
  meltini<-FstListREDinfo_noDuplicateNeighbors[which(FstListREDinfo_noDuplicateNeighbors$Pop1==TARGET),]
  meltiniSAME<-meltini[which(meltini$SameFamily=="YES"),]
  maxgeofam<-max(meltiniSAME$GEOdist, na.rm = T)
  if(maxgeofam<500){
    maxgeofam<-500
  }
  meltini<-meltini[which(meltini$GEOdist<=maxgeofam),]

  ggg<-ggplot(meltini,aes(y=Fst, x=SameFamily,color=SameFamily))
  
  ggg<-ggg+ geom_jitter(shape=16, position=position_jitter(0.2))+
    geom_violin(trim=FALSE, alpha=0.4)+
    stat_summary(fun.data="mean_sdl", 
                 geom="pointrange", color="darkmagenta") +
    ggtitle(paste0(TARGET," - ",perpopRED$glottolog.NAME[i]," - ",maxgeofam, "km radius"), 
            subtitle = paste0("W test: ", round(perpopRED$WtestGEOfilter[i]), ", p value:", round(perpopRED$WtestPvalueGEOfilter[i], digits = 4))) +
    theme_light()
  
  print(ggg)
  
}
dev.off()



#******************************************
#*#******************************************
#*#******************************************
#*#******************************************
#

####
# Figure S7 Comparison of medians on a map
########

perpopREDSHIFTMAP_FSTless<-perpopREDSHIFTMAP[which(perpopREDSHIFTMAP$WtestPvalueGEOfilter>0.01),]
perpopREDSHIFTMAP_3wilcox<-perpopREDSHIFTMAP[which(perpopREDSHIFTMAP$WtestPvalueGEOfilter<=0.01),]

base_plot + geom_point(data = perpopREDSHIFTMAP,
                       aes(x = lon.1, y = lat.1),
                       color = "black",
                       size = 0.5, shape=3)+
  geom_point(data = perpopREDSHIFTMAP,
             aes(x = lon.1, y = lat.1,color = WtestPvalueGEOfilter),
             size = 1.2, alpha=0.8)+
  scale_color_gradient2( low="darkred",mid="pink",
                         high="blue", midpoint = 0.01, na.value = "grey50", name="P value Wilcoxon Test")+
  geom_text(data=perpopREDSHIFTMAP_2wilcox, aes(x = lon.1, y = lat.1,label=PopName), size=0.8)+
  geom_point(data = perpopREDSHIFTMAP_3wilcox, aes(x = lon.1, y = lat.1), color = "darkred",  size = 1.2, alpha=0.5)+
  xlab("") + ylab("") +
  guides(alpha=FALSE, size=FALSE)


ggsave("MapWilcoxonTest_2021_.pdf", useDingbats=FALSE, width = 15, height = 8, dpi = 300, units = "in")

## FST proportion with neighbors 
# does not work on a global scale, too strong differences

base_plot + geom_point(data = perpopREDSHIFTMAP,
                       aes(x = lon.1, y = lat.1),
                       color = "black",
                       size = 0.5, shape=3)+
  geom_point(data = perpopREDSHIFTMAP[which(perpopREDSHIFTMAP$proportionFstAdjustedNeighbors<2),],
             aes(x = lon.1, y = lat.1,color = proportionFstAdjustedNeighbors),
             size = 1.2, alpha=0.8)+ 
  scale_color_viridis_c(option = "plasma")

  scale_colour_gradient2(midpoint=1, low="green", mid="yellow",
                                                            high="red",name="Proportion FST with neighbors")+
  
  scale_colour_gradientn(colours = rainbow)
  
  scale_color_viridis_c(option = "plasma")
  scale_colour_gradient2(midpoint=1, low="green", mid="gray",
                         high="orange",name="Proportion FST with neighbors")+
#  scale_color_gradient2( low="darkred",mid="pink", high="blue", midpoint = 0.01, na.value = "grey50", name="P value Wilcoxon Test")+
 # geom_text(data=perpopREDSHIFTMAP_2wilcox, aes(x = lon.1, y = lat.1,label=PopName), size=0.8)+
#  geom_point(data = perpopREDSHIFTMAP_3wilcox, aes(x = lon.1, y = lat.1), color = "darkred",  size = 1.2, alpha=0.5)+
  xlab("") + ylab("") +
  guides(alpha=FALSE, size=FALSE)




#********************************************
#********************************************
#********************************************
# ***************************************************
# ANALYSIS OF LANGUAGE ISOLATES
# ***************************************************
  #  Figure S10
  

library(ggrepel)
perpopREDISOLATE<-perpopRED[which(perpopRED$isolate=="Isolate"),]
library(patchwork)
gg<-ggplot(perpopRED)

gg<-ggplot(perpopRED) + geom_point(aes(proportionHeterozyAdjustedNeighbors,proportionFstAdjustedNeighbors), alpha=0.3)+
  geom_point(data=perpopREDISOLATE,aes(proportionHeterozyAdjustedNeighbors,proportionFstAdjustedNeighbors), color="red" )+
  theme_light()+
  geom_text_repel(data=perpopREDISOLATE,aes(proportionHeterozyAdjustedNeighbors,proportionFstAdjustedNeighbors, color=isolate, label=PopName), color="red" )+
 xlab("proportion of heterozygozity compared to neighbors")+ylab("proportion of median FST compared to neighbors")
 
#ggsave("LangIsolate_plotpoportionFSTAndmedianFST.pdf", useDingbats=FALSE)

# supplementary Isolates
ga<-ggplot(perpopRED,aes(x=isolate,y=proportionFstAdjustedNeighbors))+
 geom_jitter(shape=16, position=position_jitter(0.2))+
  geom_violin(trim=FALSE, alpha=0.4)+
  stat_summary(fun.data="mean_sdl", mult=1, 
               geom="pointrange", color="orange")+
  theme_light()+
  xlab("")+ylab("proportion of median FST")
#ggsave("LangIsolate_violinpoportionFST_2021.pdf", useDingbats=FALSE, height = 4, width = 6, units = "in")

(ga + gg) 
ggsave("LangIsolate_suppldouble_2021.pdf", useDingbats=FALSE, height = 6, width = 12, units = "in")

   #  Figure S10



#********************************************
#********************************************
#********************************************
# ***************************************************
# ANALYSIS OF REGIONAL CASE STUDIES
# ***************************************************

#  EUROPE 

coordinates<-c(35,56,-9,27)
perpopREDEUROPE<-perpopRED[which(perpopRED$lat>coordinates[1]&perpopRED$lat<coordinates[2]&perpopRED$lon>coordinates[3]&perpopRED$lon<coordinates[4]),]
map.world <- map_data(map="world")
library(ggrepel)

gg <- ggplot()
gg <- gg + theme()
gg <- gg + geom_map(data=map.world, map=map.world, aes(map_id=region), fill="white", colour="black", size=0.15)
gg<- gg+coord_quickmap(ylim=coordinates[1:2], xlim=coordinates[3:4])

FSTpropEUR<-gg + geom_point(data=perpopREDEUROPE[!is.na(perpopREDEUROPE$proportionFstAdjustedNeighbors),], 
                          
                          aes(x=lon, y=lat, color=proportionFstAdjustedNeighbors), size=3 , alpha=0.4) +
  scale_colour_gradient2(midpoint=1, low="darkgreen", mid="gray",
                         high="orange")+
  geom_label_repel(data=perpopREDEUROPE, aes(x=lon, y=lat,label=PopName, color=proportionFstAdjustedNeighbors), size=2, label.padding=0.1)+
  geom_point(data=perpopREDEUROPE, aes(x=lon, y=lat, shape=glottolog.NAME), size=1 ) +
  #ggtitle("Proportion FST with neighbors")+
  labs(shape="Language Family", colour="Poportion FST with neighbors")+
  theme(axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank(),
        legend.title = element_text(size = 8), legend.text = element_text(size = 8),
        legend.key.size = unit(0.4, 'cm'))



diffFST_EUR<- gg + geom_point(data=perpopREDEUROPE, 
                          aes(x=lon, y=lat, color=MedDiff,na.rm=T, size=MedDiffCI,
                              alpha=MedDiffCI) ) +
  scale_color_gradient2( low="yellow",mid="gold",
                         high="darkblue",  midpoint = 0, na.value = "grey50", name="Difference FST between/within")+
  geom_label_repel(data=perpopREDEUROPE, aes(x=lon, y=lat,label=PopName, color=MedDiff, na.rm= TRUE ), size=2, label.padding=0.1)+
  geom_point(data=perpopREDEUROPE, aes(x=lon, y=lat, shape=glottolog.NAME) ) +
  labs(shape="Language Family", colour="Difference between FST distribution between and within")+
  scale_alpha(range = c(0.8,0.2),name="Difference FST CI")+
  scale_size(range = c(7,2),name="Difference FST CI")+
  theme(axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank(),
        legend.title = element_text(size = 8), legend.text = element_text(size = 8),
        legend.key.size = unit(0.4, 'cm'))





# *********************************************
##### AFRICA #####
# *********************************************

# MAP OF sub saharan africa 

coordinates<-c(-35,10,-5,45)
perpopREDAFRICA<-perpopRED[which(perpopRED$lat>coordinates[1]&perpopRED$lat<coordinates[2]&perpopRED$lon>coordinates[3]&perpopRED$lon<coordinates[4]),]


gg <- ggplot()
gg <- gg + theme()
gg <- gg + geom_map(data=map.world, map=map.world, aes(map_id=region), fill="white", colour="black", size=0.15)
gg<- gg+coord_quickmap(ylim=coordinates[1:2], xlim=coordinates[3:4])

fstAfrica<-    gg + geom_point(data=perpopREDAFRICA[!is.na(perpopREDAFRICA$proportionFstAdjustedNeighbors),], 
                          
                          aes(x=lon, y=lat, color=proportionFstAdjustedNeighbors), size=3 , alpha=0.4) +
  scale_colour_gradient2(midpoint=1, low="darkgreen", mid="gray",
                         high="orange")+
  geom_label_repel(data=perpopREDAFRICA, aes(x=lon, y=lat,label=PopName, color=proportionFstAdjustedNeighbors), size=2, label.padding=0.1)+
  geom_point(data=perpopREDAFRICA, aes(x=lon, y=lat, shape=glottolog.NAME), size=1 ) +
    labs(shape="Language Family", colour="Poportion FST with neighbors")+
  scale_shape_manual(values=1:10)+
  theme(axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank(),
        legend.title = element_text(size = 8), legend.text = element_text(size = 8),
        legend.key.size = unit(0.4, 'cm'))


diffFST_AFRICA<- gg + geom_point(data=perpopREDAFRICA, 
                              aes(x=lon, y=lat, color=MedDiff,na.rm=T,size=MedDiffCI,
                                  alpha=MedDiffCI) )  +
  scale_alpha(range = c(0.8,0.2),name="Difference FST CI")+
  scale_size(range = c(7,2),name="Difference FST CI")+
  scale_color_gradient2( low="yellow",mid="gold",
                         high="darkblue",  midpoint = 0, na.value = "grey50", name="Difference FST between/within")+
  geom_label_repel(data=perpopREDAFRICA, aes(x=lon, y=lat,label=PopName, color=MedDiff, na.rm= TRUE ), size=2, label.padding=0.1)+
  geom_point(data=perpopREDAFRICA, aes(x=lon, y=lat, shape=glottolog.NAME) ) +
  labs(shape="Language Family", colour="Difference between FST distribution between and within")+
  scale_shape_manual(values=1:10)+
  theme(axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank(),
        legend.title = element_text(size = 8), legend.text = element_text(size = 8),
        legend.key.size = unit(0.4, 'cm'))




# *********************************************
##### OCEANIA #####

# MAP OF OCEANIA
library(ggrepel)

coordinates<-c(-25,24,92,170)
perpopREDOCEANIA<-perpopRED[which(perpopRED$lat>coordinates[1]&perpopRED$lat<coordinates[2]&perpopRED$lon>coordinates[3]&perpopRED$lon<coordinates[4]),]
map.world <- map_data(map="world")

gg <- ggplot()
gg <- gg + theme()
gg <- gg + geom_map(data=map.world, map=map.world, aes(map_id=region), fill="white", colour="black", size=0.15)
gg<- gg+coord_quickmap(ylim=coordinates[1:2], xlim=coordinates[3:4])

fstOCEANIA<-  gg + geom_point(data=perpopREDOCEANIA[!is.na(perpopREDOCEANIA$proportionFstAdjustedNeighbors),], 
                        
                        aes(x=lon, y=lat, color=proportionFstAdjustedNeighbors), size=3 , alpha=0.4) +
  scale_colour_gradient2(midpoint=1, low="darkgreen", mid="gray",
                         high="orange")+
  geom_label_repel(data=perpopREDOCEANIA, aes(x=lon, y=lat,label=PopName, color=proportionFstAdjustedNeighbors), size=2, label.padding=0.1)+
  geom_point(data=perpopREDOCEANIA, aes(x=lon, y=lat, shape=glottolog.NAME), size=1 ) +
  labs(shape="Language Family", colour="Poportion FST with neighbors")+
  scale_shape_manual(values=1:13)+
  theme(axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank(),
        legend.title = element_text(size = 8), legend.text = element_text(size = 8),
        legend.key.size = unit(0.4, 'cm'))


diffFST_OCEANIA<- gg + geom_point(data=perpopREDOCEANIA, 
                                  aes(x=lon, y=lat, color=MedDiff,na.rm=T,size=MedDiffCI,
                                      alpha=MedDiffCI) )  +
  scale_alpha(range = c(0.8,0.2),name="Difference FST CI")+
  scale_size(range = c(7,2),name="Difference FST CI")+
  scale_color_gradient2( low="yellow",mid="gold",
                         high="darkblue",  midpoint = 0, na.value = "grey50", name="Difference FST between/within")+
  geom_label_repel(data=perpopREDOCEANIA, aes(x=lon, y=lat,label=PopName, color=MedDiff, na.rm= TRUE ), size=2, label.padding=0.1)+
  geom_point(data=perpopREDOCEANIA, aes(x=lon, y=lat, shape=glottolog.NAME), size=1 ) +
  labs(shape="Language Family", colour="Difference between FST distribution between and within")+
  scale_shape_manual(values=1:13)+
  theme(axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank(),
        legend.title = element_text(size = 8), legend.text = element_text(size = 8),
        legend.key.size = unit(0.4, 'cm'))


# *********************************************
##### CAUCASUS #####

# MAP OF caucasus
library(ggrepel)

coordinates<-c(38,45,34,50)
perpopREDcaucasus<-perpopRED[which(perpopRED$lat>coordinates[1]&perpopRED$lat<coordinates[2]&perpopRED$lon>coordinates[3]&perpopRED$lon<coordinates[4]),]
map.world <- map_data(map="world")

gg <- ggplot()
gg <- gg + theme()
gg <- gg + geom_map(data=map.world, map=map.world, aes(map_id=region), fill="white", colour="black", size=0.15)
gg<- gg+coord_quickmap(ylim=coordinates[1:2], xlim=coordinates[3:4])


fstCaucasus<-gg + geom_point(data=perpopREDcaucasus[!is.na(perpopREDcaucasus$proportionFstAdjustedNeighbors),], 
                      
                      aes(x=lon, y=lat, color=proportionFstAdjustedNeighbors), size=3 , alpha=0.4) +
  scale_colour_gradient2(midpoint=1, low="darkgreen", mid="gray",
                         high="orange")+
  geom_label_repel(data=perpopREDcaucasus, aes(x=lon, y=lat,label=PopName, color=proportionFstAdjustedNeighbors), size=2, label.padding=0.1)+
  geom_point(data=perpopREDcaucasus, aes(x=lon, y=lat, shape=glottolog.NAME), size=1 ) +
  labs(shape="Language Family", colour="Poportion FST with neighbors")+
  scale_shape_manual(values=1:7)+
  theme(axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank(),
        legend.title = element_text(size = 8), legend.text = element_text(size = 8),
        legend.key.size = unit(0.4, 'cm'))


diffFST_CAUCASUS<- gg + geom_point(data=perpopREDcaucasus, 
                                   aes(x=lon, y=lat, color=MedDiff,na.rm=T,size=MedDiffCI,
                                       alpha=MedDiffCI) )  +
  scale_alpha(range = c(0.8,0.2),name="Difference FST CI")+
  scale_size(range = c(7,2),name="Difference FST CI")+
  scale_color_gradient2( low="yellow",mid="gold",
                         high="darkblue",  midpoint = 0, na.value = "grey50", name="Difference FST between/within")+
  geom_label_repel(data=perpopREDcaucasus, aes(x=lon, y=lat,label=PopName, color=MedDiff, na.rm= TRUE ), size=2, label.padding=0.1)+
  geom_point(data=perpopREDcaucasus, aes(x=lon, y=lat, shape=glottolog.NAME), size=1 ) +
  labs(shape="Language Family", colour="Difference between FST distribution between and within")+
  scale_shape_manual(values=1:7)+
  theme(axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank(),
        legend.title = element_text(size = 8), legend.text = element_text(size = 8),
        legend.key.size = unit(0.4, 'cm'))


##+++++++++++++++++++++++
# main figure combined
# 1: europe, 2: oceania, 3: africa, 4: caucasus
#***************************************************

library(ggpubr)


ggarrange(fstAfrica, FSTpropEUR, fstCaucasus,fstOCEANIA + rremove("x.text"), 
 #         labels = c("A", "B", "C", "D"),
          ncol = 2, nrow = 2)
ggsave("combinedMAP_FSTproportion_supple2021.pdf", useDingbats=FALSE, height = 15, width = 17)

#*Figure S7B

ggarrange(diffFST_AFRICA, diffFST_EUR, diffFST_CAUCASUS,diffFST_OCEANIA
          + rremove("x.text"), 
 #         labels = c("A", "B", "C", "D"),
          ncol = 2, nrow = 2)
ggsave("combinedMAP_FSTproportion_supple2021.pdf", useDingbats=FALSE, height = 15, width = 17)

#*Figure S10C


library(patchwork)
(fstAfrica + FSTpropEUR) | (fstCaucasus + fstOCEANIA)
ggsave("combinedMAP_FSTproportion_supple2021_nolabels.pdf", useDingbats=FALSE, height = 15, width = 17)

#***************************************************
### SINGLE POPS SPECIAL CASES
#***************************************************

# BASQUE

basqlist<-c("Spanish_Basque", "Basque", "Spanish_PaisVasco")
basq<-FstListREDinfo[c(which(FstListREDinfo$Pop1%in%basqlist)),]
basq<-as.data.frame(basq)

selectionbasq<-basq[order(basq$FstLinear),][1:30,]

colorispecial<-  MainFamilies2$COLOR[match(  names(table(selectionbasq$family2)),
                                             MainFamilies2$MainFamilies)]
colorispecial[ is.na(colorispecial)]<-c("gray50") # special color code including minor language families in grayscale

agg<-ggplot(selectionbasq)
BASQUE<- agg+ geom_label_repel( aes(GEOdist,FstLinear, color=family2,label=Pop2), size=3, label.padding=0.1)+
  geom_point(aes(GEOdist,FstLinear, color=family2))+
  ggtitle("Basque")+
  scale_color_manual(values = colorispecial)+
  theme_light()+
  labs(colour="Language Family")

#***************************************************
# MALTA  ### supplementary FIG


malta<-FstListREDinfo[c(which(FstListREDinfo$Pop1=="Maltese")),]

selectionMalta<-malta[order(malta$FstLinear),][1:30,]

colorispecial<-  MainFamilies2$COLOR[match(  names(table(selectionMalta$family2)),
                                             MainFamilies2$MainFamilies)]
#colorispecial[ is.na(colorispecial)]<-c("gray15")
# special color code including minor language families in grayscale

agg<-ggplot(selectionMalta)

MALTA<- agg+ geom_label_repel( aes(GEOdist,FstLinear, color=family2,label=Pop2), size=3, label.padding=0.1)+
  geom_point(aes(GEOdist,FstLinear, color=family2))+
  ggtitle("Maltese")+  theme_light()+
  scale_color_manual(values = colorispecial)+
  labs(colour="Language Family")

#***************************************************
# SANDAWE

SAND<-FstListREDinfo[c(which(FstListREDinfo$Pop1=="Sandawe")),]

selectionSAND<-SAND[order(SAND$FstLinear),][1:30,]

colorispecial<-  MainFamilies2$COLOR[match(  names(table(selectionSAND$family2)),
                                             MainFamilies2$MainFamilies)]
colorispecial[ is.na(colorispecial)]<-c("gray15", "gray40","gray60","gray80")
# special color code including minor language families in grayscale

agg<-ggplot(selectionSAND)

SANDAWE<- agg+ geom_label_repel( aes(GEOdist,FstLinear, color=family2,label=Pop2), size=3, label.padding=0.1)+
  geom_point(aes(GEOdist,FstLinear, color=family2))+
  ylim(0.01,0.06)+  theme_light()+
  ggtitle("Sandawe")+
  scale_color_manual(values = colorispecial)+
  labs(colour="Language Family")

#***************************************************
# Hadza

Hadza<-FstListREDinfo[c(which(FstListREDinfo$Pop1=="Hadza")),]
#Hadza<-Hadza[-which(Hadza$GEOdist<40),] #eliminate the other Hadza
selectionHadz<-Hadza[order(Hadza$FstLinear),][1:30,]

colorispecial<-  MainFamilies2$COLOR[match(  names(table(selectionHadz$family2)),
                                             MainFamilies2$MainFamilies)]
colorispecial[ is.na(colorispecial)]<-c( "gray40","gray10","gray60","gray80", "gray30")
agg<-ggplot(selectionHadz)
HADZA<- agg+ geom_label_repel( aes(GEOdist,FstLinear, color=family2,label=Pop2), size=3, label.padding=0.1)+
  geom_point(aes(GEOdist,FstLinear, color=family2))+
  ggtitle("Hadza")+ 
  scale_color_manual(values = colorispecial)+ theme_light()+
  labs(colour="Language Family")

#***************************************************
## ARMENIAN

armenian<-FstListREDinfo[c(which(FstListREDinfo$Pop1=="Armenian")),]

selectionarmenian<-armenian[order(armenian$FstLinear),][1:30,]

colorispecial<-  MainFamilies2$COLOR[match(  names(table(selectionarmenian$family2)),
                                             MainFamilies2$MainFamilies)]
colorispecial[ is.na(colorispecial)]<-c("gray80")
# special color code including minor language families in grayscale

agg<-ggplot(selectionarmenian)

ARMENIAN<- agg+ geom_label_repel( aes(GEOdist,FstLinear, color=family2,label=Pop2), size=3, label.padding=0.1)+
  geom_point(aes(GEOdist,FstLinear, color=family2))+
  ggtitle("Armenian")+
  scale_color_manual(values = colorispecial)+ theme_light()+
  labs(colour="Language Family")

#***************************************************
# Azerbaijan

azerbajan<-FstListREDinfo[c(which(FstListREDinfo$Pop1=="Azeri_Azerbajan")),]

selectionazerb<-azerbajan[order(azerbajan$FstLinear),][1:30,]

colorispecial<-  MainFamilies2$COLOR[match(  names(table(selectionazerb$family2)),
                                             MainFamilies2$MainFamilies)]
colorispecial[ is.na(colorispecial)]<-c("gray80")
# special color code including minor language families in grayscale

agg<-ggplot(selectionazerb)

AZERBAJAN<- agg+ geom_label_repel( aes(GEOdist,FstLinear, color=family2,label=Pop2), size=3, label.padding=0.1)+
  geom_point(aes(GEOdist,FstLinear, color=family2))+
  ggtitle("Azeri_Azerbajan")+
  scale_color_manual(values = colorispecial)+ theme_light()+
  labs(colour="Language Family")

#***************************************************
## combine 6 figures in one single populations
# figure S12
#***************************************************

ggarrange(SANDAWE, HADZA, BASQUE, MALTA, ARMENIAN, AZERBAJAN + rremove("x.text"), 
          # labels = c("A", "B", "C", "D"),
          ncol = 2, nrow = 3)
ggsave("combinedSinglePops_FstGeodist_supple.pdf", useDingbats=FALSE, height = 15, width = 13)

(SANDAWE+ HADZA) | (BASQUE + MALTA) | (ARMENIAN + AZERBAJAN) 





######################################################################
########################################################
##############
# language family comparisons
# geographic distance and smooth or linear regression
## testing if SAME FAMILY has an effect on Fst given GeoDist
############################
############################


# Fst is modeled as the dependent variable with SameFamily as the factor and GEOdist as the covariate

# ********************************************
# per each language family





## exlude comparisons with the enclaves
excludepops<-union(ListEnclaves, Misalligned)
coso<-union(excludepops, DRIFTONI)
write.table(coso, "listaExclude.txt")
union(ListEnclaves, DRIFTONI)
FstListREDinfo_noDuplicateNeighbors<-FstListREDinfo_noDuplicateNeighbors[-c
                                                                         (which(FstListREDinfo_noDuplicateNeighbors$Pop2%in%excludepops), 
                                                                           which(FstListREDinfo_noDuplicateNeighbors$Pop1%in%excludepops)),]

#*********************************************
### SCRIPT FROM DAMIAN FOR PLOTTING 
#*********************************************
#*
#*

FAM<-MainFamiliesTab$MainFamilies

for(f in FAM) {
  FstListREDinfo_noDuplicateNeighbors[[f]]<-(FstListREDinfo_noDuplicateNeighbors$family1==f)|(FstListREDinfo_noDuplicateNeighbors$family2==f)
  
  threshold_geo<-max(FstListREDinfo_noDuplicateNeighbors$GEOdist[FstListREDinfo_noDuplicateNeighbors$FAMILY==f],
                     na.rm = T)
     if(threshold_geo<500){
      threshold_geo<-500
    }
  
  print(f)
  print(threshold_geo)
  FstListREDinfo_noDuplicateNeighbors[[f]]<-FstListREDinfo_noDuplicateNeighbors[[f]]*(FstListREDinfo_noDuplicateNeighbors$GEOdist<=threshold_geo)
  FstListREDinfo_noDuplicateNeighbors[[f]]<-sapply(FstListREDinfo_noDuplicateNeighbors[[f]],function(x) ifelse(is.na(x),FALSE,x))
  
  }

family_distances<-pivot_longer(FstListREDinfo_noDuplicateNeighbors,
                               cols = FAM,
                               names_to = "Family_plot",
                               values_to = "Include") %>%
  filter(Include==TRUE)

family_distances %>%
  ggplot(aes(y=FstLinear,x=GEOdist))+
  geom_point(aes(fill=SameFamily),alpha=0.3,size=1, shape=21, stroke=0)+
  geom_smooth(aes(color=SameFamily),se=F)+
  geom_smooth(color="#FBB13C",linetype = "dotdash",alpha=0.4, size=1, method = lm)+
  facet_wrap(~Family_plot,
             scales="free",
             nrow=length(FAM)/4)+
  theme_minimal()+
  scale_y_sqrt()+
  scale_fill_manual(values=c("#7FDBEB","#E45B5B"))+
  scale_color_manual(values=c("#198B9F","#BA1E1E"))+
  theme(legend.position = "bottom")

#ggsave("smoothregressionPlot.pdf", useDingbats=FALSE, height = 7, width = 10)
ggsave("smoothregressionPlotwithLinear.pdf", useDingbats=FALSE, height = 7, width = 10)


for (i in 1:nrow(MainFamiliesTab)){
  TARGET<-MainFamiliesTab$MainFamilies[i]
  meltini<-FstListREDinfo_noDuplicateNeighbors[which(FstListREDinfo_noDuplicateNeighbors$family1==TARGET),]
  meltiniRED<-meltini[!is.na(meltini$GEOdist),] # eliminate the rows for which i do not have geo coordinated in one of the pops
  meltiniSAME<-meltiniREDcut[which(meltiniREDcut$SameFamily=="YES"),]
  
  maxgeofam<-max(meltiniSAME$GEOdist, na.rm = T)
  if(maxgeofam<500){
    maxgeofam<-500
  }
  meltiniREDcut<-meltiniREDcut[which(meltiniREDcut$GEOdist<=maxgeofam),]
  mod2 <- lm(FstLinear ~ GEOdist, data=meltiniREDcut) # the general simple model 
  meltiniREDcut$residuals<-mod2$residuals  # add the residuals for each pair of pops according to the general simple model
  tempWithin<-meltiniREDcut[which(meltiniREDcut$FAMILY!="DIVERSE"),] # only inside the same Family
  RatioResiduals=round(length(which(tempWithin$residuals>0))/length(which(tempWithin$residuals<0)),digits = 2)
  MainFamiliesTab$ProportionNegResid_Within[i]<-length(which(tempWithin$residuals<0))/nrow(tempWithin)
}

gi<-ggplot(MainFamiliesTab, aes(y=ProportionNegResid_Within, x=MainFamilies)) + 
geom_segment( aes(x=MainFamilies, xend=MainFamilies, y=0, yend=ProportionNegResid_Within))+
  geom_point(color="#FBB13C", size=5)+
#  geom_point(aes(color=Family), size=5)+
#  scale_color_manual(values = colorchoice)+
  theme_light()+
  labs(x="", y="Proportion of negative residuals - same family")+
  theme(axis.text.x = element_text(angle = 45, vjust=1, hjust=1)) 
ggsave("lollipop_negativeresiduals.pdf", useDingbats=FALSE, width = 5, height = 5)


# *************************
# Figure S13
# *************************

MainFamiliesTab<-cbind(MainFamilies,rep(NA,length(MainFamilies)))
MainFamiliesTab<-as.data.frame(MainFamiliesTab)
MainFamiliesTab<-MainFamiliesTab[-which(MainFamiliesTab$MainFamilies=="Tupian"),]

pdf("FamiliesFSTGeo_linearModels_2021_geocut.pdf", width = 10, height = 17)

meltiniREDTOTAL<-c()

par(mfrow=c(5, 3))

for (i in 1:nrow(MainFamiliesTab)){
  TARGET<-MainFamiliesTab$MainFamilies[i]
  meltini<-FstListREDinfo_noDuplicateNeighbors[which(FstListREDinfo_noDuplicateNeighbors$family1==TARGET),]
  
  meltiniRED<-meltini[!is.na(meltini$GEOdist),] # eliminate the rows for which i do not have geo coordinated in one of the pops
  
  meltiniRED<-meltiniRED[which(meltiniRED$FstLinear<=0.1),]   #cut the Fst larger than 0.1
  meltiniRED<-meltiniRED[which(meltiniRED$GEOdist<=10000),]   #cut with a geo radius of 10,000 km
  
  maxFSTFam<-max(meltiniRED$FstLinear[which(meltiniRED$SameFamily=="YES")])
  
  meltiniREDcut<-meltiniRED[which(meltiniRED$FstLinear<=maxFSTFam),]   #cut the Fst larger than the max Fst within Family
  

  tempWithin<-meltiniREDcut[which(meltiniREDcut$FAMILY!="DIVERSE"),] # only inside the same Family
  tempbetween<-meltiniREDcut[which(meltiniREDcut$FAMILY=="DIVERSE"),] # outside the same Family
  
  if(nrow(tempbetween)<2){
    meltiniREDcut<-meltiniRED[which(meltiniRED$FstLinear<=maxFSTFam*2),]   #cut the Fst larger than the max Fst within Family
  }
  
  meltiniSAME<-meltiniREDcut[which(meltiniREDcut$SameFamily=="YES"),]
  
  maxgeofam<-max(meltiniSAME$GEOdist, na.rm = T)
  
  if(maxgeofam<500){
    maxgeofam<-500
  }
  meltiniREDcut<-meltiniREDcut[which(meltiniREDcut$GEOdist<=maxgeofam),]
  
  
  mod1 <- lm(FstLinear ~ GEOdist * SameFamily, data=meltiniREDcut) # the general mixed model with reference level as "NO" and effect of same family on slope and intercept
  mod2 <- lm(FstLinear ~ GEOdist, data=meltiniREDcut) # the general simple model 
  meltiniREDcut$residuals<-mod2$residuals  # add the residuals for each pair of pops according to the general simple model
  meltiniREDcut$targetFamily<-TARGET
  meltiniREDTOTAL<-rbind(meltiniREDTOTAL,meltiniREDcut) ## create the sum file with all the residuals for major families
  tempWithin<-meltiniREDcut[which(meltiniREDcut$FAMILY!="DIVERSE"),] # only inside the same Family
  tempbetween<-meltiniREDcut[which(meltiniREDcut$FAMILY=="DIVERSE"),] # outside the same Family
  
  RatioResiduals=round(length(which(tempWithin$residuals>0))/length(which(tempWithin$residuals<0)),digits = 2)
  RatioResidualsREST=round(length(which(tempbetween$residuals>0))/length(which(tempbetween$residuals<0)),digits = 2)
  
  reg1 <- lm(FstLinear ~GEOdist, data=meltiniREDcut[which(meltiniREDcut$SameFamily=="YES"),]) # one model for the same family
  reg2 <- lm(FstLinear ~GEOdist, data=meltiniREDcut[which(meltiniREDcut$SameFamily=="NO"),]) # one model for the different family
  MainFamiliesTab$RatioResidualsWithin[i]<-RatioResiduals
  MainFamiliesTab$RatioResidualsBetween[i]<-RatioResidualsREST
  MainFamiliesTab$ProportionNegResid_Within[i]<-length(which(tempWithin$residuals<0))/nrow(tempWithin)
  MainFamiliesTab$ProportionNegResid_Between[i]<-length(which(tempbetween$residuals<0))/nrow(tempbetween)
  
  #  MainFamiliesTab$Intercept[i]<-round(mod1$coefficients, digits = 5)[1]
  #  MainFamiliesTab$GEOdistintercept[i]<-round(mod1$coefficients, digits = 5)[2]
  #  MainFamiliesTab$SameFamilyYESintercept[i]<-round(mod1$coefficients, digits = 5)[3]
  #  MainFamiliesTab$GEOdist_SameFamilyYES_EFFECT[i]<-round(mod1$coefficients, digits = 5)[4]
  
  plot(FstLinear ~GEOdist, data=meltiniREDcut, type='n', main=TARGET)
  points(meltiniREDcut[which(meltiniREDcut$SameFamily=="NO"),]$GEOdist,meltiniREDcut[which(meltiniREDcut$SameFamily=="NO"),]$FstLinear,col="blue", pch=1)
  points(meltiniREDcut[which(meltiniREDcut$SameFamily=="YES"),]$GEOdist,meltiniREDcut[which(meltiniREDcut$SameFamily=="YES"),]$FstLinear, col="red",pch=20)
  abline(reg1, lty=1)
  abline(reg2, lty=2)
  abline(mod2, lty=3, col="gold", lwd=3)
}

MainFamiliesTab<-MainFamiliesTab[,-2]



barplot(t(as.matrix(MainFamiliesTab[,4:5])),names.arg=MainFamiliesTab[,1],
        main = "",
        xlab = "", ylab = "proportion below regression",
        # col = MainFamilies2$COLOR[-13],
        col = c( "red","blue"),
        legend.text = c("within","between"),las=2,
        beside = TRUE) # Grouped bars

dev.off()


# write.table(perpopRED, "/Users/chiarabarbieri/Library/Mobile Documents/com~apple~CloudDocs/GeLaTo/perPopMarchMami2.txt",  row.names = F, sep = "\t", quote = F)

