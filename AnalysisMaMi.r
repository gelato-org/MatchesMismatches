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
#perpopRED<-read.table("/Users/chiarabarbieri/Library/Mobile Documents/com~apple~CloudDocs/GeLaTo/perpopREDMaMi2021.txt",   sep = "\t", header=T, as.is=T)
perpopRED<-read.table("/Users/chiarabarbieri/Library/Mobile Documents/com~apple~CloudDocs/GeLaTo/perPopMarchMami2.txt",   sep = "\t", header=T, as.is=T)

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
# Basic map
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


#Matches

ListMatches<-perpopRED$PopName [which(perpopRED$Mismatch3a=="MATCH")]

### DRIFTONI
DRIFTONI<-perpopRED[which(perpopRED$medianFSTregion>0.1&perpopRED$medianFST>0.1),]$PopName


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


### manual screening
#enclavesByMistake<-c("Yoruba"  ,"Nama", "Han-NChina" , "Bengali", "Evenk_FarEast") # these pops are not real enclaves, their linguistically unrelated pair is the one driving this genetic proximity effect
#ListEnclaves<-ListEnclaves[-which(ListEnclaves%in%(enclavesByMistake))]

perpopRED$ListEnclaves<-NA
perpopRED$ListEnclaves[which(perpopRED$PopName%in%ListEnclaves)]<-"Enclave"


#*#******************************************
#*
#* potential FIGURE 1 different distribution FST closest same family and different family
#* 

aa<-perpopRED[,c(1,44,46,16)] 
meltperpop<-melt(aa)

aa$difference<-aa$closeFstSameFamily-aa$closeFstDIFFFamily
length(which(aa$difference>0))/(404-length(which(is.na(aa$difference)))) # only for the comparisons for which i have a same family FST

[1] 0.1822917

# 18 % of the pops who have another genetic population of the same language family do have closer FST with speakers of another language family


#**************************************************
### table of cases mismatch enclaves

mismatch3a<-perpopRED[grep("MISMATCH",perpopRED$Mismatch3a),]
perpopRED$DIFFFamilyClosestFAM<-perpopRED$glottolog.NAME[match(perpopRED$DIFFFamilyClosestpop,perpopRED$PopName)]
colnameschoice<-c("PopName","glottolog.NAME", "geodistSameFamily","geodistDIFFFamily","DIFFFamilyClosestpop", "DIFFFamilyClosestFAM")
mismatch3a<-mismatch3a[,colnameschoice]
#write.table(mismatch3a[order(mismatch3a$geodistDIFFFamily),],"tableSupplMismatches3a.txt", sep="\t", row.names = F, quote=F) 



#******************************************
#*#******************************************
#*#******************************************
#*#******************************************
#*#******************************************
#*#******************************************
#*    COMPARING MEDIAN FST BETWEEN AND WITHIN
# POPULATION HEURISTIC 2
#* FST distribution between and within language families MISALIGNED
#*#******************************************
#*#******************************************


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

## exclude DRIFTONI 
## exclude drifted pops or the Fst averages will be higher than normal
DRIFTONI<-perpopRED[which(perpopRED$medianFSTregion>0.1&perpopRED$medianFST>0.1),]$PopName

FstListREDinfo_noDuplicateNeighbors<-FstListREDinfo[-which(FstListREDinfo$GEOdist<10&FstListREDinfo$glottocodeBase1==FstListREDinfo$glottocodeBase2),]
FstListREDinfo_noDuplicateNeighborsNoDrift<-FstListREDinfo_noDuplicateNeighbors[-c(which(FstListREDinfo_noDuplicateNeighbors$Pop2%in%DRIFTONI), which(FstListREDinfo_noDuplicateNeighbors$Pop1%in%DRIFTONI)),]


library(pairwiseCI)

perpopRED$MedianWithinSMALLER2<-NA

for (i in 1:nrow(perpopRED)){
  TARGET<-perpopRED$PopName[i]
  meltini<-FstListREDinfo_noDuplicateNeighborsNoDrift[which(FstListREDinfo_noDuplicateNeighborsNoDrift$Pop1==TARGET),]
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
    perpopRED$MedianWithinSMALLER2[i]<-median(meltiniSAME$FstLinear)<median(meltiniNO$FstLinear) # if TRUE, the median FST within is smaller than the FST between
  }
}

perpopRED$MedDiffCI<-perpopRED$MedDiffCIupper-perpopRED$MedDiffCIlower

ListMisaligned<-perpopRED$PopName[which(perpopRED$MedDiff<0)]


perpopRED$listSingleCases<-"ND"
perpopRED$listSingleCases[which(perpopRED$PopName%in%ListEnclaves)]<-"Enclave"
perpopRED$listSingleCases[which(perpopRED$PopName%in%ListMisaligned)]<-"Misaligned"
perpopRED$listSingleCases[which(perpopRED$PopName%in%intersect(ListMisaligned, ListEnclaves))]<-"EnclaveANDMisaligned"
perpopRED$listSingleCases[which(perpopRED$PopName%in%ListMatches)]<-"Match"
perpopRED$listSingleCases[which(perpopRED$PopName%in%intersect(ListMisaligned, ListMatches))]<-"MatchBUTMisaligned"
perpopRED$listSingleCases[which(perpopRED$PopName%in%DRIFTONI)]<-"Drifted"



### 
#*#******************************************
#*  #*#******************************************
#*    #*#******************************************
#*    
#*    




########################################################
# Figure S8 supplementary - comparison of difference of medians and CI, colored per language family
########################################################
perpopREDfamily<-perpopRED[which(perpopRED$glottolog.NAME %in% MainFamilies),]
perpopREDfamilynoDrift<-perpopREDfamily[-which(perpopREDfamily$listSingleCases=="Drifted"),]

perpopREDEnclaves<-perpopRED[which(perpopRED$Mismatch3a%in%c("MATCH", "MISMATCH","secondaryMISMATCH_FSTzeroDiffFamily")),]

perpopREDEnclaves<-perpopRED[which(perpopRED$Mismatch3a%in%c("MATCH", "MISMATCH","secondaryMISMATCH_FSTzeroDiffFamily")),]
perpopRED$enclavesAgain<-NA
perpopRED$enclavesAgain[which(perpopRED$Mismatch3a=="MATCH")]<-"Match"
perpopRED$enclavesAgain[which(perpopRED$Mismatch3a=="MISMATCH")]<-"GeneticEnclave"
perpopRED$enclavesAgain[which(perpopRED$Mismatch3a=="secondaryMISMATCH_FSTzeroDiffFamily")]<-"LinguisticEnclave"

perpopREDnoDrift<-perpopRED[-which(perpopRED$listSingleCases=="Drifted"),]
library(ggrepel)
 colorchoice2<-colorchoice
 colorchoice2[10]<-"gray20"
 colorchoice2[11:16]<-colorchoice[10:15]

 ggplot(perpopREDnoDrift,aes(MedDiff,MedDiffCI, color=MainFamilies))+
  geom_point()+
  geom_text_repel(aes(label=PopName, color=MainFamilies), size=0.5)+
scale_color_manual(values = colorchoice2)+
  theme_light()+ xlab("Difference Median FST between - within") +
  ylab("Confidence Interval median difference") 
ggsave("DifferenceMedianwithText.pdf", height = 12, width = 14, useDingbats=FALSE)



fst1<-ggplot(perpopREDnoDrift,aes(MedDiff,MedDiffCI, color=MainFamilies))+
  geom_point(alpha=0.5, size=3,show.legend = FALSE)+
  #geom_text_repel(aes(label=PopName, color=MainFamilies), size=0.5)+
  scale_color_manual(values = colorchoice2)+
  theme_light() +
xlab("Difference Median FST between - within") +
  ylab("Confidence Interval median difference") +
  geom_vline(xintercept=0)+
  scale_y_log10()


ggsave("DifferenceMedianLogScale.pdf", height = 12, width = 14, useDingbats=FALSE)


fst2<-ggplot(perpopREDfamilynoDrift ,aes(MedDiff,MedDiffCI, color=MainFamilies))+
  geom_point(alpha=0.5, size=3,show.legend = FALSE)+
  geom_text_repel(aes(label=PopName, color=MainFamilies), size=0.5, show.legend = FALSE)+
  scale_color_manual(values = colorchoice)+
  theme_light()+ xlab("Difference Median FST between - within") +
  ylab("Confidence Interval median difference") +
  geom_vline(xintercept=0)+
  scale_y_log10()+
  facet_wrap(~MainFamilies,ncol=5)
ggsave("DifferenceMedianLogScaleMajorLangFam.pdf", height = 12, width = 14, useDingbats=FALSE)

ggplot(perpopREDEnclaves,aes(MedDiff,MedDiffCI,color=Mismatch3a)) +
  geom_point(alpha=0.5, size=3,show.legend = T) +
  geom_text_repel(aes(label=PopName, color=Mismatch3a), size=0.8) +
  geom_point(perpopREDnoDrift, aes(MedDiff,MedDiffCI),alpha=0.5, size=2) +
  # scale_color_manual(values = colorchoice2)+
  theme_light() +
  xlab("Difference Median FST between - within") +
  ylab("Confidence Interval median difference") +
  geom_vline(xintercept=0)+
  theme(legend.title = element_blank())+
  scale_y_log10()


fst3<-ggplot(perpopREDnoDrift,aes(MedDiff,MedDiffCI, color=enclavesAgain))+
  geom_point(alpha=0.5, size=2)+
  geom_text_repel(data=perpopREDnoDrift[!is.na(perpopREDnoDrift$enclavesAgain),],
                  aes(MedDiff,MedDiffCI,label=PopName, color=enclavesAgain), size=1.2)+
  # scale_color_manual(values = colorchoice2)+
  theme_light() +
  xlab("Difference Median FST between - within") +
  ylab("Confidence Interval median difference") +
  geom_vline(xintercept=0)+
  theme(legend.title = element_blank())+
  scale_y_log10()

ggsave("DifferenceMedianwithTextENCLAVES.pdf", height = 12, width = 14, useDingbats=FALSE)


library(patchwork)
(fst1 | fst3) / fst2
ggsave("DifferenceMedianLogScaleGLOBALandMajorLangFam_3plots_2.pdf", height = 10, width = 12, useDingbats=FALSE)



### 
# mark the mismatches in fst distribution

perpopRED$PopName[which(perpopRED$Mismatch3a=="MISMATCH")] 

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
# supplementary S10
#**************************************************
library(reshape) 
perpopREDfamily<-perpopRED[which(perpopRED$glottolog.NAME %in% MainFamilies),]

size<-as.numeric(table(perpopREDfamily$MainFamilies))
proportionSingleCases<-table(perpopREDfamily$MainFamilies,perpopREDfamily$listSingleCases)/size

proportionSingleMELT<-melt(proportionSingleCases)
colnames(proportionSingleMELT)<-c("Family", "variable", "proportion")

gi<-ggplot(proportionSingleMELT, aes(Family, proportion))
    
    ggplot(proportionSingleMELT, aes(fill=variable, y=proportion, x=Family)) + 
      geom_bar(position="fill", stat="identity", alpha=0.8)+
      scale_fill_manual(values = c("purple", "pink", "#FC7A1E", "#8E3B46","darkblue", "69b3a2", "gray20"), name="")+
    theme(axis.text.x = element_text(angle = 45, vjust=1, hjust=1)) +xlab("") +ylab("")
  ggsave("proportionSingleCasesEachFamily_2.pdf", useDingbats=FALSE, width = 7, height = 7)

########
## the map with singular mismatches of linguistic and genetic migrants, misaligned and drifted
# Figure 1 A
########
  
library(ggrepel)

#enclavesByMistake<-c("Yoruba"  ,"Nama", "Han-NChina" , "Bengali", "Evenk_FarEast") # these pops are not real enclaves, their linguistically unrelated pair is the one driving this genetic proximity effect

perpopRED2<-perpopRED[!is.na(perpopRED$lat),]
perpopREDSHIFTMAP<-MoveDataCoord(perpopRED2)
perpopREDSHIFTMAPinterest<-perpopREDSHIFTMAP[-which(perpopREDSHIFTMAP$listSingleCases=="ND"),]
#perpopREDSHIFTMAPmismatch<-perpopREDSHIFTMAPinterest[which(perpopREDSHIFTMAPinterest$singlepops%in%c("Misalligned","Enclave")),]
perpopREDSHIFTMAPnodata<-perpopREDSHIFTMAP[which(perpopREDSHIFTMAP$Mismatch3a=="ZeroSameFamilyNeighbors"),]
#perpopREDSHIFTMAPinterest<-perpopREDSHIFTMAPinterest[-which(perpopREDSHIFTMAPinterest$PopName%in% enclavesByMistake),]
#perpopREDSHIFTMAPmismatch<-perpopREDSHIFTMAPinterest[which(perpopREDSHIFTMAPinterest$Mismatch3a%in%c("MISMATCH","secondaryMISMATCH_FSTzeroDiffFamily")),]

### variant with language family color
# and symbols for single cases

base_plot + geom_point(data = perpopREDSHIFTMAP,
                       aes(x = lon.1, y = lat.1), 
                       color="black",shape=3,size=0.5)+
  theme(legend.position="bottom")+
  geom_point(data = perpopREDSHIFTMAP[which(perpopREDSHIFTMAP$MainFamilies%in%MainFamilies),],
             aes(x = lon.1, y = lat.1, color=MainFamilies), 
             size = 2, alpha=0.6)+
  geom_point(data=perpopREDSHIFTMAPinterest, aes(x = lon.1, y = lat.1, 
                                                 shape=listSingleCases),alpha=0.5, size=2.5)+
  # geom_point(data=perpopREDSHIFTMAPnodata, aes(x = lon.1, y = lat.1),
  #            color="ghostwhite", shape= 13,alpha=0.8)+
scale_shape_manual(values = c(11,6,5,8,7,0))+
  scale_color_manual(values=colorchoice)

ggsave("WholeGelatoMap_Pacific_LittleCrossesMajorFamilies_SingleCases_rightSymbols.pdf", useDingbats=FALSE, width = 15, height = 8, dpi = 300, units = "in")



# ***************************************************
#### SUPPLEMENTARY FIGURES 
# ***************************************************
# ***************************************************
# ***************************************************



# ***************************************************
## Figure  S3

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
  scale_color_gradient(low = "blue", high = "red",name="") +
  xlab("") + ylab("") 
 # theme(legend.title = element_blank())
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
  scale_color_gradient(low = "blue", high = "red",name="") +
  xlab("") + ylab("") 
ggsave("DensityDifferentLFamilyNeighbors_noContinentFilter_map_PacificCenter_1000km_2021.pdf", useDingbats=FALSE, width = 15, height = 8, dpi = 300, units = "in")


#******************************************
# control the distribution of divergence time for each macro region
# fig S5
#******************************************
perpopREDNe<-perpopRED[which(perpopRED$USEforNe_calculation=="YES"),]
possiblepopswithNE<-perpopREDNe$PopName
meltFstREDinfoYESne<-FstListREDinfo[which(FstListREDinfo$Pop1%in%possiblepopswithNE&FstListREDinfo$Pop2%in%possiblepopswithNE),]


meltFstGlotto_infowithinREgionTMRCA<-meltFstREDinfoYESne[which(meltFstREDinfoYESne$region1==meltFstREDinfoYESne$region2),]

## exclude drifted pops or the Fst averages will be higher than normal
DRIFTONI<-perpopRED[which(perpopRED$medianFSTregion>0.1&perpopRED$medianFST>0.1),]$PopName
#DRIFTONI<-perpopRED[which(perpopRED$averageFSTregion>0.1&perpopRED$averageFST>0.1),]$PopName

meltFstGlotto_infowithinREgionTMRCA<-meltFstGlotto_infowithinREgionTMRCA[-c(which(meltFstGlotto_infowithinREgionTMRCA$Pop2%in%DRIFTONI), which(meltFstGlotto_infowithinREgionTMRCA$Pop1%in%DRIFTONI)),]

ggplot(meltFstGlotto_infowithinREgionTMRCA, aes(region1, TMRCA_doubleNe))+
  geom_jitter(shape=16, position=position_jitter(0.2))+
  geom_violin(trim=FALSE, alpha=0.4)+
  stat_summary(fun.data="mean_sdl", mult=1, 
               geom="pointrange", color="orange")+
  ylim(0,70000)+xlab("")+ylab("population divergence time, years ago")

ggsave("distribTMRCA_ContinentsViolin.pdf", useDingbats=FALSE, width=6, height = 4)




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
betweenfamilSMALLgeoplottabileSHIFTMAP2<-betweenfamilSMALLgeoplottabileSHIFTMAP2[order(betweenfamilSMALLgeoplottabileSHIFTMAP2$percentileFST, decreasing = T),]

# only different L family pairs
betweenfamilSMALLgeoplottabileSHIFTMAP3<-betweenfamilSMALLgeoplottabileSHIFTMAP2[which(betweenfamilSMALLgeoplottabileSHIFTMAP2$FAMILY=="DIVERSE"),]


fig6a<-base_plot + geom_point(data = perpopREDSHIFTMAP,
                       aes(x = lon.1, y = lat.1),
                       color = "black",
                       size = 0.5, shape=3)+
  geom_path(data=betweenfamilSMALLgeoplottabileSHIFTMAP3, aes(x = lon.1, y = lat.1, 
                                                              group=index, alpha=percentileFST, size=percentileFST, color=percentileFST))+
  scale_alpha(range = c(0.2,0.1))+
  scale_size(range = c(0.5,.1))+
  scale_colour_gradient2(midpoint=0.1, low="darkred", mid="darkorange",
                         high="white", name="Percentile FST distribution")+
  #ggtitle("weight of mismatches according to FST percentile distribution")+
  xlab("") + ylab("") +
  guides(alpha=FALSE, size=FALSE)


### Supplementary Figure S6A
ggsave("MapPairMismatchesPercentileDistribution_until10percentile_2021.pdf", useDingbats=FALSE, width = 15, height = 8, dpi = 300, units = "in")


# ***************************************************
# sensitivity test with different geo thresholds

# Figure S2
### all the close FST distances on a map
base_plot + geom_point(data = perpopREDSHIFTMAP,
                       aes(x = lon.1, y = lat.1),
                       color = "black",
                       size = 0.5, shape=3)+
  geom_path(data=betweenfamilSMALLgeoplottabileSHIFTMAP2, aes(x = lon.1, y = lat.1, 
                                                             group=index, alpha=percentileFST, color=percentileFST))+
  scale_alpha(range = c(0.2,0.01))+
  scale_size(range = c(0.5,.01))+
  scale_colour_gradient2(midpoint=0.05, low="purple", mid="yellow",
                         high="white", name="Percentile FST distribution")+
  xlab("") + ylab("") +
guides(alpha=FALSE, size=FALSE)

ggsave("MapPairFSTPercentileDistribution_2021.pdf", useDingbats=FALSE, width = 15, height = 8, dpi = 300, units = "in")



##-------------------------------------
# figure  area proportion of linguistically unrelated in close FSTs
#  supplementary Fig S6B
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

fig6B<-ggplot(plotproportion,aes(x=as.numeric(percentiles), y=percentageMismatchesinTotalSmallFST, color=as.numeric(percentiles)))+
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
FstListREDinfoDIVERSE<-FstListREDinfo[which(FstListREDinfo$FAMILY=="DIVERSE"),]
for (i in length(percentiles):1){
  percTarget<-percentiles[i]
  temp<-FstListREDinfoDIVERSE[which(FstListREDinfoDIVERSE$percentileFST<=percTarget),]
  listoni<-  unique(temp$Pop1)
  perpopRED$percentileMismatch[which(perpopRED$PopName%in%listoni)]<-percTarget
  singlepopinmismatch[i]<-length(listoni)
}



# CONNECT MISMATCHES BETWEEN MAJOR FAMILIES
# with a big circle

##### FIGURE 1C  ****************************************

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
## Figure 1D, case studies
#******************************************

# subset case studies

casestudies<-c( "Hungarian1", "Maltese",  "Azeri_Azerbajan", "Kalmyk","Tunisian","Armenian" )
#casestudies<-c( "Wayku",   "Azeri_Azerbajan", "Kalmyk","Cusco2","Tunisian","Armenian" )
#casestudies<-c( "Wayku","Cusco2")
meltoni<-NA

for (i in 1:length(casestudies)){
  TARGET<-casestudies[i]
  meltini<-FstListREDinfo_noDuplicateNeighborsNoDrift[which(FstListREDinfo_noDuplicateNeighborsNoDrift$Pop1==TARGET),]
  meltiniSAME<-meltini[which(meltini$SameFamily=="YES"),]
  maxgeofam<-max(meltiniSAME$GEOdist, na.rm = T)
    meltini<-meltini[which(meltini$GEOdist<=maxgeofam),]
    meltini$TARGET<-TARGET
    meltoni<-rbind(meltoni, meltini)
}

meltoni<-meltoni[-1,]
meltoni$TARGET<-factor(meltoni$TARGET, levels=casestudies)

  ggg<-ggplot(meltoni,aes(y=Fst, x=SameFamily,color=SameFamily))
  ggg+ geom_jitter(shape=16, position=position_jitter(0.2))+
    geom_violin(trim=FALSE, alpha=0.4)+
    stat_summary(fun.data="mean_sdl", 
                 geom="pointrange", color="gold") +
     theme_light()+
    scale_y_log10()+
  scale_color_manual(values=c("#198B9F","#BA1E1E"))+
    facet_wrap(~TARGET,nrow = 2)
  ggsave("caseStudyMisalignedVIOLIN_YESorNO.pdf", height=5,width=5, useDingbats=FALSE)

  ggsave("caseStudyMisalignedVIOLIN_YESorNO2.pdf", height=5,width=5, useDingbats=FALSE)
  
#******************************************
#*#******************************************
#*#******************************************
#*#******************************************
#



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
ggsave("combinedMAP_FSTproportion_supple2021_2.pdf", useDingbats=FALSE, height = 15, width = 19)

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


# HUNGARIAN

HUNGlist<-c("Hungarian1", "Hungarian2")
hung<-FstListREDinfo[c(which(FstListREDinfo$Pop1%in%HUNGlist)),]
hung<-as.data.frame(hung)

selectionhung<-hung[order(hung$FstLinear),][1:30,]

colorispecial<-  MainFamilies2$COLOR[match(  names(table(selectionhung$family2)),
                                             MainFamilies2$MainFamilies)]
colorispecial[ is.na(colorispecial)]<-c("gray50") # special color code including minor language families in grayscale

agg<-ggplot(selectionhung)
HUNG<- agg+ geom_label_repel( aes(GEOdist,FstLinear, color=family2,label=Pop2), size=3, label.padding=0.1)+
  geom_point(aes(GEOdist,FstLinear, color=family2))+
  ggtitle("Hungarian")+
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

#***************************************************
## combine 4 figures in one single populations, without Language Isolates
# figure S12
#***************************************************

ggarrange(HUNG, MALTA, ARMENIAN, AZERBAJAN + rremove("x.text"), 
          # labels = c("A", "B", "C", "D"),
          ncol = 2, nrow = 2)
ggsave("combinedSinglePops_FstGeodist_supple_4targetpops.pdf", useDingbats=FALSE, height = 10, width = 13)





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
#excludepops<-union(ListEnclaves, Misalligned)
#coso<-union(excludepops, DRIFTONI)
#write.table(coso, "listaExclude.txt")
#union(ListEnclaves, DRIFTONI)
FstListREDinfo_noDuplicateNeighbors<-FstListREDinfo[-which(FstListREDinfo$GEOdist<10&FstListREDinfo$glottocodeBase1==FstListREDinfo$glottocodeBase2),]

FstListREDinfo_noDuplicateNeighborsNoDrift<-FstListREDinfo_noDuplicateNeighbors[-c
                                                                         (which(FstListREDinfo_noDuplicateNeighbors$Pop2%in%DRIFTONI), 
                                                                           which(FstListREDinfo_noDuplicateNeighbors$Pop1%in%DRIFTONI)),]

#*********************************************
### FIGURE 2
#*********************************************
#*
#*

### Within/between family comparisons
# script adapted from Damian
#
# FIGURE 2
#
library(tidyverse)

large_families<-MainFamilies[-13] # exclude Tupi

pops_for_test<-FstListREDinfo_noDuplicateNeighborsNoDrift
pops_for_test$listSingleCasesPOP1<-perpopRED$listSingleCases[match(pops_for_test$Pop1, perpopRED$PopName)]
pops_for_test$listSingleCasesPOP1[which(pops_for_test$listSingleCasesPOP1=="Match")]<-"ND"
pops_for_test$listSingleCasesPOP1[which(pops_for_test$listSingleCasesPOP1=="ND")]<-NA
pops_for_test$listSingleCasesPOP1[which(pops_for_test$FAMILY=="DIVERSE")]<-NA


FAM<-MainFamilies[-13] # exclude Tupi

for(f in FAM) {
  pops_for_test[[f]]<-(pops_for_test$family1==f)|(pops_for_test$family2==f)
  
  threshold_geo<-max(pops_for_test$GEOdist[pops_for_test$FAMILY==f],
                     na.rm = T)
     if(threshold_geo<500){
      threshold_geo<-500
    }
  
  print(f)
  print(threshold_geo)
  pops_for_test[[f]]<-pops_for_test[[f]]*(pops_for_test$GEOdist<=threshold_geo)
  pops_for_test[[f]]<-sapply(pops_for_test[[f]],function(x) ifelse(is.na(x),FALSE,x))
  
  }

family_distances<-pivot_longer(pops_for_test,
                               cols = FAM,
                               names_to = "Family_plot",
                               values_to = "Include") %>%
  filter(Include==TRUE)

####  add the within comparison symbols for listSingleCasesPOP1

family_distances %>%
  ggplot(aes(y=FstLinear,x=GEOdist))+
  geom_point(aes(fill=SameFamily, size=SameFamily),shape=21, alpha=0.3, stroke=0)+
  geom_point(aes(fill=SameFamily, shape=listSingleCasesPOP1),alpha=0.3,size=3, stroke=0.4)+
  geom_smooth(aes(group=SameFamily),se=F, color="white", size=2.5)+
  geom_smooth(aes(color=SameFamily),se=F)+
  geom_smooth(color="#FBB13C",linetype = "dotdash",alpha=0.4, size=1, method = lm)+
  facet_wrap(~Family_plot,
             scales="free",
             nrow=length(FAM)/4)+
  theme_minimal()+
  scale_y_sqrt()+
  scale_fill_manual(values=c("#7FDBEB","#E45B5B"))+
  scale_color_manual(values=c("#198B9F","#BA1E1E"))+
  scale_shape_manual(values = c(25,23,7,22))+
  scale_size_manual(values = c(1,2))+
  theme(legend.position = "bottom", axis.text=element_text(size=5))+
  xlab("Geographic Distance - km") +
  ylab("Linear FST") 

#ggsave("smoothregressionPlot.pdf", useDingbats=FALSE, height = 7, width = 10)
ggsave("smoothregressionPlotwithLinear_highlightMismatches_rightSymbols.pdf", useDingbats=FALSE, height = 8, width = 12)



# write.table(perpopRED, "/Users/chiarabarbieri/Library/Mobile Documents/com~apple~CloudDocs/GeLaTo/perPopMarchMami2.txt",  row.names = F, sep = "\t", quote = F)


#*********************************************
### FIGURE 3
#*********************************************
#*
#*

### Root time divergence from genetic data
# script adapted from Damian

# elaborated table where i mark the pairs of genetic populations that pass through the root of the language family phylogeny

time_depths<-read.table("meltFstREDinfoYESneONLYfamiliesHALF_v2.txt",header = T, sep="\t", as.is=T)

famtemp<-unlist(labels(table(time_depths$FAMILY)))

time_depths$root<-sapply(time_depths$passesThroughTheRoot, 
                         function(x) ifelse(x %in% c("no","lowlevel"),"no","yes"))

### keep only connections that pass through the root

time_depths<-time_depths[which(time_depths$root=="yes"),]

ListEnclaves<-perpopRED$PopName [which(perpopRED$listSingleCases==("Enclave"))]
EnclaveANDMisaligned<-perpopRED$PopName [which(perpopRED$listSingleCases==("EnclaveANDMisaligned"))]

time_depths$listSingleCasesPOP<-NA
time_depths$listSingleCasesPOP[which(time_depths$Pop1%in%ListEnclaves)]<-"Enclave"
time_depths$listSingleCasesPOP[which(time_depths$Pop2%in%ListEnclaves)]<-"Enclave"
time_depths$listSingleCasesPOP[which(time_depths$Pop1%in%ListMisaligned)]<-"Misaligned"
time_depths$listSingleCasesPOP[which(time_depths$Pop2%in%ListMisaligned)]<-"Misaligned"
time_depths$listSingleCasesPOP[which(time_depths$Pop1%in%EnclaveANDMisaligned)]<-"EnclaveANDMisaligned"
time_depths$listSingleCasesPOP[which(time_depths$Pop2%in%EnclaveANDMisaligned)]<-"EnclaveANDMisaligned"


time_depths %>% 
  ggplot(aes(y=TMRCA_doubleNe,x=FAMILY,color=FAMILY))+
  geom_jitter(position=position_jitterdodge(dodge.width=0.5,jitter.width=0),size=2,alpha=0.3)+
  coord_flip()+
  theme_minimal()+
  geom_point(aes( shape=listSingleCasesPOP),alpha=0.6,size=3, stroke=0.9)+  # different shape for pairs with an enclave
  scale_shape_manual(values = c(25,23,22))+
  scale_color_manual(values = MainFamilies2$COLOR[which(MainFamilies2$MainFamilies%in%famtemp)])+
  theme(legend.position = "bottom",
        axis.text = element_text(size=14,color="black"),
        axis.text.x = element_text(angle = 45, vjust=1, hjust=1, size = 10),
        panel.grid.minor = element_blank())+
  guides(color = "none")+
  
  scale_y_sqrt(breaks=c(100,500,1000,2000, 3000, 4000,5000,10000,20000,40000))+
  labs(x="",y="")+
  annotate(geom="segment",y = 2897,yend=6626,x=1,xend=1,
           size=10,alpha=0.3,color="gray55")+
  annotate(geom="segment",y = 2469,yend=4546,x=2,xend=2,
           size=10,alpha=0.3,color="gray55")+
  annotate(geom="segment",y = 3314,yend=8831,x=3,xend=3,
           size=10,alpha=0.3,color="gray55")+
  annotate(geom="segment",y = 3074,yend=7213,x=4,xend=4,
           size=10,alpha=0.3,color="gray55")+
  annotate(geom="segment",y = 1236,yend=2288,x=5,xend=5,
           size=10,alpha=0.3,color="gray55")+
  annotate(geom="segment",y = 2186,yend=3499,x=6,xend=6,
           size=10,alpha=0.3,color="gray55")+
  annotate(geom="segment",y = 2446,yend=4426,x=7,xend=7,
           size=10,alpha=0.3,color="gray55")+
  annotate(geom="segment",y = 1811,yend=2594,x=8,xend=8,
           size=10,alpha=0.3,color="gray55")+
  annotate(geom="segment",y = 2068,yend=3158,x=9,xend=9,
           size=10,alpha=0.3,color="gray55")+
  annotate(geom="segment",y = 2396,yend=4199,x=10,xend=10,
           size=10,alpha=0.3,color="gray55")+
  annotate(geom="segment",y = 2495,yend=4632,x=11,xend=11,
           size=10,alpha=0.3,color="gray55")

ggsave("TimeDepth_color_withEnclaves_rightSymbols.pdf", useDingbats=FALSE, height =7, width = 10, units = "in")
