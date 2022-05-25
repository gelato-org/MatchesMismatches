#******************************************
#******************************************
### Analysis Paper "A global analysis of matches and mismatches between human genetic and linguistic histories"
# GeLaTo database
# Chiara Barbieri April 2022

#******************************************
#*
#*
#*

### read the two main files
# list of 397 populations :
perpopRED<-read.table("PerpopRED_MaMi2022.txt", header = T, sep = "\t", as.is=T)

# list of pairwise comparisons :
FstListinfo<-read.table("FstListREDinfo_MaMi2022.txt", header=T, sep="\t")


# color palette for major language families
#
Lfamil<-table(perpopRED$glottolog.NAME)
MainFamilies<-unlist(labels(Lfamil[which(Lfamil>4)])) # minimum 5 populations per Lang Family

sizes<-as.numeric(unlist((Lfamil[which(Lfamil>4)])))
perpopREDfamily<-perpopRED[which(perpopRED$glottolog.NAME %in% MainFamilies),]
colorchoice<-c( "darkorange4" ,"#93032E"  ,   "#33A02C"    , "#A6D854"  ,   "#377EB8"  ,   "#E31A1C"   ,  "#FFD92F" ,    "#FF7F00"  ,   "#666666" ,   
                "cyan4"  ,     "#BC80BD"   ,  "#FED9A6" ,    "tan3" ,       "#6A3D9A" ,    "deeppink"   )
perpopRED$MainFamilies<-NA
perpopRED$MainFamilies[which(perpopRED$glottolog.NAME%in%MainFamilies)]<-perpopRED$glottolog.NAME[which(perpopRED$glottolog.NAME%in%MainFamilies)]
# #*********
# library(RColorBrewer)
# qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
# col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
# colorchoice<-sample(col_vector, length(MainFamilies))
# #*********

MainFamilies2<-as.data.frame(MainFamilies)
MainFamilies2$COLOR<-colorchoice
MainFamilies2$names<-paste0(MainFamilies, " (",sizes, ")" )

library(ggplot2)

# color palette blue and red
darkredYES<-"#BA1E1E"
blueNO<-"#198B9F"


#******************************************
### MAP
#******************************************

library(maps)
library('maps')
library('geosphere')
library(rworldmap)


# Script from Balthasar Bickel to plot maps Pacific Centered
source("ggworld.R") # map for Pacific Center projection, modified from Balthasar Bickel
#

perpopRED$lat<-as.numeric(perpopRED$lat)
perpopRED$lon<-as.numeric(perpopRED$lon)


perpopRED2<-perpopRED[!is.na(perpopRED$lat),] # exclude populations for which i do not have geographic coordinates

perpopREDSHIFTMAP<-MoveDataCoord(perpopRED2)   # perform coordinate shift to plot Pacific centered maps

## map population location and language families
#### MAP WITH MAJOR LANGUAGE FAMILIES ASSIGNED TO COLOR CODE,  
# Basic map
#******************************************
#*
 
  # make a  base_plot with contained latitudes 
  
base_plot <- ggplot() +
  geom_polygon(data = subset(world.sh.tr.df, lat > -60 & lat < 72),
               aes(x = X, y = Y, group = group),
               fill = base_fill,
               color = base_color,
               size = .1
  ) + theme_void()

  base_plot + geom_point(data = perpopREDSHIFTMAP,
                         aes(x = lon.1, y = lat.1), 
                         color="black",shape=3,size=0.5)+ 
  theme(legend.position="bottom", legend.text=element_text(size=8), legend.title = element_blank())+
  geom_point(data = perpopREDSHIFTMAP[which(perpopREDSHIFTMAP$glottolog.NAME%in%MainFamiliesSHORT),],
             aes(x = lon.1, y = lat.1, color=glottolog.NAME), 
             size = 2, alpha=0.5)+
  scale_color_manual(values=colorchoice3, labels = MainFamilies3$names)

ggsave("WholeGelatoMap_Pacific_LittleCrossesMajorFamilies_largerFig2_2022.pdf", useDingbats=FALSE, width = 8.5, height = 5, dpi = 300, units = "in")


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
perpopRED$SamefamilyClosestpop<-NA

# how far the same family is gen closer than other families?

for (i in 1:nrow(perpopRED)){
  targetpop<-perpopRED$PopName[i]
  temp<-FstListinfo[which(FstListinfo$Pop1==targetpop),]
  
  samefamily<-temp[which(temp$FAMILY!="DIVERSE"),]
  escludiniVicini<-which(samefamily$GEOdist<10&samefamily$glottocodeBase1==samefamily$glottocodeBase2)
  if(length(escludiniVicini)!=0){
    samefamily<-samefamily[-escludiniVicini,] # exclude when there is a neighbor too close (LESS THAN 10 KM) from same exact language as DUPLICATED SAMPLE
  } 
  if(nrow(samefamily)==0){
    perpopRED$hasOtherMatches[i]<-"NO" 
  }
  else{
    perpopRED$closeFstSameFamily[i]<-sort(samefamily$FST)[1]  # the closest Fst
    perpopRED$geodistSameFamily[i]<-samefamily$GEOdist[order(samefamily$FST)][1]    # the geographic distance from the closest Fst
    perpopRED$SamefamilyClosestpop[i]<-samefamily[order(samefamily$FST),][1,]$Pop2  # the pop which closest fst in match
    
  }
  DIFFfamily<-temp[which(temp$FAMILY=="DIVERSE"),]
  perpopRED$closeFstDIFFFamily[i]<-sort(DIFFfamily$FST)[1]  # the closest Fst from a different language family
  perpopRED$geodistDIFFFamily[i]<-DIFFfamily$GEOdist[order(DIFFfamily$FST)][1]    # the geographic distance from the closest Fst of a different language family
  perpopRED$DIFFFamilyClosestpop[i]<-DIFFfamily[order(DIFFfamily$FST),][1,]$Pop2  # the pop which closest fst in mismatch
  if(length(which(sort(samefamily$FST)<perpopRED$closeFstDIFFFamily[i]))==0){
    perpopRED$SameFamilyMostDistantClosestFst[i]<-"NONE"    }
  else{
    perpopRED$SameFamilyMostDistantClosestFst[i]<-tail(samefamily$GEOdist[which(sort(samefamily$FST)<perpopRED$closeFstDIFFFamily[i])],1) # the geographic distance from the closest Fst of the same language family that is less than the closest fst from diff L family
  }
}


perpopRED$proportionFST_diff_sameFamily<- perpopRED$closeFstDIFFFamily/perpopRED$closeFstSameFamily 
perpopRED$proportionGeoDistSameDiffFamily<- perpopRED$geodistDIFFFamily/perpopRED$geodistSameFamily 
perpopRED$EnclavesMismatch<-NA
perpopRED$EnclavesMismatch[which(perpopRED$closeFstDIFFFamily==0&perpopRED$closeFstSameFamil!=0)] <- "secondaryMISMATCH_FSTzeroDiffFamily"

perpopRED$EnclavesMismatch[which(perpopRED$proportionFST_diff_sameFamily<1&perpopRED$proportionGeoDistSameDiffFamily>1)] <- "MISMATCH"
perpopRED$EnclavesMismatch[which(perpopRED$proportionFST_diff_sameFamily>1&perpopRED$proportionGeoDistSameDiffFamily<1)] <- "MATCH"

perpopRED$EnclavesMismatch[which(perpopRED$hasOtherMatches!="YES"  )]<-"ZeroSameFamilyNeighbors"

table(perpopRED$EnclavesMismatch)/nrow(perpopRED)


### table to visualize how many cases

#************************
## april 2022 
MATCH                            MISMATCH secondaryMISMATCH_FSTzeroDiffFamily             ZeroSameFamilyNeighbors 
0.130982368                         0.070528967                         0.005037783                         0.055415617 
MATCH                            MISMATCH secondaryMISMATCH_FSTzeroDiffFamily             ZeroSameFamilyNeighbors 
52                                  28                                   2                                  22 
#************************


#Matches

ListMatches<-perpopRED$PopName [which(perpopRED$EnclavesMismatch=="MATCH")]

### List of very drifted populations to flag
DRIFTONI<-perpopRED[which(perpopRED$medianFSTregion>0.1&perpopRED$medianFST>0.1),]$PopName
> DRIFTONI
[1] "Baining_Marabu"      "Chukchi"             "Eskimo_Sireniki"     "Itelmen"             "Ju_hoan_North"      
[6] "Karitiana"           "Koryak"              "Lahu"                "Nanai"               "Nganasan"           
[11] "Nganasan_UstAvam"    "Nganasan_Volochanka" "Nivh"                "Onge"                "Rennell_and_Bellona"
[16] "San"                 "She"                 "Surui"       

## genetic and linguistic enclaves 

ListEnclaves<-perpopRED$PopName [which(perpopRED$EnclavesMismatch%in%c("MISMATCH","secondaryMISMATCH_FSTzeroDiffFamily"))]

> ListEnclaves
[1] "Yoruba"            "Mengen"            "Bengali"           "Hazara"           
[5] "Kharia"            "Gui"               "Khwe"              "Nama"             
[9] "Mongola"           "Avar_Gunibsky"     "Aleut"             "Dai"              
[13] "Jew_Georgian"      "Khomani"           "Spanish_PaisVasco" "Yaquis"           
[17] "Yukagir_Forest"    "Yukagir_Tundra"    "Zapotec"           "Wayku"            
[21] "Han-NChina"        "Evenk_FarEast"     "Cocama"            "Guarani"          
[25] "Guarani_GN"        "Karitiana"         "Surui"             "Azeri_Azerbajan"  
[29] "Hungarian1"        "Hungarian2"       


perpopRED$ListEnclaves<-NA
perpopRED$ListEnclaves[which(perpopRED$PopName%in%ListEnclaves)]<-"Enclave"


#*#******************************************
#*
#*  different distribution FST closest same family and different family
#* 

aa<-perpopRED[,c(1, which(colnames(perpopRED)%in%c("glottolog.NAME","closeFstSameFamily","closeFstDIFFFamily" )))] 


aa$difference<-aa$closeFstSameFamily-aa$closeFstDIFFFamily
length(which(aa$difference>0))/(nrow(perpopRED)-length(which(is.na(aa$difference)))) # only for the comparisons for which i have a same family FST

[1] 0.1973333  ## april 2022

# ~20 % of the pops who have another genetic population of the same language family do have closer FST with speakers of another language family


# *****************************************
# ### figure 1B, target enclaves
#   *****************************************
# # target match enclave Kalmyk BedouinB
# target genetic enclave Jew_Georgian
# target linguistic enclave Hungarian
library(reshape)

targetenclaves<-c("Jew_Georgian", "Hungarian1","Himba")

targetenclavesinfo<-perpopRED[which(perpopRED$PopName%in%targetenclaves),]
targetenclavesinfo<-targetenclavesinfo[,which(colnames(targetenclavesinfo)%in%
                                      c( "PopName", "glottolog.NAME","closeFstSameFamily","geodistSameFamily"   ,
                                       "SamefamilyClosestpop" ,"closeFstDIFFFamily" , "geodistDIFFFamily" ,"DIFFFamilyClosestpop", "DIFFFamilyClosestFAM"))]

targetenclavesinfo

targetenclavesinfo[2,1]<-"Hungarian"
targetenclavesinfo[3,1]<-"Jewish Georgian"


colnames(targetenclavesinfo)[3]<-"YES"
colnames(targetenclavesinfo)[5]<-"NO"

testt<-melt(targetenclavesinfo[,-c(4,6)])
coso<-melt(targetenclavesinfo[,c(1,4,6)])
testt$GeoDist<-paste0(coso$value, " km")
testt$GeoDistTrue<-coso$value
testt$closestPOP<-testt$SamefamilyClosestpop
testt$closestPOP[which(testt$variable=="NO")]<-testt$DIFFFamilyClosestpop[which(testt$variable=="NO")]
testt$value[which(testt$value==0)]<-0.0001  # or it does not show up in the plot
testt<-testt[1:6,]
# Balthasar version
ggplot(testt, aes(x=GeoDistTrue,y= value,  color=variable))+
  geom_point(size=4) +
  scale_color_manual("Same Family", values=c(darkredYES,blueNO), position="bottom")+
  xlab("Geographic distance (km)") +
  ylab("Fst") +
  theme_minimal()+
  facet_wrap(~PopName,scales="free")
ggsave("Fig1b_toyExamplesEnclaves_2022.pdf",height=2.5,width=8, useDingbats=FALSE)

 
#******************************************
#*#******************************************
#*#******************************************
#*#******************************************
#*#******************************************
#*#******************************************
#*    COMPARING MEDIAN FST BETWEEN AND WITHIN
# POPULATION HEURISTIC case 2
#* FST distribution between and within language families MISALIGNED
#*#******************************************
#*#******************************************


## exclude neighbors who speak the same language as DUPLICATE POPS
escludiniVicini<-FstListinfo[which(FstListinfo$GEOdist<10&FstListinfo$glottocodeBase1==FstListinfo$glottocodeBase2),]
> dim(escludiniVicini)
[1] 32 25
escludiniVicini[which(escludiniVicini$case=="single"),][,c(1:3,11:12)]
Pop2               Pop1         FST       family2        region2
73128     Kove_Tamuniai               Kove 0.000000000  Austronesian        OCEANIA
85508          Mamanwa1            Mamanwa 0.050776000  Austronesian SOUTHEAST_ASIA
115332              San      Ju_hoan_North 0.001266060           Kxa         AFRICA
130035     Sulka_Watwat        Sulka_Ganai 0.000000000         Sulka        OCEANIA
14937          BedouinB           BedouinA 0.017572100  Afro-Asiatic        EURASIA
42676         Greek_WGA       Greek_Athens 0.003147220 Indo-European        EURASIA
45019         GujaratiB          GujaratiA 0.001633070 Indo-European        EURASIA
45284         GujaratiC          GujaratiB 0.003168390 Indo-European        EURASIA
45342         GujaratiC          GujaratiA 0.002983290 Indo-European        EURASIA
46884          Han_HGDP         Han-NChina 0.002480980  Sino-Tibetan        EURASIA
49935        Hungarian2         Hungarian1 0.000000000        Uralic        EURASIA
78287   Lebanese_Muslim Lebanese_Christian 0.001025050  Afro-Asiatic        EURASIA
107168      Palestinian              Druze 0.008987530  Afro-Asiatic        EURASIA
107902     Peru_Quechua             Cusco2 0.000427919      Quechuan       AMERICAS
150014   Uzbek_Tashkent              Uzbek 0.000917217        Turkic        EURASIA
150910 Vietnamese_North               Kinh 0.001434990 Austroasiatic SOUTHEAST_ASIA


## exclude DRIFTONI 
## exclude drifted pops or the Fst averages will be higher than normal
DRIFTONI<-perpopRED[which(perpopRED$medianFSTregion>0.1&perpopRED$medianFST>0.1),]$PopName

FstListREDinfo_noDuplicateNeighbors<-FstListinfo[-which(FstListinfo$GEOdist<10&FstListinfo$glottocodeBase1==FstListinfo$glottocodeBase2),]
FstListREDinfo_noDuplicateNeighborsNoDrift<-FstListREDinfo_noDuplicateNeighbors[-c(which(FstListREDinfo_noDuplicateNeighbors$Pop2%in%DRIFTONI), which(FstListREDinfo_noDuplicateNeighbors$Pop1%in%DRIFTONI)),]


library(pairwiseCI)
perpopRED$MedDiff<-NA
perpopRED$MedDiffCIlower<-NA
perpopRED$MedDiffCIupper<-NA
perpopRED$MedianWithinSMALLER2<-NA
perpopRED$MedDiffCI<-NA

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
length(ListMisaligned)/length(which(!is.na(perpopRED$MedDiff)))
#*#******************************************
[1] 0.2044728 ## 20% of populations in misalignment
#*#******************************************
> length(ListMisaligned)
[1] 64


perpopRED$listSingleCases<-"ND"
perpopRED$listSingleCases[which(perpopRED$PopName%in%ListEnclaves)]<-"Enclave"
perpopRED$listSingleCases[which(perpopRED$PopName%in%ListMisaligned)]<-"Misaligned"
perpopRED$listSingleCases[which(perpopRED$PopName%in%intersect(ListMisaligned, ListEnclaves))]<-"EnclaveANDMisaligned"
perpopRED$listSingleCases[which(perpopRED$PopName%in%ListMatches)]<-"Match"
perpopRED$listSingleCases[which(perpopRED$PopName%in%intersect(ListMisaligned, ListMatches))]<-"MatchBUTMisaligned"
perpopRED$listSingleCases[which(perpopRED$PopName%in%DRIFTONI)]<-"Drifted"




#**************************************************
### table of cases mismatch enclaves
# Supplementary Table S2

EnclavesMismatch<-perpopRED[grep("MISMATCH",perpopRED$EnclavesMismatch),]
perpopRED$DIFFFamilyClosestFAM<-perpopRED$glottolog.NAME[match(perpopRED$DIFFFamilyClosestpop,perpopRED$PopName)]
colnameschoice<-c("PopName","glottolog.NAME", "geodistSameFamily","geodistDIFFFamily","DIFFFamilyClosestpop", "DIFFFamilyClosestFAM")
EnclavesMismatch<-EnclavesMismatch[,colnameschoice]
EnclavesMismatch$MisalignedMedian<-NA
EnclavesMismatch$MisalignedMedian<-perpopRED$MedDiff[match(EnclavesMismatch$PopName,perpopRED$PopName)]
EnclavesMismatch$MisalignedMedian<-EnclavesMismatch$MisalignedMedian<0
EnclavesMismatch$MisalignedLowerCI<-NA
EnclavesMismatch$MisalignedLowerCI<-perpopRED$MedDiffCIlower[match(EnclavesMismatch$PopName,perpopRED$PopName)]
EnclavesMismatch$MisalignedLowerCI<-EnclavesMismatch$MisalignedLowerCI<0
EnclavesMismatch$MisalignedUpperCI<-NA
EnclavesMismatch$MisalignedUpperCI<-perpopRED$MedDiffCIupper[match(EnclavesMismatch$PopName,perpopRED$PopName)]
EnclavesMismatch$MisalignedUpperCI<-EnclavesMismatch$MisalignedUpperCI<0

#  write.table(EnclavesMismatch[order(EnclavesMismatch$geodistDIFFFamily),],"tableSupplEnclaves2022.txt", sep="\t", row.names = F, quote=F) 



#******************************************
## figure distribution FST SAME OR DIFFERENT FAMILY for each population
## Figure 1C, case studies
#******************************************

# subset case studies
# FstListREDinfo_noDuplicateNeighbors<-FstListinfo[-which(FstListinfo$GEOdist<10&FstListinfo$glottocodeBase1==FstListinfo$glottocodeBase2),]
# FstListREDinfo_noDuplicateNeighborsNoDrift<-FstListREDinfo_noDuplicateNeighbors[-c(which(FstListREDinfo_noDuplicateNeighbors$Pop2%in%DRIFTONI), which(FstListREDinfo_noDuplicateNeighbors$Pop1%in%DRIFTONI)),]

casestudies<-c( "Kalmyk", "Azeri_Azerbajan")

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

ggg<-ggplot(meltoni,aes(y=FST, x=SameFamily,color=SameFamily))
ggg+ geom_jitter(shape=16, position=position_jitter(0.2))+
  geom_violin(trim=FALSE, alpha=0.4)+
  stat_summary(fun.y= median, 
               geom="pointrange", color="gold") +
  theme_minimal()+
  scale_y_log10()+
  scale_color_manual(values=c(blueNO,darkredYES))+
  guides(color=FALSE)+
  facet_wrap(~TARGET,nrow = 1)
# scales="free")
ggsave("Fig1C_caseStudyMisalignedVIOLIN_YESorNO_2022.pdf", height=3,width=5.5, useDingbats=FALSE)


### 
#*#******************************************
#*  #*#******************************************
#*    #*#******************************************
#*    
#*    


######################################################################
########################################################
##############
# language family comparisons
# geographic distance and smooth or linear regression
## testing if SAME FAMILY has an effect on Fst given GeoDist
############################
############################

# 
# FstListREDinfo_noDuplicateNeighbors<-FstListinfo[-which(FstListinfo$GEOdist<10&FstListinfo$glottocodeBase1==FstListinfo$glottocodeBase2),]
# 
# FstListREDinfo_noDuplicateNeighborsNoDrift<-FstListREDinfo_noDuplicateNeighbors[-c
#                                                                                 (which(FstListREDinfo_noDuplicateNeighbors$Pop2%in%DRIFTONI), 
#                                                                                   which(FstListREDinfo_noDuplicateNeighbors$Pop1%in%DRIFTONI)),]

#*********************************************
### FIGURE 2
#*********************************************
#*
#*

### Within/between family comparisons
# original script from Damián Blasi, adapted by C. Barbieri
#
# FIGURE 2
#
library(tidyverse)

#large_families<-MainFamilies[-13] # exclude Tupi
large_families<-c(  "Atlantic-Congo" ,    "Indo-European", "Sino-Tibetan"    , "Mongolic-Khitan",
                  "Afro-Asiatic", "Turkic"   ,   "Austronesian" ,           "Uralic") # families above 8 pops


pops_for_test<-FstListREDinfo_noDuplicateNeighborsNoDrift
pops_for_test$FAMILY<-factor(pops_for_test$FAMILY, levels=large_families) # force the order when doing facet wrap in the figure

pops_for_test$listSingleCasesPOP1<-perpopRED$listSingleCases[match(pops_for_test$Pop1, perpopRED$PopName)]
pops_for_test$listSingleCasesPOP1[which(pops_for_test$listSingleCasesPOP1=="Match")]<-"ND"
pops_for_test$listSingleCasesPOP1[which(pops_for_test$listSingleCasesPOP1=="MatchBUTMisaligned")]<-"ND"
pops_for_test$listSingleCasesPOP1[which(pops_for_test$listSingleCasesPOP1=="ND")]<-NA
pops_for_test$listSingleCasesPOP1[which(pops_for_test$SameFamily=="NO")]<-NA  # do not want to plot mismatches in different lang family comparisons


FAM<-  large_families  


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

####  add the within comparison symbols for listSingleCasesPOP1, flagging mismatches
family_distances$Family_plot<-factor(family_distances$Family_plot, levels=large_families) # force the order when doing facet wrap in the figure

family_distances %>%
  
  ggplot(aes(y=FstLinear,x=GEOdist))+
  geom_point(aes(fill=SameFamily, size=SameFamily),shape=21, alpha=0.3, stroke=0)+
  geom_point(aes(fill=SameFamily, shape=listSingleCasesPOP1),alpha=0.3,size=3, stroke=0.3)+
  geom_smooth(aes(group=SameFamily),se=F, color="white", size=2.5)+
  geom_smooth(aes(color=SameFamily),se=F)+
  geom_smooth(color="#FBB13C",linetype = "dotdash",alpha=0.4, size=1, method = lm)+
  facet_wrap(~Family_plot,
             scales="free",
             nrow=length(FAM)/4)+
  theme_minimal()+
  scale_y_sqrt()+
  scale_fill_manual(values=c("#7FDBEB","#E45B5B"))+
  scale_color_manual(values=c(blueNO,darkredYES))+
  scale_shape_manual(" ",values = c(25,23,22))+
  scale_size_manual(values = c(1,2))+
  theme(legend.position = "bottom", axis.text=element_text(size=7),strip.text = element_text(size = 12))+
  xlab("Geographic Distance - km") +
  ylab("Linear FST") 

ggsave("Fig2_smoothregressionPlotwithLinear_highlightMismatches_2022.pdf", useDingbats=FALSE, height = 6, width = 10)

ggsave("Fig2_smoothregressionPlotwithLinear_highlightMismatches_2022.png",  height = 5, width = 10)




# ***************************************************
#### SUPPLEMENTARY FIGURES 
# ***************************************************
# ***************************************************
# ***************************************************



########################################################
# Figure S8 supplementary - comparison of difference of medians and CI, colored per language family
########################################################
MainFamilies<-unlist(labels(Lfamil[which(Lfamil>4)])) # minimum 5 populations per Lang Family

perpopREDfamily<-perpopRED[which(perpopRED$glottolog.NAME %in% MainFamilies),]
perpopREDfamilynoDrift<-perpopREDfamily[-which(perpopREDfamily$listSingleCases=="Drifted"),]

perpopREDEnclaves<-perpopRED[which(perpopRED$EnclavesMismatch%in%c("MATCH", "MISMATCH","secondaryMISMATCH_FSTzeroDiffFamily")),]

perpopREDEnclaves<-perpopRED[which(perpopRED$EnclavesMismatch%in%c("MATCH", "MISMATCH","secondaryMISMATCH_FSTzeroDiffFamily")),]
perpopRED$enclavesAgain<-NA
perpopRED$enclavesAgain[which(perpopRED$EnclavesMismatch=="MATCH")]<-"Match"
perpopRED$enclavesAgain[which(perpopRED$EnclavesMismatch=="MISMATCH")]<-"GeneticEnclave"
perpopRED$enclavesAgain[which(perpopRED$EnclavesMismatch=="secondaryMISMATCH_FSTzeroDiffFamily")]<-"LinguisticEnclave"

perpopREDnoDrift<-perpopRED[-which(perpopRED$listSingleCases=="Drifted"),]
library(ggrepel)

# panel A, color code language families all together 
fst1<-ggplot(perpopREDnoDrift,aes(MedDiff,MedDiffCI, color=MainFamilies))+
  geom_point(alpha=0.5, size=3,show.legend = FALSE)+
  #geom_text_repel(aes(label=PopName, color=MainFamilies), size=0.5)+
  scale_color_manual(values = colorchoice)+
  theme_light() +
xlab("Difference Median FST between - within") +
  ylab("Confidence Interval median difference") +
  geom_vline(xintercept=0)+
  scale_y_log10()


#ggsave("DifferenceMedianLogScale.pdf", height = 12, width = 14, useDingbats=FALSE)

# panel C, separate language families facet wrap

fst2<-ggplot(perpopREDfamilynoDrift ,aes(MedDiff,MedDiffCI, color=MainFamilies))+
  geom_point(alpha=0.5, size=3,show.legend = FALSE)+
  geom_text_repel(aes(label=PopName, color=MainFamilies), size=0.5, show.legend = FALSE)+
  scale_color_manual(values = colorchoice)+
  theme_light()+ xlab("Difference Median FST between - within") +
  ylab("Confidence Interval median difference") +
  geom_vline(xintercept=0)+
  scale_y_log10()+
  facet_wrap(~MainFamilies,ncol=6)
#ggsave("DifferenceMedianLogScaleMajorLangFam.pdf", height = 12, width = 14, useDingbats=FALSE)

# panel B, color code previous enclaves and matches 

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

#ggsave("DifferenceMedianwithTextENCLAVES.pdf", height = 12, width = 14, useDingbats=FALSE)


library(patchwork)
patchwork<-(fst1 | fst3) / fst2
patchwork + plot_annotation(tag_levels = 'A')
ggsave("DifferenceMedianLogScaleGLOBALandMajorLangFam_3plots_2022.pdf", height = 10, width = 12, useDingbats=FALSE)
ggsave("DifferenceMedianLogScaleGLOBALandMajorLangFam_3plots_2022.png", height = 10, width = 12, dpi = 150)


########
## the map with singular mismatches of linguistic and genetic migrants, misaligned and drifted
# Figure S10
########
  
library(ggrepel)


perpopRED2<-perpopRED[!is.na(perpopRED$lat),]
perpopREDSHIFTMAP<-MoveDataCoord(perpopRED2)
perpopREDSHIFTMAPinterest<-perpopREDSHIFTMAP[-which(perpopREDSHIFTMAP$listSingleCases=="ND"),]
#perpopREDSHIFTMAPmismatch<-perpopREDSHIFTMAPinterest[which(perpopREDSHIFTMAPinterest$singlepops%in%c("Misalligned","Enclave")),]
perpopREDSHIFTMAPnodata<-perpopREDSHIFTMAP[which(perpopREDSHIFTMAP$EnclavesMismatch=="ZeroSameFamilyNeighbors"),]
#perpopREDSHIFTMAPinterest<-perpopREDSHIFTMAPinterest[-which(perpopREDSHIFTMAPinterest$PopName%in% enclavesByMistake),]
#perpopREDSHIFTMAPmismatch<-perpopREDSHIFTMAPinterest[which(perpopREDSHIFTMAPinterest$EnclavesMismatch%in%c("MISMATCH","secondaryMISMATCH_FSTzeroDiffFamily")),]

### variant with language family color
# and symbols for single cases

base_plot + geom_point(data = perpopREDSHIFTMAP,
                       aes(x = lon.1, y = lat.1), 
                       color="black",shape=3,size=0.5)+
  theme(legend.position="bottom",legend.title = element_blank())+
  geom_point(data = perpopREDSHIFTMAP[which(perpopREDSHIFTMAP$MainFamilies%in%MainFamilies),],
             aes(x = lon.1, y = lat.1, color=MainFamilies), 
             size = 2, alpha=0.6)+
  geom_point(data=perpopREDSHIFTMAPinterest, aes(x = lon.1, y = lat.1, 
                                                 shape=listSingleCases),alpha=0.5, size=2.5)+
  # geom_point(data=perpopREDSHIFTMAPnodata, aes(x = lon.1, y = lat.1),
  #            color="ghostwhite", shape= 13,alpha=0.8)+
scale_shape_manual(values = c(11,6,5,8,7,0))+
  scale_color_manual(values=colorchoice, labels = MainFamilies2$names)

ggsave("Fig_s10_WholeGelatoMap_Pacific_LittleCrossesMajorFamilies_SingleCases_rightSymbols_2022.pdf", useDingbats=FALSE, width = 14, height = 8, dpi = 300, units = "in")
ggsave("Fig_s10_WholeGelatoMap_Pacific_LittleCrossesMajorFamilies_SingleCases_rightSymbols_2022.png",  width = 14, height = 8, dpi = 300, units = "in")



# ***************************************************
## Figure  S3

library(plyr)

threshold<-c(500,1000,1500,2000,2500,3000)
blockneighborthresholdGeo<-c()
for (i in 1:length(threshold)){
  testino2<-ddply(FstListinfo, "Pop1", function(x) length(which(x$GEOdist<threshold[i]&x$FAMILY=="DIVERSE"))) # no continental filter

  colnames(testino2)<-c("PopName","GeoCloseUnrelated1")
  testinomerged<-testino2
  testinomerged$GeoCloseUnrelated<-testinomerged$GeoCloseUnrelated1/2
  testinomerged$threshold<-threshold[i]
  
  blockneighborthresholdGeo<-rbind(blockneighborthresholdGeo,testinomerged)
}


blockneighborthresholdGeoINFO<-merge(perpopRED,blockneighborthresholdGeo)



# frequencies at least one neighbor from diff family
freqOneneighbor<-c()
for (i in 1:length(threshold)){
  thresholdtemp<-threshold[i]
popsNumberofNeighbors<-ddply(FstListinfo, "Pop1", function(x) length(which(x$GEOdist<thresholdtemp&x$FAMILY=="DIVERSE"))) # no continental filter
freqOneneighbor[i]<-length(which(popsNumberofNeighbors$V1>0))/nrow(perpopRED)
}
valuestop<-round(freqOneneighbor*100)
[1]  57  85  93  98  99 100  # add manually on the violin plot figure
valuestop<-paste0(valuestop, "%")

# violin plots with number of populations from a different language family
# fig S3B


 blockneighborthresholdGeoINFO$threshold<-as.character(blockneighborthresholdGeoINFO$threshold)
blockneighborthresholdGeoINFO$threshold[which(blockneighborthresholdGeoINFO$threshold=="500")]<-" 500"
ga<-ggplot(blockneighborthresholdGeoINFO,aes(x=threshold,y=GeoCloseUnrelated1))
S3B <-ga+ geom_jitter(shape=16, position=position_jitter(0.2))+
  geom_violin(trim=FALSE, alpha=0.4)+
  stat_summary(fun.data="mean_sdl", 
               geom="pointrange", color="orange")+
  theme_minimal()+
  theme(axis.title.y = element_text(size = 8)) +
  ylim(0,150)+xlab("Geographic radius in km")+
ylab("Number of populations from a different language family") +
annotate("text", x = 1:6, y = 150, label = valuestop)


# ***************************************************
# evaluate the incidence of neighbors from a different L famiy for each population

# count how many populations from different L families in a radius of 1000 km in the same continent

popsNumberofNeighbors<-ddply(FstListinfo, "Pop1", function(x) length(which(x$GEOdist<1000&x$FAMILY=="DIVERSE"))) # no continental filter
perpopRED$nNeighborsDiffFamily<-popsNumberofNeighbors
colnames(popsNumberofNeighbors)<-c("PopName","GeoCloseUnrelated")
perpopRED3<-merge(perpopRED,popsNumberofNeighbors)
perpopRED4<-perpopRED3[-which(perpopRED3$GeoCloseUnrelated==0),]
perpopRED4_SHIFTMAP<-MoveDataCoord(perpopRED4)

# plot on a map the density of unrelated pop geographically close
# fig S3A


S3A <- base_plot + geom_point(data = perpopREDSHIFTMAP,
                       aes(x = lon.1, y = lat.1),
                       color = "black",
                       size = 0.5, shape=3)+
  geom_point(data=perpopRED4_SHIFTMAP,  aes(x = lon.1, y = lat.1,color=as.numeric(GeoCloseUnrelated)), alpha=0.5, size=2)+
  # ggtitle("number of neighbors from a different language family")+
  scale_color_gradient(low = "blue", high = "red",name="") +
  xlab("") + ylab("") 
 # theme(legend.title = element_blank())
# ggsave("DensityPopulationsWithDifferentLFamilyNeighbors_noContinentFilter_map_PacificCenter_1000km_2021_.pdf", useDingbats=FALSE, width = 15, height = 8, dpi = 300, units = "in")

# count how many LANGUAGE FAMILIES  in a radius of 1000 km in the same continent

popsNumberofNeighborFamilies<-ddply(FstListinfo, "Pop1", function(x) length(unique(x[which(x$GEOdist<1000&x$FAMILY=="DIVERSE"),]$family2))) # no continental filter
perpopRED$nDiffFamilies1000km<-popsNumberofNeighbors

colnames(popsNumberofNeighborFamilies)<-c("PopName","GeoCloseUnrelatedFamilies")
perpopRED3<-merge(perpopRED,popsNumberofNeighborFamilies)
perpopRED4<-perpopRED3[-which(perpopRED3$GeoCloseUnrelatedFamilies==0),]
perpopRED4_SHIFTMAP<-MoveDataCoord(perpopRED4)

# plot on a map the density of unrelated Language Families geographically close
# fig S3C


S3C <- base_plot + geom_point(data = perpopREDSHIFTMAP,
                       aes(x = lon.1, y = lat.1),
                       color = "black",
                       size = 0.5, shape=3)+
  geom_point(data=perpopRED4_SHIFTMAP,  aes(x = lon.1, y = lat.1,color=as.numeric(GeoCloseUnrelatedFamilies)), alpha=0.5, size=2)+
  # ggtitle("number of neighbors from a different language family")+
  scale_color_gradient(low = "blue", high = "red",name="") +
  xlab("") + ylab("") 
# ggsave("DensityDifferentLFamilyNeighbors_noContinentFilter_map_PacificCenter_1000km_2021.pdf", useDingbats=FALSE, width = 15, height = 8, dpi = 300, units = "in")

## assemble
library(patchwork)

layout <- "
AAAA
AAAA
#BB#
CCCC
CCCC
"
S3A + S3B + S3C + 
  plot_layout(design = layout)+ plot_annotation(tag_levels = 'A')


ggsave("Fig_s03_evaluateIncidenceLanguageFamilies2022.pdf", height = 12, width = 10, useDingbats=FALSE)
ggsave("Fig_s03_evaluateIncidenceLanguageFamilies2022.png", height = 12, width = 10, dpi = 150)


#******************************************
# control the distribution of divergence time for each macro region
# fig S5
#******************************************
perpopREDNe<-perpopRED[which(perpopRED$USEforNe_calculation=="YES"),]
possiblepopswithNE<-perpopREDNe$PopName
meltFstREDinfoYESne<-FstListinfo[which(FstListinfo$Pop1%in%possiblepopswithNE),]
meltFstREDinfoYESne<-meltFstREDinfoYESne[which(meltFstREDinfoYESne$Pop2%in%possiblepopswithNE),]


meltFstGlotto_infowithinREgionTMRCA<-meltFstREDinfoYESne[which(meltFstREDinfoYESne$region1==meltFstREDinfoYESne$region2),]
meltFstGlotto_infowithinREgionTMRCA<-meltFstGlotto_infowithinREgionTMRCA[which(meltFstGlotto_infowithinREgionTMRCA$case=="single" ),]


## exclude drifted pops or the Fst averages will be higher than normal
#DRIFTONI<-perpopRED[which(perpopRED$medianFSTregion>0.1&perpopRED$medianFST>0.1),]$PopName

meltFstGlotto_infowithinREgionTMRCA<-meltFstGlotto_infowithinREgionTMRCA[-c(which(meltFstGlotto_infowithinREgionTMRCA$Pop2%in%DRIFTONI)),]
meltFstGlotto_infowithinREgionTMRCA<-meltFstGlotto_infowithinREgionTMRCA[-c(which(meltFstGlotto_infowithinREgionTMRCA$Pop1%in%DRIFTONI)),]
meltFstGlotto_infowithinREgionTMRCA$TMRCA_doubleNe<-as.numeric(meltFstGlotto_infowithinREgionTMRCA$TMRCA_doubleNe)

ggplot(meltFstGlotto_infowithinREgionTMRCA, aes(region1, TMRCA_doubleNe))+
  geom_jitter(shape=16, position=position_jitter(0.2))+
  geom_violin(trim=FALSE, alpha=0.4)+
  stat_summary(fun.data="mean_sdl", 
               geom="pointrange", color="orange")+
  ylim(0,70000)+xlab("")+ylab("population divergence time, years ago")+
  theme_light() 
  

ggsave("FigS05_distribTMRCA_ContinentsViolin2022.pdf", useDingbats=FALSE, width=6, height = 4)




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
FstListGlotto_infowithinREgion<-FstListinfo[which(FstListinfo$region1==FstListinfo$region2),]

## exclude drifted pops or the Fst averages will be higher than normal
#DRIFTONI<-perpopRED[which(perpopRED$medianFSTregion>0.1&perpopRED$medianFST>0.1),]$PopName

FstListGlotto_infowithinREgionNoDrif<-FstListGlotto_infowithinREgion[-c(which(FstListGlotto_infowithinREgion$Pop2%in%DRIFTONI), which(FstListGlotto_infowithinREgion$Pop1%in%DRIFTONI)),]
FstListGlottoIBD_infoNoDrift<-FstListinfo[-c(which(FstListinfo$Pop2%in%DRIFTONI), which(FstListinfo$Pop1%in%DRIFTONI)),]

# set a series of percentiles up to 0.5
percentiles<-seq(from = .002, to = .5, by = .002)
#percentilesCoarse<-seq(from = .005, to = .5, by = .005)

# calculate thresholds per continents

library(plyr)
medianGEO<-ddply(FstListGlotto_infowithinREgionNoDrif, "region1", function(x) median(x$FST))
valuespercentileRegions<-data.frame(row.names =medianGEO$region1)
for (i in 1:length(percentiles)){
  percTarget<-percentiles[i]
  temp<-ddply(FstListGlotto_infowithinREgionNoDrif, "region1", function(x) quantile(x$FST, percTarget))
  valuespercentileRegions[,i]<-temp[,2]
  #colnames(valuespercentileRegions[i])<-colnames(temp)[2]
}
colnames(valuespercentileRegions)<-percentiles
medianGEO<-cbind(medianGEO,valuespercentileRegions)
colnames(medianGEO)[2]<-"medianRegion"

aggiuntaWorld<-c()  # the list of percentile threshold on global distribution

for (i in 1:length(percentiles)){
  percTarget<-percentiles[i]
  aggiuntaWorld[i]<-quantile(FstListGlottoIBD_infoNoDrift$FST, percTarget)
}
aggiuntaWorldline<-c("WORLD", median(FstListGlottoIBD_infoNoDrift$FST),aggiuntaWorld) 
medianGEO<-rbind(medianGEO,aggiuntaWorldline)


# assign the lowest percentile to each pair fst

FstListinfo$percentileFST<-NA

for (i in 1:nrow(FstListinfo)){
  if(FstListinfo$region1[i]!=FstListinfo$region2[i] ){
    FstListinfo$percentileFST[i]<- percentiles[which(FstListinfo$FST[i]<aggiuntaWorld)][1]
  }
  else {
    reftemp<-  valuespercentileRegions[which( medianGEO$region1== FstListinfo$region2[i]),] # list of percentiles corresponding to the continent of the two populations
    FstListinfo$percentileFST[i]<- percentiles[which(FstListinfo$FST[i]<reftemp)][1]
    
  }
}

# prepare plottable file for pairwise connections

betweenfamilSMALLgeoplottabile1<-FstListinfo
betweenfamilSMALLgeoplottabile1$index<-c(1:nrow(FstListinfo))
betweenfamilSMALLgeoplottabile1<-betweenfamilSMALLgeoplottabile1[,-which(colnames(betweenfamilSMALLgeoplottabile1)%in%c("lat2","lon2"))]

betweenfamilSMALLgeoplottabile2<-FstListinfo
betweenfamilSMALLgeoplottabile2$index<-c(1:nrow(FstListinfo))
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
ggsave("MapPairMismatchesPercentileDistribution_until10percentile_2022_.pdf", useDingbats=FALSE, width = 15, height = 8, dpi = 300, units = "in")


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

ggsave("MapPairFSTPercentileDistribution_2022.pdf", useDingbats=FALSE, width = 15, height = 8, dpi = 300, units = "in")



##-------------------------------------
# figure  area proportion of linguistically unrelated in close FSTs
#  supplementary Fig S6B
##-------------------------------------



#FstListREDinfo_diverse<-FstListREDinfo[which(FstListREDinfo$FAMILY=="DIVERSE"),]
FstListGlotto_infowithinREgion<-FstListinfo[which(FstListinfo$region1==FstListinfo$region2),]
FstListREDinfo_BetweenRegion<-FstListinfo[which(FstListinfo$region1!=FstListinfo$region2),]



numberTotalSmallFST<-c()
numberTotalSmallFSTDIVERSE<-c()
percentageMismatchesinTotalSmallFST<-c()

meltFstGlottoIBD_infoWithPercentageMismatch<-data.frame(row.names = FstListinfo$popslistemp)

for (j in 1:length(percentiles)){
  mismatchlist<-c()
  percTarget<-percentiles[j]
  for (i in 1:5){  # the five macro continents
    mismatchlist1<-FstListGlotto_infowithinREgion[which(FstListGlotto_infowithinREgion$region1==medianGEO[i,1]),] # select continent
    thresholdtemp<-as.numeric(medianGEO[i,which(colnames(medianGEO)==percTarget)]) # select percentile value fst
    mismatchlist1<-mismatchlist1[which(mismatchlist1$FST<= thresholdtemp),]
    mismatchlist<-rbind(mismatchlist,mismatchlist1)
  }
  thresholdtempworld<-as.numeric(medianGEO[6,which(colnames(medianGEO)==percTarget)])
  mismatchlistBEtweenregion<-FstListREDinfo_BetweenRegion[which(FstListREDinfo_BetweenRegion$FST<= thresholdtempworld),] # use the target quantile all over the world Fst WORLD
  TheMismatches<-rbind(mismatchlist, mismatchlistBEtweenregion)
  numberTotalSmallFST[j]<-nrow(TheMismatches) # number of pairs within FST percentile
  numberTotalSmallFSTDIVERSE[j]<-length(which(TheMismatches$FAMILY=="DIVERSE"))  # number of pairs in mismatch
  
  percentageMismatchesinTotalSmallFST[j]<-length(which(TheMismatches$FAMILY=="DIVERSE"))/nrow(TheMismatches)
  meltFstGlottoIBD_infoWithPercentageMismatch[,j]<-NA
  meltFstGlottoIBD_infoWithPercentageMismatch[TheMismatches$popslistemp,j]<-percTarget
}
colnames(meltFstGlottoIBD_infoWithPercentageMismatch)<-percentiles

listperpercentileNumberCases<-rbind(numberTotalSmallFST,numberTotalSmallFSTDIVERSE,percentageMismatchesinTotalSmallFST)
colnames(listperpercentileNumberCases)<-percentiles

# count the cases below each percentile
f<-function(x, output){sort(x)[1]}
FstListinfo$percentileFSTdoublecheck<-apply(meltFstGlottoIBD_infoWithPercentageMismatch,
                                               1, f)

plot(FstListinfo$percentileFSTdoublecheck, FstListinfo$percentileFST) # check if the two slightly different method give the same percentiles

plotproportion<-data.frame(percentiles,nrow(FstListinfo)-numberTotalSmallFST,numberTotalSmallFST-numberTotalSmallFSTDIVERSE,numberTotalSmallFSTDIVERSE)
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

ggsave("ProportionaMismatchesSensitivityThresholdFST_2022_better.pdf", useDingbats=FALSE, width = 4, height = 2)



### for each population, i annotate the percentile FST threshold in which they appear in a mismatch pair
perpopRED$percentileMismatch<-NA
singlepopinmismatch<-c()
FstListREDinfoDIVERSE<-FstListinfo[which(FstListinfo$FAMILY=="DIVERSE"),]
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

meltFstGlottoIBD_infoLessPercentile2<-FstListinfo[which(FstListinfo$percentileFST<0.1),] # the closest 0.1 percentile


meltFstREDinfoMAJORFAM<-meltFstGlottoIBD_infoLessPercentile2[which(meltFstGlottoIBD_infoLessPercentile2$FAMILY=="DIVERSE"),]
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
meltFstREDinfo222<-FstListinfo[which(FstListinfo$family1%in%MainFamilies),]
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

colfunc <- colorRampPalette(c("darkred", "darkorange"))
listcolors<-colfunc(length((listvalues)))
listvalues<-cbind(listvalues,listcolors)

matFSTpercentilemeltcolor<- listvalues[,2][match(matFSTpercentilemelt[,3],listvalues[,1])]

matpercentilecolor<-matrix(matFSTpercentilemeltcolor,
                           nrow(matAdjustSizeROUND),
                           ncol(matAdjustSizeROUND))


pdf("circlize_Mismatches_families_adjustSampleSize_below10percentile2022.pdf", width = 8, height = 8)
ff<-chordDiagramFromMatrix(matAdjustSizeROUND, directional = 0, 
                           transparency = 0.1, symmetric=T, order = orderspecial, 
                           col=matpercentilecolor,grid.col = grid.col)
dev.off()




#******************************************
#*#******************************************
#*#******************************************
#*#******************************************
#
  ##
  # Figure S7A map 
  ########

  perpopREDSHIFTMAP_Mismatched<-perpopREDSHIFTMAP[which(perpopREDSHIFTMAP$MedDiff<0),] # fst between is smaller than the fst between (mismatch)
  #perpopREDSHIFTMAP_3wilcox<-perpopREDSHIFTMAP[which(perpopREDSHIFTMAP$WtestPvalueGEOfilter<=0.01&perpopREDSHIFTMAP$WtestMeanPvalueGEOfilter_resample <=0.01),] # significant
perpopREDSHIFTMAP$MedDiff<-as.numeric(perpopREDSHIFTMAP$MedDiff)
 S7A<- base_plot + geom_point(data = perpopREDSHIFTMAP,
                         aes(x = lon.1, y = lat.1),
                         color = "black",
                         size = 0.5, shape=3)+
    geom_point(data = perpopREDSHIFTMAP,
               aes(x = lon.1, y = lat.1,color = MedDiff, alpha=0.5-MedDiffCI),
               size = 2)+
    scale_color_gradient2( low="yellow",mid="gold",
                           high="darkblue",  midpoint = 0,  name="Difference FST between/within")+
    # scale_color_gradient2( low="darkred",mid="pink",
    #                        high="blue", midpoint = 0.00, na.value = "grey50", name="")+
    geom_text(data=perpopREDSHIFTMAP_Mismatched, aes(x = lon.1, y = lat.1,label=PopName), size=0.8)+
    # geom_point(data = perpopREDSHIFTMAP_Mismatched, aes(x = lon.1, y = lat.1), color = "darkred",  size = 1.2, alpha=0.5)+
    xlab("") + ylab("") +
    theme(legend.position="bottom")+
    guides(alpha=FALSE, size=FALSE)
  
   ggsave("MapDifferentFSTmedianBetweenWithin_newColors2022.pdf", useDingbats=FALSE, width = 12, height = 8, dpi = 300, units = "in")
  
  


#********************************************
#********************************************
#********************************************
# ***************************************************
# ANALYSIS OF REGIONAL CASE STUDIES
# ***************************************************
   base_fill = 'gray70'
   base_color = 'black'
   textsize<-4
   ThresholdMedDiff<-0.0
   
   
   # *********************************************
   ##### EUROPE #####
   # *********************************************
   
coordinates<-c(35,56,-9,27)
perpopREDEUROPE<-perpopRED[which(perpopRED$lat>coordinates[1]&perpopRED$lat<coordinates[2]&perpopRED$lon>coordinates[3]&perpopRED$lon<coordinates[4]),]
map.world <- map_data(map="world")
library(ggrepel)

gg <- ggplot()
gg <- gg + theme_light()
gg <- gg + geom_map(data=map.world, map=map.world, aes(map_id=region), fill=base_fill, colour="black", size=0.15)

gg<- gg+coord_quickmap(ylim=coordinates[1:2], xlim=coordinates[3:4])

diffFST_EUR<- gg + geom_point(data=perpopREDEUROPE, 
                          aes(x=lon, y=lat, color=MedDiff,na.rm=T, size=MedDiffCI,
                              alpha=MedDiffCI) ) +
  scale_color_gradient2( low="yellow",mid="gold",
                         high="darkblue",  midpoint = 0, na.value = "grey50", name="Difference FST between/within")+
  geom_label_repel(data=perpopREDEUROPE[which(perpopREDEUROPE$MedDiff<ThresholdMedDiff),], 
                   aes(x=lon, y=lat,label=PopName, color=MedDiff, na.rm= TRUE ), size=textsize, label.padding=0.1)+
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
gg <- gg + theme_light()
gg <- gg + geom_map(data=map.world, map=map.world, aes(map_id=region), fill=base_fill, colour="black", size=0.15)
gg<- gg+coord_quickmap(ylim=coordinates[1:2], xlim=coordinates[3:4])

diffFST_AFRICA<- gg + geom_point(data=perpopREDAFRICA, 
                              aes(x=lon, y=lat, color=MedDiff,na.rm=T,size=MedDiffCI,
                                  alpha=MedDiffCI) )  +
  scale_alpha(range = c(0.8,0.2),name="Difference FST CI")+
  scale_size(range = c(7,2),name="Difference FST CI")+
  scale_color_gradient2( low="yellow",mid="gold",
                         high="darkblue",  midpoint = 0, na.value = "grey50", name="Difference FST between/within")+
  geom_label_repel(data=perpopREDAFRICA[which(perpopREDAFRICA$MedDiff<ThresholdMedDiff),], aes(x=lon, y=lat,label=PopName, color=MedDiff, na.rm= TRUE ), size=textsize, label.padding=0.1)+
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
gg <- gg + theme_light()
gg <- gg + geom_map(data=map.world, map=map.world, aes(map_id=region), fill=base_fill, colour="black", size=0.15)
gg<- gg+coord_quickmap(ylim=coordinates[1:2], xlim=coordinates[3:4])

diffFST_OCEANIA<- gg + geom_point(data=perpopREDOCEANIA, 
                                  aes(x=lon, y=lat, color=MedDiff,na.rm=T,size=MedDiffCI,
                                      alpha=MedDiffCI) )  +
  scale_alpha(range = c(0.8,0.2),name="Difference FST CI")+
  scale_size(range = c(7,2),name="Difference FST CI")+
  scale_color_gradient2( low="yellow",mid="gold",
                         high="darkblue",  midpoint = 0, na.value = "grey50", name="Difference FST between/within")+
  geom_label_repel(data=perpopREDOCEANIA[which(perpopREDOCEANIA$MedDiff<ThresholdMedDiff),], aes(x=lon, y=lat,label=PopName, color=MedDiff, na.rm= TRUE ), size=textsize, label.padding=0.1)+
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
gg <- gg + theme_light()
gg <- gg + geom_map(data=map.world, map=map.world, aes(map_id=region), fill=base_fill, colour="black", size=0.15)
gg<- gg+coord_quickmap(ylim=coordinates[1:2], xlim=coordinates[3:4])


diffFST_CAUCASUS<- gg + geom_point(data=perpopREDcaucasus, 
                                   aes(x=lon, y=lat, color=MedDiff,na.rm=T,size=MedDiffCI,
                                       alpha=MedDiffCI) )  +
  scale_alpha(range = c(0.8,0.2),name="Difference FST CI")+
  scale_size(range = c(7,2),name="Difference FST CI")+
  scale_color_gradient2( low="yellow",mid="gold",
                         high="darkblue",  midpoint = 0, na.value = "grey50", name="Difference FST between/within")+
  geom_label_repel(data=perpopREDcaucasus[which(perpopREDcaucasus$MedDiff<ThresholdMedDiff),], aes(x=lon, y=lat,label=PopName, color=MedDiff, na.rm= TRUE ), size=textsize, label.padding=0.1)+
  geom_point(data=perpopREDcaucasus, aes(x=lon, y=lat, shape=glottolog.NAME), size=1 ) +
  labs(shape="Language Family", colour="Difference between FST distribution between and within")+
  scale_shape_manual(values=1:7)+
  theme(axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank(),
        legend.title = element_text(size = 8), legend.text = element_text(size = 8),
        legend.key.size = unit(0.4, 'cm'))




# *********************************************
# main figure combined
# 1: europe, 2: oceania, 3: africa, 4: caucasus
# *********************************************


# Figure S7
library(patchwork)

patchwork<-S7A / ( diffFST_AFRICA | diffFST_EUR ) / ( diffFST_CAUCASUS | diffFST_OCEANIA)
patchwork + plot_annotation(tag_levels = 'A')
ggsave("Fig_S07_combinedMAP_FSTproportion2022.pdf", height = 18, width = 15, useDingbats=FALSE)
ggsave("Fig_S07_combinedMAP_FSTproportion2022.png", height = 22, width = 18)

#***************************************************
### SINGLE POPS SPECIAL CASES
#***************************************************


# HUNGARIAN

HUNGlist<-c("Hungarian1", "Hungarian2")
hung<-FstListinfo[c(which(FstListinfo$Pop1%in%HUNGlist)),]
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


malta<-FstListinfo[c(which(FstListinfo$Pop1=="Maltese")),]

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
## ARMENIAN

armenian<-FstListinfo[c(which(FstListinfo$Pop1=="Armenian")),]

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

azerbajan<-FstListinfo[c(which(FstListinfo$Pop1=="Azeri_Azerbajan")),]

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
## combine 4 figures in one single populations, without Language Isolates
# figure  S9
#***************************************************
library(ggpubr)

ggarrange(HUNG, MALTA, ARMENIAN, AZERBAJAN + rremove("x.text"), 
          # labels = c("A", "B", "C", "D"),
          ncol = 2, nrow = 2)
ggsave("combinedSinglePops_FstGeodist_supple_4targetpops_rightColor_rightGeo2022.pdf", useDingbats=FALSE, height = 10, width = 13)




#*********************************************
### FIGURE 3
#*********************************************
#*
#*

### Root time divergence from genetic data
# original script from Damián Blasi, adapted by Barbieri

# elaborated table where i mark the pairs of genetic populations that pass through the root of the language family phylogeny
library(plyr)

library(tidyverse)

MainFamilies<-unlist(labels(Lfamil[which(Lfamil>5)])) # minimum 10 populations per Lang Family

perpopREDNe<-perpopRED[which(perpopRED$USEforNe_calculation=="YES"),] # a subset of populations for which i can calculate effective population size
possiblepopswithNE<-perpopREDNe$PopName

meltFstREDinfoYESne<-FstListREDinfo_noDuplicateNeighborsNoDrift[which(FstListREDinfo_noDuplicateNeighborsNoDrift$Pop1%in%possiblepopswithNE
                                                                      &FstListREDinfo_noDuplicateNeighborsNoDrift$Pop2%in%possiblepopswithNE),]

meltFstREDinfoYESneSameFamily<-meltFstREDinfoYESne[which(meltFstREDinfoYESne$FAMILY!="DIVERSE"),] # comparisons within families only
meltFstREDinfoYESneSameFamilyHALF<-meltFstREDinfoYESneSameFamily[which(meltFstREDinfoYESneSameFamily$case=="single"),] # the original file is double side matrix

ListEnclaves<-perpopRED$PopName [which(perpopRED$listSingleCases==("Enclave"))]
EnclaveANDMisaligned<-perpopRED$PopName [which(perpopRED$listSingleCases==("EnclaveANDMisaligned"))]
ListMisaligned<-perpopRED$PopName[which(perpopRED$MedDiff<0)]
ListAligned<-perpopRED$PopName[which(perpopRED$MedDiff>0)]


# flag comparisons with mismatches
meltFstREDinfoYESneSameFamilyHALF$listSingleCasesPOP<-"no mismatches"
meltFstREDinfoYESneSameFamilyHALF$listSingleCasesPOP[which(meltFstREDinfoYESneSameFamilyHALF$Pop1%in%ListEnclaves)]<-"Enclave"
meltFstREDinfoYESneSameFamilyHALF$listSingleCasesPOP[which(meltFstREDinfoYESneSameFamilyHALF$Pop2%in%ListEnclaves)]<-"Enclave"
meltFstREDinfoYESneSameFamilyHALF$listSingleCasesPOP[which(meltFstREDinfoYESneSameFamilyHALF$Pop1%in%ListMisaligned)]<-"Misaligned"
meltFstREDinfoYESneSameFamilyHALF$listSingleCasesPOP[which(meltFstREDinfoYESneSameFamilyHALF$Pop2%in%ListMisaligned)]<-"Misaligned"
meltFstREDinfoYESneSameFamilyHALF$listSingleCasesPOP[which(meltFstREDinfoYESneSameFamilyHALF$Pop1%in%EnclaveANDMisaligned)]<-"EnclaveANDMisaligned"
meltFstREDinfoYESneSameFamilyHALF$listSingleCasesPOP[which(meltFstREDinfoYESneSameFamilyHALF$Pop2%in%EnclaveANDMisaligned)]<-"EnclaveANDMisaligned"

time_depths<-meltFstREDinfoYESneSameFamilyHALF
time_depths<-time_depths[which(time_depths$root=="yes"),] ### keep only connections that pass through the root (which i manually screened)

time_depths<-time_depths[which(time_depths$FAMILY%in%MainFamilies),]
famtemp<-unlist(labels(table(time_depths$FAMILY)))[which(table(time_depths$FAMILY)>5)]
famtempnames<-c("Afro-Asiatic"  ,  "Atlantic-Congo (Bantu)"   ,   "Austronesian","Indo-European", "Khoe", "Daghestanian" , "Quechuan" ,
                "Turkic", "Uralic")

 

time_depths<-time_depths[which(time_depths$FAMILY%in%famtemp),]  ### keep only Main Families, min n comparison > 3

dim(time_depths)
[1] 230  35

dodg <- 0.2 # to scale the time reconstruction from historical linguistics below Generalized Bayesian Dating
typedosh<-13
coloretemp<-"gray40"
sizetemp<-4


time_depths %>% 
  ggplot(aes(y=TMRCA_doubleNe,x=FAMILY,shape=listSingleCasesPOP,  color=FAMILY))+
  #ggplot(time_depths, aes(y=TMRCA_doubleNe,x=FAMILY,shape=listSingleCasesPOP))+
  # geom_jitter(alpha=0.7, width = 0.15, size=4)+
  geom_point(alpha=0.7, width = 0.15, size=4, color=darkredYES)+
  coord_flip()+
  theme_minimal()+
 # geom_point(aes( shape=listSingleCasesPOP),alpha=0.6,size=3, stroke=0.9)+  # different shape for pairs with an enclave
  scale_shape_manual("", values = c(25,23,22, 19))+

 # scale_color_manual(values = MainFamilies2$COLOR[which(MainFamilies2$MainFamilies%in%famtemp)])+
  theme(legend.position = "bottom",
        axis.text = element_text(size=14,color="black"),
        # axis.text.x = element_text(angle = 45, vjust=1, hjust=1, size = 10),
         axis.text.x = element_text( size = 8),
           panel.grid.minor = element_blank())+
  scale_x_discrete(limits=rev, labels=rev(famtempnames))+
  scale_y_sqrt(breaks=c(100,500,1000,2000, 3000, 4000,5000,10000,20000,40000))+
  labs(x="",y="")+ guides(color = FALSE, size = FALSE) +
  #manually annotate time from Generalized Bayesian Dating 
  annotate(geom="segment",y = 2495,yend=4632,x=1,xend=1,
           size=10,alpha=0.3,color="gray55")+  # Uralic
  annotate(geom="segment",y = 2396,yend=4199,x=2,xend=2,
           size=10,alpha=0.3,color="gray55")+
  annotate(geom="segment",y = 1811,yend=2594,x=3,xend=3,
           size=10,alpha=0.3,color="gray55")+
  annotate(geom="segment",y = 2446,yend=4426,x=4,xend=4,
           size=10,alpha=0.3,color="gray55")+
  annotate(geom="segment",y = 1236,yend=2288,x=5,xend=5,
           size=10,alpha=0.3,color="gray55")+
  annotate(geom="segment",y = 3074,yend=7213,x=6,xend=6,
           size=10,alpha=0.3,color="gray55")+
  annotate(geom="segment",y = 3314,yend=8831,x=7,xend=7,
           size=10,alpha=0.3,color="gray55")+
  annotate(geom="segment",y = 2469,yend=4546,x=8,xend=8,
           size=10,alpha=0.3,color="gray55")+
  annotate(geom="segment",y = 2897,yend=6626,x=9,xend=9,
           size=10,alpha=0.3,color="gray55") + # Afro-Asiatic
#manually annotate time from Historical linguistic evidence
annotate(geom="segment", linetype=typedosh, y = 4000,yend=5000,x=1-dodg,xend=1-dodg,
         size=sizetemp,alpha=0.8,color=coloretemp)+ #Uralic
  annotate(geom="segment", linetype=typedosh, y = 6000,yend=7000,x=1-dodg,xend=1-dodg,
           size=sizetemp,alpha=0.8,color=coloretemp)+ # uralic again
  annotate(geom="segment", linetype=typedosh, y = 2000,yend=2500,x=2-dodg,xend=2-dodg,
           size=sizetemp,alpha=0.8,color=coloretemp)+ # Turkic
  annotate(geom="segment", linetype=typedosh, y = 1000,yend=2000,x=3-dodg,xend=3-dodg,
           size=sizetemp,alpha=0.8,color=coloretemp)+ # Quechua
  annotate(geom="segment", linetype=typedosh, y = 2000,yend=4000,x=5-dodg,xend=5-dodg,
           size=sizetemp,alpha=0.8,color=coloretemp)+ # Khoe
  annotate(geom="segment", linetype=typedosh, y = 5500,yend=6500,x=6-dodg,xend=6-dodg,
           size=sizetemp,alpha=0.8,color=coloretemp)+ # IE
  annotate(geom="segment", linetype=typedosh, y = 7500,yend=8500,x=6-dodg,xend=6-dodg,
           size=sizetemp,alpha=0.8,color=coloretemp)+ # IE again
  annotate(geom="segment", linetype=typedosh, y = 5000,yend=6000,x=7-dodg,xend=7-dodg,
           size=sizetemp,alpha=0.8,color=coloretemp)+ # Austronesian
    annotate(geom="segment",linetype=typedosh,y = 4000,yend=5000,x=8-dodg,xend=8-dodg,
           size=sizetemp,alpha=0.8,color=coloretemp)+ # Bantu
  annotate(geom="segment",linetype=typedosh,y = 4000,yend=11000,x=9-dodg,xend=9-dodg,
           size=sizetemp,alpha=0.8,color=coloretemp) # Afro-Asiatic


ggsave("Fig3_TimeDepth_red_withEnclaves_reducedSet_withLegend_noJitter2022.pdf", useDingbats=FALSE, height =12, width = 18, units = "cm")
ggsave("Fig3_TimeDepth_red_withEnclaves_reducedSet_withLegend_noJitter2022.png",  height =6, width = 9, units = "in")


# ---------------------------------------------------------
#### the distribution of divergence time for all language families, plus drifted 
#### Figure supplementary S11
# jitter
# ---------------------------------------------------------
# 
# perpopREDNe<-perpopRED[which(perpopRED$USEforNe_calculation=="YES"),]
# possiblepopswithNE<-perpopREDNe$PopName
# 
# meltFstREDinfoYESne<-FstListREDinfo[which(FstListREDinfo$Pop1%in%possiblepopswithNE&FstListREDinfo$Pop2%in%possiblepopswithNE),]

# included Drifted populations
meltFstREDinfoYESne<-FstListREDinfo_noDuplicateNeighbors[which(FstListREDinfo_noDuplicateNeighbors$Pop1%in%possiblepopswithNE
                                                                      &FstListREDinfo_noDuplicateNeighbors$Pop2%in%possiblepopswithNE),]
meltFstREDinfoYESneSameFamily<-meltFstREDinfoYESne[which(meltFstREDinfoYESne$FAMILY!="DIVERSE"),] # comparisons within families only
meltFstREDinfoYESneSameFamilyHALF<-meltFstREDinfoYESneSameFamily[which(meltFstREDinfoYESneSameFamily$case=="single"),] # the original file is double side matrix

ListEnclaves<-perpopRED$PopName [which(perpopRED$listSingleCases==("Enclave"))]
EnclaveANDMisaligned<-perpopRED$PopName [which(perpopRED$listSingleCases==("EnclaveANDMisaligned"))]
ListMisaligned<-perpopRED$PopName[which(perpopRED$MedDiff<0)]


# flag comparisons with mismatches
meltFstREDinfoYESneSameFamilyHALF$listSingleCasesPOP<-"no mismatches"
meltFstREDinfoYESneSameFamilyHALF$listSingleCasesPOP[which(meltFstREDinfoYESneSameFamilyHALF$Pop1%in%ListEnclaves)]<-"Enclave"
meltFstREDinfoYESneSameFamilyHALF$listSingleCasesPOP[which(meltFstREDinfoYESneSameFamilyHALF$Pop2%in%ListEnclaves)]<-"Enclave"
meltFstREDinfoYESneSameFamilyHALF$listSingleCasesPOP[which(meltFstREDinfoYESneSameFamilyHALF$Pop1%in%ListMisaligned)]<-"Misaligned"
meltFstREDinfoYESneSameFamilyHALF$listSingleCasesPOP[which(meltFstREDinfoYESneSameFamilyHALF$Pop2%in%ListMisaligned)]<-"Misaligned"
meltFstREDinfoYESneSameFamilyHALF$listSingleCasesPOP[which(meltFstREDinfoYESneSameFamilyHALF$Pop1%in%EnclaveANDMisaligned)]<-"EnclaveANDMisaligned"
meltFstREDinfoYESneSameFamilyHALF$listSingleCasesPOP[which(meltFstREDinfoYESneSameFamilyHALF$Pop2%in%EnclaveANDMisaligned)]<-"EnclaveANDMisaligned"
meltFstREDinfoYESneSameFamilyHALF$listSingleCasesPOP[which(meltFstREDinfoYESneSameFamilyHALF$Pop1%in%DRIFTONI)]<-"Drifted"
meltFstREDinfoYESneSameFamilyHALF$listSingleCasesPOP[which(meltFstREDinfoYESneSameFamilyHALF$Pop2%in%DRIFTONI)]<-"Drifted"

time_depths<-meltFstREDinfoYESneSameFamilyHALF
time_depths<-time_depths[which(time_depths$root=="yes"),] ### keep only connections that pass through the root (which i manually screened)

time_depths<-as.data.frame(time_depths)

famtemp<-labels(which(table(time_depths$FAMILY)>0))
famtempnames2<-c("Afro-Asiatic"  ,  "Atlantic-Congo (Bantu)"   ,   "Austronesian","Chukotko-Kamchatkan", "Eskimo-Aleut","Indo-European", "Khoe", 
                 "Kx'a" ,"Daghestanian" , "Quechuan" ,"Tungusic"   ,         "Tupian"   ,
                "Turkic", "Uralic")

COLOR<-rep("Gray50",length(levels(time_depths$FAMILY)) )
MainFamilies3<-cbind(famtemp,COLOR)
colnames(MainFamilies3)[1]<-"MainFamilies"
MainFamilies3<-as.data.frame(MainFamilies3)
MainFamilies3$COLOR<-MainFamilies2$COLOR[ match(MainFamilies3$MainFamilies,MainFamilies2$MainFamilies)]
MainFamilies3$COLOR[is.na(MainFamilies3$COLOR)]<-"Gray50"
rownames(MainFamilies3)<-as.character(MainFamilies3$MainFamilies)
MainFamilies3<-MainFamilies3[famtemp,]

time_depths %>% 
  # ggplot(aes(y=TMRCA_doubleNe,x=FAMILY,color=FAMILY))+
  ggplot(aes(y=TMRCA_doubleNe,x=FAMILY,shape=listSingleCasesPOP, color=FAMILY))+
  geom_jitter(alpha=0.7, width = 0.2)+
  coord_flip()+
  theme_minimal()+
  #geom_point(aes( shape=listSingleCasesPOP),alpha=0.6,size=3, stroke=0.9)+  # different shape for pairs with an enclave
  scale_shape_manual(values = c(11,25,23,22,19))+
  scale_color_manual(values = MainFamilies3$COLOR)+
  theme(legend.position = "bottom",legend.title = element_blank(),
        axis.text = element_text(size=14,color="black"),
        axis.text.x = element_text(angle = 45, vjust=1, hjust=1, size = 10),
        panel.grid.minor = element_blank())+
  guides(color = "none")+
  scale_x_discrete(limits=rev, labels=rev(famtempnames2))+
  scale_y_sqrt(breaks=c(100,500,1000,2000, 3000, 4000,5000,10000,20000,40000))+
  labs(x="",y="")

ggsave("FigS11_TimeDepth_color_withEnclavesandDrifted_rightSymbols_allFamilies2022.png",  height =6, width = 8, units = "in")
ggsave("FigS11_TimeDepth_color_withEnclavesandDrifted_rightSymbols_allFamilies2022.pdf", useDingbats=FALSE, height =6, width = 8, units = "in")


#*************************************************
### RESAMPLING SENSITIVITY TEST
# Figure S11

Lfamil<-table(perpopRED$glottolog.NAME)
MainFamiliesSHORT<-unlist(labels(Lfamil[which(Lfamil>7)])) # minimum 8 populations per Lang Family

MINIMUMPOP<-8


# exclude Drifted
# perpopREDnoDrift<-perpopRED[-which(perpopRED$listSingleCases=="Drifted"),]
# LfamilnoDrift<-table(perpopREDnoDrift$glottolog.NAME)[which(table(perpopREDnoDrift$glottolog.NAME)>7)]
# MainFamiliesnoDrift<-unlist(labels(LfamilnoDrift[which(LfamilnoDrift>7)]))

# all language families
LfamiliesTot<-table(perpopRED$glottolog.NAME)
LfamiliesTotNAMES<-unlist(labels(LfamiliesTot))


#****************************
## SUBSAMPLE the language families to a  maximum n size = 8
#****************************

nRepetitions<-100
EnclavesResampling<-matrix(NA, nRepetitions, 3)
enclaveTAG<-c("MATCH",  "MISMATCH" , "ZeroSameFamilyNeighbors")
colnames(EnclavesResampling)<-enclaveTAG

differenceSUB<-NA # percentage of the pops who have another genetic population of the same language family do have closer FST with speakers of another language family

ProportionMisaligned<- NA
ProportionAligned<- NA

for (k in 1:nRepetitions){
  
  subsetpops<-c()
  
  
  for (i in 1:length(LfamiliesTotNAMES)){
    TARGETFAM<-LfamiliesTotNAMES[i]
    if (LfamiliesTot[i]<8){
      subsetpopsADD<-perpopRED$PopName[which(perpopRED$glottolog.NAME==TARGETFAM)]
    } else {
      subsetpopsADD<-perpopRED$PopName[sample(which(perpopRED$glottolog.NAME==TARGETFAM),MINIMUMPOP)]
    }
    subsetpops<-c(subsetpops,subsetpopsADD)
  }
  
  perpopSUBSAMPLE<-perpopRED[which(perpopRED$PopName%in%subsetpops),]
  
  # LfamilSUB<-table(perpopSUBSAMPLE$glottolog.NAME)
  
  FstLisinfoSUBSAMPLE<-FstListinfo[which(FstListinfo$Pop1%in%subsetpops&FstListinfo$Pop2%in%subsetpops),]
  SUBDRIFTONI<-perpopSUBSAMPLE[which(perpopSUBSAMPLE$medianFSTregion>0.1&perpopSUBSAMPLE$medianFST>0.1),]$PopName
  
  #***********************
  #* HERE the session to count enclaves
  #*
  ## enclaves
  #***********************
  #*
  
  
  SUBhasOtherMatches<-NA
  SUBcloseFstSameFamily<-NA  
  SUBgeodistSameFamily<-NA  
  SUBcloseFstDIFFFamily<-NA  
  SUBgeodistDIFFFamily<-NA  
  
  # how far the same family is gen closer than other families?
  
  for (i in 1:nrow(perpopSUBSAMPLE)){
    targetpop<-perpopSUBSAMPLE$PopName[i]
    temp<-FstLisinfoSUBSAMPLE[which(FstLisinfoSUBSAMPLE$Pop1==targetpop),]
    
    samefamily<-temp[which(temp$FAMILY!="DIVERSE"),]
    escludiniVicini<-which(samefamily$GEOdist<10&samefamily$glottocodeBase1==samefamily$glottocodeBase2)
    if(length(escludiniVicini)!=0){
      samefamily<-samefamily[-escludiniVicini,] # exclude when there is a neighbor too close (LESS THAN 10 KM) from same exact language as DUPLICATED SAMPLE
    } 
    if(nrow(samefamily)==0){
      SUBhasOtherMatches[i]<-"NO" 
    }
    else{
      SUBcloseFstSameFamily[i]<-sort(samefamily$FST)[1]  # the closest Fst
      SUBgeodistSameFamily[i]<-samefamily$GEOdist[order(samefamily$FST)][1]    # the geographic distance from the closest Fst
      perpopRED$SamefamilyClosestpop[i]<-samefamily[order(samefamily$FST),][1,]$Pop2  # the pop which closest fst in match
      
    }
    DIFFfamily<-temp[which(temp$FAMILY=="DIVERSE"),]
    SUBcloseFstDIFFFamily[i]<-sort(DIFFfamily$FST)[1]  # the closest Fst from a different language family
    SUBgeodistDIFFFamily[i]<-DIFFfamily$GEOdist[order(DIFFfamily$FST)][1]    # the geographic distance from the closest Fst of a different language family
  }
  
  
  SUBproportionFST_diff_sameFamily<- SUBcloseFstDIFFFamily/SUBcloseFstSameFamily
  SUBproportionGeoDistSameDiffFamily<- SUBgeodistDIFFFamily/SUBgeodistSameFamily 
  SUBEnclavesMismatch<-NA
  #SUBEnclavesMismatch[which(SUBcloseFstDIFFFamily==0&SUBcloseFstSameFamily!=0)] <- "secondaryMISMATCH_FSTzeroDiffFamily"
  
  SUBEnclavesMismatch[which(SUBproportionFST_diff_sameFamily<1&SUBproportionGeoDistSameDiffFamily>1)] <- "MISMATCH"
  SUBEnclavesMismatch[which(SUBproportionFST_diff_sameFamily>1&SUBproportionGeoDistSameDiffFamily<1)] <- "MATCH"
  
  SUBEnclavesMismatch[which(SUBhasOtherMatches=="NO"  )]<-"ZeroSameFamilyNeighbors" ## note that this is constant as all of them are in families with less than 8 samples
  
  EnclavesResampling[k,]<-(table(SUBEnclavesMismatch)/nrow(perpopSUBSAMPLE))[enclaveTAG]
  
  
  
  
  
  difference<-SUBcloseFstSameFamily-SUBcloseFstDIFFFamily
  differenceSUB[k]<-length(which(difference>0))/(nrow(perpopSUBSAMPLE)-length(which(is.na(difference)))) # only for the comparisons for which i have a same family FST
  
  
  #***********************
  #* HERE the session to count MISALIGNED
  
  #***********************
  #*
  
  SUBDRIFTONI<-perpopSUBSAMPLE[which(perpopSUBSAMPLE$medianFSTregion>0.1&perpopSUBSAMPLE$medianFST>0.1),]$PopName
  
  SUBFstListREDinfo_noDuplicateNeighbors<-FstLisinfoSUBSAMPLE[-which(FstLisinfoSUBSAMPLE$GEOdist<10&FstLisinfoSUBSAMPLE$glottocodeBase1==FstLisinfoSUBSAMPLE$glottocodeBase2),]
  SUBFstListREDinfo_noDuplicateNeighborsNoDrift<-SUBFstListREDinfo_noDuplicateNeighbors[-c(which(SUBFstListREDinfo_noDuplicateNeighbors$Pop2%in%SUBDRIFTONI), which(SUBFstListREDinfo_noDuplicateNeighbors$Pop1%in%SUBDRIFTONI)),]
  
  SUBMedDiff<-NA
  
  for (i in 1:nrow(perpopSUBSAMPLE)){
    TARGET<-perpopRED$PopName[i]
    
    targetpop<-perpopSUBSAMPLE$PopName[i]
    # temp<-FstLisinfoSUBSAMPLE[which(FstLisinfoSUBSAMPLE$Pop1==targetpop),]
    
    meltini<-SUBFstListREDinfo_noDuplicateNeighborsNoDrift[which(SUBFstListREDinfo_noDuplicateNeighborsNoDrift$Pop1==targetpop),]
    meltiniSAME<-meltini[which(meltini$SameFamily=="YES"),]
    maxgeofam<-max(meltiniSAME$GEOdist, na.rm = T)
    if(maxgeofam<500){
      maxgeofam<-500
    }
    meltini<-meltini[which(meltini$GEOdist<=maxgeofam),]
    meltiniNO<-meltini[which(meltini$SameFamily!="YES"),]
    
    if(length(which(meltini$SameFamily=="YES"))>2&
       length(which(meltini$SameFamily=="NO"))>2){   # need to have at least 3 comparisons between and within family
      
      SUBMedDiff[i]=median(meltiniNO$FstLinear)-median(meltiniSAME$FstLinear)
    }  
  }
  
  ProportionMisaligned[k]<- length(which(SUBMedDiff<0))/length(which(!is.na(SUBMedDiff)))
  ProportionAligned[k]<- length(which(SUBMedDiff>0))/length(which(!is.na(SUBMedDiff)))
  
}


FinalSchemeResampling<-as.data.frame(cbind(differenceSUB,EnclavesResampling[,1:2],ProportionAligned,ProportionMisaligned))
FinalSchemeResampling$repetition<-c(1:nRepetitions)
MELTFinalSchemeResampling <- melt(FinalSchemeResampling, id.vars='repetition', 
                                  measure.vars=c("differenceSUB","MATCH",  "MISMATCH", "ProportionAligned", "ProportionMisaligned"))

valuesToCompare<-c(
  length(which(aa$difference>0))/(nrow(perpopRED)-length(which(is.na(aa$difference)))), # proportion of the pops who  do have closer FST with speakers of another language family
  (table(perpopRED$EnclavesMismatch)/nrow(perpopRED))[1:2], # proportion enclaves match and mismatch MaMi total
  length(ListAligned)/length(which(!is.na(perpopRED$MedDiff))),  # proportion aligned MaMi total
  length(ListMisaligned)/length(which(!is.na(perpopRED$MedDiff))))  # proportion misaligned MaMi total

valuesToCompare<-as.numeric(valuesToCompare)

ggg<-ggplot(data=MELTFinalSchemeResampling,aes(x=variable, y=value)) 

ggg+ geom_jitter(shape=16, position=position_jitter(0.2))+
  geom_violin(trim=FALSE, alpha=0.4)+
  stat_summary(fun.y= median, 
               geom="pointrange", color="gold") +
  annotate("point", x = 1:5, y = valuesToCompare, color="purple", size=5, alpha=0.7) +
  labs(y="",x = "")+
  scale_x_discrete(labels=c("Closest to\nLinguistically unrelated", "Match", "Mismatch", "Proportion\nAligned",
                            "Proportion\nMisaligned"))+
  theme_minimal()
ggsave("FIG_SX_RESAMPLING.pdf", height=4,width=5.5, useDingbats=FALSE)




#************************************************************
##### script patch
 ## dbmem, spatial autocorrelation analysis
# script written by Alexandros Sotiropoulos alexandros.sotiropoulos@botinst.uzh.ch 
# October 2021
#************************************************************
# Figure S4


library(ade4)
library(sp)
library(vcfR) 
library(adegenet) 
library(adegraphics) 
library(pegas) 
library(StAMPP) 
library(lattice) 
library(gplots) 
library(ape) 
library(ggmap)
library(MASS)

library(memgene)
library("maps")
library(geodist)

library(rgdal)


coords_chiara <-    # tab-separated file for all pops with lon and lat
aa.D.ind_world_chiara <-      # tab-separated file for all pops matrix Fst distances, same number of rows as coords_chiara
aa.D.ind_world_chiara2 <- do.call(rbind, lapply(aa.D.ind_world_chiara, as.numeric))
typeof(aa.D.ind_world_chiara2)
typeof(coords_chiara)

dm_world <- mgQuick(aa.D.ind_world_chiara2, coords_chiara, longlat = TRUE, truncation = NULL, transformation = NULL, forwardPerm = 10000, forwardAlpha = 0.05,
                    finalPerm = 10000, doPlot = 3, verbose = TRUE) #excluding JPN-AUS

write.table(x = dm_world$memgene, file = "dbmem_all_chiara_100perm")
coords_chiara <- read.csv ("listpops_Lon_Lat_REDUCED.txt", sep ="\t")     # tab-separated file for all pops 
aa.D.ind_world_chiara <- read.csv ("fstShapeMatrixMaMi_REDUCED.txt", sep ="\t",header =T, row.names = 1)     # tab-separated file for all pops 
aa.D.ind_world_chiara2 <- do.call(rbind, lapply(aa.D.ind_world_chiara, as.numeric))
typeof(aa.D.ind_world_chiara2)
typeof(coords_chiara)

dm_world <- mgQuick(aa.D.ind_world_chiara2, coords_chiara, longlat = TRUE, truncation = NULL, transformation = NULL, forwardPerm = 100, forwardAlpha = 0.05,
                    finalPerm = 10, doPlot = 3, verbose = TRUE) #excluding JPN-AUS

dm_world_morePerm <- mgQuick(aa.D.ind_world_chiara2, coords_chiara, longlat = TRUE, truncation = NULL, transformation = NULL, forwardPerm = 1000, forwardAlpha = 0.05,
                             finalPerm = 1000, doPlot = 3, verbose = TRUE) #excluding JPN-AUS
dm_world_evenmorePerm <- mgQuick(aa.D.ind_world_chiara2, coords_chiara, longlat = TRUE, truncation = NULL, transformation = NULL, forwardPerm = 10000, forwardAlpha = 0.05,
                                 finalPerm = 10000, doPlot = 3, verbose = TRUE) #excluding JPN-AUS

map('world',boundary=T, fill=T, col="grey")
mgMap(coords_chiara, dm_world$memgene[, 1], add.plot = TRUE, legend = T)
mgMap(coords_chiara, dm_world$memgene[, 2], add.plot = TRUE, legend = T)
mgMap(coords_chiara, dm_world$memgene[, 3], add.plot = TRUE, legend = T)



