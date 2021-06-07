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
#******************************************
#******************************************
#******************************************
# ANALYSIS PERCENTILES
#******************************************
# OVERVIEW OF MISMATCHES WITH CLOSE FST DISTANCES
#******************************************
#******************************************

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
# figure 1 area proportion of linguistically unrelated in close FSTs
# now supplementary
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

##### FIGURE 1C ****************************************

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

### tabl

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

enclavesByMistake<-c("Yoruba"  ,"Nama", "Han-NChina" ) # these pops are not real enclaves, their linguistically unrelated pair is the one driving this genetic proximity effect
ListEnclaves<-ListEnclaves[-which(ListEnclaves%in%(enclavesByMistake))]


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

colorchoice2<-colorchoice
colorchoice2[10]<-"gray20"
colorchoice2[11:16]<-colorchoice[10:15]

ga<-ggplot(meltperpop, aes(variable,value , color=MainFamilies))
ga+ 
  #geom_violin(trim=FALSE, alpha=0.4)+
  #  stat_summary(fun.data="mean_sdl", mult=1, geom="pointrange", color="gray20", alpha=0.5)+
  geom_boxplot(color="gray20", alpha=0.7, outlier.shape = NA)+
  geom_jitter(shape=16, position=position_jitter(0.2))+
  scale_color_manual(values = colorchoice2)+
  theme_light()
ggsave("distributCloseFSTsameFamDiffFam.pdf", useDingbats=FALSE, height = 5, width = 4)

#*#******************************************
#* 
#* 
#* 




### 
#*#******************************************
#*  #*#******************************************
#*  now the bulk of the analysis with the Mismatch Indicators
#*    #*#******************************************
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


## exlude comparisons with the enclaves

FstListREDinfo_noDuplicateNeighbors<-FstListREDinfo_noDuplicateNeighbors[-c
                                              (which(FstListREDinfo_noDuplicateNeighbors$Pop2%in%ListEnclaves), 
                                                which(FstListREDinfo_noDuplicateNeighbors$Pop1%in%ListEnclaves)),]

#******************************************
#*#******************************************
#*#******************************************
#*#******************************************
#*#******************************************
#*#******************************************
#*    WILCOXON TEST
# Mismatch Indicator 1
#* FST distribution between and within language families
#*#******************************************
#*#******************************************
#*
#*
#*
#*# now with geo filter Wilcoxon test
# consider a geographic range equal to the one of the most distant from the same family
# if it is too small, expand until 500 km


## figure distribution and Wilcox test for each population



perpopRED$WtestGEOfilter<-NA
perpopRED$WtestPvalueGEOfilter<-NA


for (i in 1:nrow(perpopRED)){
  TARGET<-perpopRED$PopName[i]
  meltini<-FstListREDinfo_noDuplicateNeighbors[which(FstListREDinfo_noDuplicateNeighbors$Pop1==TARGET),]
  meltiniSAME<-meltini[which(meltini$SameFamily=="YES"),]
  maxgeofam<-max(meltiniSAME$GEOdist, na.rm = T)
  if(maxgeofam<500){
    maxgeofam<-500
  }
  meltini<-meltini[which(meltini$GEOdist<maxgeofam),]
  
  if(length(which(meltini$SameFamily=="YES"))>2){
    wtest<-wilcox.test(Fst ~ SameFamily, exact= F, data=meltini, alternative = "greater")  # Wilcoxon
    perpopRED$WtestGEOfilter[i]<-wtest$statistic
    perpopRED$WtestPvalueGEOfilter[i]<-wtest$p.value
    
  }
}



> length(which(perpopRED$WtestPvalueGEOfilter<0.01))
[1] 212

length(which(!is.na(perpopRED$WtestPvalueGEOfilter)))
[1] 330

> 212/330
[1] 0.6424242

# 64 % of the time the distribution Fst within and between is significantly different
# percentage calculated over the populations for which it is possible to calculate Wilcoxon test (enough between language family comparisons)


#******************************************
## figure distribution and Wilcox test for each Major Language Family
##supplementary

meltinitony<-c()
wtest<-c()
for (i in 1:length(MainFamilies)){
  TARGET<-MainFamilies[i]
  meltini<-FstListREDinfo_noDuplicateNeighbors[which(FstListREDinfo_noDuplicateNeighbors$family1==TARGET),]
  meltiniSAME<-meltini[which(meltini$SameFamily=="YES"),]
  maxgeofam<-max(meltiniSAME$GEOdist, na.rm = T)
  if(maxgeofam<500){
    maxgeofam<-500
  }
  meltini<-meltini[which(meltini$GEOdist<maxgeofam),]
  #wtest<-wilcox.test(Fst ~ SameFamily, exact=F, data=meltini, alternative = "greater")
  meltini$familytarget<-TARGET
  meltini$maxgeofam<-maxgeofam
  wtest[i]<-wilcox.test(Fst ~ SameFamily, exact=F, data=meltini, alternative = "greater")$p.value
  meltini$wtestFamily<-wtest[i]
  meltinitony<-rbind(meltinitony,meltini)
}
round(wtest,digits = 4)

#### supplementaryfigure !!!
# Figure S7

ggg<-ggplot(meltinitony,aes(y=FstLinear, x=SameFamily,shape=SameFamily, fill=maxgeofam))

ggg+ geom_jitter(shape=16, position=position_jitter(0.2), color="gray80", size=0.8)+
  geom_violin(trim=FALSE, alpha=0.4)+
  stat_summary(fun.data="mean_sdl", 
               geom="pointrange", color="darkblue") +
  scale_fill_gradient( low="darkred",
                       high="cyan", name="max radius comparison in km")+
   # ggtitle( subtitle = paste0("W test  p value:", round(meltinitony$wtestFamily, digits = 4)))+
  theme_light()+
  facet_wrap(~familytarget,ncol=3,scales = "free_y")
ggsave("family_plotFST_density_YESorNOGEOfilter_2021.pdf", useDingbats=FALSE, height = 8, width = 8)


#******************************************
## figure distribution and Wilcox test for each population
##supplementary Figure S9

pdf("EACHPOP_violin_YESorNO_GEOradius_2021.pdf")  


#listplot<-list()

for (i in 1:nrow(perpopRED)){
  TARGET<-perpopRED$PopName[i]
  meltini<-FstListREDinfo_noDuplicateNeighbors[which(FstListREDinfo_noDuplicateNeighbors$Pop1==TARGET),]
  meltiniSAME<-meltini[which(meltini$SameFamily=="YES"),]
  maxgeofam<-max(meltiniSAME$GEOdist, na.rm = T)
  if(maxgeofam<500){
    maxgeofam<-500
  }
  meltini<-meltini[which(meltini$GEOdist<maxgeofam),]
  # wtest<-wilcox.test(Fst ~ SameFamily, data=meltini)
  
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


######################################################################
########################################################
##############
# MISMATCH INDICATOR 3b geographic distance and linear regression
## testing if SAME FAMILY has an effect on Fst given GeoDist
############################
############################


# Fst is modeled as the dependent variable with SameFamily as the factor and GEOdist as the covariate

# ********************************************
# per each language family

MainFamiliesTab<-cbind(MainFamilies,rep(NA,length(MainFamilies)))
MainFamiliesTab<-as.data.frame(MainFamiliesTab)

pdf("FamiliesFSTGeo_linearModels_2021.pdf", width = 10, height = 17)

meltiniREDTOTAL<-c()

par(mfrow=c(5, 3))

for (i in 1:length(MainFamilies)){
  TARGET<-MainFamilies[i]
  meltini<-FstListREDinfo[which(FstListREDinfo$family1==TARGET),]
  meltiniRED<-meltini[which(meltini$FstLinear<=max(meltini$FstLinear[which(meltini$SameFamily=="YES")])),]   #cut the Fst larger than the max Fst within Family
  
  meltiniRED<-meltiniRED[!is.na(meltiniRED$GEOdist),] # eliminate the rows for which i do not have geo coordinated in one of the pops
  mod1 <- lm(FstLinear ~ GEOdist * SameFamily, data=meltiniRED) # the general mixed model with reference level as "NO" and effect of same family on slope and intercept
  mod2 <- lm(FstLinear ~ GEOdist, data=meltiniRED) # the general simple model 
  meltiniRED$residuals<-mod2$residuals  # add the residuals for each pair of pops acording to the general simple model
  meltiniRED$targetFamily<-TARGET
  meltiniREDTOTAL<-rbind(meltiniREDTOTAL,meltiniRED) ## create the sum file with all the residuals for major families
  
  tempWithin<-meltiniRED[which(meltiniRED$FAMILY!="DIVERSE"),] # only inside the same Family
  tempbetween<-meltiniRED[which(meltiniRED$FAMILY=="DIVERSE"),] # outside the same Family
  RatioResiduals=round(length(which(tempWithin$residuals>0))/length(which(tempWithin$residuals<0)),digits = 2)
  RatioResidualsREST=round(length(which(tempbetween$residuals>0))/length(which(tempbetween$residuals<0)),digits = 2)
  
  reg1 <- lm(FstLinear ~GEOdist, data=meltiniRED[which(meltiniRED$SameFamily=="YES"),]) # one model for the same family
  reg2 <- lm(FstLinear ~GEOdist, data=meltiniRED[which(meltiniRED$SameFamily=="NO"),]) # one model for the different family
  MainFamiliesTab$RatioResiduals[i]<-RatioResiduals
  MainFamiliesTab$RatioResidualsREST[i]<-RatioResidualsREST
  MainFamiliesTab$Intercept[i]<-round(mod1$coefficients, digits = 5)[1]
  MainFamiliesTab$GEOdistintercept[i]<-round(mod1$coefficients, digits = 5)[2]
  MainFamiliesTab$SameFamilyYESintercept[i]<-round(mod1$coefficients, digits = 5)[3]
  MainFamiliesTab$GEOdist_SameFamilyYES_EFFECT[i]<-round(mod1$coefficients, digits = 5)[4]
  
  plot(FstLinear ~GEOdist, data=meltiniRED, type='n', main=TARGET)
  points(meltiniRED[which(meltiniRED$SameFamily=="NO"),]$GEOdist,meltiniRED[which(meltiniRED$SameFamily=="NO"),]$FstLinear,col="red", pch=1)
  points(meltiniRED[which(meltiniRED$SameFamily=="YES"),]$GEOdist,meltiniRED[which(meltiniRED$SameFamily=="YES"),]$FstLinear, col="blue",pch=20)
  # text(meltiniRED$GEOdist,meltiniRED$FstLinear, labels = meltiniRED$popslistemp, cex=0.5)
  abline(reg1, lty=1)
  abline(reg2, lty=2)
  abline(mod2, lty=3, col="lightblue", lwd=3)
}

dev.off()

#*********************************
# NOW FOR EACH POPULATION
perpopREDGEO<-perpopRED[-which(is.na(perpopRED$lat)),]
FstListREDinfoGEO<-FstListREDinfo[-which(is.na(FstListREDinfo$GEOdist)),]

perpopRED$ratioresidualsWithin<-NA
perpopRED$ratioresidualsBetween<-NA

pdf("EACHPOP_FSTGeo_linearModels_2021.pdf")

par(mfrow=c(5,3), mar = c(2, 2, 2, 2))

for (i in 1:nrow(perpopREDGEO)){ 
  TARGET<-perpopREDGEO$PopName[i]
  meltini<-FstListREDinfoGEO[which(FstListREDinfoGEO$Pop1==TARGET),]
  meltiniRED<-meltini[which(meltini$FstLinear<=max(meltini$FstLinear[which(meltini$SameFamily=="YES")])*2),]   #cut the Fst larger than TWICE the max Fst within Family, but also outlier FST >0.1
  meltiniRED<-meltiniRED[which(meltiniRED$FstLinear<=0.1),]   #cut the Fst larger than TWICE the max Fst within Family, but also outlier FST >0.1
  
  if(length(table(meltiniRED$SameFamily))>1){
    if(table(meltiniRED$SameFamily)[[1]]>2&table(meltiniRED$SameFamily)[[2]]>2){ # there have to be enough data points to perform the analysis, some populations will be excluded
      mod1 <- lm(FstLinear ~ GEOdist * SameFamily, data=meltiniRED) # the general mixed model with reference level as "NO" and effect of same family on slope and intercept
      mod2 <- lm(FstLinear ~ GEOdist, data=meltiniRED) # the general simple model 
      meltiniRED$residuals<-mod2$residuals  # add the residuals for each pair of pops acording to the general simple model
      tempWithin<-meltiniRED[which(meltiniRED$FAMILY!="DIVERSE"),] # only inside the same Family
      tempbetween<-meltiniRED[which(meltiniRED$FAMILY=="DIVERSE"),] # outside the same Family
      RatioResiduals=round(length(which(tempWithin$residuals>0))/length(which(tempWithin$residuals<0)),digits = 2)
      RatioResidualsREST=round(length(which(tempbetween$residuals>0))/length(which(tempbetween$residuals<0)),digits = 2)
      perpopRED$ratioresidualsWithin[i]<-RatioResiduals
      perpopRED$ratioresidualsBetween[i]<-RatioResidualsREST
      reg1 <- lm(FstLinear ~GEOdist, data=meltiniRED[which(meltiniRED$SameFamily=="YES"),])
      reg2 <- lm(FstLinear ~GEOdist, data=meltiniRED[which(meltiniRED$SameFamily=="NO"),])
      
      plot(FstLinear ~GEOdist, data=meltiniRED, type='n', main=paste0(TARGET," - ",perpopREDGEO$glottolog.NAME[i]),
           cex.main=0.5,cex.lab=0.5,cex.axis=0.5)
      points(meltiniRED[which(meltiniRED$SameFamily=="NO"),]$GEOdist,meltiniRED[which(meltiniRED$SameFamily=="NO"),]$FstLinear,col="red", pch=1)
      points(meltiniRED[which(meltiniRED$SameFamily=="YES"),]$GEOdist,meltiniRED[which(meltiniRED$SameFamily=="YES"),]$FstLinear, col="blue",pch=20)
      
      abline(reg1, lty=1)
      abline(reg2, lty=2)
      abline(mod2, lty=3, col="lightblue", lwd=3)
     }
  }
}
dev.off()

# write.table(perpopRED, "/Users/chiarabarbieri/Library/Mobile Documents/com~apple~CloudDocs/GeLaTo/perPopMarchMami2.txt",  row.names = F, sep = "\t", quote = F)


#************************************************************
# calculate Mismatch Indicator variables for each family
#### Figure 2 panel below
#************************************************************
library(reshape)
size<-as.numeric(table(perpopREDfamily$MainFamilies))

perLANGminiplot<-perpopREDfamily[,c("MainFamilies", "percentileMismatch","WtestPvalueGEOfilter","ratioresidualsWithin" )]

MELTperLANGminiplot<-melt(perLANGminiplot)
levels(MELTperLANGminiplot$variable)<-c("MI1_GeneticallyClose", "MI2_FSTdistributions", "MI3b_GeographicNeighborsIBD")

ga<-ggplot(MELTperLANGminiplot,aes(x=MainFamilies,y=value, color=MainFamilies))
ga+ 
  #geom_violin(trim=FALSE, alpha=0.4)+
#  stat_summary(fun.data="mean_sdl", mult=1, geom="pointrange", color="gray20", alpha=0.5)+
  geom_boxplot(color="gray20", alpha=0.7, outlier.shape = NA)+
  geom_jitter(shape=16, position=position_jitter(0.2))+
  scale_color_manual(values = colorchoice)+
  theme_light()+
  theme(axis.text.x = element_text(angle = 45, vjust=1, hjust=1),legend.position = "none") +  
  facet_wrap(~variable,ncol=1, scales = "free_y")

ggsave("MismatchIndicatorJitter_perfamily.pdf", useDingbats=FALSE, height = 6, width = 8, units = "in")


# ---------------------------------------------------------
#### the distribution of divergence time for each language family 
#### Figure 2 panel below
# three different plots, did not decide yet
# ---------------------------------------------------------

perpopREDNe<-perpopRED[which(perpopRED$USEforNe_calculation=="YES"),]
possiblepopswithNE<-perpopREDNe$PopName

meltFstREDinfoYESne<-FstListREDinfo[which(FstListREDinfo$Pop1%in%possiblepopswithNE&FstListREDinfo$Pop2%in%possiblepopswithNE),]

meltFstREDinfoYESneSameFamily<-meltFstREDinfoYESne[which(meltFstREDinfoYESne$FAMILY!="DIVERSE"),]
meltFstREDinfoYESneONLYfamilies<-meltFstREDinfoYESneSameFamily[which(meltFstREDinfoYESneSameFamily$family1%in%MainFamilies),]
#DRIFTONI<-perpopRED[which(perpopRED$medianFSTregion>0.1&perpopRED$medianFST>0.1),]$PopName
meltFstREDinfoYESneONLYfamilies<-meltFstREDinfoYESneONLYfamilies[-which(meltFstREDinfoYESneONLYfamilies$Pop1%in%DRIFTONI),]
meltFstREDinfoYESneONLYfamilies<-meltFstREDinfoYESneONLYfamilies[-which(meltFstREDinfoYESneONLYfamilies$Pop2%in%DRIFTONI),]

lb <- function(x) mean(x) - sd(x)
ub <- function(x) mean(x) + sd(x)

sumld <- meltFstREDinfoYESneONLYfamilies %>% 
#  select(-FAMILY) %>% 
  group_by(FAMILY) %>% 
  summarise_all(funs(mean, median, lower = lb, upper = ub))

library(tidyverse)
source("https://gist.githubusercontent.com/benmarwick/2a1bb0133ff568cbe28d/raw/fb53bd97121f7f9ce947837ef1a4c65a73bffb3f/geom_flat_violin.R")

g <-  ggplot(data = meltFstREDinfoYESneONLYfamilies, 
         aes(x = FAMILY, y = TMRCA_doubleNe/1000, fill = FAMILY)) +
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
  geom_point(aes(y = TMRCA_doubleNe/1000, color = FAMILY), 
             position = position_jitter(width = .15), size = .5, alpha = 0.8) +
  geom_boxplot(width = .1, outlier.shape = NA, alpha = 0.8) +
 #  geom_point(data = sumld, aes(x = FAMILY, y = TMRCA_doubleNe_median/1000), 
  #           position = position_nudge(x = 0.3), size = 2.5) +
  #geom_errorbar(data = sumld, aes(ymin = TMRCA_doubleNe_95_lower/1000, ymax = TMRCA_doubleNe_95_upper/1000, y = TMRCA_doubleNe_mean/1000), position = position_nudge(x = 0.3), width = 0) +
 # expand_limits(x = 5.25) +
  guides(fill = FALSE) +
  guides(color = FALSE) +
  scale_color_manual(values = MainFamilies2$COLOR[which(MainFamilies2$MainFamilies%in%meltFstREDinfoYESneONLYfamilies$family1)])+
scale_fill_manual(values = MainFamilies2$COLOR[which(MainFamilies2$MainFamilies%in%meltFstREDinfoYESneONLYfamilies$family1)])+
coord_flip() + # flip or not
  ylim(0,20)+   # cut at 20,000 years ago
  labs(y="Divergence time, k years ago", x="") +
  theme_bw() 
g
ggsave("DensityTMRCAlimitY_noDrifted_MajorFamilies2021_rainclouds.pdf", useDingbats=FALSE, width = 7, height = 5)



ggplot(meltFstREDinfoYESneONLYfamilies, aes(TMRCA_doubleNe/1000 , fill=FAMILY))+
  geom_density(alpha=0.4)+ 
  facet_wrap( ~ FAMILY, ncol=3,scales = "free_y")+
  xlim(0,20)+
  #  ggtitle("TMRCA density major language families, exclude pop drifted")+
  labs(x="Divergence time, k years ago", y="") +
  theme(axis.text.y=element_blank())+
  scale_fill_manual(values = MainFamilies2$COLOR[which(MainFamilies2$MainFamilies%in%meltFstREDinfoYESneONLYfamilies$family1)])

ggsave("DensityTMRCAlimitY_noDrifted_MajorFamilies2021.pdf", useDingbats=FALSE, width = 7, height = 5)



ga<-ggplot(meltFstREDinfoYESneONLYfamilies, aes(FAMILY,TMRCA_doubleNe/1000 , color=FAMILY))
ga+ 
  #geom_violin(trim=FALSE, alpha=0.4)+
  #  stat_summary(fun.data="mean_sdl", mult=1, geom="pointrange", color="gray20", alpha=0.5)+
  geom_boxplot(color="gray20", alpha=0.7, outlier.shape = NA)+
  geom_jitter(shape=16, position=position_jitter(0.2))+
  scale_color_manual(values = colorchoice)+
  theme(axis.text.x = element_text(angle = 45, vjust=1, hjust=1),legend.position = "none")+
  ylim(0,20)
  
ggsave("TMRCAlimitY_noDriftedJitter_perfamily.pdf", useDingbats=FALSE, height = 3, width = 8, units = "in")





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



####
# Figure S8 Wilcox map - Mismatch Indicator 2
########

perpopREDSHIFTMAP_2wilcox<-perpopREDSHIFTMAP[which(perpopREDSHIFTMAP$WtestPvalueGEOfilter>0.01),]
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



########
## the map with singular mismatches of linguistic and genetic migrants - Mismatch Indicator 3a
# Figure 1 potential
########
library(ggrepel)

perpopRED2<-perpopRED[!is.na(perpopRED$lat),]
perpopREDSHIFTMAP<-MoveDataCoord(perpopRED2)
perpopREDSHIFTMAPinterest<-perpopREDSHIFTMAP[!is.na(perpopREDSHIFTMAP$Mismatch3a),]
perpopREDSHIFTMAPmismatch<-perpopREDSHIFTMAPinterest[which(perpopREDSHIFTMAPinterest$Mismatch3a%in%c("MISMATCH","secondaryMISMATCH_FSTzeroDiffFamily")),]
#perpopREDSHIFTMAPNoMismatch<-perpopREDSHIFTMAP[which(perpopREDSHIFTMAP$hasOtherMatches=="NO"),]

base_plot + geom_point(data = perpopREDSHIFTMAP,
                       aes(x = lon.1, y = lat.1),
                       color = "black",
                       size = 0.5, shape=3)+
  geom_point(data=perpopREDSHIFTMAPinterest, aes(x = lon.1, y = lat.1, 
                                         color=Mismatch3a),alpha=0.8)+
  geom_text_repel(data=perpopREDSHIFTMAPmismatch, aes(x = lon.1, y = lat.1, 
                                               label=PopName), size=2)+
  xlab("") + ylab("") +
  scale_colour_manual(values = c("blue", "red", "orange", "yellow"), name="", 
                      labels=c("Match", "Mismatch Genetic Enclave","Mismatch Linguistic Enclave", "no neighbors from the same Language Family"))

ggsave("MapMismatchIndicator3a_.pdf", useDingbats=FALSE, width = 15, height = 8, dpi = 300, units = "in")

#**************************************************
### table of cases

mismatch3a<-perpopRED[grep("MISMATCH",perpopRED$Mismatch3a),]
perpopRED$DIFFFamilyClosestFAM<-perpopRED$glottolog.NAME[match(perpopRED$DIFFFamilyClosestpop,perpopRED$PopName)]
colnameschoice<-c("PopName","glottolog.NAME", "geodistSameFamily","geodistDIFFFamily","DIFFFamilyClosestpop", "DIFFFamilyClosestFAM")
mismatch3a<-mismatch3a[,colnameschoice]
#write.table(mismatch3a[order(mismatch3a$geodistDIFFFamily),],"tableSupplMismatches3a.txt", sep="\t", row.names = F, quote=F) 

#**************************************************
### lollipop proportion 3a
size<-as.numeric(table(perpopREDfamily$MainFamilies))

proportion3a<-table(perpopREDfamily$MainFamilies,perpopREDfamily$Mismatchsingular)/size
proportion3aMELT<-melt(proportion3a)
colnames(proportion3aMELT)<-c("Family", "variable", "proportion")

gi<-ggplot(proportion3aMELT, aes(Family, proportion))

gi+ geom_segment( aes(x=Family, xend=Family, y=0, yend=proportion))+
  geom_point(aes(color=Family), size=5)+
scale_color_manual(values = colorchoice)+
  theme_light()+
  theme(axis.text.x = element_text(angle = 45, vjust=1, hjust=1)) +
  facet_wrap(~variable,ncol=1)
ggsave("lollipop_3a.pdf", useDingbats=FALSE, width = 8, height = 5)



#********************************************
#********************************************
#********************************************
# ***************************************************
# ANALYSIS OF LANGUAGE ISOLATES
# ***************************************************

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
  scale_colour_gradient2(midpoint=1, low="green", mid="gray",
                         high="orange")+
  geom_label_repel(data=perpopREDEUROPE, aes(x=lon, y=lat,label=PopName, color=proportionFstAdjustedNeighbors), size=2, label.padding=0.1)+
  geom_point(data=perpopREDEUROPE, aes(x=lon, y=lat, shape=glottolog.NAME), size=1 ) +
  #ggtitle("Proportion FST with neighbors")+
  labs(shape="Language Family", colour="Poportion FST with neighbors")+
  theme(axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank())

perpopREDEUROPESIGNF<-perpopREDEUROPE[which(perpopREDEUROPE$WtestPvalueGEOfilter<0.01),]

wilcEUR<- gg + geom_point(data=perpopREDEUROPE, 
                        aes(x=lon, y=lat, color=WtestPvalueGEOfilter,na.rm=T), size=3 , alpha=0.4) +
  scale_color_gradient2( low="darkred",mid="pink",
                         high="blue",  midpoint = 0.2, na.value = "grey50", name="P value Wilcoxon Test")+
  geom_label_repel(data=perpopREDEUROPE, aes(x=lon, y=lat,label=PopName, color=WtestPvalueGEOfilter, na.rm= TRUE ), size=2, label.padding=0.1)+
  geom_point(data = perpopREDEUROPESIGNF, aes(x = lon, y = lat), color = "darkred",  size = 2, alpha=0.5)+
  geom_point(data=perpopREDEUROPE, aes(x=lon, y=lat, shape=glottolog.NAME), size=1 ) +
  labs(shape="Language Family", colour="P value Wilcoxon Test")+
  theme(axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank())

# *********************************************
##### AFRICA #####
# *********************************************

# MAP OF sub saharan africa 

coordinates<-c(-36,21,-20,52)
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
  geom_point(data=perpopREDAFRICA, aes(x=lon, y=lat, shape=MainFamilies), size=1 ) +
    labs(shape="Language Family", colour="Poportion FST with neighbors")+
  #scale_shape_manual(values=1:9)+
  theme(axis.title=element_blank(), axis.text=element_blank(), axis.ticks=element_blank())


perpopREDAFRICASIGNF<-perpopREDAFRICA[which(perpopREDAFRICA$WtestPvalueGEOfilter<0.01),]

WilcAfrica<- gg + geom_point(data=perpopREDAFRICA, 
                        aes(x=lon, y=lat, color=WtestPvalueGEOfilter,na.rm=T), size=3 , alpha=0.4) +
  scale_color_gradient2( low="darkred",mid="pink",
                         high="blue",  midpoint = 0.1, na.value = "grey50", name="P value Wilcoxon Test")+
  geom_label_repel(data=perpopREDAFRICA, aes(x=lon, y=lat,label=PopName, color=WtestPvalueGEOfilter ), size=2, label.padding=0.1)+
  geom_point(data = perpopREDAFRICASIGNF, aes(x = lon, y = lat), color = "darkred",  size = 2, alpha=0.5)+
  geom_point(data=perpopREDAFRICA, aes(x=lon, y=lat, shape=MainFamilies), size=1 ) +
  labs(shape="Language Family", colour="P value Wilcoxon Test")+
 #scale_shape_manual(values=1:9) +
  theme(axis.title=element_blank(), axis.text=element_blank(), axis.ticks=element_blank())

# *********************************************
##### OCEANIA #####

# MAP OF OCEANIA

coordinates<-c(-25,24,92,180)
perpopREDOCEANIA<-perpopRED[which(perpopRED$lat>coordinates[1]&perpopRED$lat<coordinates[2]&perpopRED$lon>coordinates[3]&perpopRED$lon<coordinates[4]),]
map.world <- map_data(map="world")
library(ggrepel)

gg <- ggplot()
gg <- gg + theme()
gg <- gg + geom_map(data=map.world, map=map.world, aes(map_id=region), fill="white", colour="black", size=0.15)
gg<- gg+coord_quickmap(ylim=coordinates[1:2], xlim=coordinates[3:4])

fstOCEANIA<-  gg + geom_point(data=perpopREDOCEANIA[!is.na(perpopREDOCEANIA$proportionFstAdjustedNeighbors),], 
                        
                        aes(x=lon, y=lat, color=proportionFstAdjustedNeighbors), size=3 , alpha=0.4) +
  scale_colour_gradient2(midpoint=1, low="darkgreen", mid="gray",
                         high="orange")+
  geom_label_repel(data=perpopREDOCEANIA, aes(x=lon, y=lat,label=PopName, color=proportionFstAdjustedNeighbors), size=2, label.padding=0.1)+
  geom_point(data=perpopREDOCEANIA, aes(x=lon, y=lat, shape=MainFamilies), size=1 ) +
  labs(shape="Language Family", colour="Poportion FST with neighbors")+
  #scale_shape_manual(values=1:8)+
  theme(axis.title=element_blank(), axis.text=element_blank(), axis.ticks=element_blank())

perpopREDOCEANIASIGNF<-perpopREDOCEANIA[which(perpopREDOCEANIA$WtestPvalueGEOfilter<0.01),]

wilcOCEANIA<- gg + geom_point(data=perpopREDOCEANIA, 
                        aes(x=lon, y=lat, color=WtestPvalueGEOfilter,na.rm=T), size=3 , alpha=0.4) +
  scale_color_gradient2( low="darkred",mid="pink",
                         high="blue",  midpoint = 0.05, na.value = "grey50", name="P value Wilcoxon Test")+
  geom_label_repel(data=perpopREDOCEANIA, aes(x=lon, y=lat,label=PopName, color=WtestPvalueGEOfilter, na.rm= TRUE ), size=2, label.padding=0.1)+
  geom_point(data = perpopREDOCEANIASIGNF, aes(x = lon, y = lat), color = "darkred",  size = 2, alpha=0.5)+
  geom_point(data=perpopREDOCEANIA, aes(x=lon, y=lat, shape=MainFamilies), size=1 ) +
  labs(shape="Language Family", colour="P value Wilcoxon Test")+
  theme(axis.title=element_blank(), axis.text=element_blank(), axis.ticks=element_blank())

# *********************************************
##### CAUCASUS #####

# MAP OF caucasus

coordinates<-c(38,48,34,54)
perpopREDcaucasus<-perpopRED[which(perpopRED$lat>coordinates[1]&perpopRED$lat<coordinates[2]&perpopRED$lon>coordinates[3]&perpopRED$lon<coordinates[4]),]
map.world <- map_data(map="world")
library(ggrepel)

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
  theme(axis.title=element_blank(), axis.text=element_blank(), axis.ticks=element_blank())

perpopREDcaucasusSIGNF<-perpopREDcaucasus[which(perpopREDcaucasus$WtestPvalueGEOfilter<0.01),]

wilcCaucasus<- gg + geom_point(data=perpopREDcaucasus, 
                        aes(x=lon, y=lat, color=WtestPvalueGEOfilter,na.rm=T), size=3 , alpha=0.4) +
  scale_color_gradient2( low="darkred",mid="pink",
                         high="blue",  midpoint = 0.05, na.value = "grey50", name="P value Wilcoxon Test")+
  geom_label_repel(data=perpopREDcaucasus, aes(x=lon, y=lat,label=PopName, color=WtestPvalueGEOfilter, na.rm= TRUE ), size=2, label.padding=0.1)+
  geom_point(data = perpopREDcaucasusSIGNF, aes(x = lon, y = lat), color = "darkred",  size = 2, alpha=0.5)+
  geom_point(data=perpopREDcaucasus, aes(x=lon, y=lat, shape=glottolog.NAME), size=1 ) +
  labs(shape="Language Family", colour="P value Wilcoxon Test")+
  scale_shape_manual(values=1:7)+
  theme(axis.title=element_blank(), axis.text=element_blank(), axis.ticks=element_blank())

##+++++++++++++++++++++++
# main figure combined
# 1: europe, 2: oceania, 3: africa, 4: caucasus

library(ggpubr)
ggarrange(WilcAfrica, wilcEUR, wilcCaucasus,wilcOCEANIA + rremove("x.text"), 
          labels = c("A", "B", "C", "D"),
          ncol = 2, nrow = 2)
ggsave("combinedMAP_wilcox_Figure4_2021.pdf", useDingbats=FALSE, height = 15, width = 17)

ggarrange(fstAfrica, FSTpropEUR, fstCaucasus,fstOCEANIA + rremove("x.text"), 
          labels = c("A", "B", "C", "D"),
          ncol = 2, nrow = 2)
ggsave("combinedMAP_FSTproportion_supple2021.pdf", useDingbats=FALSE, height = 15, width = 17)


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

ggarrange(SANDAWE, HADZA, BASQUE, MALTA, ARMENIAN, AZERBAJAN + rremove("x.text"), 
          # labels = c("A", "B", "C", "D"),
          ncol = 2, nrow = 3)
ggsave("combinedSinglePops_FstGeodist_supple.pdf", useDingbats=FALSE, height = 15, width = 13)

(SANDAWE+ HADZA) | (BASQUE + MALTA) | (ARMENIAN + AZERBAJAN) 


