## summarize and analyze results from GAPIT. ##
load("170206phenotype_analysis.RData")
# names of traits that were analyzed
dimnames(allOneSiteBLUPs)[[2]]   
dimnames(allMultiSiteBLUPs)[[2]]
# Set of SNPs
load("170308forGAPIT.RData") 
SNPnames <- myGM$Name

# load one CSV of P-values to get the format
testP <- read.csv("GAPIT/output/GAPIT..HU.Ccirc.year.3.GWAS.Results.csv")
head(testP)

# matrix to hold FDR-corrected P-values for single-site traits
oneSiteFDRP <- matrix(NA, nrow = length(SNPnames), ncol = dim(Y.single.site)[2] - 1,
                      dimnames = list(SNPnames, dimnames(Y.single.site)[[2]][-1]))
# fill in matrix  ### will need to load new data
for(trait in dimnames(oneSiteFDRP)[[2]]){
  thistab <- read.csv(paste("GAPIT/output/GAPIT..",trait,".GWAS.Results.csv", sep = ""), stringsAsFactors = FALSE)
  oneSiteFDRP[thistab$SNP,trait] <- thistab$FDR_Adjusted_P.values
}
# what is significant?
oneSiteSig <- oneSiteFDRP < 0.05
oneSiteTotalSig <- rowSums(oneSiteSig, na.rm = TRUE) # number significant by SNP
hist(oneSiteTotalSig[oneSiteTotalSig > 0])
table(oneSiteTotalSig)
   # one site*trait*year: 326 sig SNPs
   # two site*trait*year:  24 sig SNPs
   # three              :   1
   # four               :   2


# same for multi-site traits
multiSiteFDRP <- matrix(NA, nrow = length(SNPnames), ncol = dim(allMultiSiteBLUPs)[2],
                        dimnames = list(SNPnames, dimnames(allMultiSiteBLUPs)[[2]]))
for(trait in dimnames(multiSiteFDRP)[[2]]){
  if(file.exists(paste("GAPIT/output/GAPIT..",trait,".GWAS.Results.csv", sep = ""))){
    thistab <- read.csv(paste("GAPIT/output/GAPIT..",trait,".GWAS.Results.csv", sep = ""),
                        stringsAsFactors = FALSE)
  } else {
    thistab <- data.frame(FDR_Adjusted_P.values = NA)
  }
  multiSiteFDRP[thistab$SNP,trait] <- thistab$FDR_Adjusted_P.values
}
# missing ZJU.KNU.CmL.year.2 and ZJU.KNU.CmN.year.2, which had zero heritability
# multi-site significance
multiSiteSig <- multiSiteFDRP < 0.05
multiSiteTotalSig <- rowSums(multiSiteSig, na.rm = TRUE)
table(multiSiteTotalSig)
  # one:   137
  # two:    29
  # three:  43
  # four:    2
  # ten:     1
  # fifteen! 1
for(i in which(multiSiteTotalSig > 3)){
  print(SNPnames[i])
  print(dimnames(multiSiteSig)[[2]][which(multiSiteSig[i,])])
}


## start setup for export
# columns for each location
siteColsSingle <- lapply(c("HU", "NEF", "ZJU"),
                         function(x) grep(x, dimnames(oneSiteSig)[[2]]))
names(siteColsSingle) <- c("HU", "NEF", "ZJU")
siteColsMulti <- list(HU.NEF = 1:29, ZJU.KNU = 30:53, HU.NEF.CHA.CSU.KNU = 54:79,
                      HU.NEF.CHA.CSU = 80:108, HU.NEF.CHA.CSU.KNU.ZJU = 109:132)
# columns for each trait
traitNames <- c("Biomass.yield.dry.weight","Ccirc", "Bcirc", "Ccirc.Bcirc", "CmL", "CmN", 
                "IntL", "CmDW", "CmVol", "CmDens", "CmD_I1", "CmD_LI", "Total.culms",
                "Prop.repro.culms", "Culm.per.area")


# for single location trait BLUPs, compare number significant SNPs to heritability
sigSNPvsHerit <- data.frame(traitBLUP = dimnames(oneSiteSig)[[2]],
                            sigSNPs = colSums(oneSiteSig, na.rm = TRUE),
                            heritability = rep(NA, dim(oneSiteSig)[2]))
for(loc in dimnames(myHerit)[[1]]){
  for(trait in dimnames(myHerit)[[2]]){
    ltc <- paste(loc, trait, sep = ".")
    if(ltc %in% sigSNPvsHerit$traitBLUP){
      sigSNPvsHerit[match(ltc, sigSNPvsHerit$traitBLUP),"heritability"] <- 
        myHerit[loc, trait]
    }
  }
}
plot(sigSNPvsHerit$heritability, sigSNPvsHerit$sigSNPs,
     xlab = "Heritability", ylab = "Significant SNPs identified",
     main = "Single site trait BLUPs")
sigSNPvsHerit[which(sigSNPvsHerit$sigSNPs > 200),]

# look at skewness and kurtosis of trait BLUPs; does that explain number of significant SNPs?
require(moments)
sigSNPvsHerit$skewness <- rep(NA, dim(sigSNPvsHerit)[1])
sigSNPvsHerit$kurtosis <- rep(NA, dim(sigSNPvsHerit)[1])
for(i in 1:dim(sigSNPvsHerit)[1]){
  tb <- as.character(sigSNPvsHerit$traitBLUP[i])
  sigSNPvsHerit$skewness[i] <- skewness(Y.single.site[[tb]], na.rm = TRUE)
  sigSNPvsHerit$kurtosis[i] <- kurtosis(Y.single.site[[tb]], na.rm = TRUE)
}
plot(sigSNPvsHerit$skewness, sigSNPvsHerit$sigSNPs,
     xlab = "Skewness", ylab = "Significant SNPs identified",
     main = "Single site trait BLUPs", xlim = c(-6,6))
manysig <- which(sigSNPvsHerit$sigSNPs > 75)
text(sigSNPvsHerit$skewness[manysig], sigSNPvsHerit$sigSNPs[manysig],
     sigSNPvsHerit$traitBLUP[manysig])
plot(sigSNPvsHerit$kurtosis, sigSNPvsHerit$sigSNPs,
     xlab = "Kurtosis", ylab = "Significant SNPs identified",
     main = "Single site trait BLUPs", xlim = c(-1, 55))
text(sigSNPvsHerit$kurtosis[manysig], sigSNPvsHerit$sigSNPs[manysig],
     sigSNPvsHerit$traitBLUP[manysig])


# look at number of individuals with data vs. number of significant SNPs identified
sigSNPvsHerit$nInd <- rep(NA, dim(sigSNPvsHerit)[1])
for(i in 1:dim(sigSNPvsHerit)[1]){
  thistrait <- as.character(sigSNPvsHerit$traitBLUP[i])
  sigSNPvsHerit$nInd[i] <- sum(!is.na(Y.single.site[[thistrait]]))
}
plot(sigSNPvsHerit$nInd, sigSNPvsHerit$sigSNPs, xlab = "Individuals with data",
     ylab = "Significant SNPs identified",
     main = "Single-trait BLUPs") # below 200 invididuals -> all traits above 300 sig SNPs
text(sigSNPvsHerit$nInd[manysig], sigSNPvsHerit$sigSNPs[manysig],
     sigSNPvsHerit$traitBLUP[manysig])
plot(sigSNPvsHerit$nInd, sigSNPvsHerit$skewness, xlab = "Individuals with data",
     ylab = "Skewness") # no relationship between  nInd and skewness
plot(sigSNPvsHerit$nInd, sigSNPvsHerit$kurtosis, xlab = "Individuals with data",
     ylab = "Kurtosis") # no relationship between nInd and kurtosis

# do the same for multi-site BLUPs (tallying # sig SNPs, trait statistics)
sigSNPvsHeritMulti <- data.frame(traitBLUP = dimnames(multiSiteSig)[[2]],
                                 sigSNPs = colSums(multiSiteSig, na.rm = TRUE),
                                 heritability = rep(NA, dim(multiSiteSig)[2]))
for(loc in dimnames(comboHerit)[[2]]){
  for(trait in dimnames(comboHerit)[[1]]){
    ltc <- paste(loc, trait, sep = ".")
    if(ltc %in% sigSNPvsHerit$traitBLUP){
      sigSNPvsHeritMulti[match(ltc, sigSNPvsHeritMulti$traitBLUP),"heritability"] <- 
        comboHerit[trait, loc]
    }
  }
}
sigSNPvsHeritMulti$skewness <- rep(NA, dim(sigSNPvsHeritMulti)[1])
sigSNPvsHeritMulti$kurtosis <- rep(NA, dim(sigSNPvsHeritMulti)[1])
sigSNPvsHeritMulti$nInd <- rep(NA, dim(sigSNPvsHeritMulti)[1])
for(i in 1:dim(sigSNPvsHeritMulti)[1]){
  tb <- as.character(sigSNPvsHeritMulti$traitBLUP[i])
  sigSNPvsHeritMulti$skewness[i] <- skewness(Y.multi.site[[tb]], na.rm = TRUE)
  sigSNPvsHeritMulti$kurtosis[i] <- kurtosis(Y.multi.site[[tb]], na.rm = TRUE)
  sigSNPvsHeritMulti$nInd[i] <- sum(!is.na(Y.multi.site[[tb]]))
}
hist(sigSNPvsHeritMulti$nInd)
plot(sigSNPvsHeritMulti$skewness, sigSNPvsHeritMulti$sigSNPs)

# look at strange SNPs
plot(Y.single.site$ZJU.Bcirc.year.2, Y.single.site$ZJU.Biomass.yield.dry.weight.year.2)
Y.single.site$Taxa[which(Y.single.site$ZJU.Bcirc.year.2 > 20)] # PMS-496
plot(Y.single.site$ZJU.Bcirc.year.2, Y.single.site$ZJU.Bcirc.year.3)
plot(Y.single.site$ZJU.Bcirc.year.2, Y.single.site$NEF.Bcirc.year.3)
  ## ZJU: huge plant, probably error
hist(Y.single.site$NEF.Bcirc.year.3, breaks = 20)
Y.single.site$Taxa[which(Y.single.site$NEF.Bcirc.year.2 > 8)]
plot(Y.single.site$NEF.Bcirc.year.2, Y.single.site$NEF.Bcirc.year.3)
plot(Y.single.site$HU.Bcirc.year.3, Y.single.site$NEF.Bcirc.year.3)
  ## NEF: some outliers just one year, some both

## filter which traits to keep for GWAS analysis
# more than 200 individuals, skewness not exceeding 2 or -2
singleTraitsKeep <- row.names(sigSNPvsHerit)[sigSNPvsHerit$nInd > 200]
plot(sigSNPvsHerit[singleTraitsKeep,"skewness"], sigSNPvsHerit[singleTraitsKeep,"sigSNPs"])
singleTraitsKeep <- singleTraitsKeep[abs(sigSNPvsHerit[singleTraitsKeep,"skewness"]) <= 2]
plot(sigSNPvsHerit[singleTraitsKeep,"kurtosis"], sigSNPvsHerit[singleTraitsKeep,"sigSNPs"])
singleTraitsKeep[sigSNPvsHerit[singleTraitsKeep,"sigSNPs"] > 20]

multiTraitsKeep <- row.names(sigSNPvsHeritMulti)[abs(sigSNPvsHeritMulti$skewness) <= 2]
multiTraitsKeep <- multiTraitsKeep[!is.na(multiTraitsKeep)]

#save(SNPnames, oneSiteFDRP, oneSiteSig, multiSiteFDRP, multiSiteSig, siteColsSingle, siteColsMulti, 
#     traitNames, sigSNPvsHerit, sigSNPvsHeritMulti, singleTraitsKeep, multiTraitsKeep,
#     file = "170320GWASresults.RData")
load("170320GWASresults.RData")

### set up chromosome plots ###
# chromosome lengths
Sbicolor3chrlen <- c(80884392, 77742459, 74386277, 68658214,
                     71854669, 61277060, 65505356, 62686529,
                     59416394, 61233695)
# get alignment and annotation information
myalign <- read.csv("SNP annotations/161104markers.csv", header = TRUE, stringsAsFactors = FALSE)
alignedRows <- which(!is.na(myalign$Sbicolor3SNPpos))
alignedRows <- alignedRows[!alignedRows %in% grep("super", myalign$Sbicolor3Chr)]
# get SNP densities
alignDensities <- list()
length(alignDensities) <- 10
chrNumByRows <- as.integer(substring(myalign$Sbicolor3Chr[alignedRows], 4, 5))
for(chr in 1:10){
  alignDensities[[chr]] <- density(myalign$Sbicolor3SNPpos[alignedRows][chrNumByRows == chr]/1e6,
                                   from = 0, to = Sbicolor3chrlen[chr]/1e6)
}

traitNamesPrint <- c("Dry biomass yield (g per plant)",
                     "Compressed circumference (cm)", "Basal circumference (cm)", 
                     "Compressed circumference/basal circumference", "Culm length (cm)",
                     "Culm node number", "Internode length (cm)", "Culm dry weight (g)",
                     "Culm volume (cm3)", "Culm density (g cm-3)", 
                     "Diameter of basal internode (mm)", "Diameter of topmost internode (mm)",
                     "Total number of culms", "Proportion of reproductive culms",
                     "Culms per footprint (cm-2)")
names(traitNamesPrint) <- traitNames
traitSym <- c("Yld", "CC", "BC", "CC/BC", "CmL", "CmNdN", "IntL", "CmDW",
              "CmV", "CmDW/V", "DBI", "DTI", "TCmN", "RCmN/TCmN", "TCmN/A")
names(traitSym) <- traitNames

manhattan <- function(colHeader, drawTitle = TRUE, drawYlab = TRUE){ # colHeader is column name from oneSiteFDRP or multiSiteFDRP
  isSingle <- colHeader %in% singleTraitsKeep
  isMulti <- colHeader %in% multiTraitsKeep
  if(!isSingle && !isMulti){ # if we don't have P-values, make an empty plot
    plot(NA, xlim = c(0,1), ylim = c(0,1), axes = FALSE, xlab = "", ylab = "")
  } else {
    par(mgp = c(1.5, 0.5, 0))
    sigline <- -log10(0.05) # position of line indicating significance
    # index by SNPs output by TagDigger
    indexByMyalign <- match(myalign$Name, gsub(".", "-", dimnames(oneSiteFDRP)[[1]], fixed = TRUE))
    # get P-values and log transform
    if(isSingle){
      theseLogP <- -log10(oneSiteFDRP[indexByMyalign,colHeader])
    }
    if(isMulti){
      theseLogP <- -log10(multiSiteFDRP[indexByMyalign,colHeader])
    }
    # transpose SNP positions into positions for plotting
    transposeByChr <- sapply(1:11, function(x) ifelse(x ==1, 0, sum(Sbicolor3chrlen[1:(x-1)])))
    names(transposeByChr) <- c(paste("Chr0", 1:9, sep = ""), "Chr10", "unaligned")
    transposeByChr["unaligned"] <- transposeByChr["unaligned"] + 1e7
    newChr <- myalign$Sbicolor3Chr
    newChr[newChr == ""] <- "unaligned"
    newChr[grep("super", newChr)] <- "unaligned"
    newPos <- myalign$Sbicolor3SNPpos
    newPos[newChr == "unaligned"] <- 1:sum(newChr == "unaligned")
    newPos <- newPos + transposeByChr[newChr]
    # make colors
    chrcol <- c(rainbow(12)[c(1,7,2,8)], "yellow3", rainbow(12)[c(9,4,10,5,11,6)])
    names(chrcol) <- names(transposeByChr)
    # get max log P to display
    Ymax <- max(c(sigline+0.1, theseLogP+0.1), na.rm = TRUE)
    # make scatter plot
    par(bg = "white")
    plot(newPos, theseLogP, pch = 16, col = chrcol[newChr], axes = FALSE, ylim = c(0, Ymax),
         ylab = "", xlab = "SNP position")#ifelse(drawTitle, "SNP position", ""))
    axis(2)
    if(drawYlab){
      title(ylab = expression(italic(-log[10])(italic(P)[FDR])))
    }
    abline(h = sigline, lty = 2)
    # overlay SNP densities
    maxDens <- max(sapply(alignDensities, function(x) max(x$y)))
    densScale <- Ymax/3/maxDens # use as scaling factor for y-axis of density plot
    for(chr in 1:10){
      points(alignDensities[[chr]]$x * 1e6 + transposeByChr[[chr]],
             alignDensities[[chr]]$y * densScale, type = 'l')
      # also add chromosome labels while we are looping
      #if(drawTitle){
        mtext(names(transposeByChr)[chr], side = 1, line = 0, at = mean(transposeByChr[c(chr, chr+1)]),
              cex = par("cex"))
      #}
    }
    #if(drawTitle){
      mtext("UA", side = 1, line = 0, at = transposeByChr["unaligned"], cex = par("cex"))
    #}
    # make title for whole plot
    if(drawTitle){
      trtNm <- traitNames[sapply(traitNames, function(x) length(grep(x, colHeader))>0)]
      if(length(trtNm) > 1) trtNm <- trtNm[3]
#      trtSym <- traitSym[trtNm]
      site <- strsplit(colHeader, paste(".",trtNm, sep = ""))[[1]][1]
      site <- gsub("CHA.CSU", "CSU.UI", site)
      site <- gsub("ZJU.KNU", "KNU.ZJU", site)
      site <- gsub("CHA", "UI", site)
      site <- gsub(".", "+", site, fixed = TRUE)
      yr <- substring(colHeader, nchar(colHeader), nchar(colHeader))
      title(main = paste(traitNamesPrint[trtNm], site, "year", yr))
    }
  }
}

manhattan("NEF.Ccirc.year.3")
manhattan("NEF.Biomass.yield.dry.weight.year.3")

pdf("Fig1_manhattan_190326.pdf", width = 6.5, height = 3, pointsize = 10)
par(mar = c(2.7, 2.7, 0.5, 0.5))
manhattan("HU.NEF.CHA.CSU.KNU.ZJU.Biomass.yield.dry.weight.year.3",
          drawTitle = FALSE)
dev.off()

# draw array of 60 manhattan plots
#tiff("170320manhattan_plots.tiff", width = 7*300, height = 9*300, pointsize = 10, res = 300, 
#     compression = "lzw")
par(mfrow = c(15, 4), mar = c(1.1, 3.1, 1.1, 1.1))#, oma = c(0,3.1,4.1,0))
for(tr in traitNames){
  for(site in c("HU.NEF.CHA.CSU.KNU", "ZJU")){
    for(yr in c(2,3)){
      colHeader <- paste(site, tr, "year", yr, sep = ".")
      if(tr == "Ccirc.Bcirc"){
        colHeader <- gsub("year.", "year", colHeader, fixed = TRUE)
      }
      manhattan(colHeader, drawTitle = FALSE, drawYlab = site == "HU.NEF.CHA.CSU.KNU" && yr == 2 && tr == "Culm.per.area")
#      if(site == "HU.NEF.CHA.CSU.KNU" && yr == 2){
#        mtext(traitSym[tr], side = 2, line = 0.5, at = mean(par("yaxp")[1:2]), outer = TRUE)
#      }
#      if(tr == "Biomass.yield.dry.weight"){
#        mtext()
#      }
    }
  }
}
dev.off()

# draw all manhattan plots for supplemental info
plotsDrawn <- 0
manhattanPerPage <- 4
for(tr in traitNames){
  cat(tr, sep = "\n")
  for(site in c("HU", "NEF", "ZJU", "HU.NEF", "HU.NEF.CHA.CSU", "HU.NEF.CHA.CSU.KNU",
                "HU.NEF.CHA.CSU.KNU.ZJU", "ZJU.KNU")){
    cat(site, sep = "\n")
    for(yr in 2:4){
      colHeader <- paste(site, tr, "year", yr, sep = ".")
      cat(colHeader, sep = "\n")
      if(tr == "Ccirc.Bcirc"){
        colHeader <- gsub("year.", "year", colHeader, fixed = TRUE)
      }
      if(!colHeader %in% c(singleTraitsKeep, multiTraitsKeep)){
        next
      } 
      if(plotsDrawn %% manhattanPerPage == 0){ # start a new file if necessary
#        tiff(paste("manhattan plot images/170320manhattan", plotsDrawn/manhattanPerPage, ".tiff", sep = ""),
#             width = 7*300, height = 9*300, res = 300, pointsize = 12, compression = "lzw")
        par(mfrow = c(manhattanPerPage, 1), mar = c(3.1, 3.1, 4.1, 2.1))
      }
      manhattan(colHeader, drawTitle = TRUE) # draw the manhattan plot

      plotsDrawn <- plotsDrawn + 1 # increment number of plots drawn
      if(plotsDrawn %% manhattanPerPage == 0){ # close file if done
        dev.off()
      }
    }
  }
}
dev.off()

# tally of individuals included in GWAS
IndWBLUP <- dimnames(myBLUP)[[1]][apply(myBLUP, 1, function(x) !all(is.na(x)))] # 569
sum(rowMeans(is.na(myA.EM.Msi$imputed[IndWBLUP,])) < 1)
sum(IndWBLUP %in% dimnames(myA.EM.Msi$imputed)[[1]]) # 568 (no KMS018)
write.csv(data.frame(ind = dimnames(myBLUP)[[1]],
                     inGWAS = ifelse(dimnames(myBLUP)[[1]] %in% IndWBLUP &
                                       dimnames(myBLUP)[[1]] %in% dimnames(myA.EM.Msi$imputed)[[1]],
                     "Y", "N")), file = "170112list of in GWAS.csv")


tdMarkers <- read.csv("SNP annotations/161201markers.csv", stringsAsFactors = FALSE)

# look at SNPs id'd multiple times, using only the traits we are keeping
multiHitSNP <- dimnames(oneSiteSig)[[1]][rowSums(cbind(oneSiteSig[,singleTraitsKeep], 
                                                       multiSiteSig[,multiTraitsKeep]), na.rm = TRUE) > 1]
multiHitSNP <- multiHitSNP[!is.na(multiHitSNP)]
length(multiHitSNP) # 92 SNPs

for(snp in multiHitSNP){
  cat(paste(c(snp,singleTraitsKeep[which(oneSiteSig[snp,singleTraitsKeep])],
           multiTraitsKeep[which(multiSiteSig[snp,multiTraitsKeep])]), collapse = ","),
       sep = "\n", file = "170320multihitSNPs.csv", append = TRUE)
}

# number of sig yield SNPs for year 2 and year 3?
na.omit(dimnames(oneSiteSig)[[1]][rowSums(cbind(oneSiteSig[,grep("Biomass.yield.dry.weight.year.3", singleTraitsKeep, value = TRUE)],
      multiSiteSig[,grep("Biomass.yield.dry.weight.year.3", multiTraitsKeep, value = TRUE)]), na.rm = TRUE) > 0])

na.omit(dimnames(oneSiteSig)[[1]][rowSums(cbind(oneSiteSig[,grep("Biomass.yield.dry.weight.year.2", singleTraitsKeep, value = TRUE)],
                                        multiSiteSig[,grep("Biomass.yield.dry.weight.year.2", multiTraitsKeep, value = TRUE)]), na.rm = TRUE) > 0])


# sorted through them in Excel to find ones that really have different traits/sites/years
SNPsforTable <- read.table("170320SNPsForTable.txt", header = FALSE, stringsAsFactors = FALSE,
                           sep = "\t")[[1]]
#write.csv(tdMarkers[match(gsub(".", "-", SNPsforTable, fixed = TRUE), tdMarkers$Original.marker.name),],
#          file = "170320interesting_SNP_align.csv")

# number of genes covered by all of the SNPs
sum(myalign$Sbicolor.gene.name != "") # 20,611 SNPs
geneNameLength <- sapply(myalign$Sbicolor.gene.name, nchar)
hist(geneNameLength)
sum(geneNameLength == 16) # 17921
sum(geneNameLength == 33) #  2604
sum(geneNameLength == 50) #    68
sum(geneNameLength == 67) #     1
table(geneNameLength)
myalign$Sbicolor.gene.name[geneNameLength == 13] # 16 genes
myalign$Sbicolor.gene.name[geneNameLength == 27] # 1 SNP, 2 genes

17921 + 2*2604 + 3 *68 + 4 + 16 + 2 # 23,355 genes

unlistGenes <- unlist(strsplit(myalign$Sbicolor.gene.name, split = ";"))
length(unique(unlistGenes)) # 14,062 genes
geneTab <- table(unlistGenes)
names(geneTab)[geneTab > 1]

14062/(2*34211)

# number of genes with SNPs actually in them
myalign2 <- read.csv("SNP annotations/170111markers.csv", stringsAsFactors = FALSE)
genes <- myalign2$Sbicolor.gene.name[match(tdMarkers$Original.marker.name,
                                           myalign2$Name)]
geneNameLength <- sapply(genes, nchar)
hist(geneNameLength)
table(geneNameLength)
sum(geneNameLength == 13) #    11
sum(geneNameLength == 16) # 16499
sum(geneNameLength == 33) #   314
sum(geneNameLength == 50) #     1
sum(geneNameLength > 0) # 16,825 SNPs in genes

unlistGenes <- unlist(strsplit(genes, split = ";"))
length(unique(unlistGenes)) # 11,013 genes

# look up significance of SNPs from genes ID'd in Slavov et al.
slavovSNP <- c('PstI.TP15894','PstI.TP353084','NsiI.TP613270','NsiI.TP740457', 'NsiI.TP441912',	'PstI.TP880089',	
               'PstI.TP935379',	'NsiI.TP55659', 'PstI.TP434991', 'PstI.TP1292743')
sum(oneSiteSig[slavovSNP,singleTraitsKeep], na.rm = TRUE)
sum(multiSiteSig[slavovSNP,multiTraitsKeep], na.rm = TRUE)

# Get correlation of SNPs with population structure
myQ$Sichuan <- 1 - rowSums(myQ[,-1])
Qcorr <- matrix(NA, nrow = 46177, ncol = 8,
                dimnames = list(SNPnames, names(myQ)[-1]))
for(snp in 1:46177){
  if(snp %% 1000 == 0){
    print(snp)
  }
  for(pop in 1:8){
    Qcorr[snp,pop] <- cor(myGD[[snp+1]], myQ[[pop+1]], use = "pairwise.complete.obs")
  }
}
hist(Qcorr)
head(Qcorr)
which(Qcorr > 0.8 | Qcorr < -0.8, arr.ind = TRUE)

# Correlation of traits with population structure
TraitCorr <- matrix(NA, nrow = 149+219, ncol = 8, 
                    dimnames = list(c(names(Y.single.site)[-1], names(Y.multi.site)[-1]),
                                    names(myQ)[-1]))
for(tr in 1:219){
  for(pop in 1:8){
    TraitCorr[tr, pop] <- cor(Y.single.site[[tr+1]], myQ[[pop+1]], use = "pairwise.complete.obs")
  }
}
for(tr in 1:149){
  for(pop in 1:8){
    TraitCorr[tr+219, pop] <- cor(Y.multi.site[[tr+1]], myQ[[pop+1]], use = "pairwise.complete.obs")
  }
}
hist(TraitCorr)
#write.csv(TraitCorr, file = "161201TraitQcorrelations.csv")

# read in table of SNPs processed by TagDigger
tdMarkers <- read.csv("SNP annotations/161201markers.csv", stringsAsFactors = FALSE)
# check if excluded ones have significant associations
dim(tdMarkers)
length(SNPnames) # 74 thrown out by TagDigger
tdTrash <- SNPnames[!gsub(".", "-", SNPnames, fixed = TRUE) %in% tdMarkers$Original.marker.name]
dim(oneSiteSig)
sum(oneSiteSig[tdTrash,], na.rm = TRUE)
which(oneSiteSig[tdTrash,], arr.ind = TRUE)
dimnames(oneSiteSig)[[2]][which(oneSiteSig["NsiI.TP510394",])] # NEF.Bcirc.year.3
sum(is.na(oneSiteSig[tdTrash,]))
# get a vector to match rows from the TagDigger output with SNPs in our set
identical(SNPnames, dimnames(oneSiteFDRP)[[1]])
TDrowMatch <- match(tdMarkers$Original.marker.name, gsub(".", "-", SNPnames, fixed = TRUE))
#write.csv(cbind(Qcorr[TDrowMatch,], oneSiteFDRP[TDrowMatch, singleTraitsKeep],
#      multiSiteFDRP[TDrowMatch, multiTraitsKeep]), file = "170320FDR_P_values_and_Q_correlation.csv")

# get major and minor alleles for each SNP
load("../../Seq/SNP set 160106/160322filteredSNPs.RData")
load("160324EMimputedSNP_Msi.RData")
dimnames(myA.EM.Msi$A)[[1]]
SNP_alleles <- lapply(strsplit(tdMarkers$Tag.sequence, split = ""),
                         function(x) {
                           sort(x[c(which(x == "[") + 1, which(x == "]") - 1)])
                         }
)
names(SNP_alleles) <- tdMarkers$Marker.name
tdFreq <- colMeans(datacomb3[dimnames(myA.EM.Msi$A)[[1]],tdMarkers$Original.marker.name], na.rm = TRUE)/2 # Msi only
hist(tdFreq)
whichAlMinor <- ifelse(tdFreq > 0.5, 1,2)
hist(whichAlMinor)
tdMarkers$Minor.allele <- sapply(1:dim(tdMarkers)[1], function(i){
  SNP_alleles[[i]][whichAlMinor[i]]
})
table(tdMarkers$Minor.allele)
# minor allele frequency
tdMarkers$Minor.allele.freq <- tdFreq
tdMarkers$Minor.allele.freq[tdFreq > 0.5] <- 1 - tdFreq[tdFreq > 0.5]
hist(tdMarkers$Minor.allele.freq)

#write.csv(tdMarkers[,c("Marker.name","Minor.allele", "Minor.allele.freq")], file = "170130minoralleles.csv")

# export one sheet of P-values per trait ## need to load alEffects and domEffects functions at end of this file
traitSym <- c("Yld", "CC", "BC", "CC/BC", "CmL", "CmNdN", "IntL", "CmDW",
              "CmV", "CmDW/V", "DBI", "DTI", "TCmN", "RCmN/TCmN", "TCmN/A")
names(traitSym) <- traitNames
sigAny <- rowSums(oneSiteSig[TDrowMatch,singleTraitsKeep], na.rm = TRUE) +
  rowSums(multiSiteSig[TDrowMatch,multiTraitsKeep], na.rm = TRUE) > 0 # 298 SNPs significant for any retained trait
outmatSigAny <- matrix(NA, nrow = sum(sigAny), ncol = length(singleTraitsKeep) + length(multiTraitsKeep),
                       dimnames = list(tdMarkers$Marker.name[sigAny], 
                                       paste("V", 1:(length(singleTraitsKeep) + length(multiTraitsKeep)))))
currcolSigAny <- 1
for(tn in traitNames){
  st <- grep(tn, singleTraitsKeep, value = TRUE)
  mt <- grep(tn, multiTraitsKeep, value = TRUE)
  if(tn %in% c("Ccirc", "Bcirc")){
    st <- st[-grep("Ccirc.Bcirc", st)]
    mt <- mt[-grep("Ccirc.Bcirc", mt)]
  }
  outmat <- matrix(nrow = length(TDrowMatch), ncol = length(st)+length(mt),
                   dimnames = list(tdMarkers$Marker.name, paste("V", 1:(length(st)+length(mt)))))
  currcol <- 1 # track which column we are on
  for(site in c("HU", "NEF", "ZJU")){ # single sites
    st2 <- grep(site, st, value = TRUE)
    yrs <- substring(st2, nchar(st2), nchar(st2))
    nYr <- length(yrs)
    if(nYr == 0) next
    dimnames(outmat)[[2]][currcol:(currcol+nYr-1)] <- paste(traitSym[tn], site, "year", yrs)
    outmat[,currcol:(currcol+nYr-1)] <- oneSiteFDRP[TDrowMatch,st2]
    currcol <- currcol + nYr
  }
  # multi-sites
  myrs <- substring(mt, nchar(mt), nchar(mt))
  msites <- sapply(mt, function(x){
    thisstring <- strsplit(x, tn)[[1]][1]
    thisstring <- substring(thisstring, 1, nchar(thisstring)-1)
    thisstring <- gsub("CHA.CSU", "CSU+UI", thisstring, fixed = TRUE)
    thisstring <- gsub(".", "+", thisstring, fixed = TRUE)
    if(thisstring == "ZJU+KNU") thisstring <- "KNU+ZJU"
    return(thisstring)
    })
  dimnames(outmat)[[2]][currcol:dim(outmat)[2]] <- paste(traitSym[tn], msites, "year", myrs)
  outmat[,currcol:dim(outmat)[2]] <- multiSiteFDRP[TDrowMatch, mt]
  # rows with any significant SNPs
  sigrows <- rowSums(outmat < 0.05, na.rm = TRUE) > 0
  
#  write.csv(cbind(tdMarkers[sigrows,], outmat[sigrows,]), 
#            file = paste("170320", gsub("/", " ", traitSym[tn]), "FDR P-values subset.csv"),
#            na = "")
#  write.csv(cbind(tdMarkers[,], outmat[,]), 
#            file = paste("170320", gsub("/", " ", traitSym[tn]), "FDR P-values.csv"),
#            na = "")
  
  origNames <- tdMarkers$Original.marker.name[sigrows]
  # allele effect magnitudes, in both positive and negative direction
  theseyears <- sapply(c(st,mt), function(x){
    yr <- substring(x, nchar(x), nchar(x))
  })
  thesesites <- sapply(c(st,mt), function(x){
    site <- substring(x, 1, nchar(x) - 7 + (tn == "Ccirc.Bcirc") - nchar(tn) - 1)
  })
  thisorder <- order(factor(thesesites, levels = c("HU", "NEF", "ZJU", "HU.NEF", "ZJU.KNU", "HU.NEF.CHA.CSU.KNU", 
                                                 "HU.NEF.CHA.CSU", "HU.NEF.CHA.CSU.KNU.ZJU")), theseyears)
  theseyears <- theseyears[thisorder]
  thesesites <- thesesites[thisorder]
  maglist <- mapply(alEffects, site = thesesites, year = theseyears,
                    MoreArgs = list(trait = tn, SNPs = gsub("-", ".", origNames)), SIMPLIFY = FALSE)
  domlist <- mapply(domEffects, site = thesesites, year = theseyears,
                    MoreArgs = list(trait = tn, SNPs = gsub("-", ".", origNames)), SIMPLIFY = FALSE)
  # matrix to fill in allele effect magnitudes in the right direction
  effectMat <- matrix(NA, nrow = sum(sigrows), ncol = dim(outmat)[2],
                      dimnames = list(tdMarkers$Marker.name[sigrows], dimnames(outmat)[[2]]))
  domMat <- effectMat # for filling in dominance deviations
  
  # minor alleles, frequencies, and effects
  thesePhen <- cbind(Y.single.site[,st], Y.multi.site[,mt])
  
  associationSign <- rep("mixed", length(origNames))
  for(i in 1:length(origNames)){ # loop through significant markers for this trait
    nm <- origNames[i]
    theseGen <- myGD[,gsub("-", ".", nm)]
    if(whichAlMinor[sigrows][i] == 1) theseGen <- 2 - theseGen # flip so minor allele is positive
    # get just BLUPs with significant association for this marker
    whichPhenSig <- outmat[sigrows,,drop = FALSE][i,] < 0.05
    whichPhenSig[is.na(whichPhenSig)] <- FALSE
#    thisPhen2 <- thesePhen[,whichPhenSig, drop = FALSE]
    thisCov <- sapply(1:dim(thesePhen)[2], function(x) cov(theseGen, thesePhen[,x], use = "complete.obs"))
    if(all(thisCov[whichPhenSig] >= 0)) associationSign[i] <- "+"
    if(all(thisCov[whichPhenSig] <= 0)) associationSign[i] <- "-"
    if(all(thisCov[whichPhenSig] == 0)) associationSign[i] <- "0"
    # fill in effectMat
    for(j in 1:dim(outmat)[2]){
      if(thisCov[j] > 0){
        effectMat[i,j] <- maglist[[j]][i,"Pos.allele.effect"]
      } else {
        effectMat[i,j] <- -1 * maglist[[j]][i, "Neg.allele.effect"]
      }
      domMat[i,j] <- domlist[[j]][i]
    }
  }
  

  outmatSigAny[,(1:dim(outmat)[2])+currcolSigAny-1] <- outmat[sigAny,]
  dimnames(outmatSigAny)[[2]][(1:dim(outmat)[2])+currcolSigAny-1] <- dimnames(outmat)[[2]]
  currcolSigAny <- currcolSigAny + dim(outmat)[2]
  
  dimnames(outmat)[[2]] <- paste("P value", dimnames(outmat)[[2]])
  dimnames(effectMat)[[2]] <- paste("Allelic effect", dimnames(effectMat)[[2]])
  dimnames(domMat)[[2]] <- paste("Dominance deviation", dimnames(domMat)[[2]])
#  write.csv(cbind(tdMarkers[sigrows,,drop = FALSE],
#                  associationSign, outmat[sigrows,, drop = FALSE], effectMat, domMat), 
#            file = paste("170912allele effects and p values ", gsub("/", " ", traitSym[tn]), ".csv", sep = ""),
#            na = "")
}
#write.csv(cbind(tdMarkers[sigAny,], outmatSigAny), file = "170320sig SNPs only FDR P-values.csv",
#          na = "")

# make a table of number of significant hits by trait, site and year
sigTally <- data.frame(Trait = c(rep("Yld", 3), rep(traitSym[-1], each = 2)),
                       Year = c(2:4, rep(c(2,3), times = length(traitNames)-1)),
                       HU = rep(NA, 31), NEF = rep(NA, 31), ZJU = rep(NA,31),
                       HU.NEF = rep(NA,31), HU.NEF.CHA.CSU = rep(NA, 31), 
                       HU.NEF.CHA.CSU.KNU = rep(NA, 31), 
                       HU.NEF.CHA.CSU.KNU.ZJU = rep(NA, 31), ZJU.KNU = rep(NA, 31))
for(tr in traitNames){
  for(site in c("HU", "NEF", "ZJU", "HU.NEF", "HU.NEF.CHA.CSU", "HU.NEF.CHA.CSU.KNU",
                "HU.NEF.CHA.CSU.KNU.ZJU", "ZJU.KNU")){
    isSingle <- site %in% c("HU", "NEF", "ZJU")
    for(yr in c(2,3,4)){
      thisName <- paste(site, tr, "year", yr, sep = ".")
      if(tr == "Ccirc.Bcirc"){
        thisName <- gsub("year.", "year", thisName, fixed = TRUE)
      }
      if((isSingle && !thisName %in% singleTraitsKeep) || 
         (!isSingle && !thisName %in% multiTraitsKeep)){
        next
      }
      thisrow <- which(sigTally$Trait == traitSym[tr] & sigTally$Year == yr)
      if(isSingle){
        sigTally[thisrow, site] <- sum(oneSiteSig[,thisName], na.rm = TRUE)
      } else {
        sigTally[thisrow, site] <- sum(multiSiteSig[,thisName], na.rm = TRUE)
      }
    }
  }
}
#write.csv(sigTally, file = "170320 tally of significant SNPs.csv", na = "", row.names = FALSE)


## Get SNP genotypes for SNPs of interest, for supplementary table  #edit from here
intSNPs <- read.table("170410SNPsfortable.txt", header = FALSE, stringsAsFactors = FALSE)[[1]]
PstIHM <- read.table("../../Seq/SNP set 160106/PstI HapMap/HapMap.hmp.txt",
                     sep = "\t", header = TRUE, row.names = 1)
samCol <- 11:605
PstIintSNPs <- grep("PstI", intSNPs, value = TRUE)

PstIintSNPs <- gsub("PstI.", "", PstIintSNPs)
PstGen <- PstIHM[PstIintSNPs, samCol]
rm(PstIHM)
NsiIintSNPs <- grep("NsiI", intSNPs, value = TRUE)
NsiIintSNPs <- gsub("NsiI.", "", NsiIintSNPs)
NsiIHM <- read.table("../../Seq/SNP set 160106/NsiI HapMap/HapMap.hmp.txt",
                     sep = "\t", header = TRUE, row.names = 1)
NsiIGen <- NsiIHM[NsiIintSNPs, samCol]
rm(NsiIHM)
rm(NsiIGen)
# some of the NsiI genotypes had been merged after running UNEAK
# try getting things from final SNP set
load("../../Seq/SNP set 160106/160322filteredSNPs.RData")
dim(datacomb3) # 46,177 SNPs, 594 ind
dimnames(datacomb3)[[1]]
dim(allOneSiteBLUPs)
dimnames(allOneSiteBLUPs)[[1]][!dimnames(allOneSiteBLUPs)[[1]] %in% dimnames(datacomb3)[[1]]]
  # KMS018 omitted (not done for NsiI)
str(datacomb3[,gsub(".", "-", intSNPs, fixed = TRUE)]) # numeric SNPs, unimputed, subsetted for interesting ones

# get universal names
intSNP_u <- tdMarkers$Marker.name[match(gsub(".", "-", intSNPs, fixed = TRUE),
                                        tdMarkers$Original.marker.name)]
# get alleles
intSNP_alleles <- lapply(strsplit(tdMarkers$Tag.sequence[match(gsub(".", "-", intSNPs, fixed = TRUE),
                             tdMarkers$Original.marker.name)], split = ""),
       function(x) {
         sort(x[c(which(x == "[") + 1, which(x == "]") - 1)])
       }
)
intSNP_gennames <- lapply(intSNP_alleles, function(x){
  paste(x[c(1,1,2)], x[c(1,2,2)], sep = "")
})
# get genotypes as characters
intSNP_chargen <- matrix("", nrow = dim(datacomb3)[1], ncol = length(intSNPs),
                         dimnames = list(dimnames(datacomb3)[[1]], intSNP_u))
for(i in 1:length(intSNPs)){
  intSNP_chargen[,i] <- intSNP_gennames[[i]][datacomb3[,gsub(".", "-", intSNPs[i], fixed = TRUE)]+1]
}
intSNP_chargen[is.na(intSNP_chargen)] <- ""

identical(match(dimnames(datacomb3)[[1]], dimnames(allOneSiteBLUPs)[[1]]),
          sort(match(dimnames(datacomb3)[[1]], dimnames(allOneSiteBLUPs)[[1]]))) # same order, just missing KMS018

#write.csv(intSNP_chargen, file = "170410interesting_SNP_genotypes.csv")
#check
NsiIHM["TP71276",samCol]

# look at culm density SNP, out of curiosity
plot(datacomb3[,"NsiI-TP427434"], 
     allMultiSiteBLUPs[,"HU.NEF.CHA.CSU.CmDens.year.3"][match(dimnames(datacomb3)[[1]], 
                                                          dimnames(allMultiSiteBLUPs)[[1]])])
abline(lm(allMultiSiteBLUPs[,"HU.NEF.CHA.CSU.CmDens.year.3"][match(dimnames(datacomb3)[[1]], 
                                                                   dimnames(allMultiSiteBLUPs)[[1]])] ~
            datacomb3[,"NsiI-TP427434"]), col = "red")
boxplot(allMultiSiteBLUPs[,"HU.NEF.CHA.CSU.CmDens.year.3"][match(dimnames(datacomb3)[[1]], 
                                                                 dimnames(allMultiSiteBLUPs)[[1]])] ~
          datacomb3[,"NsiI-TP427434"])

# did we really export all interesting SNPs?
dim(multiSiteSig)
dimnames(multiSiteSig)[[2]][multiSiteSig["PstI.TP250403",]]
dimnames(oneSiteSig)[[2]][oneSiteSig["PstI.TP250403",]]

# make function to double-check information in multi-hit SNPs table
checkSig <- function(mrkrnum){ # mrkrnum is six-digit number for UIMiscanthus number
  mrkrname <- sprintf("UIMiscanthus%06d", mrkrnum)
  oldname <- tdMarkers$Original.marker.name[match(mrkrname, tdMarkers$Marker.name)]
  oldname <- gsub("-", ".", oldname)
  for(tr in traitNames){
    for(site in c('HU', 'NEF', 'ZJU', 'HU.NEF', 'HU.NEF.CHA.CSU', 'HU.NEF.CHA.CSU.KNU', 'HU.NEF.CHA.CSU.KNU.ZJU',
                  'ZJU.KNU')){
      for(yr in c(2,3)){
        if(tr == "Ccirc.Bcirc"){
          trName <- paste(site, ".", tr, ".year", yr, sep = "")
        } else {
          trName <- paste(site, tr, "year", yr, sep = ".")
        }
        if((trName %in% singleTraitsKeep && !is.na(oneSiteSig[oldname,trName]) && oneSiteSig[oldname,trName]) ||
           (trName %in% multiTraitsKeep && !is.na(multiSiteSig[oldname,trName]) && multiSiteSig[oldname,trName])){
          cat(paste(traitSym[tr], site, yr), sep = "\n")
        }
      }
    }
  }
}

checkSig(022102)

## Estimate allele effects ####
# function for reverse Box-Cox transformation
bcBack <- function(x1, lda){
  if(lda == 0){
    x <-  exp(x1)
  } else {
    x <- (lda * x1 + 1)^(1/lda)
  }
  return(x)
}
# function to estimate allele effects
alEffects <- function(trait, site, year, SNPs){
  if(trait == "Ccirc.Bcirc"){
    trName <- paste(site, ".", trait, ".year", year, sep = "")
  } else {
    trName <- paste(site, trait, "year", year, sep = ".")
  }
  if(!trName %in% c(singleTraitsKeep, multiTraitsKeep)){
    stop(paste("Trait", trName, "not part of GWAS"))
  }
  trNameShort <- substring(trName, nchar(site) + 2, nchar(trName))
  # get median values
  if(trName %in% singleTraitsKeep){
    medVal <- median(mydata[[trNameShort]][mydata$Loc == site], na.rm = TRUE)
  }
  if(trName %in% multiTraitsKeep){
    sites <- strsplit(site, ".", fixed = TRUE)[[1]]
    medVal <- median(mydata[[trNameShort]][mydata$Loc %in% sites], na.rm = TRUE)
  }
  # search for allele effects file
  filename <- paste("GAPIT.", trName, "Allelic_Effect_Estimates.csv", sep = ".")
  if(file.exists(paste("GAPIT/output/", filename, sep = ""))){
    filename <- paste("GAPIT/output/", filename, sep = "")
  }
  # import allele effects
  efftab <- read.csv(filename, stringsAsFactors = FALSE, row.names = 1)
  theseeff <- abs(efftab[SNPs, "Allelic.Effect.Estimate"])
  # get lambda for BC transformation
  thislambda <- BClambdas[match(trNameShort, names(mydata)[traits])]
  # transform the median
  medBC <- bc(medVal, thislambda)
  ## get transformed median value for this site and year
  plusBack <- bcBack(medBC + theseeff, thislambda) - medVal
  minusBack <- medVal - bcBack(medBC - theseeff, thislambda)
  
  return(data.frame(row.names = SNPs, Pos.allele.effect = plusBack,
                    Neg.allele.effect = minusBack, 
                    BC.effect = theseeff))
}

alEffects("Biomass.yield.dry.weight", "HU.NEF.CHA.CSU.KNU.ZJU", 3, 
          dimnames(multiSiteSig)[[1]][multiSiteSig[,"HU.NEF.CHA.CSU.KNU.ZJU.Biomass.yield.dry.weight.year.3"]])
alEffects("Biomass.yield.dry.weight", "ZJU", 3, 
          dimnames(multiSiteSig)[[1]][multiSiteSig[,"HU.NEF.CHA.CSU.KNU.ZJU.Biomass.yield.dry.weight.year.3"]])

alEffects("Biomass.yield.dry.weight", "NEF", 3, 
          paste("PstI.TP", c(237577,395258,582235,645986,669168), sep = ""))

# function to estimate dominance effects
domEffects <- function(trait, site, year, SNPs){
  if(trait == "Ccirc.Bcirc"){
    trName <- paste(site, ".", trait, ".year", year, sep = "")
  } else {
    trName <- paste(site, trait, "year", year, sep = ".")
  }
  if(!trName %in% c(singleTraitsKeep, multiTraitsKeep)){
    stop(paste("Trait", trName, "not part of GWAS"))
  }
  trNameShort <- substring(trName, nchar(site) + 2, nchar(trName))
  if(trName %in% names(Y.single.site)){
    thisBLUP <- Y.single.site[[trName]]
  } else {
    thisBLUP <- Y.multi.site[[trName]]
  }
  # get lambda for BC transformation
  thislambda <- BClambdas[match(trNameShort, names(mydata)[traits])]
  # get median values
  if(trName %in% singleTraitsKeep){
    medVal <- median(mydata[[trNameShort]][mydata$Loc == site], na.rm = TRUE)
  }
  if(trName %in% multiTraitsKeep){
    sites <- strsplit(site, ".", fixed = TRUE)[[1]]
    medVal <- median(mydata[[trNameShort]][mydata$Loc %in% sites], na.rm = TRUE)
  }
  # transform the median
  medBC <- bc(medVal, thislambda)
  
  outdom <- numeric(length(SNPs))
  names(outdom) <- SNPs
  for(snp in SNPs){
    thesegen <- myGD[[snp]]
    zeros <- sapply(thesegen, function(x) isTRUE(all.equal(x, 0)))
    ones <- sapply(thesegen, function(x) isTRUE(all.equal(x, 1)))
    twos <- sapply(thesegen, function(x) isTRUE(all.equal(x, 2)))
    # observed values for homozygotes
    if(mean(thesegen) > 1){ # 0 is rare allele
      phenHomRare <- mean(thisBLUP[zeros], na.rm = TRUE) + medBC
      phenHomCommon <- mean(thisBLUP[twos], na.rm = TRUE) + medBC
    } else { # 2 is the rare allele
      phenHomRare <- mean(thisBLUP[twos], na.rm = TRUE) + medBC
      phenHomCommon <- mean(thisBLUP[zeros], na.rm = TRUE) + medBC
    }
    # observed and expected values for heterozygotes
    phenHet <- mean(thisBLUP[ones], na.rm = TRUE) + medBC
    alEff <- (phenHomRare - phenHomCommon)/2
    expPhenHet <- phenHomCommon + alEff #(phenHomRare + phenHomCommon)/2
#    print(snp)
#    print(c(phenHomRare, phenHomCommon, phenHet, expPhenHet))
    # skip if we don't have all three genotypic categories
    if(any(is.na(c(phenHomRare, phenHomCommon, phenHet)))){
      outdom[snp] <- NA
      next
    }
    # back-transform to original units
    phenHetBack <- bcBack(phenHet, thislambda)
    expPhenHetBack <- bcBack(expPhenHet, thislambda)
    # record difference between observed and expected values for hets
    outdom[snp] <- phenHetBack - expPhenHetBack
  }
  return(outdom)
}

# test out domEffects; get sixsiteyieldSNPs vector from code below
domEffects("Biomass.yield.dry.weight", "HU.NEF.CHA.CSU.KNU.ZJU", 3, sixsiteyieldSNPs)

## look at unusual allele effects estimates ####
plot(myGD$PstI.TP88378, Y.single.site$ZJU.Biomass.yield.dry.weight.year.3,
     xlim = c(0,2), xlab = "Genotype (including imputed)", ylab = "Phenotype") # imputed ones driving association
grid()
abline(lm(Y.single.site$ZJU.Biomass.yield.dry.weight.year.3 ~ myGD$PstI.TP88378), col = "red")
hist(myGD$PstI.TP88378)
plot(myGD$PstI.TP1237941, Y.single.site$ZJU.Biomass.yield.dry.weight.year.3) # one het driving association
grid()
hist(myGD$PstI.TP1237941)
plot(myGD$NsiI.TP380536, Y.single.site$ZJU.Ccirc.year.2)

mean(myGD$PstI.TP88378[!is.na(Y.single.site$ZJU.Biomass.yield.dry.weight.year.3)])

plot(myGD$NsiI.TP616411, Y.multi.site$HU.NEF.CHA.CSU.KNU.ZJU.Biomass.yield.dry.weight.year.3,
     xlim = c(0,2))
plot(myGD$NsiI.TP616411, Y.multi.site$HU.NEF.CHA.CSU.KNU.ZJU.CmN.year.3,
     xlim = c(0,2))

plot(myGD$NsiI.TP360284, Y.multi.site$ZJU.KNU.Biomass.yield.dry.weight.year.2,
     xlim = c(0,2))

plot(myGD$PstI.TP213479, Y.multi.site$HU.NEF.CHA.CSU.KNU.ZJU.CmD_I1.year.3) # big culm diameter hit, ACA9 gene
plot(myGD$PstI.TP213479, Y.multi.site$HU.NEF.CHA.CSU.KNU.ZJU.Culm.per.area.year.3)
# individuals with minor allele of this gene
minAl092590 <- as.character(myGD$Taxa[which(myGD$PstI.TP213479 < 1.1)]) # hits multiple genetic groups
myQ[match(minAl092590, as.character(myQ$Taxa)),]

# SNPs significant for yield at the 6 combined sites in year 3, and which taxa where driving
sixsiteyieldSNPs <- dimnames(multiSiteSig)[[1]][which(multiSiteSig[,"HU.NEF.CHA.CSU.KNU.ZJU.Biomass.yield.dry.weight.year.3"])]
sixsiteyieldSNPs
mydapc <- read.csv("161213dapcgroups.csv", stringsAsFactors = FALSE, row.names = 1)
plotcol <- c("yellow", "blue", "red", "orange", "darkgreen", "purple",
             "cyan", "yellow3", "yellow4")[mydapc[as.character(Y.multi.site$Taxa),1]]
pdf("170630yieldSNPsPlots.pdf")
for(snp in sixsiteyieldSNPs){
  plot(myGD[[snp]], Y.multi.site$HU.NEF.CHA.CSU.KNU.ZJU.Biomass.yield.dry.weight.year.3,
       xlab = "Genotype", ylab = "Year 3 yield BLUP across all sites", main = snp, col = plotcol)
  if(mean(myGD[[snp]], na.rm = TRUE) < 1){
    indRareAl <- which(myGD[[snp]] >= 1)
    thispos <- 2
  } else {
    indRareAl <- which(myGD[[snp]] <= 1)
    thispos <- 4
  }
  text(myGD[[snp]][indRareAl], Y.multi.site$HU.NEF.CHA.CSU.KNU.ZJU.Biomass.yield.dry.weight.year.3[indRareAl],
       Y.multi.site$Taxa[indRareAl], pos = thispos)
}
dev.off()
  # yes, these look real, are driven by a variety of different lines (although mostly PMS)