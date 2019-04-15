# estimate Q_ST
library(lme4)
library(qqman)

# get data for random effects models, from before
load("170206phenotype_analysis.RData")

identical(rownames(allOneSiteBLUPs), rownames(mydapc))

# is variance among BLUPs the same as genetic variance? (no)
#var(allOneSiteBLUPs[mydapc$DAPC_cluster == 5, "NEF.Biomass.yield.dry.weight.year.3"], na.rm = TRUE) # 37.89

#mySubset <- which(mydata$dapcWithHyb == 5 & mydata$Loc == "NEF" & Year3OK)
#myPhen <- bc(mydata$Biomass.yield.dry.weight.year.3[mySubset], BClambdas[2])
#myRep <- as.factor(mydata$Rep[mySubset])
#myGen <- as.factor(mydata$Sample.name[mySubset])

#myModel <- lmer(myPhen ~ (1|myRep) + (1|myGen))
#summary(myModel) # genetic variance is 35.71

# Qst function ####
# subset must be logical
# group1 and 2 are group numbers to test
# phen is a vector of the phenotypes
Qst <- function(phen, group1, group2, subset, groups = mydata$dapcWithHyb,
                reps = mydata$Rep, locs = mydata$Loc, genotypes = mydata$Sample.name){
  # subsets for DAPC groups
  subset1 <- subset & groups == group1 & !is.na(phen)
  subset2 <- subset & groups == group2 & !is.na(phen)
  subset1[is.na(subset1)] <- FALSE
  subset2[is.na(subset2)] <- FALSE
  # loci with data for both groups
  locsCommon <- intersect(unique(locs[subset1]), unique(locs[subset2]))
  # quit if this was only measured at one site
  if(length(locsCommon) <= 1) return(NA)
  # further subset by locations that have both genetic groups, and convert to numeric index
  subset1 <- which(subset1 & locs %in% locsCommon)
  subset2 <- which(subset2 & locs %in% locsCommon)
  # get combined subset with both genetic groups
  subsetA <- c(subset1, subset2)
  
  # make vectors for model
  phenA <- phen[subsetA]
  repsA <- factor(reps[subsetA])
  locsA <- factor(locs[subsetA])
  gensA <- factor(genotypes[subsetA])
  groupsA <- factor(groups[subsetA])
  # run model
  modelA <- lmer(phenA ~ (1|locsA/repsA) + (1|groupsA:locsA) + (1|groupsA/gensA) + (1|gensA:locsA))
  # extract genetic variance
  varGroup <- as.data.frame(VarCorr(modelA))$vcov[6]
  varGeno <- as.data.frame(VarCorr(modelA))$vcov[2]

  # return Qst
  return(varGroup/(varGroup + varGeno))
}

# see http://www.g3journal.org/content/ggg/2/11/1427.full.pdf for method

# test Qst function ####
test <- Qst(bc(mydata$Biomass.yield.dry.weight.year.3, 0.4), 
    group1 = 3, group2 = 5, subset = Year3OK & temperateSubset)
test <- Qst(bc(mydata$Biomass.yield.dry.weight.year.3, 0.4), 
            group1 = 2, group2 = 5, subset = Year3OK & temperateSubset)
Qst(bc(mydata$Biomass.yield.dry.weight.year.3, 0.4), 
    group1 = 2, group2 = 3, subset = Year3OK & temperateSubset)
Qst(bc(mydata$Biomass.yield.dry.weight.year.3, 0.4), 
    group1 = 2, group2 = 3, subset = Year3OK)

# Run pairwise Qst on everything ####
genGroupNums <- c(1:6, 8, 9)
geneticGroups <- as.character(mydapc$DAPC_group[match(genGroupNums, mydapc$DAPC_cluster)])

qstAll <- qstTemperate <- array(NA_real_, dim = c(8, 8, length(traits)),
                          dimnames = list(geneticGroups, geneticGroups, names(mydata)[traits]))
names(BClambdas) <- names(mydata)[traits]

temperateSubset <- mydata$Loc != "ZJU"

for(tr in names(mydata)[traits]){
  print(tr)
  thisphen <- bc(mydata[[tr]], BClambdas[tr]) # Box-Cox transform
  if(match(tr, names(mydata)) %in% year2traits){
    yrSub <- Year2OK 
  } else {
    yrSub <- Year3OK
  }
  
  for(g1 in 1:7){
    for(g2 in (g1+1):8){
      qstAll[g1, g2, tr] <- Qst(thisphen, genGroupNums[g1], genGroupNums[g2],
                                subset = yrSub)
      qstTemperate[g1, g2, tr] <- Qst(thisphen, genGroupNums[g1], genGroupNums[g2],
                                      subset = yrSub & temperateSubset)
    }
  }
}

# Export ####
#outconTemperate <- file("180411qst_temperate.csv", open = "wt")
for(i in 1:dim(qstTemperate)[3]){
  cat(paste(dimnames(qstTemperate)[[3]][i], ",", sep = ""), sep = "\n", file = outconTemperate)
  write.csv(qstTemperate[,,i], file = outconTemperate)
}
close(outconTemperate)

#outconAll <- file("180411qst_all.csv", open = "wt")
for(i in 1:dim(qstAll)[3]){
  cat(paste(dimnames(qstAll)[[3]][i], ",", sep = ""), sep = "\n", file = outconAll)
  write.csv(qstAll[,,i], file = outconAll)
}
close(outconAll)

# Try BLUPs with effect of genetic group removed ####
mydata$Rep.fact <- factor(mydata$Rep)
mydata$DAPC.fact <- factor(mydata$dapcWithHyb)
yldModel <- lmer(Biomass.yield.dry.weight.year.3 ~ (1|Loc/Rep.fact) + (1|DAPC.fact:Loc) + (1|DAPC.fact/Sample.name) + (1|Sample.name:Loc),
                 data = mydata, subset = Year3OK & temperateSubset & MsiRows)
summary(yldModel)

str(ranef(yldModel))
newYldBLUPs <- ranef(yldModel)$`Sample.name:DAPC.fact`
newYldBLUPs$Sample.name <- gsub(":[[:digit:]]$", "", rownames(newYldBLUPs))
newYldBLUPs$DAPC <- as.integer(substring(rownames(newYldBLUPs), nchar(rownames(newYldBLUPs)), nchar(rownames(newYldBLUPs))))
head(newYldBLUPs)
newYldBLUPs$OldBLUP <- allMultiSiteBLUPs[newYldBLUPs$Sample.name, "HU.NEF.CHA.CSU.KNU.Biomass.yield.dry.weight.year.3"]

dapccol <- c("yellow", "blue", "red", "orange", "darkgreen", "purple",
             "cyan", "yellow3", "yellow4")
plot(newYldBLUPs[["OldBLUP"]], newYldBLUPs[[1]], col = dapccol[newYldBLUPs$DAPC],
     xlab = "Original BLUP", ylab = "BLUP with genetic group effect removed", pch = 16)

# genomic prediction from the two versions
library(rrBLUP)
load("160324EMimputedSNP_Msi.RData")
thisKfold <- 10
newYldBLUPs <- newYldBLUPs[newYldBLUPs$Sample.name != "KMS018",]

predValuesNew <- predValuesOld <- matrix(NA_real_, nrow = nrow(newYldBLUPs), ncol = 100,
                                         dimnames = list(newYldBLUPs$Sample.name, NULL))

for(j in 1:100){
  print(j)
  scrambleInd <- sample(nrow(newYldBLUPs)) # random individual order for this iteration
  genPredOutputNew <- list()
  length(genPredOutputNew) <- thisKfold # you can save some processing time by pre-specifying the length of a vector or list
  genPredOutputOld <- genPredOutputNew
  nPerRep <- floor(nrow(newYldBLUPs)/thisKfold)
  
  for(i in 1:thisKfold){
    # identify individuals for training and prediction sets
    firstind <- (i-1) * nPerRep + 1
    if(i == thisKfold){
      lastind <- nrow(newYldBLUPs)
    } else {
      lastind <- i * nPerRep
    }
    train <- scrambleInd[-(firstind:lastind)]
    pred <- scrambleInd[firstind:lastind]
    
    # phenotypes for training set only
    thisphenNew <- thisphenOld <- rep(NA, nrow(newYldBLUPs))
    thisphenNew[train] <- newYldBLUPs[train,1]
    thisphenOld[train] <- newYldBLUPs[train,"OldBLUP"]
    genPredOutputNew[[i]] <- kin.blup(data = data.frame(pheno = thisphenNew, geno = newYldBLUPs$Sample.name),
                                   geno = "geno", pheno = "pheno", GAUSS = FALSE,
                                   K = myA.EM.Msi$A[newYldBLUPs$Sample.name, newYldBLUPs$Sample.name])
    genPredOutputOld[[i]] <- kin.blup(data = data.frame(pheno = thisphenOld, geno = newYldBLUPs$Sample.name),
                                                        geno = "geno", pheno = "pheno", GAUSS = FALSE,
                                                        K = myA.EM.Msi$A[newYldBLUPs$Sample.name, newYldBLUPs$Sample.name])
                                   
    predValuesNew[pred,j] <- genPredOutputNew[[i]]$g[pred]
    predValuesOld[pred,j] <- genPredOutputOld[[i]]$g[pred]
  }
}

# evaluate results
newYldBLUPs$PredNew <- rowMeans(predValuesNew)
newYldBLUPs$PredOld <- rowMeans(predValuesOld)
#save(newYldBLUPs, file = "180411yield_predictions.RData")
load("180411yield_predictions.RData")

# tiff("180413yield_pred_fig.tiff", width = 6.5*300, height = 4*300, res = 300, pointsize = 9,
#      compression = "lzw")
layout(matrix(c(1,3,2,3), nrow = 2, ncol = 2), heights = c(3.5, 0.5))
par(mgp = c(2, 0.5, 0), mar = c(3.1, 3.1, 3.1, 1.1))
plot(newYldBLUPs$OldBLUP, newYldBLUPs$PredOld, bg = dapccol[newYldBLUPs$DAPC],
     xlab = "Genotypic BLUP", ylab = "Predicted breeding value of genotype",
     pch = 21, cex = 1.2)
text(-12, 4.5, expression(italic(R)^2 == 0.47))
title(main = "A", adj = 0)
plot(newYldBLUPs[[1]], newYldBLUPs$PredNew, bg = dapccol[newYldBLUPs$DAPC],
     xlab = "Genotype-within-genetic-group BLUP", ylab = "Predicted breeding value within genetic group",
     pch = 21, cex = 1.2)
text(-500, 350, expression(italic(R)^2 == 0.36))
title(main = "B", adj = 0)
# color legend
par(mar = rep(0, 4))
plot(NA, axes = FALSE, col = "white", xlim = c(0,16), ylim = c(0, 9))
legtext <- c("S Japan", "N Japan", "Korea/N China", "Sichuan", "Yangtze-Qinling",
             "SE China/tropical", "Ornamentals", "US naturalized")
thisy <- 6 - ((((1:8)+1) %% 2) * 3)
spc <- 1.25
xcoord <- strwidth("Genetic groups:") * spc
xcoord <- c(xcoord, xcoord[1] + strwidth("N Japan  ") * spc)
xcoord <- c(xcoord, xcoord[2] + strwidth("Korea/N China") * spc)
xcoord <- c(xcoord, xcoord[3] + strwidth("SE China/tropical") * spc)
thisx <- rep(xcoord, each = 2)
points(thisx, thisy, pch = 21, cex = 2, bg = dapccol[-7])
text(thisx, thisy, legtext, pos = 4)
text(0, 6, "Genetic groups:", pos = 4)

dev.off()

cor(newYldBLUPs[[1]], newYldBLUPs$PredNew) # 0.355
cor(newYldBLUPs$OldBLUP, newYldBLUPs$PredOld) # 0.478

# how much of BLUP variation vs. GEBV variation is explained by genetic groups
summary(lm(newYldBLUPs$OldBLUP ~ newYldBLUPs$DAPC)) # R-squared of 0.022
summary(lm(newYldBLUPs$PredOld ~ newYldBLUPs$DAPC)) # R-squared of 0.056
# when variance due to genetic group is removed
summary(lm(newYldBLUPs[[1]] ~ newYldBLUPs$DAPC))    # R-squared of -0.002, non-significant
summary(lm(newYldBLUPs$PredNew ~ newYldBLUPs$DAPC)) # R-squared of 0.002, non-significant

# get genomic prediction accuracies for all traits shown in main manuscript, using genotype-within-genetic-group BLUPs
year3traits <- traits[c(2,5,8,11,14,17,20,23,25,28,30,33,36,39,41)]
names(mydata)[year3traits]
myaccuracies <- matrix(NA_real_, nrow = length(year3traits), ncol = 8,
                       dimnames = list(names(mydata)[year3traits], 
                                       c("HU", "NEF", "CSU", "CHA", "KNU", "ZJU", "HU.NEF.CHA.CSU.KNU", "HU.NEF.CHA.CSU.KNU.ZJU")))

# function to get mean prediction accuracy across all reps
predacc <- function(blup, predmat){
  mean(sapply(1:ncol(predmat), 
              function(i) cor(blup, predmat[,i])))
}

for(tr in names(mydata)[year3traits]){
  print(tr)
  for(sc in colnames(myaccuracies)){
    print(sc)
    # load and skip compuations if already done
    thisfile <- paste("180411newblups/", sc, ".", tr, ".RData", sep = "")
    if(file.exists(thisfile)){
      load(thisfile)
      myaccuracies[tr, sc] <- predacc(myblups[[1]], predvalues[,i])
      next
    }
    sites <- strsplit(sc, "\\.")[[1]] # locations for this dataset
    # subset of rows in my data to use
    thissubset <- Year3OK & MsiRows & mydata$Loc %in% sites & !is.na(mydata[[tr]])
    # skip if trait missing
    if(sum(thissubset) == 0) next
    # vectors for model
    phen <- mydata[[tr]][thissubset]
    groups <- factor(mydata$DAPC.fact[thissubset])
    locs <- factor(mydata$Loc[thissubset])
    reps <- factor(mydata$Rep.fact[thissubset])
    geno <- factor(mydata$Sample.name[thissubset])
    # skip if not all the locs are there
    if(length(levels(locs)) != length(sites)) next
    # model
    if(length(sites) > 1){
      mymodel <- lmer(phen ~ (1|locs/reps) + (1|groups:locs) + (1|groups/geno) + (1|geno:locs))
    } else {
      mymodel <- lmer(phen ~ (1|reps) + (1|groups/geno))
    }
    
    # extract BLUPs
    myblups <- ranef(mymodel)$`geno:groups`
    myblups$Sample.name <- gsub(":[[:digit:]]$", "", rownames(myblups))
    myblups <- myblups[myblups$Sample.name != "KMS018",]
    
    # 100 reps of 10-fold cross validation
    predvalues <- matrix(NA_real_, nrow = nrow(myblups), ncol = 100,
                         dimnames = list(myblups$Sample.name, NULL))
    for(j in 1:100){
      print(j)
      scrambleInd <- sample(nrow(myblups)) # random individual order for this iteration
      genPredOutput <- list()
      length(genPredOutput) <- thisKfold 
      nPerRep <- floor(nrow(myblups)/thisKfold)
      
      for(i in 1:thisKfold){
        # identify individuals for training and prediction sets
        firstind <- (i-1) * nPerRep + 1
        if(i == thisKfold){
          lastind <- nrow(myblups)
        } else {
          lastind <- i * nPerRep
        }
        train <- scrambleInd[-(firstind:lastind)]
        pred <- scrambleInd[firstind:lastind]
        
        # phenotypes for training set only
        thisphen <- rep(NA, nrow(myblups))
        thisphen[train] <- myblups[train,1]
        genPredOutput[[i]] <- kin.blup(data = data.frame(pheno = thisphen, geno = myblups$Sample.name),
                                          geno = "geno", pheno = "pheno", GAUSS = FALSE,
                                          K = myA.EM.Msi$A[myblups$Sample.name, myblups$Sample.name])

        predvalues[pred,j] <- genPredOutput[[i]]$g[pred]
      }
    }
    save(myblups, predvalues, file = thisfile)
    # prediction accuracy
    myaccuracies[tr, sc] <- predacc(myblups[[1]], predvalues)
  }
}
#save(myaccuracies, file = "180411newblups/accuracies190221.RData")
#write.csv(myaccuracies, "190221pred_accuracies_within_group.csv")

# look at relationship between BLUP and GEBV for Bcirc at NEF and ZJU
load("180411newblups/NEF.Bcirc.year.3.RData")
plot(myblups[[1]], rowMeans(predvalues))
plot(myblups[[1]], rowMeans(predvalues), xlim = c(-150,150), ylim = c(-2, 2))

load("180411newblups/ZJU.Bcirc.year.3.RData")
plot(myblups[[1]], rowMeans(predvalues))

#### export all BLUPs and GEBVs to spreadsheet ####
# set up trait, site, and file names
blupfiles <- list.files("180411newblups", pattern = "year\\.?3\\.RData$",
                        full.names = TRUE)
traitSym <- c("Yld", "CC", "BC", "CC/BC", "CmL", "CmNdN", "IntL", "CmDW",
              "CmV", "CmDW/V", "DBI", "DTI", "TCmN", "RCmN/TCmN", "TCmN/A")
traitNames <- c("Biomass.yield.dry.weight", "Ccirc", "Bcirc", "Ccirc.Bcirc",
                "CmL", "CmN", "IntL", "CmDW", "CmVol", "CmDens", "CmD_I1", "CmD_LI", "Total.culms",
                "Prop.repro.culms", "Culm.per.area")
site_order <- c("HU", "NEF", "CSU", "CHA", "KNU", "ZJU",
                "HU.NEF.CHA.CSU.KNU", "HU.NEF.CHA.CSU.KNU.ZJU")

trait_out_tab <- data.frame(site = rep(site_order, times = length(traitNames)),
                            traitName = rep(traitNames, each = length(site_order)),
                            traitSym = rep(traitSym, each = length(site_order)),
                            stringsAsFactors = FALSE)
trait_out_tab$fileName <- paste("180411newblups/", trait_out_tab$site, ".",
                                trait_out_tab$traitName, ".year.3.RData", sep = "")
trait_out_tab$fileName[trait_out_tab$traitSym == "CC/BC"] <- 
  gsub("year\\.3", "year3", trait_out_tab$fileName[trait_out_tab$traitSym == "CC/BC"])
trait_out_tab$colHeader <- paste(trait_out_tab$traitSym,
                                 gsub("CHA", "UI", trait_out_tab$site), "year 3")

trait_out_tab <- trait_out_tab[trait_out_tab$fileName %in% blupfiles, ]

# matrices of BLUPs and GEBVs
blupout <- matrix(NA_real_, nrow = nrow(allOneSiteBLUPs), ncol = nrow(trait_out_tab),
                  dimnames = list(rownames(allOneSiteBLUPs),
                                  trait_out_tab$colHeader))
gebvout <- blupout
gebv_se_out <- blupout

# loop through BLUPs and GEBVs and fill in
for(i in 1:nrow(trait_out_tab)){
  load(trait_out_tab$fileName[i])
  theseind <- row.names(predvalues)
  stopifnot(identical(theseind, myblups$Sample.name))
  
  blupout[theseind, trait_out_tab$colHeader[i]] <- myblups[[1]]
  gebvout[theseind, trait_out_tab$colHeader[i]] <- rowMeans(predvalues)
  gebv_se_out[theseind, trait_out_tab$colHeader[i]] <-
    apply(predvalues, 1, sd)/sqrt(ncol(predvalues))
}

blupout[1:20,1:5]
gebvout[1:20,1:5]
plot(blupout[,1], gebvout[,1])
plot(blupout[,50], gebvout[,50])
plot(blupout[,90], gebvout[,90])
plot(blupout[,95], gebvout[,95])

#write.csv(blupout, "190218eqn34blups.csv", na = "")

gebv_comb_out <- matrix(NA_real_, nrow = nrow(gebvout), ncol = 2 * ncol(gebvout),
                        dimnames = list(rownames(gebvout),
                                        paste(rep(colnames(gebvout), each = 2),
                                              c("Breeding value", "Standard error"))))
gebv_comb_out[,(1:ncol(gebvout)) * 2 - 1] <- gebvout
gebv_comb_out[,(1:ncol(gebvout)) * 2] <- gebv_se_out

#write.csv(gebv_comb_out, "190218eqn34gebvs.csv", na = "")

#### Genomic prediction within genetic groups ####
# lump groups to increase numbers within groups
# 11/15 don't lump Sichuan w/ SE China
mydapc$consolidated_groups <- as.character(mydapc$DAPC_group)
mydapc$consolidated_groups[mydapc$consolidated_groups == "Ornamental"] <- "S Japan"
mydapc$consolidated_groups[mydapc$consolidated_groups == "US naturalized"] <- "S Japan"
table(mydapc$consolidated_groups)
cons_groups <- c("S Japan", "N Japan", "Korea/N China", "Yangtze-Qinling", "SE China/tropical")
identical(rownames(mydapc), rownames(allOneSiteBLUPs))

# set up matrix for holding accuracies
siteCombos <- c("HU", "NEF", "CSU", "CHA", "KNU", "ZJU",
                "HU.NEF.CHA.CSU.KNU", "HU.NEF.CHA.CSU.KNU.ZJU")
accuracies <- matrix(NA_real_, nrow = length(year3traits), 
                     ncol = length(siteCombos) * length(cons_groups),
                     dimnames = list(names(mydata[year3traits]),
                                     paste(rep(siteCombos, each = length(cons_groups)),
                                       rep(cons_groups, times = length(siteCombos)))))
# set up matrix to hold predicted values
predvalues_mean <- matrix(NA_real_, nrow = nrow(allOneSiteBLUPs),
                     ncol = length(siteCombos) * length(year3traits),
                     dimnames = list(rownames(allOneSiteBLUPs),
                                     paste(rep(siteCombos, each = length(year3traits)),
                                           rep(names(mydata)[year3traits], 
                                               times = length(siteCombos)),
                                           sep = ".")))
predvalues_se <- predvalues_mean

# loop for genomic prediction
thisKfold <- 10
for(tr in names(mydata)[year3traits]){
  print(tr)
  for(sc in siteCombos){
    print(sc)
    # file for saving output
    thisfile <- paste("181115select_within_group/", sc, ".", tr, ".RData", sep = "")

    # extract BLUPs
    blupname <- paste(sc, tr, sep = ".")
    if(blupname %in% colnames(allOneSiteBLUPs)){
      theseblups <- allOneSiteBLUPs[,blupname]
    } else if(blupname %in% colnames(allMultiSiteBLUPs)){
      theseblups <- allMultiSiteBLUPs[,blupname]
    } else {
      next
    }
    
    # 100 reps of 10-fold cross validation
    predvalues <- matrix(NA_real_, nrow = nrow(predvalues_mean), ncol = 100,
                         dimnames = list(rownames(predvalues_mean), NULL))
    # loop through groups
    for(gp in cons_groups){
      # BLUPs just within this group
      grpblup <- theseblups[mydapc$consolidated_groups == gp]
      grpblup <- grpblup[!is.na(grpblup) & names(grpblup) != "KMS018"]
      if(length(grpblup) < 50) next # skip if too few individuals with data
      print(gp)
      # skip calculations if already done
      if(file.exists(thisfile)){
        load(thisfile)
        accuracies[tr, paste(sc, gp)] <- predacc(grpblup, predvalues[names(grpblup),])
        next
      }
      # loop through reps
      for(j in 1:100){
        print(j)
        scrambleInd <- sample(names(grpblup)) # random individual order for this iteration
        genPredOutput <- list()
        length(genPredOutput) <- thisKfold 
        nPerRep <- floor(length(grpblup)/thisKfold)
        
        for(i in 1:thisKfold){
          # identify individuals for training and prediction sets
          firstind <- (i-1) * nPerRep + 1
          if(i == thisKfold){
            lastind <- length(grpblup)
          } else {
            lastind <- i * nPerRep
          }
          train <- scrambleInd[-(firstind:lastind)]
          pred <- scrambleInd[firstind:lastind]
          
          # phenotypes for training set only
          thisphen <- rep(NA, length(grpblup))
          names(thisphen) <- names(grpblup)
          thisphen[train] <- grpblup[train]
          genPredOutput[[i]] <- kin.blup(data = data.frame(pheno = thisphen, geno = names(thisphen)),
                                         geno = "geno", pheno = "pheno", GAUSS = FALSE,
                                         K = myA.EM.Msi$A[names(thisphen), names(thisphen)])
          
          predvalues[pred,j] <- genPredOutput[[i]]$g[pred]
        } # end of loop through k-fold
      } # end of loop through 100 reps
      
      # prediction accuracy
      accuracies[tr, paste(sc, gp)] <- predacc(grpblup, predvalues[names(grpblup),])
      
    } # end of loop through genetic groups
    
#    save(predvalues, file = thisfile)
    predvalues_mean[,paste(sc, tr, sep = ".")] <- rowMeans(predvalues)
    predvalues_se[,paste(sc, tr, sep = ".")] <- apply(predvalues, 1, sd)/sqrt(ncol(predvalues))
    
  } # end of loop through sites
} # end of loop through traits

#save(predvalues_mean, predvalues_se, accuracies, 
#     file = "181115select_within_group/accuracies190221.RData")

#write.csv(accuracies, file = "181115select_within_group/190221accuracies.csv")

# how many individuals for each trait?
colSums(!is.na(predvalues_mean))
for(sc in siteCombos){
  for(gp in cons_groups){
    n <- sum(!is.na(predvalues_mean[mydapc$consolidated_groups == gp,
                                    paste(sc, "Ccirc.year.3", sep = ".")]))
    cat(paste(sc, gp, n), sep = "\n")
  }
}


# S Japan - is prediction based on wild vs. ornamental?
SJapan <- which(mydapc$consolidated_groups == "S Japan")
pdf("180918_SJapan_northern_sites_yield_accuracy.pdf")
plot(allMultiSiteBLUPs[SJapan,"HU.NEF.CHA.CSU.KNU.Biomass.yield.dry.weight.year.3"],
     predvalues_mean[SJapan, "HU.NEF.CHA.CSU.KNU.Biomass.yield.dry.weight.year.3"],
     pch = 21, bg = c(`S Japan` = "yellow", Ornamental = "yellow3", 
                      `US naturalized` = "yellow4")[as.character(mydapc$DAPC_group[SJapan])],
     xlab = "Yield BLUPs", ylab = "Yield GEBVs")
legend(-19, 6, c("S Japan", "Ornamental", "US naturalized"),
       pch = 21, pt.bg = c("yellow", "yellow3", "yellow4"))
dev.off()

## Export GEBVs to spreadsheet
traitSym <- c("Yld", "CC", "BC", "CC/BC", "CmL", "CmNdN", "IntL", "CmDW",
              "CmV", "CmDW/V", "DBI", "DTI", "TCmN", "RCmN/TCmN", "TCmN/A")
traitNames <- c("Biomass.yield.dry.weight", "Ccirc", "Bcirc", "Ccirc.Bcirc",
                "CmL", "CmN", "IntL", "CmDW", "CmVol", "CmDens", "CmD_I1", "CmD_LI", "Total.culms",
                "Prop.repro.culms", "Culm.per.area")
site_order <- c("HU", "NEF", "CSU", "CHA", "KNU", "ZJU",
                "HU.NEF.CHA.CSU.KNU", "HU.NEF.CHA.CSU.KNU.ZJU")
site_names <- c("HU", "NEF", "CSU", "UI", "KNU", "ZJU",
                "HU+NEF+CSU+UI+KNU", "HU+NEF+CSU+UI+KNU+ZJU")

gebvout <- matrix(nrow = nrow(allMultiSiteBLUPs), ncol = 218,
                  dimnames = list(rownames(allMultiSiteBLUPs), 1:218))
currcol <- 1
for(i in 1:length(traitNames)){
  tr <- names(mydata)[year3traits][i]
  for(j in 1:length(site_order)){
    sc <- site_order[j]
    thisfile <- paste("181115select_within_group/", sc, ".", tr, ".RData", sep = "")
    if(!file.exists(thisfile)) next
    load(thisfile)
    outmat <- cbind(rowMeans(predvalues), apply(predvalues, 1, sd)/sqrt(ncol(predvalues)))
    colnames(outmat) <- paste(traitSym[i], site_names[j], "Year 3",
                              c("Breeding value", "Standard error"))
    gebvout[rownames(predvalues), currcol + (0:1)] <- outmat
    colnames(gebvout)[currcol + (0:1)] <- colnames(outmat)
    currcol <- currcol + 2
  }
}
#write.csv(gebvout, file = "181115GEBVs_within_group.csv")

#### Regress BLUPs on PC axes and get residuals ####
# do PCA
mypca <- prcomp(myA.EM.Msi$imputed, scale. = TRUE)
plot(mypca$sdev[1:100]) # maybe first 9?
plot(mypca$x[,1], mypca$x[,2])
myPCs <- data.frame(mypca$x[,1:20], factor(mydapc[rownames(mypca$x),"DAPC_group"]))
names(myPCs)[21] <- "DAPC_group"
levels(myPCs$DAPC_group)

# look at PC axes to understand which are useful
library(ggplot2)
library(dplyr)

ggplot(myPCs, aes(x = PC1, y = PC9, col = DAPC_group)) +
  geom_point()

rownames(myPCs)[myPCs$PC5 > 25 & myPCs$PC1 < -25]
rownames(myPCs)[myPCs$PC5 > 125 & myPCs$PC1 < -25]

# PC1 = Japan v. mainland
# PC2 = mainland N v. S
# PC3 = Japan N v. S, orn. v. S
# PC4 = Yangtze-Qinling v. other mainland groups
# PC5 = separate a few ornamentals (Msa hybrids)
# PC6 = Sichuan v. SE China/tropical, S. Japan v. ornamentals
# PC7 = variation w/in orn and naturalized
# PC8 = S Japan v. ornamentals
# PC9 = variation within ornamentals

# how many PCs to use?
testlm <- lm(allMultiSiteBLUPs[rownames(myPCs),
                               "HU.NEF.CHA.CSU.KNU.Biomass.yield.dry.weight.year.3"] ~ 
               as.matrix(myPCs[,1:5]))
summary(testlm) # PC2 highly sig, PC1 and 5 slightly

testlm2 <- lm(allMultiSiteBLUPs[rownames(myPCs),
                               "HU.NEF.CHA.CSU.KNU.Biomass.yield.dry.weight.year.3"] ~ 
               as.matrix(myPCs[,1:9]))
summary(testlm2) # PC 8 and 9 highly sig, raises R2 a lot with respect to just first five

testlm3 <- lm(allMultiSiteBLUPs[rownames(myPCs),
                                "HU.NEF.CHA.CSU.KNU.Biomass.yield.dry.weight.year.3"] ~ 
                as.matrix(myPCs[,1:15]))
summary(testlm3) # PC14 highly sig - variation within ornamentals

# get R-squared value by number of PCs included
myR2 <- numeric(20)
for(i in 1:20){
  myR2[i] <- summary(lm(allMultiSiteBLUPs[rownames(myPCs),
                                          "HU.NEF.CHA.CSU.KNU.Biomass.yield.dry.weight.year.3"] ~ 
                          as.matrix(myPCs[,1:i])))$r.squared
}
plot(myR2) # first 9 looks good.  9 has twice the R2 of 8
grid()

# plots of pcs vs yield
myPCs$yield <- allMultiSiteBLUPs[rownames(myPCs),
                                 "HU.NEF.CHA.CSU.KNU.Biomass.yield.dry.weight.year.3"]
myPCs$culms <- allMultiSiteBLUPs[rownames(myPCs),
                                 "HU.NEF.CHA.CSU.KNU.Total.culms.year.3"]
ggplot(myPCs, aes(x = PC2, y = yield)) +
  geom_point() +
  geom_smooth(method = lm)

# how are the ornamentals ordered on PC9?  Obvious pattern?
filter(myPCs, DAPC_group == "Ornamental") %>%
ggplot(aes(x = PC9, y = yield)) +
  geom_point() +
  geom_smooth(method = lm)

filter(myPCs, DAPC_group == "Ornamental") %>%
  ggplot(aes(x = PC9, y = culms)) +
  geom_point() +
  geom_smooth(method = lm) # not related to number of culms

orn <- row.names(myPCs)[myPCs$DAPC_group == "Ornamental"]
myPCs[orn[order(myPCs[orn,"PC9"])],c("PC9", "yield")]

# proportion of marker variance for first 9?
sum(mypca$sdev[1:9]^2)/sum(mypca$sdev^2) # 25.4% of marker variation
sum(mypca$sdev[1:5]^2)/sum(mypca$sdev^2) # 20.4% for first five

## make residuals
AllMultiSiteResiduals <- matrix(NA_real_, nrow = nrow(allMultiSiteBLUPs),
                                ncol = ncol(allMultiSiteBLUPs),
                                dimnames = dimnames(allMultiSiteBLUPs))
for(i in 1:ncol(allMultiSiteBLUPs)){
  thislm <- lm(allMultiSiteBLUPs[rownames(myPCs),i] ~ as.matrix(myPCs[,1:9]))
  AllMultiSiteResiduals[names(thislm$residuals),i] <-
    thislm$residuals
}
AllOneSiteResiduals <- matrix(NA_real_, nrow = nrow(allOneSiteBLUPs),
                                ncol = ncol(allOneSiteBLUPs),
                                dimnames = dimnames(allOneSiteBLUPs))
for(i in 1:ncol(allOneSiteBLUPs)){
  thislm <- lm(allOneSiteBLUPs[rownames(myPCs),i] ~ as.matrix(myPCs[,1:9]))
  AllOneSiteResiduals[names(thislm$residuals),i] <-
    thislm$residuals
}

#save(AllMultiSiteResiduals, AllOneSiteResiduals, file = "180927residuals_9PCs.RData")

## Do genomic prediction with residuals
library(rrBLUP)
thisKfold <- 10
year3traits <- traits[c(2,5,8,11,14,17,20,23,25,28,30,33,36,39,41)]
names(mydata)[year3traits]
residuals_accuracies <- matrix(NA_real_, nrow = length(year3traits), ncol = 8,
                       dimnames = list(names(mydata)[year3traits], 
                                       c("HU", "NEF", "CSU", "CHA", "KNU", "ZJU", "HU.NEF.CHA.CSU.KNU", "HU.NEF.CHA.CSU.KNU.ZJU")))
#for(tr in names(mydata)[year3traits]){
for(tr in names(mydata)[year3traits][-1]){
  print(tr)
  for(sc in colnames(residuals_accuracies)){
#  for(sc in colnames(residuals_accuracies)[7:8]){
    print(sc)
    sc.tr <- paste(sc, tr, sep = ".")
    thisfile <- paste("180927residuals_accuracy_9PCs/", sc.tr, ".RData", sep = "")
    sites <- strsplit(sc, "\\.")[[1]] # locations for this dataset
    
    # residuals to use
    if(sc.tr %in% colnames(AllOneSiteResiduals)){
      thisres <- AllOneSiteResiduals[,sc.tr]
    } else if(sc.tr %in% colnames(AllMultiSiteResiduals)){
      thisres <- AllMultiSiteResiduals[,sc.tr]
    } else {
      next
    }
    thisres <- thisres[!is.na(thisres)]
    
    # 100 reps of 10-fold cross validation
    predvalues <- matrix(NA_real_, nrow = length(thisres), ncol = 100,
                         dimnames = list(names(thisres), NULL))
    for(j in 1:100){
      print(j)
      scrambleInd <- sample(length(thisres)) # random individual order for this iteration
      genPredOutput <- list()
      length(genPredOutput) <- thisKfold 
      nPerRep <- floor(length(thisres)/thisKfold)
      
      for(i in 1:thisKfold){
        # identify individuals for training and prediction sets
        firstind <- (i-1) * nPerRep + 1
        if(i == thisKfold){
          lastind <- length(thisres)
        } else {
          lastind <- i * nPerRep
        }
        train <- scrambleInd[-(firstind:lastind)]
        pred <- scrambleInd[firstind:lastind]
        
        # phenotypes for training set only
        thisphen <- rep(NA, length(thisres))
        thisphen[train] <- thisres[train]
        genPredOutput[[i]] <- kin.blup(data = data.frame(pheno = thisphen, geno = names(thisres)),
                                       geno = "geno", pheno = "pheno", GAUSS = FALSE,
                                       K = myA.EM.Msi$A[names(thisres), names(thisres)])
        
        predvalues[pred,j] <- genPredOutput[[i]]$g[pred]
      }
    }
    save(predvalues, file = thisfile)
    # prediction accuracy
    residuals_accuracies[tr, sc] <- cor(thisres, rowMeans(predvalues)) # incorrect but not using this
  }
}
#save(residuals_accuracies, file = "180927residuals_accuracy_9PCs/accuracy.RData")
#write.csv(residuals_accuracies, file = "180927residuals_accuracy_9PCs/accuracy.csv")
hist(residuals_accuracies) # centered below zero
# combined sites still moderate for culm related traits, everything else pretty bad.
# culm length, node number, internode length, diameter of basal internode, culm volume

#### Get residuals of BLUPs after fitting to DAPC group ID ####
## make residuals
AllMultiSiteResiduals <- matrix(NA_real_, nrow = nrow(allMultiSiteBLUPs),
                                ncol = ncol(allMultiSiteBLUPs),
                                dimnames = dimnames(allMultiSiteBLUPs))
for(i in 1:ncol(allMultiSiteBLUPs)){
  thislm <- lm(allMultiSiteBLUPs[,i] ~ mydapc$DAPC_group)
  AllMultiSiteResiduals[names(thislm$residuals),i] <-
    thislm$residuals
}
AllOneSiteResiduals <- matrix(NA_real_, nrow = nrow(allOneSiteBLUPs),
                              ncol = ncol(allOneSiteBLUPs),
                              dimnames = dimnames(allOneSiteBLUPs))
for(i in 1:ncol(allOneSiteBLUPs)){
  thislm <- lm(allOneSiteBLUPs[,i] ~ mydapc$DAPC_group)
  AllOneSiteResiduals[names(thislm$residuals),i] <-
    thislm$residuals
}

#save(AllMultiSiteResiduals, AllOneSiteResiduals, file = "180928residuals_DAPC.RData")
load("180928residuals_DAPC.RData")

## Do genomic prediction with residuals
library(rrBLUP)
thisKfold <- 10
year3traits <- traits[c(2,5,8,11,14,17,20,23,25,28,30,33,36,39,41)]
names(mydata)[year3traits]
residuals_accuracies <- matrix(NA_real_, nrow = length(year3traits), ncol = 8,
                               dimnames = list(names(mydata)[year3traits], 
                                               c("HU", "NEF", "CSU", "CHA", "KNU", "ZJU", "HU.NEF.CHA.CSU.KNU", "HU.NEF.CHA.CSU.KNU.ZJU")))
for(tr in names(mydata)[year3traits]){
  print(tr)
  for(sc in colnames(residuals_accuracies)){
    print(sc)
    sc.tr <- paste(sc, tr, sep = ".")
    thisfile <- paste("180928residuals_accuracy_DAPC/", sc.tr, ".RData", sep = "")
    sites <- strsplit(sc, "\\.")[[1]] # locations for this dataset

    # residuals to use
    if(sc.tr %in% colnames(AllOneSiteResiduals)){
      thisres <- AllOneSiteResiduals[,sc.tr]
    } else if(sc.tr %in% colnames(AllMultiSiteResiduals)){
      thisres <- AllMultiSiteResiduals[,sc.tr]
    } else {
      next
    }
    thisres <- thisres[!is.na(thisres)]
    thisres <- thisres[names(thisres) != "KMS018"]
    
    # skip analysis if already done
    if(file.exists(thisfile)){
      load(thisfile)
      residuals_accuracies[tr, sc] <- predacc(thisres, predvalues)
      next
    }
    
    # 100 reps of 10-fold cross validation
    predvalues <- matrix(NA_real_, nrow = length(thisres), ncol = 100,
                         dimnames = list(names(thisres), NULL))
    for(j in 1:100){
      print(j)
      scrambleInd <- sample(length(thisres)) # random individual order for this iteration
      genPredOutput <- list()
      length(genPredOutput) <- thisKfold 
      nPerRep <- floor(length(thisres)/thisKfold)
      
      for(i in 1:thisKfold){
        # identify individuals for training and prediction sets
        firstind <- (i-1) * nPerRep + 1
        if(i == thisKfold){
          lastind <- length(thisres)
        } else {
          lastind <- i * nPerRep
        }
        train <- scrambleInd[-(firstind:lastind)]
        pred <- scrambleInd[firstind:lastind]
        
        # phenotypes for training set only
        thisphen <- rep(NA, length(thisres))
        thisphen[train] <- thisres[train]
        genPredOutput[[i]] <- kin.blup(data = data.frame(pheno = thisphen, geno = names(thisres)),
                                       geno = "geno", pheno = "pheno", GAUSS = FALSE,
                                       K = myA.EM.Msi$A[names(thisres), names(thisres)])
        
        predvalues[pred,j] <- genPredOutput[[i]]$g[pred]
      }
    }
    save(predvalues, file = thisfile)
    # prediction accuracy
    residuals_accuracies[tr, sc] <- predacc(thisres, predvalues)
  }
}
#save(residuals_accuracies, file = "180928residuals_accuracy_DAPC/accuracy190221.RData")
#write.csv(residuals_accuracies, file = "180928residuals_accuracy_DAPC/accuracy190221.csv")

## Export residuals to spreadsheet
load("180928residuals_DAPC.RData")
names(mydata)[year3traits]

traitSym <- c("Yld", "CC", "BC", "CC/BC", "CmL", "CmNdN", "IntL", "CmDW",
              "CmV", "CmDW/V", "DBI", "DTI", "TCmN", "RCmN/TCmN", "TCmN/A")
traitNames <- c("Biomass.yield.dry.weight", "Ccirc", "Bcirc", "Ccirc.Bcirc",
                "CmL", "CmN", "IntL", "CmDW", "CmVol", "CmDens", "CmD_I1", "CmD_LI", "Total.culms",
                "Prop.repro.culms", "Culm.per.area")
site_order <- c("HU", "NEF", "CSU", "CHA", "KNU", "ZJU",
                "HU.NEF.CHA.CSU.KNU", "HU.NEF.CHA.CSU.KNU.ZJU")
site_names <- c("HU", "NEF", "CSU", "UI", "KNU", "ZJU",
                "HU+NEF+CSU+UI+KNU", "HU+NEF+CSU+UI+KNU+ZJU")
# ordered matrix to export
resout <- matrix(nrow = nrow(AllMultiSiteResiduals), ncol = 0,
                  dimnames = list(rownames(AllMultiSiteResiduals), character(0)))
for(tr in traitNames){
  # column names in existing matrices
  thiscol <- paste(site_order, tr, 
                   ifelse(tr == "Ccirc.Bcirc", "year3","year.3"),
                   sep= ".")
  thiscolS <- thiscol[1:6]
  thiscolS <- thiscolS[thiscolS %in% colnames(AllOneSiteResiduals)]
  thiscolM <- thiscol[7:8]
  thiscolM <- thiscolM[thiscolM %in% colnames(AllMultiSiteResiduals)]
  # output column names
  outcol <- paste(traitSym[match(tr, traitNames)], site_names, "year 3")
  # output matrix
  outmat <- cbind(AllOneSiteResiduals[,thiscolS], AllMultiSiteResiduals[,thiscolM])
  colnames(outmat) <- outcol[match(c(thiscolS, thiscolM), thiscol)]
  resout <- cbind(resout, outmat)
}
#write.csv(resout, file = "190218residuals_DAPC_values.csv")

## export GEBVs from residuals to spreadsheet
gebvout <- matrix(nrow = nrow(AllMultiSiteResiduals), ncol = 218,
                  dimnames = list(rownames(AllMultiSiteResiduals), 1:218))
currcol <- 1
for(i in 1:length(traitNames)){
  tr <- names(mydata)[year3traits][i]
  for(j in 1:length(site_order)){
    sc <- site_order[j]
    thisfile <- paste("180928residuals_accuracy_DAPC/", sc, ".", tr, ".RData", sep = "")
    if(!file.exists(thisfile)) next
    load(thisfile)
    outmat <- cbind(rowMeans(predvalues), apply(predvalues, 1, sd)/sqrt(ncol(predvalues)))
    colnames(outmat) <- paste(traitSym[i], site_names[j], "Year 3",
                              c("Breeding value", "Standard error"))
    gebvout[rownames(predvalues), currcol + (0:1)] <- outmat
    colnames(gebvout)[currcol + (0:1)] <- colnames(outmat)
    currcol <- currcol + 2
  }
}
#write.csv(gebvout, file = "181003residual_DAPC_GEBVs.csv")

## Figure comparing genomic prediction with BLUPs and residuals ####
newYldBLUPs <- read.csv("181116values_for_genomic_prediction_fig.csv", row.names = 1,
                        stringsAsFactors = FALSE)
names(dapccol) <- c("S Japan", "N Japan", "Korea/N China","Sichuan","Yangtze-Qinling",
                    "SE China/tropical", "Msa", "Ornamental", "US naturalized")
dapcpch <- c(22, 22, 21, 21, 24, 25, 21, 24, 25)
newYldBLUPs$DAPC.num <- match(newYldBLUPs$DAPC.group, names(dapccol))

pdf("Fig2_yield_pred_190326.pdf", width = 6.5, height = 6.5, pointsize = 10)
layout(matrix(c(1,3,5,2,4,5), nrow = 3, ncol = 2), heights = c(3, 3, 0.5))
par(mgp = c(2, 0.5, 0), mar = c(3.1, 3.1, 3.1, 1.1))

plot(newYldBLUPs$BLUP_eqn2, newYldBLUPs$GEBV2, bg = dapccol[newYldBLUPs$DAPC.num],
     xlab = "Genotypic BLUP (Eqn. 2)", ylab = "GEBV",
     pch = dapcpch[newYldBLUPs$DAPC.num], cex = 1.2)
text(-12, 4.5, expression(italic(R) == 0.47))
title(main = "(a)", adj = 0, font = 2)

plot(newYldBLUPs$BLUP_eqn2, newYldBLUPs$GEBV2_grp, bg = dapccol[newYldBLUPs$DAPC.num],
     xlab = "Genotypic BLUP (Eqn. 2)", ylab = "GEBV within genetic group",
     pch = dapcpch[newYldBLUPs$DAPC.num], cex = 1.2)
text(-12, 2.2, expression(italic(R) == 0.45))
title(main = "(b)", adj = 0, font = 2)

plot(newYldBLUPs$BLUP_eqn4, newYldBLUPs$GEBV4, bg = dapccol[newYldBLUPs$DAPC.num],
     xlab = "Genotype-within-genetic group BLUP (Eqn. 4)", ylab = "GEBV",
     pch = dapcpch[newYldBLUPs$DAPC.num], cex = 1.2)
text(-604, 361, expression(italic(R) == 0.35))
title(main = "(c)", adj = 0, font = 2)

plot(newYldBLUPs$residual_eqn5, newYldBLUPs$GEBV5, bg = dapccol[newYldBLUPs$DAPC.num],
     xlab = "Residual of BLUP fitted on genetic group (Eqn. 5)", 
     ylab = "GEBV",
     pch = dapcpch[newYldBLUPs$DAPC.num], cex = 1.2)
text(-10, 4, expression(italic(R) == 0.31))
title(main = "(d)", adj = 0, font = 2)
# color legend
par(mar = rep(0, 4))
plot(NA, axes = FALSE, col = "white", xlim = c(0,16), ylim = c(0, 9))
legtext <- c("S Japan", "N Japan", "Korea/N China", "Sichuan", "Yangtze-Qinling",
             "SE China/tropical", "Ornamentals", "US naturalized")
thisy <- 6 - ((((1:8)+1) %% 2) * 3)
spc <- 1.25
xcoord <- strwidth("Genetic groups:") * spc
xcoord <- c(xcoord, xcoord[1] + strwidth("N Japan  ") * spc)
xcoord <- c(xcoord, xcoord[2] + strwidth("Korea/N China") * spc)
xcoord <- c(xcoord, xcoord[3] + strwidth("SE China/tropical") * spc)
thisx <- rep(xcoord, each = 2)
points(thisx, thisy, pch = dapcpch[-7], cex = 2, bg = dapccol[-7])
text(thisx, thisy, legtext, pos = 4)
text(0, 6, "Genetic groups:", pos = 4)

dev.off()

# how much of correlation is from ornamentals?
cor(newYldBLUPs$residual, newYldBLUPs$GEBV2, use = "complete.obs") # 0.32
notorn <- which(newYldBLUPs$DAPC.group != "Ornamental")
cor(newYldBLUPs$residual[notorn], newYldBLUPs$GEBV2[notorn], use = "complete.obs") # 0.22

cor(newYldBLUPs$BLUP, newYldBLUPs$GEBV1, use = "complete.obs") # 0.48
cor(newYldBLUPs$BLUP[notorn], newYldBLUPs$GEBV1[notorn], use = "complete.obs") # 0.46
