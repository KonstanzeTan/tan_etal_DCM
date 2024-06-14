#' -----------------------------------------------------------------------------
#' Generating background set of CpGs for enrichment analysis (eQTM and TFBS)
#' 
#' @description Script generates a background consisting of 1000 permuted sets of 
#' CpGs, each equal in length to the number of DCM sentinel CpG set (n=194).
#' A sliding-window approach is implemented to select EPIC array CpGs of 
#' equivalent methylation levels and variablity to each sentinel CpG
#'
#' @author Konstanze Tan <konstanz001@e.ntu.edu.sg> (adapted from: McAllan et al.
#' 2023 (doi: 10.1038/s41467-023-38439-z)
#'
#' @date Tues Mar 7 2023
#' -----------------------------------------------------------------------------
# load sentinels

sentinels<-read.delim("sentinels_194_ID", header=F)
sentinels<-sentinels$V1

# ------------------------------------------------------------------------------
# load mean and sd methylation values for all EPIC CpG sites analysed in EWAS
# mean and sd calculated separately in DCM (n1) and Normal (n2) groups; then combined
## combined mean formula: (n1*mean1+n2*mean2)/(n1+n2)
## combined sd formula: sqrt(((n1-1)*s1^2 + (n2-1)*s2^2) / (n1+n2-2))

msd<-readRDS("msd_ewas_838624cpgs_magnet_formatted.rds")

# ------------------------------------------------------------------------------
# generate background table of all CpGs excluding i. sentinel CpGs and 
# ii. CpGs within 5kb of a sentinel 

# subset dataframe of mean sd for sentinel CpGs 
hits=data.frame(msd[as.character(sentinels),]) 
# store sentinel CpG ID as a numeric vector
hits$CG=as.character(rownames(hits)) 

# Identify CpGs within 5kb of sentinel CpGs
for (c in 1:length(hits$CG)) {
  print(c)
  cg <- hits$CG[c]
  cg.coord <- msd[cg, ] # Retrieve coordinate of the sentinel
  
  # Retrieve coordinates of all CpGs within 5kb of the sentinel
  cis.cg.coord <- msd[msd$chr == cg.coord$chr & abs(msd$pos - cg.coord$pos) < 5000, ]
  
  # Retrieve CpG identifiers
  if (c == 1) {
    cis.cgs <- cis.cg.coord$CpG
  } else {
    cis.cgs <- c(cis.cgs, cis.cg.coord$CpG)
  }
}

cis.cgs <- unique(cis.cgs)
# Exclude sentinels from background dataframe
backg <- data.frame(msd[!(rownames(msd) %in% hits$CG), ])
backg <- data.frame(msd[!(rownames(msd) %in% cis.cgs), ])

# Set permutation parameters
parameters <- expand.grid(seq(0.0025, 0.025, 0.0025), seq(0.025, 0.25, 0.025), stringsAsFactors = FALSE)
colnames(parameters) <- c("SD", "Mean")

# Create empty matrix to be filled by permutated CpGs
matches <- matrix(nrow = nrow(hits), ncol = 1000)
colnames(matches) <- as.character(seq(from = 1, to = 1000, 1))
rownames(matches) <- hits$CG

# Create matrix to store match parameters (lowest threshold at which >1000 permuted CpG sites obtained)
matchParameters <- matrix(nrow = nrow(hits), ncol = 3)
colnames(matchParameters) <- c("mean_threshold", "sd_threshold", "n_matches")
rownames(matchParameters) <- hits$CG

# Identify matched CpGs for each sentinel
for (h in 1:nrow(hits)) {
  cpg.match <- NA
  len <- NA
  hit <- hits[h, ]
  print(h)
  
  # Calculate number of matches across the full range of mean SD thresholds
  for (para in 1:nrow(parameters)) {
    n <- length(sample(rownames(backg[abs(backg$meanMeth - hit$meanMeth) < parameters[para, 2] & abs(backg$sdMeth - hit$sdMeth) < parameters[para, 1], ])))
    if (para == 1) {
      sel <- c(parameters[para, 2], parameters[para, 1], n)
    } else {
      sel <- rbind(sel, c(parameters[para, 2], parameters[para, 1], n))
    }
  }
  
  colnames(sel) <- c("mean", "sd", "n_matches") # Number of rows = number of parameters
  
  # Select combinations of mean SD thresholds which give n > 1000 for a particular sentinel
  sel <- sel[which(sel[, 3] > 1000), ]
  meanthres <- sel[1, 1]
  sdthres <- sel[1, 2]
  nthresh <- sel[1, 3]
  
  # Identify matches for each hit, filling by columns without replacement
  for (p in 1:1000) {
    set.seed(p)
    print(p)
    if (exists('cgm')) { rm(cgm) }
    cgms <- sample(rownames(backg[abs(backg$meanMeth - hit$meanMeth) < meanthres & abs(backg$sdMeth - hit$sdMeth) < sdthres, ]))
    cgm <- setdiff(cgms, cpg.match)[1]
    if (!is.na(cgm)) {
      cpg.match <- c(cpg.match, cgm)
    } else {
      cpg.match <- c(cpg.match, as.character(rownames(hit))) # Temporary solution
      print('no matching marker found')
    }
  }
  
  matches[h, ] <- na.omit(cpg.match)
  matchParameters[h, ] <- c(meanthres, sdthres, nthresh)
}

# Save results
save(matches, matchParameters, file = "1000_permutations_194_sentinels.RData")

# ------------------------------------------------------------------------------
print("SessionInfo:")
# ------------------------------------------------------------------------------
sessionInfo()
