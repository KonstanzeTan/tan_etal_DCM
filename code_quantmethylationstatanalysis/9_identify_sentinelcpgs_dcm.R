#' -----------------------------------------------------------------------------
#' Assign CpGs to unique genomic loci 
#' 
#' @description This script assigns candidate sentinel CpGs to unique genomic 
#' loci (>1Mb apart), identifies the lead signal at each locus and highlights 
#' loci wth more than one candidate CpG for conditional analysis to identify
#' independent signals.
#'
#' @author Konstanze Tan <konstanz001@e.ntu.edu.sg>
#'
#' @date Tues Mar 6 2024
#' -----------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# load required libraries
library(dplyr)

# ------------------------------------------------------------------------------
# Assign CpGs to genomic loci

# load candidate sentinel CpGs (disc Bonferroni P<0.05, rep const. dir.)
## columns: CpG, chr, pos, EWAS P
check_sentinel<-readRDS("candidate_196_sentinel.rds") 

# arrange candidate sentinel CpGs by location 
check_sentinel<-check_sentinel %>%
  group_by(chr) %>%
  arrange(pos, .by_group = TRUE)

check_sentinel$loci_ID<-0

for(cg in 1:length(check_sentinel$CpG)){
  cpg=check_sentinel[cg,]
  cpgChr = cpg$chr
  cpgPos = cpg$pos
  if(cg==1) {check_sentinel[cg,5]=1}
  if(cg>1)  {
    preceedingCpgs = check_sentinel[cg-1,] 
    check_sentinel[cg,5]=ifelse(preceedingCpgs$chr == cpgChr & 
                                  abs(preceedingCpgs$pos-cpgPos)<500000,
                                preceedingCpgs$loci_ID,preceedingCpgs$loci_ID+1)
  }
}

# ------------------------------------------------------------------------------
# Identify lead signal in locus 

# arrange candidate sentinel CpGs by EWAS P
check_sentinel<-check_sentinel %>%
  group_by(chr) %>%
  arrange(disc_P_baconcorr, .by_group = TRUE)

for(cg in 1:length(check_sentinel$CpG)){
  cpg=check_sentinel[cg,]
  cpgChr = cpg$chr
  cpgPos = cpg$pos
  if(cg==1) {lead="Y"}
  if(cg>1)  {
    preceedingCpgs = check_sentinel[1:cg-1,] 
    keep=ifelse(length(which(preceedingCpgs$chr == cpgChr & 
                               abs(preceedingCpgs$pos-cpgPos)<500000))==0,T,NA)
    lead=c(lead,keep)
  }
}

check_sentinel$lead<-lead

# ------------------------------------------------------------------------------
# format data for conditional analysis

# identify genomic loci containing >1 candidate sentinel CpG 
loci_freq<-as.data.frame(table(check_sentinel$loci_ID))
names(loci_freq)<-c("loci_ID", "Freq")
multiple_candidates_loci_ID<-as.numeric(loci_freq[loci_freq$Freq>1,]$loci_ID)

# subset check_sentinel df for loci with >1 CpG
multiple_candidates_loci<- check_sentinel[check_sentinel$loci_ID 
                                          %in% multiple_candidates_loci_ID,]

# split into lead and nonlead dataframes
lead_cpgs <- multiple_candidates_loci[which(multiple_candidates_loci$lead 
                                            == TRUE), ]
nonlead_cpgs <- multiple_candidates_loci[is.na(multiple_candidates_loci$lead), ]
names(lead_cpgs)[1]<-c("lead_cpg")
names(nonlead_cpgs)[1]<-c("nonlead_cpg")

# merge to create combined df specifying pairs for conditional analysis 
## df columns: loci_id, lead_cpg, non-lead cpg
cond_analysis_pairs <- merge(lead_cpgs[,c("loci_ID", "lead_cpg")], 
                             nonlead_cpgs[,c("loci_ID", "nonlead_cpg")], 
                             by="loci_ID", all.x=TRUE) 

# save output 
saveRDS(check_sentinel, "candidate_196_sentinels_leadstatus_loci.rds")
saveRDS(cond_analysis_pairs, "cond_analysis_pairs.rds")

# ------------------------------------------------------------------------------
sessionInfo()
