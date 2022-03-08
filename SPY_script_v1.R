require(data.table)
library(dplyr)

# set base values

pG <- -2
pA <- -2
pC <- 2
pT <- 3

setwd("~/RNAseq_data_et_analyses/SPY_script")
getwd()

#import files with intronic sequences and tab delimited nucleotides matrix
SPY_seq <- as.data.frame(fread("test_SPY_seq"))
SPY_tab <- as.data.frame(fread("test_SPY_tab"))

SPY_tab$exon <- NULL
seq_size <- as.numeric(ncol(SPY_tab))

# find AG exclusion zone size (agez)
SPY_seq$agez <- sapply(1:nrow(SPY_seq), function(i) nchar(gsub(".*AG","",substring(SPY_seq[i,"intron_sequence"], first =1, last = (nchar(SPY_seq[i,"intron_sequence"])-2))))+1)

# find SPY score
# create a matrix for pyrimidines
SPY_tab_M <- data.frame()
SPY_tab_M <- (SPY_tab=="C" | SPY_tab=="T")

# find pyrimidine stretches (PY) of 9 bases
SPY_tab_M9 <- data.frame()
for (i in 1:nrow(SPY_tab_M)) for (j in (9:(seq_size-4)))
{SPY_tab_M9[i,j-8] <- (SPY_tab_M[i,j-8] & SPY_tab_M[i,j]) &
  (SPY_tab_M[i,j-7] | (SPY_tab_M[i,j-6] & SPY_tab_M[i,j-5] & SPY_tab_M[i,j-4])) &
  (SPY_tab_M[i,j-6] | (SPY_tab_M[i,j-7] & SPY_tab_M[i,j-5] & SPY_tab_M[i,j-4])) &
  (SPY_tab_M[i,j-5] | (SPY_tab_M[i,j-7] & SPY_tab_M[i,j-6] & SPY_tab_M[i,j-4]) | (SPY_tab_M[i,j-6] & SPY_tab_M[i,j-4] & SPY_tab_M[i,j-3] & SPY_tab_M[i,j-2])) &
  (SPY_tab_M[i,j-4] | (SPY_tab_M[i,j-7] & SPY_tab_M[i,j-6] & SPY_tab_M[i,j-5] & SPY_tab_M[i,j-3]) | (SPY_tab_M[i,j-6] & SPY_tab_M[i,j-5] & SPY_tab_M[i,j-3] & SPY_tab_M[i,j-2]) | (SPY_tab_M[i,j-5] & SPY_tab_M[i,j-3] & SPY_tab_M[i,j-2] & SPY_tab_M[i,j-1])) &
  (SPY_tab_M[i,j-3] | (SPY_tab_M[i,j-6] & SPY_tab_M[i,j-5] & SPY_tab_M[i,j-4] & SPY_tab_M[i,j-2]) | (SPY_tab_M[i,j-4] & SPY_tab_M[i,j-2] & SPY_tab_M[i,j-1])) &
  (SPY_tab_M[i,j-2] | (SPY_tab_M[i,j-4] & SPY_tab_M[i,j-3] & SPY_tab_M[i,j-1])) &
  (SPY_tab_M[i,j-1] | (SPY_tab_M[i,j-4] & SPY_tab_M[i,j-3] & SPY_tab_M[i,j-2]))
}

# find PY of 8 bases
SPY_tab_M8 <- data.frame()
for (i in 1:nrow(SPY_tab_M)) for (j in (8:(seq_size-4)))
{SPY_tab_M8[i,j-7] <-  (SPY_tab_M[i,j-7] & SPY_tab_M[i,j]) & 
  (SPY_tab_M[i,j-6] | (SPY_tab_M[i,j-5] & SPY_tab_M[i,j-4] & SPY_tab_M[i,j-3])) &
  (SPY_tab_M[i,j-5] | (SPY_tab_M[i,j-6] & SPY_tab_M[i,j-4] & SPY_tab_M[i,j-3])) &
  (SPY_tab_M[i,j-4] | (SPY_tab_M[i,j-6] & SPY_tab_M[i,j-5] & SPY_tab_M[i,j-3]) | (SPY_tab_M[i,j-5] & SPY_tab_M[i,j-3] & SPY_tab_M[i,j-2] & SPY_tab_M[i,j-1])) &
  (SPY_tab_M[i,j-3] | (SPY_tab_M[i,j-6] & SPY_tab_M[i,j-5] & SPY_tab_M[i,j-4] & SPY_tab_M[i,j-2]) | (SPY_tab_M[i,j-4] & SPY_tab_M[i,j-2] & SPY_tab_M[i,j-1])) &
  (SPY_tab_M[i,j-2] | (SPY_tab_M[i,j-4] & SPY_tab_M[i,j-3] & SPY_tab_M[i,j-1])) &
  (SPY_tab_M[i,j-1] | (SPY_tab_M[i,j-4] & SPY_tab_M[i,j-3] & SPY_tab_M[i,j-2]))
}

# calculate score for each PY of 9 bases
SPY_tab_M9v <- data.frame()
for (i in 1:nrow(SPY_tab_M)) for (j in (9:(seq_size-4)))
{v <- 0
for (k in ((j-8):j)) (v <- v + (SPY_tab[i,k]=="G")*pG+(SPY_tab[i,k]=="A")*pA+(SPY_tab[i,k]=="C")*pC+(SPY_tab[i,k]=="T")*pT)
SPY_tab_M9v[i,j-8] <- v*SPY_tab_M9[i,j-8]
}

# calculate score for each PY of 8 bases
SPY_tab_M8v <- data.frame()
for (i in 1:nrow(SPY_tab_M)) for (j in 8:(seq_size-4))
{v <- 0
for (k in ((j-7):j)) (v <- v + (SPY_tab[i,k]=="G")*pG+(SPY_tab[i,k]=="A")*pA+(SPY_tab[i,k]=="C")*pC+(SPY_tab[i,k]=="T")*pT)
SPY_tab_M8v[i,j-7] <- v*SPY_tab_M8[i,j-7]
}

# find recursively the SPY value within the seq_size window
SPY_tab_SPY_M <- matrix(NA,nrow=nrow(SPY_tab_M),ncol=seq_size-11)
SPY_tab_SPY_M[,seq_size-11]<- SPY_tab_M8v[,seq_size-11]
for (i in 1:nrow(SPY_tab_M)) for (k in 1:7) SPY_tab_SPY_M[i,seq_size-11-k] <- max(SPY_tab_SPY_M[i,seq_size-11-k+1], SPY_tab_M9v[i,seq_size-11-k], SPY_tab_M8v[i,seq_size-11-k] )
for (i in 1:nrow(SPY_tab_M)) for (k in 8:8) SPY_tab_SPY_M[i,seq_size-11-k] <- max(SPY_tab_SPY_M[i,seq_size-11-k+1], SPY_tab_M9v[i,seq_size-11-k], SPY_tab_SPY_M[i,seq_size-11]+SPY_tab_M8v[i,seq_size-11-k])
for (i in 1:nrow(SPY_tab_M)) for (k in 9:(seq_size-12)) SPY_tab_SPY_M[i,seq_size-11-k] <- max(SPY_tab_SPY_M[i,seq_size-11-k+1], SPY_tab_SPY_M[i,seq_size-11-k+9]+SPY_tab_M9v[i,seq_size-11-k], SPY_tab_SPY_M[i,seq_size-11-k+8]+SPY_tab_M8v[i,seq_size-11-k])

SPY_seq[,"SPY"] <- SPY_tab_SPY_M[,1]

# find  the SPY value within the agez
for (i in 1:nrow(SPY_tab))
  SPY_seq[i,"SPY_agez"] <- SPY_tab_SPY_M[i,min(seq_size-11,seq_size-SPY_seq[i,"agez"])]

#end


