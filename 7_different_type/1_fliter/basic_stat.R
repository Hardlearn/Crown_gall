library(tidyverse)
dup_Am <- Am_hog[duplicated(Am_hog[,c("gene")]) | 
                      duplicated(Am_hog[,c("gene")], fromLast = TRUE),]
same_dup_am <- dup_Am[duplicated(dup_Am[,c('gene','phylostrata')]) |
                           duplicated(dup_Am[,c('gene',"phylostrata")], fromLast = TRUE),]
diff_Am_hog <- anti_join(dup_Am, same_dup_am)

dup_Ej <- Ej_hog[duplicated(Ej_hog[,c("gene")]) |
                      duplicated(Ej_hog[,c("gene")], fromLast = TRUE),]
same_dup_Ej <- dup_Ej[duplicated(dup_Ej[,c('gene','phylostrata')]) |
                           duplicated(dup_Ej[,c('gene','phylostrata')], fromLast = TRUE),]
diff_Ej_hog <- anti_join(dup_Ej, same_dup_Ej)

dup_Go <- Go_hog[duplicated(Go_hog[,c("gene")]) |
                      duplicated(Go_hog[,c("gene")], fromLast = TRUE),]
same_dup_Go <- dup_Go[duplicated(dup_Go[,c('gene','phylostrata')]) |
                           duplicated(dup_Go[,c('gene','phylostrata')], fromLast = TRUE),]
diff_Go_hog <- anti_join(dup_Go, same_dup_Go)


## protein specific
setwd("F:/project/scp_files/scp_files")
Am_oma <- read.table("Am_oma_root.txt")
colnames(Am_oma) <- c("transcript","protein","hog")
Ej_oma <- read.table("Ej_oma_root.txt")
colnames(Ej_oma) <- c('transcript',"protein","hog")
Go_oma <- read.table("Go_oma_root.txt")
colnames(Go_oma) <- c("transcript","protein","hog")

#Am
isoform_Am_hog <- Am_oma[duplicated(Am_oma[,c("transcript")]) | 
                                         duplicated(Am_oma[,c("transcript")], fromLast = TRUE),]
isoform_same_Am_hog <- isoform_Am_hog[duplicated(isoform_Am_hog[,c("transcript","hog")]) | 
                                           duplicated(isoform_Am_hog[,c("transcript","hog")], fromLast = TRUE),]
isoform_diif_Am_hog <- anti_join(isoform_Am_hog, isoform_same_Am_hog)

#Ej
isoform_Ej_hog <- Ej_oma[duplicated(Ej_oma[,c("transcript")]) | 
                              duplicated(Ej_oma[,c("transcript")], fromLast = TRUE),]
isoform_same_Ej_hog <- isoform_Ej_hog[duplicated(isoform_Ej_hog[,c("transcript","hog")]) | 
                                           duplicated(isoform_Ej_hog[,c("transcript","hog")], fromLast = TRUE),]
isoform_diif_Ej_hog <- anti_join(isoform_Ej_hog, isoform_same_Ej_hog)

#Go
isoform_Go_hog <- Go_oma[duplicated(Go_oma[,c("transcript")]) | 
                              duplicated(Go_oma[,c("transcript")], fromLast = TRUE),]
isoform_same_Go_hog <- isoform_Go_hog[duplicated(isoform_Go_hog[,c("transcript","hog")]) | 
                                           duplicated(isoform_Go_hog[,c("transcript","hog")], fromLast = TRUE),]
isoform_diif_Go_hog <- anti_join(isoform_Go_hog, isoform_same_Go_hog)

## isoform_exp_level
Am_TMM.mat <- read.table("./dl/Am/Am_salmon.isoform.TMM.EXPR.matrix.txt", row.names = NULL)
Ej_TMM.mat <- read.table("./dl/Ej/Ej_salmon.isoform.TMM.EXPR.matrix.txt", row.names = NULL)
Go_TMM.mat <- read.table("./dl/Go/Go_salmon.isoform.TMM.EXPR.matrix.txt", row.names = NULL)
colnames(Am_TMM.mat)[1] <- "transcript"
colnames(Ej_TMM.mat)[1] <- "transcript"
colnames(Go_TMM.mat)[1] <- "transcript"

# isoform to gene correct
Am_TMM_diff.mat <- merge(Am_TMM.mat, diff_Am_hog, by = "transcript")
Ej_TMM_diff.mat <- merge(Ej_TMM.mat, diff_Ej_hog, by = "transcript")
Go_TMM_diff.mat <- merge(Go_TMM.mat, diff_Go_hog, by = "transcript")

## mean_value
Am_TMM_diff.mat$mean <- (Am_TMM_diff.mat$Am_gall_rep1_TranscriptAbundance + 
     Am_TMM_diff.mat$Am_gall_rep2_TranscriptAbundance + 
     Am_TMM_diff.mat$Am_gall_rep3_TranscriptAbundance + 
     Am_TMM_diff.mat$Am_leaf_rep1_TranscriptAbundance + 
     Am_TMM_diff.mat$Am_leaf_rep2_TranscriptAbundance +
     Am_TMM_diff.mat$Am_leaf_rep3_TranscriptAbundance) / 6

Ej_TMM_diff.mat$mean <- (Ej_TMM_diff.mat$Ej_gall_rep1_TranscriptAbundance + 
                         Ej_TMM_diff.mat$Ej_gall_rep2_TranscriptAbundance + 
                         Ej_TMM_diff.mat$Ej_gall_rep3_TranscriptAbundance + 
                         Ej_TMM_diff.mat$Ej_leaf_rep1_TranscriptAbundance + 
                         Ej_TMM_diff.mat$Ej_leaf_rep2_TranscriptAbundance + 
                         Ej_TMM_diff.mat$Ej_leaf_rep3_TranscriptAbundance) / 6

Go_TMM_diff.mat$mean <- (Go_TMM_diff.mat$Go_mature_gall_rep1_TranscriptAbundance +
                              Go_TMM_diff.mat$Go_mature_gall_rep2_TranscriptAbundance +
                              Go_TMM_diff.mat$Go_mature_gall_rep3_TranscriptAbundance + 
                              Go_TMM_diff.mat$Go_mature_leaf_rep1_TranscriptAbundance + 
                              Go_TMM_diff.mat$Go_mature_leaf_rep2_TranscriptAbundance +
                              Go_TMM_diff.mat$Go_mature_leaf_rep3_TranscriptAbundance +
                              Go_TMM_diff.mat$Go_young_gall_rep1_TranscriptAbundance + 
                              Go_TMM_diff.mat$Go_young_gall_rep2_TranscriptAbundance +
                              Go_TMM_diff.mat$Go_young_gall_rep3_TranscriptAbundance +
                              Go_TMM_diff.mat$Go_young_leaf_rep1_TranscriptAbundance +
                              Go_TMM_diff.mat$Go_young_leaf_rep2_TranscriptAbundance +
                              Go_TMM_diff.mat$Go_young_leaf_rep3_TranscriptAbundance) / 12

write.table(Am_TMM_diff.mat, "Am_TMM_diff.mat.txt")
write.table(Ej_TMM_diff.mat, "Ej_TMM_diff.mat.txt")
write.table(Go_TMM_diff.mat, "Go_TMM_diff.mat.txt")

## 去除重复
Am_correct_hog <- read.table("Am_correct_hog.txt")
Ej_correct_hog <- read.table("Ej_correct_hog.txt")
Go_correct_hog <- read.table("Go_correct_hog.txt")
colnames(Am_correct_hog) <- c("gene","hog_correct","hog_id")
colnames(Ej_correct_hog) <- c("gene","hog_correct","hog_id")
colnames(Go_correct_hog) <- c("gene","hog_correct","hog_id")

library(dplyr)

A <- Am_correct_hog %>% distinct(gene, .keep_all = TRUE)
E <- Ej_correct_hog %>% distinct(gene, .keep_all = TRUE)
G <- Go_correct_hog %>% distinct(gene, .keep_all = TRUE)

write.table(A, "Am_hog_corrected.txt")
write.table(E, "Ej_hog_corrected.txt")
write.table(G, "Go_hog_corrected.txt")
