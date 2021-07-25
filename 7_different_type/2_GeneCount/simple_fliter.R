## 1. loading data

# fliter hog data
Am_hog_corrected <- read.table("./1_fliter/Am_hog_corrected.txt")
colnames(Am_hog_corrected)[2] <- "phylostrata_cor"
colnames(Am_hog_corrected)[3] <- "hog_cor"

Ej_hog_corrected <- read.table("./1_fliter/Ej_hog_corrected.txt")
colnames(Ej_hog_corrected)[2] <- "phylostrata_cor"
colnames(Ej_hog_corrected)[3] <- "hog_cor"

Go_hog_corrected <- read.table("./1_fliter/Go_hog_corrected.txt")
colnames(Go_hog_corrected)[2] <- "phylostrata_cor"
colnames(Go_hog_corrected)[3] <- "hog_cor"

# raw hog data
Am_hog_raw <- read.table("./0_database/Am_oma_result.1.txt")
colnames(Am_hog_raw) <- c("transcript", "hog", "id", "phylostrata")

Ej_hog_raw <- read.table("./0_database/Ej_oma_result.1.txt")
colnames(Ej_hog_raw) <- c("transcript", "hog", "id", "phylostrata")

Go_hog_raw <- read.table("./0_database/Go_oma_result.1.txt")
colnames(Go_hog_raw) <- c("transcript", "hog", "id", "phylostrata")

# transcript map gene 
Am.trans2gene <- read.table("./0_database/Am.gene_trans_map.txt")
colnames(Am.trans2gene) <- c('gene', 'transcript')
Ej.trans2gene <- read.table("./0_database/Ej.gene_trans_map.txt")
colnames(Ej.trans2gene) <- c('gene', 'transcript')
Go.trans2gene <- read.table("./0_database/Go.gene_trans_map.txt")
colnames(Go.trans2gene) <- c('gene', 'transcript')

## 2. Merge hog_cor cols and transcript map 
Am_hog_corrected <- merge(Am_hog_corrected, Am.trans2gene, by = "gene")
Ej_hog_corrected <- merge(Ej_hog_corrected, Ej.trans2gene, by = "gene")
Go_hog_corrected <- merge(Go_hog_corrected, Go.trans2gene, by = "gene")

## 3. Merge raw data and transcript map
Am_hog_raw <- merge(Am_hog_raw, Am.trans2gene, by = "transcript")
Ej_hog_raw <- merge(Ej_hog_raw, Am.trans2gene, by = "transcript")
Go_hog_raw <- merge(Go_hog_raw, Am.trans2gene, by = "transcript")

write.table(Am_hog_raw, "./1_fliter/Am_hog_raw.txt")
write.table(Ej_hog_raw, "./1_fliter/Ej_hog_raw.txt")
write.table(Go_hog_raw, "./1_fliter/Go_hog_raw.txt")
