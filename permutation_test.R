library(dplyr)
ego_result <- read.csv("./result/ego_result_indi.csv",header = FALSE)
colnames(ego_result) <- c("GOID","GOTERM", "gene_id")
data <- merge(ego_result, anno_UCMC, by = "gene_id")
data <- data %>% distinct(gene_id, GOID, age, .keep_all = TRUE)

propotion.UCMC <- function(data){
     data_MC = data[data$age=="MC",]
     GO_MC = as.data.frame(table(data_MC$GOID))
     colnames(GO_MC)[1] = "GOID"
     colnames(GO_MC)[2] = "countMC"
     
     data_UC = data[data$age=="UC",]
     GO_UC = as.data.frame(table(data_UC$GOID))
     colnames(GO_UC)[1] = "GOID"
     colnames(GO_UC)[2] = "countUC"
     
     df = merge(GO_MC,GO_UC,by="GOID",all = TRUE)
     df = replace(df, is.na(df),0)
     df$total = df$countMC + df$countUC
     df$proportion.MC = df$countMC/df$total
     df$proportion.UC = df$countUC/df$total
     return(df)
}

df <- propotion.UCMC(data)

random <- function(data){
     GO_ID = data$GOID
     data$GOID = sample(x=GO_ID,size=length(GO_ID),replace = F)
     return(data)
}

data1 = random(data)

df1 <- propotion.UCMC(data1)

# initial
countvector.MC = as.data.frame(integer(dim(df)[1]))
colnames(countvector.MC)[1] = "countMC"
countvector.UC = as.data.frame(integer(dim(df)[1]))
colnames(countvector.UC)[1] = "countUC"


# Permutation

Permutation <- function(data,n){
     cnt <- 1
     repeat{
          data1 = random(data)
          df1 <- propotion.UCMC(data1)
          loc = match(df$GOID,df1$GOID)
          df1 = df1[c(loc),]
          
          countvector.MC[(df1$proportion.MC >= df$proportion.MC) == TRUE,1] <- countvector.MC[(df1$proportion.MC >= df$proportion.MC) == TRUE,1] + 1
          countvector.UC[(df1$proportion.UC >= df$proportion.UC) == TRUE,1] <- countvector.UC[(df1$proportion.UC >= df$proportion.UC) == TRUE,1] + 1
          
          cnt <- cnt + 1
          if(cnt>n){
               break
          }
     }
     
     countdata <- cbind(countvector.MC,countvector.UC)
     return(countdata)

}

countdata = Permutation(data,10000)

# Pvalue
df$MC.Pvalue <- countdata$countMC / 10000
df$UC.Pvalue <- countdata$countUC / 10000

# GO term UC or MC
df$age = NA
df[df$proportion.UC > df$proportion.MC & df$UC.Pvalue < 0.05,]$age = "UC"
df[df$proportion.MC > df$proportion.UC & df$MC.Pvalue < 0.05,]$age = "MC"

write.csv(df,"./result/permutation_result.csv")


