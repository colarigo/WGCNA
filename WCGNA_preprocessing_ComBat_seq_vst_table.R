library(tidyverse)
library(data.table)
library(DESeq2)
library(purrr)
library(GGally)
library(xlsx)
library(tidyr)
library(sva)
library(dplyr)

### Rename count files to correspond to info
dd <- list.files(path = getwd(), pattern= "*.txt", full.names = F)
df_count_matrix <- tibble(file_name = dd) %>% 
                    mutate(file_cont = map(file_name,fread,data.table = F))  %>%
                    unnest(cols = c(file_cont)) %>% 
                    mutate(file_name = gsub(pattern="*.txt",replacement="",file_name))  %>% 
                    spread(key = file_name , value = V2) %>%
                    dplyr::slice(6:nrow(.)) 

#--- filter genes if read count is less than 10 in all columns
df_count_matrix_filtered =as.data.frame(df_count_matrix)
row.names(df_count_matrix_filtered) <- df_count_matrix_filtered$V1 
df_count_matrix_filtered <- subset(df_count_matrix_filtered, select=-c(V1))
df_count_matrix_filtered <- df_count_matrix_filtered %>% filter_all(any_vars(. >= 10))
#df_count_matrix_filtered <- df_count_matrix_filtered %>% filter_all(any_vars(. >= 5*ncol(df_count_matrix_filtered)/2))


#Remove low expressed genes
df_count_matrix_filtered = as.data.frame(df_count_matrix) %>% mutate(Median=apply(df_count_matrix[,-1],1,median)) %>% filter(Median >=10) %>% select(!Median)
row.names(df_count_matrix_filtered) <- df_count_matrix_filtered$V1 
df_count_matrix_filtered <- subset(df_count_matrix_filtered, select=-c(V1))

# filter genes using template
template = read.csv("path/to/template.txt")
template_row = template$Gene
df_count_matrix_filtered = as.data.frame(df_count_matrix %>% filter(df_count_matrix$V1 %in% template_row))
row.names(df_count_matrix_filtered) <- df_count_matrix_filtered$V1 
df_count_matrix_filtered <- subset(df_count_matrix_filtered, select=-c(V1))


#Check for NAs
NAs = any(is.na(df_count_matrix_filtered))
if(NAs){cat(c('Found',sum(is.na(df_count_matrix_filtered)),'Missing Data Values/n'),sep=' ')}
df_count_matrix_filtered[is.na(df_count_matrix_filtered)] <- 0

#Optional save matrix
#write.csv(df_count_matrix_filtered, file = "TB_count_matrix.csv")

 #--- add conditions
condition <- as.data.frame(read.csv(list.files(pattern = "*_batch.csv|*_info.csv"),h = T))
condition=subset(condition, condition$Sample %in% colnames(df_count_matrix_filtered))
#Get same order as condition
# check order
df_count_matrix_filtered <- df_count_matrix_filtered[,condition$Sample]

batch = as.factor(condition$Time_point)
#group = as.factor(condition$Resistance)
adjusted_counts <-ComBat_seq(as.matrix(df_count_matrix_filtered), batch = batch, group = NULL, covar_mod=NULL, full_mod=TRUE, 
                       shrink=FALSE, shrink.disp=FALSE, gene.subset.n=NULL)


colData <- as.data.frame(colnames(adjusted_counts))
colData$condition <- condition$Plasmopara[match(colData[,1], condition$Sample)]
colnames(colData) <- c("colData", "condition")
colData$condition=factor(colData$condition, levels=c("0", "1"))

dds <- DESeqDataSetFromMatrix(countData = adjusted_counts,
                                        colData = colData,design = ~ condition)
vsd <- vst(dds, blind=TRUE)
write.csv(assay(vsd), file = paste0(gsub("_.*", "", dd[1]),"_ComBat_seq_vst.txt"))
write.csv(assay(vsd)[,condition$Sample[!gsub("\\d-\\d{2}h_Pv", "", condition$Sample) %in% "C"]], file = paste0(gsub("_.*", "", dd[1]),"inoc_ComBat_seq_vst.txt"))

