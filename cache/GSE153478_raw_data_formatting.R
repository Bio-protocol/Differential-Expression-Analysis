#----------------------------
####    GSE153478 raw data formatting
#----------------------------
setwd('~/path/to/folder/GSE153478/')
library(dplyr)
library(readr)

# download GSE153478 raw data
for (i in 9:44) {
  input<- paste('GSM46449', formatC(i, width = 2, flag = '0'), sep = '')
  getGEOSuppFiles(input, makeDirectory = F, baseDir = 'GSE153478',
                  fetch_files = TRUE, filter_regex = NULL)
}

setwd('~/path/to/file/GSE153478_series_matrix.txt')
gse <- getGEO(filename="GSE153478_series_matrix.txt")
acc <- pData(gse)[,grepl("geo_accession",names(pData(gse)))]
id <- pData(gse)[,grepl("title",names(pData(gse)))]

pdata<- data.frame(cbind(acc, id))

#The data.frame was cleaned by replacing long character strings with shorter ones, 
#and turning the minute variable into a factor with the correctly ordered levels.

time<- do.call(rbind,strsplit(pdata$id,'_'))[,2]
t<- sub('hr', '', time)
t[which(t=='control')]<- 0

pdataclean <- data.frame(strain=ifelse(grepl("CASP1",pdata$id),"double","single"),
                         hour=t,
                         replicate=paste0("r",sub("(.*)_rep","",pdata$id)),
                         id= pdata$id,
                         row.names=pdata$acc)

pdataclean$strain <- relevel(factor(pdataclean$strain), "single")
pdataclean$hour <- factor(pdataclean$hour, levels=c("0","2","4","8","16","24"))
str(pdataclean)
coldata <- DataFrame(pdataclean)


setwd('~/path/to/folder/GSE153478/')
file_list<- list.files()
file_list

df <- list.files(path=".", full.names = TRUE) %>% 
  lapply(read.table) %>% 
  bind_cols

df2<- df[,seq(2, ncol(df), by=2)] ## Select even rows
colnames(df2)<- do.call(rbind,strsplit(file_list,'_'))[,1]

count<- df2[-1, ]
#rownames(count)<- df[,1][-1]

reads<- data.frame(lapply(count, as.numeric))
rownames(reads)<- df[,1][-1]
head(reads)
class(reads)
str(reads)

library("annotate")
metadata <- pmid2MIAME("33686223")
metadata@url <- "http://www.ncbi.nlm.nih.gov/pubmed/33686223"

library("SummarizedExperiment")
arab_time <- SummarizedExperiment(assays=SimpleList(counts=as.matrix(reads)),
                                  colData=coldata,
                                  metadata=list(metadata))

head(counts)
class(counts)
str(counts)

setwd('~/bigdata/DEG/')
save(arab_time, file="arab_time.RData")
load(file = 'arab_time.RData')

saveRDS(arab_time, file = "arab_time.rds")

dim(assays(arab_time)$counts)
assays(arab_time)$counts[1:5,1:5]
