#Iman Safaee
#System Biology Project
#Gene expression in Normal and Tumor kidney cancer

#install and library needed packages
install.packages("BiocManager")
BiocManager :: install("oligo")
BiocManager :: install("Biobase", force = TRUE)
BiocManager :: install("Biostring", force = TRUE)
BiocManager :: install("affy", force = TRUE)

#bringing needed packages to R
library(affy)
library(oligo)

#giving the path of our data to R
celpath <- "F:/Daneshghah/term 7/System biology/Project"

#list our data
list <- list.files(celpath,full.names = TRUE)
list

#reading our data (oligo)
data <- read.celfiles(list)
data

#reading our data (affy)
data.affy <- ReadAffy(celfile.path = celpath)
data.affy

#showing all perfect math
pms.affy <- affy:: pm(data.affy)
pms.oligo <- oligo:: pm(data)

#our data is complex
str(data)

#separating phenodata from our data
ph <- data@phenoData
ph
#or use this function
phenoData(data)

#our phenodata
str(ph)
ph@data

#changing column to Normal and Tuomor
ph@data[,1] <- c("Normal1","Tumor1","Normal2","Tumor2","Normal3","Tumor3","Normal4"
                 ,"Tumor4","Normal5","Tumor5")
ph@data

#extracting the CDF of our chip
cdfName(data.affy)

#extracting probe sets
f.names <- featureNames(data)
head(f.names)

f.names.affy <- featureNames(data.affy)
head(f.names.affy)

#extracting probe names
p.names <- probeNames(data)
head(p.names)


p.names.affy <- affy::probeNames(data.affy)
head(p.names.affy)

#probe IDs of a probe set
affy:: pm(data.affy,"1007_s_at")


#chip image of one probe set
oligo::image(data[,1])

#saving image
setwd("F:/Daneshghah/term 7/System biology/Project")
jpeg("sample.jpeg")
oligo::image(data[,1])
dev.off()

pdf("sample.pdf")
oligo::image(data[,1])
dev.off()

#saving 10 image together
for(i in 1:10){
  name = paste("image","i",".jpg", sep = "")
  jpeg(name)
  oligo::image(data[,i],main = ph@data$sample[i])
  dev.off()
}

#histogram
hist(data[,1:10], lwd = 2, lty = 1, which = 'pm', col = 'red',
     ylab= 'Density', xlab = 'log2 intensities', main = 'Histogram of raw data' )

#showing difference between Normal and Tumor Cells
color <- c('green','black','green','black','green','black','green','black'
           ,'green','black')
hist(data[,1:10], lwd = 4, lty = 2, which = 'pm', col = color,
     ylab= 'Density', xlab = 'log2 intensities', main = 'Histogram of raw data' )

#Ready to Analyze :)

#Quality control before normalization
#Histogram and Pseudo image

#boxplot
boxplot(data, which = 'pm', col = 'dark red', names = ph@data$sample)

#saving boxplot image
boxplot <- "boxplot.jpg"
jpeg(boxplot)
boxplot(data, which = 'pm', col = color , names = ph@data$sample)
dev.off()

#MA Plot

MAplot(data,which=1)#oligo
ma.plot(data.affy,which=1)#affy

for(i in 1:10)
{
  name = paste("MA Plot", "i", ".jpg", sep = "")
  jpeg(name)
  oligo::MAplot(data, which = i)
  dev.off()
}


#Normalization
#1 background correction
#2 quantile normalization
#3 summerization

data.rma <- oligo:: rma(data)
data.matrix = exprs(data.rma)
data.matrix

#Quality control after normalization
data.matrix["1007_s_at",]
pm(data,"1007_s_at")

#boxplot
name = "boxplotnorm.jpg" 
jpeg(name)
boxplot(data.matrix,col=color,names=ph@data$sample)
dev.off()

#histogram
hist(data.rma,col=color,lty=1)

#MA Plot
for (i in 1:10) 
{ 
  name = paste("MAplotnorm",i,".jpg",sep="") 
  jpeg(name) 
  oligo::MAplot(data.rma,which=i) 
  dev.off() 
}

#Analyzing our data after normalization 

#install Limma package
BiocManager::install("limma") 
library(limma)

#Defining two Tumor an Normal value
Cancer <- c(0,1,0,1,0,1,0,1,0,1)
Normal <- c(1,0,1,0,1,0,1,0,1,0)  

#bind them together
Design <- cbind(Cancer,Normal)  


fit <- lmFit (data.rma, design = Design)
fit <- eBayes(fit) 
contrast.matrix <- makeContrasts(Cancer-Normal, levels=Design) 
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2) 


results <- decideTests(fit2,adjust.method="BH", p.value=0.05, lfc=1) 
sum(results > 0,na.rm = T) 
sum(results < 0,na.rm = T) 


write.table(results, "cancer-normal p-value.txt", sep="\t") 
file=paste("LIMMA cancer-normal p value 0.01 lfc=0.txt", ".txt",sep = "_") 
write.fit(fit2, results, file, adjust="BH", sep="\t") 

#Genes with high expression
upgene <- results@.Data[results@.Data>0,]
upgene
upgene <- as.data.frame(upgene)

class(results)
str(results)


#Genes with low expression
downgene <- results@.Data[results@.Data<0,]
downgene
downgene <- as.data.frame(downgene)

#binding them
colnames(downgene) <- "upgene"
DEG <- rbind(upgene,downgene)

#export
write.csv(DEG, "Diffrenetially Expressed Genes.csv")

#changing prob IDs to Genes names
BiocManager::install("biomaRt")
library(biomaRt)

listMarts()
ensembl=useMart("ensembl")


ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl") 

listDatasets(ensembl)

#filter
filters = listFilters(ensembl) 
filters[1:5,]
filters
#attributes
attributes = listAttributes(ensembl) 
attributes[1:5,] 
attributes
#value
affid<-rownames(DEG)

anot <-  getBM(attributes=c('affy_hg_u133a','external_gene_name',
'ensembl_gene_id','description'), filters = 'affy_hg_u133a'
,values = affid,mart = ensembl)
anot

#The End


