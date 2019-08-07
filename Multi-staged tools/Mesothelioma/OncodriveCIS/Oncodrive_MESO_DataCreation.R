### Data generation for Mesothelioma dataset for OncoDrive-CIS analysis
library(reshape2)
rm(list =ls())
load("/home/anita/Benchmarking/two_omics/MesotheliomaCompleteDataAnalysis/MesotheliomaRawDataset.Rdata")
rm(meso.me)

## Gene expression
topIQR <- sort(apply(meso.ge,1,IQR), decreasing = T)
topIQR2 <- topIQR[topIQR>0.9] #
ge.table <- meso.ge[attr(topIQR2, "names"),]
ge.table <- cbind(rownames(ge.table), ge.table)
colnames(ge.table)[1] <- "GeneID"
ge.table[1:5,1:5]


# Copy number data
meso.cnv <- meso.cnv[attr(topIQR2, "names"),]
# Changing -1 to -2 and 1 to 2 as required by the tool to include heterozygous deletion
table(meso.cnv[,2])
meso.cnv[meso.cnv==-1] <- -2
meso.cnv[meso.cnv==1] <- 2
table(meso.cnv[,2])
#changing the cnv table from wide to long
meso.cnv <- cbind(rownames(meso.cnv),meso.cnv)
colnames(meso.cnv)[1] <- "GeneID"
cnv.table <- melt(meso.cnv) 
head(cnv.table)


## Creating sample file: this file contains a list samples and a column containing 0 or 1, for normal and tumour sample, respectively
# This dataset contains no normal samples
sample.table <- data.frame(colnames(meso.cnv))
sample.table$Code <- 1
sample.table <- sample.table[-1,]

setwd("/home/anita/Benchmarking/two_omics/MesotheliomaCompleteDataAnalysis/OncodriveCIS/")
write.table(ge.table, file= "GeneExp_MESO.tsv", sep = "\t", row.names = F, quote = F)
write.table(cnv.table, file = "CNV_MESO.tsv", sep = "\t", row.names = F, quote = F)
write.table(sample.table, file = "SampleInfo_MESO.tsv", sep = "\t", row.names = F, quote = F)


# Please remove column header in CNV and sample file before running Oncodrive-CIS 
# sed -i -e "1d" SampleInfo_MESO.tsv

# ./oncodrivecis.py --expression=/home/anita/Benchmarking/two_omics/MesotheliomaCompleteDataAnalysis/OncodriveCIS/GeneExp_MESO.tsv --cnv=/home/anita/Benchmarking/two_omics/MesotheliomaCompleteDataAnalysis/OncodriveCIS/CNV_MESO.tsv --samples=/home/anita/Benchmarking/two_omics/MesotheliomaCompleteDataAnalysis/OncodriveCIS/SampleInfo_MESO.tsv --output=CompMESO.Output -p fff
# Parsing the cnv file /home/anita/Benchmarking/two_omics/MesotheliomaCompleteDataAnalysis/OncodriveCIS/CNV_MESO.tsv...
# All the samples from the file /home/anita/Benchmarking/two_omics/MesotheliomaCompleteDataAnalysis/OncodriveCIS/SampleInfo_MESO.tsv has info of expression and cnv.
# There are 0 normal samples and 87 tumor samples
# Folder has been created: CompMESO.Output
# Opening expression file and calculating EIS scores for CNV = 2...
# Traceback (most recent call last):
#   File "./oncodrivecis.py", line 108, in <module>
#   main()
# File "./oncodrivecis.py", line 99, in main
# oncodriveCIS(options, cnv_dict)
# File "./oncodrivecis.py", line 60, in oncodriveCIS
# amp_eisN_dict, amp_eisT_dict = calculate_eis(options, cnv_dict, '2')
# File "/home/anita/Benchmarking/two_omics/MesotheliomaCompleteDataAnalysis/OncodriveCIS/lib.py", line 173, in calculate_eis
# eisT_subdict = get_eis(exp_header_s, tumor_dip_positions, tumor_alt_positions, exp_l_s)
# File "/home/anita/Benchmarking/two_omics/MesotheliomaCompleteDataAnalysis/OncodriveCIS/lib.py", line 69, in get_eis
# eis = eis_calculation(float(exp_l_s[tumor_alt_position]), ref_median, ref_IQ, tumor_alt_IQ)
# File "/home/anita/Benchmarking/two_omics/MesotheliomaCompleteDataAnalysis/OncodriveCIS/lib.py", line 52, in eis_calculation
# eis = (exp_v - ref_median) / (abs(ref_IQ) + abs(tumor_alt_IQ))
# ZeroDivisionError: float division by zero
