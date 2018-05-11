############################################################################
#             Methylation Calculator
# 
# Written by: Leslie Oldfield 
# Date created: December 13th, 2016
# Last edited: December 23rd, 2016
############################################################################

## Set variables
tumour.directory <- "/directory/of/bismark/methylationextractor/output"
tumour.samples <- list.files(tumour.directory, pattern= '.bismark.cov')
area.of.interest <- "chr3:37034616-37036342" # hg19 MLH1 region of interest (promoter region)


# Read through all bismark.cov files for all samples and calculate percentage methylated within area.of.interest
for (i in 1:length(tumour.samples)) {
  df <- read.table(tumour.samples[i], header= FALSE, stringsAsFactors = FALSE) #reads in data
  colnames(df) <- c("chromosome", "start_position", "end_position", "methylation_percentage", "count_methylated", "count_unmethylated")
  regions.subset <- subset(df, chromosome == paste(strsplit(area.of.interest, ":")[[1]][1]) & start_position >= 
                             as.numeric(paste(strsplit((strsplit(area.of.interest, ":")[[1]][2]), "-")[[1]][1])) & start_position <=
                             as.numeric(paste(strsplit((strsplit(area.of.interest, ":")[[1]][2]), "-")[[1]][2])))
  percentage.methylated <- sum(regions.subset$methylation_percentage) / nrow(regions.subset) #calculates percentage methylation
  sample.name <- paste((strsplit(tumour.samples[i], ".")[[1]][1]))   
  write.table(percentage.methylated, paste(tumour.samples[i], "mlh1methylation", "txt", sep= "."), quote= FALSE, col.names = FALSE, row.names = FALSE)
  
}