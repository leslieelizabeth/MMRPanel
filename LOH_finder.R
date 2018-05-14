##############################################################
#####              Detecting Somatic LOH                 #####
#####                                                    #####
#####          Last edited: March 11th, 2017             #####
##############################################################

###############################################################
## Note: check, is the T or N colum first? T isn't always second in the VCF file
## Change col.sel to 1 or 2 depending on the location of the T column in Mutect
## 2= tumour is in column after normal, 1= tumour is before normal sample in column

library(VariantAnnotation)

#### Read in a file

file <- "/path/to/VCF/file.vcf"

#### Load functions
getMutectAf <- function(sample, col.sel=2){ #change col.sel depending on if T or N is first 
  print("** Notice: MuTect mode is only set up for single sample run at the moment")
  geno.gt <- geno(sample)$GT[,col.sel]   # Genotype call
  geno.dp <- geno(sample)$DP[,col.sel]   # Total depth
  geno.ad <- geno(sample)$AD[,col.sel]   # Depth of A and B allele
  
  #Calculates the Allelic Fraction
  allele.cnt.df <- do.call("rbind", geno.ad)
  af.df <- allele.cnt.df / geno.dp  #convert to allelic fraction
  af.df <- cbind(af.df, allele.cnt.df, 
                 geno.dp, gsub(":.+", "", rownames(af.df)))
  colnames(af.df) <- c("a", "b", "raw_a", "raw_b", "depth", "chr")
  
  return(af.df)
}

# Removes factors and adds genomic positions to AF dataframe
formatAfDf <- function(af.df){
  # Double check the formats and factors to ensure numeric and not character
  for(each.col in c("a", "b", "raw_a", "raw_b")){
    af.df[,each.col] <- as.numeric(as.character(af.df[,each.col]))
  }
  af.df <- as.data.frame(af.df)
  
  # Add genomic position
  m <- regexpr(":.+?_",rownames(af.df))
  chr.pos <- gsub("[:_]", "", regmatches(rownames(af.df), m))
  af.df[,'pos'] <- chr.pos
  return(af.df)
}

########## Extracting allelic fractions
sample <- readVcf(file, "hg19")
af.df <- getMutectAf(sample, col.sel=2) #change to 1 if tumour is first
af.df <- formatAfDf(af.df)
na.rows <- apply(is.na(af.df), 1, any)
af.df <- af.df[!na.rows,]

shorter.name <- basename(file)
names <- gsub(".processed.bam.vcf", "",shorter.name)
sample.name <- paste(names) ##do the basename thing 
write.table(af.df, paste(sample.name, "allelicfractions", "txt", sep = "."), quote= FALSE, col.names = TRUE, row.names = TRUE, sep= "\t")


########## Plotting for the MMR genes

## Extract allelic fractions for the positions in the MMR genes
msh2.epcam <- af.df[af.df$chr == 'chr2' & af.df$pos >= 47596287 & af.df$pos <= 47710367, ]
msh2 <- af.df[af.df$chr == 'chr2' & af.df$pos >= 47630206 & af.df$pos <= 47710367, ]
msh6 <- af.df[af.df$chr == 'chr2' & af.df$pos >= 48010221 & af.df$pos <= 48034092, ]
mlh1 <- af.df[af.df$chr == 'chr3' & af.df$pos >= 37034841 & af.df$pos <= 37092337, ]
pms2 <- af.df[af.df$chr == 'chr7' & af.df$pos >= 6012870 & af.df$pos <= 6048737, ]
coverage.limit <- c(1:5) 
msh2.epcam.cov <- msh2.epcam[!msh2.epcam$depth %in% coverage.limit, ]
msh2.cov <- msh2[!msh2$depth %in% coverage.limit, ]
msh6.cov <- msh6[!msh6$depth %in% coverage.limit, ]
mlh1.cov <- mlh1[!mlh1$depth %in% coverage.limit, ]
pms2.cov <- pms2[!pms2$depth %in% coverage.limit, ]

## Plot each figure
# msh2 + epcam 
png(paste(sample.name, "msh2.epcam", "png", sep= "."))
par(mfrow=c(2,1))
plot(as.numeric(as.character(msh2.epcam.cov$pos)),as.numeric(as.character(msh2.epcam.cov$b)), col= "red", xlab= "Position in MSH2_EPCAM", 
     ylab= "Allelic Fraction",ylim= c(0, 1.0), main= "EPCAM+MSH2, Coverage Adjusted AF")
points(as.numeric(as.character(msh2.epcam.cov$pos)), as.numeric(as.character(msh2.epcam.cov$a)), col= "blue")
abline(h=0.5,col='black',lwd=3)
plot(density(as.numeric(as.character(msh2.epcam.cov$b))), col="blue", xlab= "Allelic Fraction", 
     ylab= "Density", main= "EPCAM+MSH2, Coverage Adjusted DP")
lines(x=density(as.numeric(as.character(msh2.epcam.cov$a)))$x, y=density(as.numeric(as.character(msh2.epcam.cov$a)))$y, col="red")
dev.off() 

# msh2 only
png(paste(sample.name, "msh2", "png", sep= "."))
par(mfrow=c(2,1))
plot(as.numeric(as.character(msh2.cov$pos)),as.numeric(as.character(msh2.cov$b)), col= "red", xlab= "Position in MSH2", 
     ylab= "Allelic Fraction",ylim= c(0, 1.0), main= "MSH2, Coverage Adjusted AF")
points(as.numeric(as.character(msh2.cov$pos)), as.numeric(as.character(msh2.cov$a)), col= "blue")
abline(h=0.5,col='black',lwd=3)
plot(density(as.numeric(as.character(msh2.cov$b))), col="blue", xlab= "Allelic Fraction", 
     ylab= "Density", main= "MSH2, Coverage Adjusted DP")
lines(x=density(as.numeric(as.character(msh2.cov$a)))$x, y=density(as.numeric(as.character(msh2.cov$a)))$y, col="red")
dev.off() 

# msh6 
png(paste(sample.name, "msh6", "png", sep= "."))
par(mfrow=c(2,1))
plot(as.numeric(as.character(msh6.cov$pos)),as.numeric(as.character(msh6.cov$b)), col= "red", xlab= "Position in MSH6", 
     ylab= "Allelic Fraction",ylim= c(0, 1.0), main= "MSH6, Coverage Adjusted AF")
points(as.numeric(as.character(msh6.cov$pos)), as.numeric(as.character(msh6.cov$a)), col= "blue")
abline(h=0.5,col='black',lwd=3)
plot(density(as.numeric(as.character(msh6.cov$b))), col="blue", xlab= "Allelic Fraction", 
     ylab= "Density", main= "MSH6, Coverage Adjusted DP")
lines(x=density(as.numeric(as.character(msh6.cov$a)))$x, y=density(as.numeric(as.character(msh6.cov$a)))$y, col="red")
dev.off() 

# pms2
png(paste(sample.name, "pms2", "png", sep= "."))
par(mfrow=c(2,1))
plot(as.numeric(as.character(pms2.cov$pos)),as.numeric(as.character(pms2.cov$b)), col= "red", xlab= "Position in PMS2", 
     ylab= "Allelic Fraction",ylim= c(0, 1.0), main= "PMS2, Coverage Adjusted AF")
points(as.numeric(as.character(pms2.cov$pos)), as.numeric(as.character(pms2.cov$a)), col= "blue")
abline(h=0.5,col='black',lwd=3)
plot(density(as.numeric(as.character(pms2.cov$b))), col="blue", xlab= "Allelic Fraction", 
     ylab= "Density", main= "PMS2, Coverage Adjusted DP")
lines(x=density(as.numeric(as.character(pms2.cov$a)))$x, y=density(as.numeric(as.character(pms2.cov$a)))$y, col="red")
dev.off() 

# mlh1 
png(paste(sample.name, "mlh1", "png", sep= "."))
par(mfrow=c(2,1))
plot(as.numeric(as.character(mlh1.cov$pos)),as.numeric(as.character(mlh1.cov$b)), col= "red", xlab= "Position in MLH1", 
     ylab= "Allelic Fraction",ylim= c(0, 1.0), main= "MLH1, Coverage Adjusted AF")
points(as.numeric(as.character(mlh1.cov$pos)), as.numeric(as.character(mlh1.cov$a)), col= "blue")
abline(h=0.5,col='black',lwd=3)
plot(density(as.numeric(as.character(mlh1.cov$b))), col="blue", xlab= "Allelic Fraction", 
     ylab= "Density", main= "MLH1, Coverage Adjusted DP")
lines(x=density(as.numeric(as.character(mlh1.cov$a)))$x, y=density(as.numeric(as.character(mlh1.cov$a)))$y, col="red")
dev.off() 



