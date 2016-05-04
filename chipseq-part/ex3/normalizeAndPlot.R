# Normalization, comparing N1 and T1 and making plots in R

# Import the package preprocessCore (used for quantile normalization)
library(preprocessCore)

# Read the files
N1cCHIP = read.table("N1_CHIP_countR.txt", sep="\t")
T1cCHIP = read.table("T1_CHIP_countR.txt", sep="\t")
T1cINPUT = read.table("T1_INPUT_countR.txt", sep="\t")
N1cINPUT = read.table("N1_INPUT_countR.txt", sep="\t")

# Set column names
colnames(N1cCHIP) = c("peak", "N1cCHIPcount")
colnames(T1cCHIP) = c("peak", "T1cCHIPcount")
colnames(T1cINPUT) = c("peak", "T1cINPUTcount")
colnames(N1cINPUT) = c("peak", "N1cINPUTcount")

# Merge the four count tables into one, and setting missing values (no reads) to 0
readcounts = merge(N1cCHIP, N1cINPUT, by="peak", ALL=TRUE)
readcounts = merge(readcounts,T1cCHIP, by="peak", ALL=TRUE)
readcounts = merge(readcounts,T1cINPUT, by="peak", ALL=TRUE)
readcounts[is.na(readcounts)] = 0

# Add four columns with size normalization counts
#################################################################### The numbers must be dynamically
readcounts$N1cCHIPsizenorm = readcounts$N1cCHIPcount/0.836576
readcounts$N1cINPUTsizenorm = readcounts$N1cINPUTcount/0.924283
readcounts$T1cCHIPsizenorm = readcounts$T1cCHIPcount/1.065605
readcounts$T1cINPUTsizenorm = readcounts$T1cINPUTcount/1.024771

# Get a table with the CHIP counts for normalization
fornorm = cbind(readcounts$N1cCHIPcount, readcounts$T1cCHIPcount)

# Normalize
normcounts = normalize.quantiles(as.matrix(fornorm), copy=TRUE)

# Add two columns with the normalized counts in the readcount table
readcounts$N1cCHIPnorm = normcounts[,1]
readcounts$T1cCHIPnorm = normcounts[,2]

# Add columns with the logfold change values for each of the three CHIP-counts
readcounts$logfoldChangeRaw = log2( (readcounts$N1cCHIPcount+1)/(readcounts$T1cCHIPcount+1) )
readcounts$logfoldChangeSnorm = log2( (readcounts$N1cCHIPsizenorm+1)/(readcounts$T1cCHIPsizenorm+1) )
readcounts$logfoldChangeQnorm = log2( (readcounts$N1cCHIPnorm+1)/(readcounts$T1cCHIPnorm+1) )

#Sort the readcount table based on logfoldChange for the quantile normalized counts and write it to file
ordered_readcounts = readcounts[order(readcounts$logfoldChangeQnorm),]
write.table(ordered_readcounts, "ordered_readcounts.txt", sep="\t", quote=FALSE, row.names=FALSE)

#Plot CHIP vs INPUT for N1 using size normalized counts
pdf("readcount_plot_N1c.pdf")
plot(log2(readcounts$N1cCHIPsizenorm+1), log2(readcounts$N1cINPUTsizenorm+1), main="Size normalized readcounts N1c", xlab="CHIP", ylab="INPUT", pch=16, xlim=c(1,12), ylim=c(1,12))
abline(0,1,col="red")
dev.off()

#Plot CHIP vs INPUT for T1 using size normalized counts
pdf("readcount_plot_T1c.pdf")
plot(log2(readcounts$T1cCHIPsizenorm+1), log2(readcounts$T1cINPUTsizenorm+1), main="Size normalized readcounts T1c", xlab="CHIP", ylab="INPUT", pch=16, xlim=c(1,12), ylim=c(1,12))
abline(0,1,col="red")
dev.off()

#Plot the quantile normalized counts of N1c vs T1c
pdf("readcount_plot_CHIP_N1c_T1c.pdf")
plot(log2(readcounts$N1cCHIPnorm+1), log2(readcounts$T1cCHIPnorm+1), main="N1c vs T1c", xlab="N1c log2(normalized counts)", ylab="T1c log2(normalised counts)", pch=16, xlim=c(1,13), ylim=c(1,13))
abline(0,1,col="red")
dev.off()

#Plot histograms of logfoldChange for raw count, size normalized count and quantile normalized #counts
pdf("histogram_logfoldChange_quantile_normalized_counts.pdf")
hist(readcounts$logfoldChangeQnorm, breaks=60, col="yellow", main="log2( (N1cCHIPnorm+1)/(T1cCHIPnorm+1) )")
abline(v=0, col="red", lwd=2)
dev.off()

pdf("histogram_logfoldChange_size_normalized_counts.pdf")
hist(readcounts$logfoldChangeSnorm, breaks=60, col="yellow", main="log2( (N1cCHIPsizenorm+1)/(T1cCHIPsizenorm+1) )")
abline(v=0, col="red", lwd=2)
dev.off()

pdf("histogram_logfoldChange_raw_counts.pdf")
hist(readcounts$logfoldChangeRaw, breaks=60, col="yellow", main="log2( (N1cCHIPcounts+1)/(T1cCHIPcounts+1) )")
abline(v=0, col="red", lwd=2)
dev.off()
