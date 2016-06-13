# Read IBAQ data
fn <- c("../interactome/160511_proteinGroups_IBAQ_interactome.tsv",
        "../interactome/160511_proteinGroups_IBAQ_totalproteome.tsv");
d <- sapply(fn, function(x) read.delim(x, sep = "\t", header = TRUE, stringsAsFactors = FALSE));

# Check column names of tables
sapply(d, colnames);

# Extract IBAQ values
IBAQ <- list(interactome = d[[1]][, c("Protein.IDs","log10.Median.CCL")],
             proteome = d[[2]][, c("Protein.IDs", "Log10.Median.IBAQ")]);

# Get rid of "Err:502" entries and turn into
# floats
IBAQ <- lapply(IBAQ, function(x)
               cbind.data.frame(x[, 1], as.numeric(gsub("Err:502", "NA", x[, 2]))));

# Get rid of NAs
IBAQ <- lapply(IBAQ, function(x) x[complete.cases(x), ]);

# Match protein ideas
commonProteins <- intersect(IBAQ[[1]][, 1], IBAQ[[2]][, 1]);
data <- cbind.data.frame(
    IBAQ[[1]][match(commonProteins, IBAQ[[1]][, 1]), 2],
    IBAQ[[2]][match(commonProteins, IBAQ[[2]][, 1]), 2]);
rownames(data) <- commonProteins;
colnames(data) <- names(IBAQ);

# Venn diagram
library(limma);
allProteins <- union(IBAQ[[1]][, 1], IBAQ[[2]][, 1]);
vdata <- cbind.data.frame(ifelse(allProteins %in% IBAQ[[1]][, 1], 1, 0),
                          ifelse(allProteins %in% IBAQ[[2]][, 1], 1, 0));
rownames(vdata) <- allProteins;
colnames(vdata) <- names(IBAQ);
vdata <- vennCounts(vdata);
pdf("venn.pdf");
vennDiagram(vdata);
dev.off();

# Scatterplot
pdf("scatterplot.pdf");
plot(data[, 1], data[, 2],
     xlim = c(0, 10), ylim = c(0, 10),
     xlab = sprintf("log10 IBAQ (%s)", colnames(data)[1]),
     ylab = sprintf("log10 IBAQ (%s)", colnames(data)[2]));
abline(a = 0, b = 1, col = "red", lty = 3);
cor <- cor.test(data[, 1], data[, 2]);
text(8, 1, sprintf("r = %4.3f", cor$estimate));
dev.off();

# Histograms
pdf("histogram_proteinsFromInteractome.pdf", width = 14);
par(mfrow = c(1, 2));
hist(data[, 1], breaks = 20,
     xlab = "log10 IBAQ",
     ylab = "Abundance",
     xlim = c(0, 10),
     main = sprintf("IBAQ from %s data\nmean(log10 IBAQ) = %4.3f, N = %i",
         colnames(data)[1], mean(data[, 1]), length(data[, 1])),
     font.main = 1);
hist(data[, 2], breaks = 20,
     xlab = "log10 IBAQ",
     ylab = "Abundance",
     xlim = c(0, 10),
     main = sprintf("IBAQ from %s data\nmean(log10 IBAQ) = %4.3f, N = %i",
         colnames(data)[2], mean(data[, 2]), length(data[, 2])),
     font.main = 1);
ttest1 <- t.test(data[, 1], data[, 2]);
ttest2 <- t.test(data[, 1] - mean(data[, 1]), data[, 2] - mean(data[, 2]));
dev.off();

# Density plot
pdf("densityplot_interactome_proteome.pdf", width = 14);
par(mfrow = c(1, 2));
plot(density(data[, 1]), col = "red",
     xlab = "log10 IBAQ",
     xlim = c(0, 10),
     ylim = c(0, 0.5),
     lwd = 2,
     main = colnames(data)[1],
     font.main = 1);
lines(density(IBAQ[[1]][, 2]), col = "black", lwd = 2);
legend("topleft",
       c(sprintf("all proteins in interactome (N = %i)",
                 length(IBAQ[[1]][, 2])),
         sprintf("proteins in interactome and proteome (N = %i)",
                 length(data[, 1]))),
       col = c("black", "red"), lwd = c(2, 2), bty = "n");
plot(density(data[, 2]), col = "red",
     xlab = "log10 IBAQ",
     xlim = c(0, 10),
     ylim = c(0, 0.5),
     lwd = 2,
     main = colnames(data)[2],
     font.main = 1);
lines(density(IBAQ[[2]][, 2]), col = "black", lwd = 2);
legend("topleft",
       c(sprintf("all proteins in proteome (N = %i)",
                 length(IBAQ[[2]][, 2])),
         sprintf("proteins in proteome and interactome (N = %i)",
                 length(data[, 2]))),
       col = c("black", "red"), lwd = c(2, 2), bty = "n");
dev.off();


# Read protein IDs for different groups
# Need to do manually, because of different file formats
groupIDs <- list();
groupIDs[[length(groupIDs) + 1]] <- as.vector(IBAQ[[2]][, 1]);
groupIDs[[length(groupIDs) + 1]] <- as.vector(read.table("../interactome/GO_RNAbinding.tsv",
                                               header = FALSE)[, 2]);
groupIDs[[length(groupIDs) + 1]] <- as.vector(read.delim("../interactome/ss_candidates.tsv",
                                               header = TRUE)[, 3]);
groupIDs[[length(groupIDs) + 1]] <- as.vector(read.delim("../interactome/ns_candidates.tsv",
                                               header = TRUE)[, 3]);
groupIDs[[length(groupIDs) + 1]] <- as.vector(read.delim("../interactome/unknown_RBD.tsv",
                                               header = TRUE)[, 2]);
names(groupIDs) <- c(names(IBAQ)[2],
                     "GO:RNA binding (from proteome)",
                     "ss candidates (interactome)",
                     "ns candidates (interactome)",
                     "unknown RBD (interactome)");


# Density plot of IBAQ values based on proteome data,
# split up into different groups
library(RColorBrewer);
pdf("densityplot_proteome.pdf", width = 7);
par(mfrow = c(1, 1));
plot(density(IBAQ[[2]][, 2]), col = "black",
     xlab = "log10 IBAQ",
     xlim = c(0, 10),
     ylim = c(0, 0.6),
     lwd = 2,
     main = "log10 IBAQ values from proteome data for groups as indicated",
     font.main = 1);
alpha <- 255;
col <- c(
    rgb(95, 178, 51, alpha = alpha, maxColorValue = 255),
    rgb(106, 127, 147, alpha = alpha, maxColorValue = 255),
    rgb(245, 114, 6, alpha = alpha, maxColorValue = 255),
    rgb(235, 15, 19, alpha = alpha, maxColorValue = 255),
    rgb(143, 47, 139, alpha = alpha, maxColorValue = 255),
    rgb(19, 150, 219, alpha = alpha, maxColorValue = 255));
col <- col[1:(length(groupIDs) - 1)];
for (i in 2:length(groupIDs)) {
    lines(density(IBAQ[[2]][which(IBAQ[[2]][, 1] %in% groupIDs[[i]]), 2]),
          col = col[i - 1], type = "l", lwd = 2);
}
legend("topleft",
       sprintf("%s (N = %i)", names(groupIDs), sapply(groupIDs, length)),
       col = c("black", col),
       lwd = rep(2, length(groupIDs)),
           bty = "n");
dev.off();
