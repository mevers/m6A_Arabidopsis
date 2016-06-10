# Read in data
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
     main = sprintf("%s\nmean(log10 IBAQ) = %4.3f", colnames(data)[1], mean(data[, 1])),
     font.main = 1);
hist(data[, 2], breaks = 20,
     xlab = "log10 IBAQ",
     ylab = "Abundance",
     xlim = c(0, 10),
     main = sprintf("%s\nmean(log10 IBAQ) = %4.3f", colnames(data)[2], mean(data[, 2])),
     font.main = 1);
ttest1 <- t.test(data[, 1], data[, 2]);
ttest2 <- t.test(data[, 1] - mean(data[, 1]), data[, 2] - mean(data[, 2]));
dev.off();


pdf("densityplot.pdf", width = 14);
par(mfrow = c(1, 2));
plot(density(data[, 1]), col = "black",
     xlab = "log10 IBAQ",
     xlim = c(0, 10),
     ylim = c(0, 0.5),
     main = colnames(data)[1],
     font.main = 1);
lines(density(IBAQ[[1]][, 2]), col = "red");
legend("topleft", c("all proteins in interactome", "proteins in interactome and proteome "),
       col = c("red", "black"), lwd = c(1, 1), bty = "n");
plot(density(data[, 2]), col = "black",
     xlab = "log10 IBAQ",
     xlim = c(0, 10),
     ylim = c(0, 0.5),
     main = colnames(data)[2],
     font.main = 1);
lines(density(IBAQ[[2]][, 2]), col = "red");
legend("topleft", c("all proteins in proteome", "proteins in proteome and interactome"),
       col = c("red", "black"), lwd = c(1, 1), bty = "n");
dev.off();
