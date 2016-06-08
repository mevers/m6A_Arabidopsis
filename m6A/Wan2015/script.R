library(AnnotationDbi);
library(GenomicRanges);
library(GenomicFeatures);

# Load transcriptome
txdb <- loadDb("../../RNAModR/txdb_TAIR10.sqlite");
tx <- exonsBy(txdb, by = "tx", use.names = TRUE);

# Read transcript coordinates
fn <- c("../../m6A/Wan2015/File12_flower.summit.bed",
		"../../m6A/Wan2015/File12_leaf.summit.bed",
		"../../m6A/Wan2015/File12_root.summit.bed");
gr <- GRanges();
for (i in 1:length(fn)) {
	d <- read.table(fn[i]); 
	gr <- append(gr, GRanges(d[, 1], 
				             IRanges(d[, 2], d[, 3])));
}

# Get transcript sequences
library("BSgenome.Athaliana.TAIR.TAIR9");
genome <- get("BSgenome.Athaliana.TAIR.TAIR9");
txSeq <- extractTranscriptSeqs(genome, tx);

# Filter entries
commonTxNames <- intersect(as.character(seqnames(gr)), 
						   names(txSeq));
txSeq <- txSeq[which(names(txSeq) %in% commonTxNames)];
gr <- gr[which(seqnames(gr) %in% commonTxNames)];

# m6A motifs (Wan et al.)
# https://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0839-2
motif <- c("AAACT", "AAACA", "AAACC",
		   "GAACT", "GAACA", "GAACC",
		   "AGACT", "AGACA", "AGACC",
		   "GGACT", "GGACA", "GGACC");

# Get gr sequences and match motifs
# Note: This process is very slow.
gr2 <- gr;
pb <- txtProgressBar(min = 0, max = length(gr), width = 80, style = 3);
for (i in 1:length(gr)) {
	setTxtProgressBar(pb, i);
	txName <- as.character(seqnames(gr[i]));
	seq <- substr(txSeq[[txName]], start(gr[i]), end(gr[i]));
	m <- sapply(motif, function(x) matchPattern(x, seq));
	m <- unlist(IRangesList(sapply(m, ranges)));
	m <- m[order(abs(start(m) - 0.5 * length(seq)))];
	if (length(m) > 0) {
		start(gr2)[i] <- start(gr)[i] + start(m)[1] + 1;
		end(gr2)[i] <- start(gr2)[i];
	} else {
#		print(sprintf("i = %i", i));		
	}
}
close(pb);
gr3 <- gr2[which(width(gr2) == 1)];

# gr3 transcript to genomic coordinates
map <- mapFromTranscripts(gr3, tx);
rtracklayer::export(map, "m6A_Wan2015_TAIR10.bed");