# Load the library
library(RNAModR);

# Read in genomic loci of sites from BED file
bedFile <- "m6Aseq_m6A_Luo2014_TAIR10.bed";
sites <- ReadBED(bedFile, collapseRange = TRUE);

# Map sites to transcriptome
posSites <- SmartMap(sites, id = "m6A", refGenome = "TAIR10");

# Generate null distribution
negSites <- GenerateNull(posSites, method = "ntAbund", nt = "A");

# Restrict sites to sites in 5'UTR, CDS, 3'UTR
posSites <- FilterTxLoc(posSites, c("5'UTR", "CDS", "3'UTR"));
negSites <- FilterTxLoc(negSites, c("5'UTR", "CDS", "3'UTR"));

# Plot sequence logos
pdf("figures/seqLogo.pdf", width = 15, height = 10);
PlotSeqLogo(posSites);
PlotSeqLogo(negSites);
dev.off();

# Plot distribution of sites across tx sections
pdf("figures/sectionAbundance.pdf", height = 5, width = 15);
par(mfrow = c(1,2));
PlotSectionDistribution(posSites);
PlotSectionDistribution(negSites);
dev.off();

# Plot enrichment of sites across tx sections
pdf("figures/sectionEnrichment.pdf");
PlotSectionEnrichment(posSites, negSites);
dev.off();

# Plot spatial abundance within tx sections
# Note: This will take a long time.
pdf("figures/spatialAbundance.pdf", width = 15, height = 10);
par(mfrow = c(2,1));
PlotSpatialDistribution(posSites, filter = c("CDS", "3'UTR"), absolute = TRUE, binWidth = 20);
PlotSpatialDistribution(negSites, filter = c("CDS", "3'UTR"), absolute = TRUE, binWidth = 20);
dev.off();

# Plot spatial enrichment within tx sections
pdf("figures/spatialEnrichment.pdf", width = 15, height = 10);
PlotSpatialEnrichment(posSites, negSites, filter = c("CDS", "3'UTR"));
dev.off();

# Get all PAS motif loci from transcriptome
motif <- c("AATAAA", "ATTAAA", "AGTAAA",
           "TATAAA", "AAGAAA", "AATACA",
           "AATATA", "CATAAA", "AATGAA",
           "GATAAA", "ACTAAA", "AATAGA");
motifSites <- GetMotifLoc(motif, filter = c("CDS", "3'UTR"), refGenome = "TAIR10");
WriteTxLocToBED(motifSites, "annotation/PASmotifs.bed");

# Plot enrichment of relative distances m6A/non-m6A to PAS motifs in tx
pdf("figures/relDistEnrichment_sites_PASmotif.pdf", width = 15, height = 10);
PlotRelDistEnrichment(posSites, negSites, motifSites, flank = 200, binWidth = 10);
dev.off();

# Venn diagram of overlaps between m6A and PAS motifs
pdf("figures/overlap_sites_PASmotif.pdf", width = 9, height = 5);
PlotOverlap(FilterTxLoc(posSites, filter = c("CDS", "3'UTR")), motifSites);
dev.off();

# Get sense PAC loci from Wu et al. and map
# sites to transcriptome
sites <- ReadBED("PAC_sense_Wu2011_TAIR9.bed");
PACsites <- SmartMap(sites, id = "PAC_sense", refGenome = "TAIR10");
PACsites <- FilterTxLoc(PACsites, filter = c("CDS", "3'UTR"));

# Plot distribution of PAC sites across transcript sections
pdf("figures/sectionAbundance_PACsites.pdf", height = 5, width = 7.5)
PlotSectionDistribution(PACsites);
dev.off();

# Plot enrichment of relative distances m6A/non-m6A to PAC sites in tx
pdf("figures/relDistEnrichment_PACsites.pdf", width = 15, height = 10);
PlotRelDistEnrichment(posSites, negSites, PACsites, flank = 200, binWidth = 10);
dev.off();

# Venn diagram of overlaps between m6A and PAC sites
pdf("figures/overlap_sites_PACsites.pdf", width = 9, height = 5);
PlotOverlap(FilterTxLoc(posSites, filter = c("CDS", "3'UTR")), PACsites);
dev.off();

# Plot distribution of relative distances m6A/non-m6A to PAC sites
pdf("figures/relDistAbundance_m6A_PACsites.pdf", width = 15, height = 10);
PlotRelDistDistribution(posSites, PACsites, flank = 200, binWidth = 10);
dev.off();
pdf("figures/relDistAbundance_null-m6A_PACsites.pdf", width = 15, height = 10);
PlotRelDistDistribution(negSites, PACsites, flank = 200, binWidth = 10, doBootstrap = FALSE);
dev.off();

# Get cleavage sites in all protein coding genes from Sherstnev et al. and
# map sites to transcriptome
sites <- ReadBED("CleavageSites_proteinCoding_Sherstnev2012_TAIR10.bed");
Csites <- SmartMap(sites, id = "CleavageSites", refGenome = "TAIR10");
Csites <- FilterTxLoc(Csites, filter = c("CDS", "3'UTR"));

# Plot distribution of cleavage sites across transcript sections
pdf("figures/sectionAbundance_Csites.pdf", height = 5, width = 7.5)
PlotSectionDistribution(Csites);
dev.off();

# Plot enrichment of relative distances m6A/non-m6A to cleavage sites in tx
pdf("figures/relDistEnrichment_sites_Csites.pdf", width = 15, height = 10);
PlotRelDistEnrichment(posSites, negSites, Csites, flank = 200, binWidth = 10);
dev.off();

# Venn diagram of overlaps between m6A and cleavage sites
pdf("figures/overlap_sites_Csites.pdf", width = 9, height = 5);
PlotOverlap(FilterTxLoc(posSites, filter = c("CDS", "3'UTR")), Csites);
dev.off();

# YTH motifs
motif <- c("GAATAC", "TAATAC",
           "GCCTAC", "TCCTAC",
           "GGGTAC", "TGGTAC",
           "GAATAC", "TAATAC",
           "GCCTAC", "TCCTAC",
           "GGGTAC", "TGGTAC");
motifSites <- GetMotifLoc(motif, filter = c("CDS", "3'UTR"), refGenome = "TAIR10");
WriteTxLocToBED(motifSites, "annotation/PASmotifs.bed");

PlotRelDistEnrichment(posSites, negSites, motifSites, flank = 200, binWidth = 10);
