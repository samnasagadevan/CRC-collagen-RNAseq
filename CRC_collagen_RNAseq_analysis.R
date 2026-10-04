################################################################################
#
#  Collagen-rich matrix in colorectal cancer: bulk RNA-seq analysis
#
#  Manuscript : [title / DOI]
#  Contact    : [name, affiliation, email]
#
#  This script reproduces the RNA-seq panels of the manuscript.
#
#  PART A - Patient-derived organoids (PDOs)
#           Four PDO lines (PDO23, PDO25, PDO30, PDO33) grown in
#           Matrigel (MG) or Matrigel + collagen 1 (MGCOL), either untreated
#           (control) or treated with 5-fluorouracil (5FU).
#
#             Contrast 1: MGCOL control vs MG control  -> Figure 3b, 3c, 3e, 3f, 3g, 6n
#             Contrast 2: MGCOL 5FU     vs MG 5FU      -> Figure 7e, 7f, 7g
#
#  PART B - TCGA colon adenocarcinoma (TCGA-COAD)
#           Patients grouped by COL1A1/COL1A2 expression (Low / Intermediate /
#           High); differential expression High vs Low.
#                                                      -> Figure 6a, 6n
#
#  The methods and parameters are those used to make the published figures.
#  Gene set enrichment uses random permutations, so NES and p-values can
#  differ from the paper in the last digit (the original runs had no seed).
#
#  HOW TO RUN
#    1. Put the input files in the data/ folder (see data/README.md).
#    2. From the project folder:   Rscript CRC_collagen_RNAseq_analysis.R
#    3. Figures and tables are written to results/.
#
#  PACKAGES (install once)
#    install.packages(c("BiocManager", "ggplot2", "pheatmap"))
#    BiocManager::install(c("DESeq2", "sva", "edgeR", "limma",
#                           "clusterProfiler", "org.Hs.eg.db", "AnnotationDbi",
#                           "fgsea", "EnhancedVolcano"))
#
#  CONTENTS
#    Section 0  Settings
#    Section 1  Gene sets
#    Section 2  Packages and helper functions
#    Section 3  PART A - load PDO data
#    Section 4  PART A - batch correction and quality-control PCA
#    Section 5  PART A - differential expression (DESeq2)
#    Section 6  PART A - per-contrast results, plots and GSEA
#    Section 7  PART A - signature bar plots (Figure 6n, PDOs)
#    Section 8  PART B - load TCGA-COAD data
#    Section 9  PART B - COL1A1/COL1A2 patient groups (Figure 6a)
#    Section 10 PART B - differential expression (limma-voom)
#    Section 11 PART B - signature bar plots (Figure 6n, TCGA)
#    Section 12 Session information
#
################################################################################



################################################################################
# Section 0. Settings
#
# Everything a user may need to change is in this section.
################################################################################

# ---- Input and output folders ------------------------------------------------

DATA_DIR    <- "data"
RESULTS_DIR <- "results"


# ---- Input files -------------------------------------------------------------

# PDO raw counts: first column = Ensembl gene ID, then one column per sample
PDO_COUNTS_FILE <- file.path(DATA_DIR, "Combined_4PDO_rawdata.csv")

# PDO sample sheet: one row per sample (see data/sample_metadata_template.csv)
PDO_METADATA_FILE <- file.path(DATA_DIR, "sample_metadata.csv")

# TCGA-COAD raw counts: column 'Entrez_Gene_Id', then one column per patient
TCGA_COUNTS_FILE <- file.path(DATA_DIR, "TCGA_COAD_read_counts.csv")

# Intestinal marker genes (Figure 6n): one gene symbol per row, column 'Genes'
INTESTINAL_MARKERS_FILE <- file.path(DATA_DIR, "Int_epi_genelist.csv")


# ---- PDO experimental design -------------------------------------------------

# The four groups expected in the 'Group' column of the sample sheet.
# 'treated' means treated with 5FU.
EXPECTED_GROUPS <- c("MGcontrol", "MGCOLcontrol", "MGtreated", "MGCOLtreated")

# Contrasts to test. Each compares two levels of 'Group':
#   numerator   = group expected to go up when log2 fold change is positive
#   denominator = reference group
#   label       = name used in plot titles
#   volcano     = axis limits and fold-change line for the volcano plot
CONTRASTS <- list(

  MGCOL_vs_MG_control = list(
    numerator   = "MGCOLcontrol",
    denominator = "MGcontrol",
    label       = "MGCOL vs MG",
    volcano     = list(xlim = c(-3, 3), ylim = c(0, 3), fc_cutoff = 0.2)   # Figure 3c
  ),

  MGCOL_vs_MG_5FU = list(
    numerator   = "MGCOLtreated",
    denominator = "MGtreated",
    label       = "MGCOL 5FU vs MG 5FU",
    volcano     = list(xlim = c(-3, 3), ylim = c(0, 7), fc_cutoff = 0.5)   # Figure 7f
  )
)

# Contrast whose fold changes are shown next to TCGA in Figure 6n
SIGNATURE_CONTRAST <- "MGCOL_vs_MG_control"


# ---- Statistical thresholds --------------------------------------------------

PADJ_CUTOFF <- 0.05   # adjusted p-value: DEG tables, volcano plots, KEGG GSEA
LFC_CUTOFF  <- 1      # |log2 fold change| for the "significant DEGs" tables

KEGG_PERMUTATIONS  <- 1000    # permutations for KEGG GSEA (as in the paper)
FGSEA_PERMUTATIONS <- 10000   # permutations for gene-set GSEA (as in the paper)

TREAT_LFC <- 0.5      # TCGA: limma treat() tests |log2 fold change| > 0.5

SEED <- 1234          # random seed set before every permutation-based step



################################################################################
# Section 1. Gene sets
#
# Gene symbols are kept exactly as used for the published analysis, so that
# results can be reproduced. Some are older names (e.g. CTGF, now CCN2);
# genes not found in the data are ignored by fgsea.
################################################################################

# ---- 1a. Gene sets tested by GSEA on the PDO data (fgsea) --------------------
#
# Each entry has:
#   title  = plot title
#   source = where the gene list comes from
#   genes  = gene symbols
#
# The first two are shown in Figure 3f and 3g. They are subsets of the
# MSigDB Hallmark gene sets (42 of 200 EMT genes; 51 of 200 TNFA genes).

FGSEA_GENE_SETS <- list(

  Hallmark_EMT = list(
    title  = "Hallmark_Epithelial Mesenchymal Transition",
    source = "MSigDB Hallmark EMT (42-gene subset) - Figure 3f",
    genes  = c(
      "GADD45A", "FAS", "TNFRSF11B", "SERPINE1", "COL5A2", "PMEPA1", "MMP14", "LAMA3",
      "BASP1", "APLP1", "GPC1", "CADM1", "SDC1", "CD44", "DKK1", "SGCB",
      "CXCL8", "CXCL1", "QSOX1", "ABI3BP", "RGS4", "JUN", "IGFBP2", "SDC4",
      "PLOD3", "INHBA", "LAMC1", "PLAUR", "ECM1", "ADAM12", "COL16A1", "BMP1",
      "RHOB", "ITGA2", "DCN", "MMP2", "SERPINH1", "CALD1", "AREG", "GLIPR1",
      "EDIL3", "DAB2"
    )
  ),

  Hallmark_TNFA_NFKB = list(
    title  = "Hallmark_TNFa signalling - NFkB",
    source = "MSigDB Hallmark TNFA signaling via NFKB (51-gene subset) - Figure 3g",
    genes  = c(
      "PLK2", "CDKN1A", "GADD45A", "BTG2", "NINJ1", "SERPINE1", "IER5", "LAMB3",
      "DRAM1", "PMEPA1", "PTPRE", "ABCA1", "JAG1", "SQSTM1", "TNFSF9", "PTGS2",
      "CEBPD", "CCL20", "BIRC3", "CSF1", "FOSL1", "CD44", "MCL1", "BTG3",
      "CXCL3", "PDE4B", "TNIP2", "ZC3H12A", "KDM6B", "NFKBIA", "LIF", "CXCL1",
      "JUN", "SDC4", "YRDC", "PFKFB3", "RNF19B", "FOSL2", "RELA", "CFLAR",
      "MAP3K8", "INHBA", "PLAUR", "TANK", "DUSP4", "SNN", "PLPP3", "BTG1",
      "TNIP1", "IL23A", "RHOB"
    )
  ),

  Hallmark_Inflammatory_Response = list(
    title  = "Hallmark_Inflammatory Response",
    source = "MSigDB Hallmark inflammatory response",
    genes  = c(
      "ABCA1", "ABI1", "ACVR1B", "ACVR2A", "ADM", "ADORA2B", "ADRM1", "AHR", "APLNR", "AQP9",
      "ATP2A2", "ATP2B1", "ATP2C1", "AXL", "BDKRB1", "BEST1", "BST2", "BTG2", "C3AR1", "C5AR1",
      "CALCRL", "CCL17", "CCL2", "CCL20", "CCL22", "CCL24", "CCL5", "CCL7", "CCR7", "CCRL2",
      "CD14", "CD40", "CD48", "CD55", "CD69", "CD70", "CD82", "CDKN1A", "CHST2", "CLEC5A",
      "CMKLR1", "CSF1", "CSF3", "CSF3R", "CX3CL1", "CXCL10", "CXCL11", "CXCL6", "CXCL9", "CXCR6",
      "CYBB", "DCBLD2", "EBI3", "EDN1", "EIF2AK2", "EMP3", "ADGRE1", "EREG", "F3", "FFAR2",
      "FPR1", "FZD5", "GABBR1", "GCH1", "GNA15", "GNAI3", "GP1BA", "GPC3", "GPR132", "GPR183",
      "HAS2", "HBEGF", "HIF1A", "HPN", "HRH1", "ICAM1", "ICAM4", "ICOSLG", "IFITM1", "IFNAR1",
      "IFNGR2", "IL10", "IL10RA", "IL12B", "IL15", "IL15RA", "IL18", "IL18R1", "IL18RAP", "IL1A",
      "IL1B", "IL1R1", "IL2RB", "IL4R", "IL6", "IL7R", "CXCL8", "INHBA", "IRAK2", "IRF1",
      "IRF7", "ITGA5", "ITGB3", "ITGB8", "KCNA3", "KCNJ2", "KCNMB2", "KIF1B", "KLF6", "LAMP3",
      "LCK", "LCP2", "LDLR", "LIF", "LPAR1", "LTA", "LY6E", "LYN", "MARCO", "MEFV",
      "MEP1A", "MET", "MMP14", "MSR1", "MXD1", "MYC", "NAMPT", "NDP", "NFKB1", "NFKBIA",
      "NLRP3", "NMI", "NMUR1", "NOD2", "NPFFR2", "OLR1", "OPRK1", "OSM", "OSMR", "P2RX4",
      "P2RX7", "P2RY2", "PCDH7", "PDE4B", "PDPN", "PIK3R5", "PLAUR", "PROK2", "PSEN1", "PTAFR",
      "PTGER2", "PTGER4", "PTGIR", "PTPRE", "PVR", "RAF1", "RASGRP1", "RELA", "RGS1", "RGS16",
      "RHOG", "RIPK2", "RNF144B", "ROS1", "RTP4", "SCARF1", "SCN1B", "SELE", "SELL", "SELENOS",
      "SEMA4D", "SERPINE1", "SGMS2", "SLAMF1", "SLC11A2", "SLC1A2", "SLC28A2", "SLC31A1",
      "SLC31A2", "SLC4A4", "SLC7A1", "SLC7A2", "SPHK1", "SRI", "STAB1", "TACR1", "TACR3",
      "TAPBP", "TIMP1", "TLR1", "TLR2", "TLR3", "TNFAIP6", "TNFRSF1B", "TNFRSF9", "TNFSF10",
      "TNFSF15", "TNFSF9", "TPBG", "VIP"
    )
  ),

  Mesenchymal_phenotype = list(
    title  = "Mesenchymal phenotype",
    source = "Tuan et al. 2014; Hussey et al. 2012",
    genes  = c(
      "GAS1", "CXCL12", "GLYR1", "FHL1", "FERMT2", "C1S", "FYN", "WIPF1", "SERPING1", "SERPINF1",
      "VCAM1", "TCF4", "SRPX", "DPT", "CALD1", "PTGIS", "CD163", "C1R", "FXYD6", "IGF1",
      "NAP1L3", "MRC1", "QKI", "MS4A4A", "ANK2", "LY96", "ZFPM2", "CSRP2", "EFEMP1", "RARRES2",
      "PTPRC", "PLEKHO1", "RGS2", "F13A1", "JAM2", "CHRDL1", "AP1S2", "MYLK", "DDR2", "DSE",
      "SACS", "GLIPR1", "CXCL13", "FLRT2", "AKT3", "COL6A2", "DPYSL3", "PDZRN3", "ZEB2", "CCL2",
      "MAFB", "SFRP1", "SYNE3", "MFAP4", "MAF", "UCHL1", "TUBB6", "HEG1", "KCNJ8", "AKAP12",
      "EVI2A", "COL14A1", "ECM2", "PLN", "OLFML3", "STON1", "SLIT2", "BICC1", "SOBP", "CLIC4",
      "ENPP2", "SAMSN1", "TPM2", "ASPN", "IGFBP5", "MOXD1", "AKAP2", "SLC2A3", "OLFML2B", "ANGPTL2",
      "PCOLCE", "COLEC12", "CTSK", "IL10RA", "C1orf54", "CEP170", "TNS1", "CLEC2B", "JAM3", "SEPT6",
      "GREM1", "ZCCHC24", "CRYAB", "SFRP4", "RUNX1T1", "FGL2", "MS4A6A", "PTRF", "GIMAP4", "TWIST1",
      "GFPT2", "LHFP", "CXCR4", "SPARC", "VSIG4", "GPM6B", "TRPC1", "SNAI2", "GUCY1B3", "PLXNC1",
      "FLI1", "MYH10", "CSF2RB", "TNC", "COL5A2", "GNG11", "CAV1", "SDC2", "PTGDS", "NR3C1",
      "SYNM", "FAP", "NUAK1", "WWTR1", "MPDZ", "EFEMP2", "GIMAP6", "KIAA1462", "CCL8", "COL15A1",
      "CHN1", "CRISPLD2", "PDGFC", "GEM", "ISLR", "GZMK", "SPARCL1", "BNC2", "BGN", "MEOX2",
      "ITM2A", "IFFO1"
    )
  ),

  TGFb_Harmonizome = list(
    title  = "TGFbeta signalling",
    source = "Harmonizome",
    genes  = c(
      "ACVR1", "ACVR1B", "ACVR1C", "ACVR2A", "ACVR2B", "ACVRL1", "AMHR2", "ATF2",
      "BAMBI", "BMP1", "BMP10", "BMP15", "BMP2", "BMP3", "BMP4", "BMP5",
      "BMP6", "BMP7", "BMP8A", "BMP8B", "BMPR1A", "BMPR1B", "BMPR2", "CITED1",
      "CITED2", "CREBBP", "DCP1B", "EP300", "FKBP1A", "FOSL1", "FOXH1", "GDF1",
      "GDF10", "GDF11", "GDF15", "GDF2", "GDF3", "GDF5", "GDF6", "GDF7", "GDF9",
      "GDNF", "HRAS", "INHBA", "INHBB", "INHBC", "INHBE", "JUN", "JUNB", "JUND",
      "LEFTY1", "LEFTY2", "MAP3K7", "MAP3K7CL", "MAPK1", "MAPK10", "MAPK11",
      "MAPK12", "MAPK13", "MAPK14", "MAPK3", "MAPK8", "MAPK9", "MSTN", "NODAL",
      "NRAS", "RRAS", "SKI", "SKIL", "SMAD1", "SMAD2", "SMAD3", "SMAD4", "SMAD6",
      "SMAD7", "SMAD9", "SMURF1", "SMURF2", "SNIP1", "TAB1", "TGFB2", "TGFB3",
      "TGFBR1", "TGFBR2", "TLL1", "TLL2", "ZFYVE9"
    )
  ),

  TGFb_Smad_WP5382 = list(
    title  = "TGFB / Smad signaling (WP5382)",
    source = "WikiPathways WP5382",
    genes  = c(
      "YWHAZ", "YWHAE", "YWHAB", "YWHAH", "YWHAQ", "YWHAG", "YWHAS", "CNN2", "CTGF",
      "COL1A1", "COL3A1", "FN1", "SERPINE1", "POSTN", "SPARC", "SMAD2", "SMAD3",
      "SMAD4", "WWTR1", "TGFB1", "TGFBR3", "TGFBR2", "TGFBR1"
    )
  ),

  TGFb_EMT_WP3859 = list(
    title  = "TGF-beta signaling for epithelial-mesenchymal transition (WP3859)",
    source = "WikiPathways WP3859",
    genes  = c(
      "AKT2", "RUNX2", "CDH16", "CDH1", "TNC", "MAPK3", "SMAD3", "SMAD4", "ID1", "MAPK1",
      "FN1", "SNAI1", "CDH6", "VIM", "SMAD2", "TGFB1", "SNAI2", "AKT1", "AKT3", "CDH2"
    )
  ),

  FOROUTAN_TGFB_EMT_DN = list(
    title  = "FOROUTAN_TGFB_EMT_DN",
    source = "MSigDB C2",
    genes  = c(
      "TSPAN1", "MPZL2", "SPRY1", "CITED2", "VAV3", "CEBPD", "AGR2", "SLC27A2", "TBC1D8", "PLAAT3",
      "CP", "ADORA2B", "CXADR", "CYB5A", "CYP1B1", "TMEM30B", "DEFB1", "SERPINB1", "ELF3", "EMP1",
      "EPAS1", "ERBB3", "EREG", "ALDH1A3", "ALDH3A2", "PEG10", "GSE1", "SYNE2", "NUP210", "RAB38",
      "RAB26", "TMT1A", "EHF", "SMPDL3B", "GLDC", "SLCO4A1", "ANK3", "CFH", "FOXA2", "HPGD",
      "BIRC3", "AQP3", "IMPA2", "INHBB", "JAG2", "AREG", "KRT15", "KRT19", "LAMA5", "LCN2",
      "ABLIM1", "LY6E", "EPCAM", "MBP", "KITLG", "MITF", "MMP7", "MUC1", "CEACAM6", "HOOK1",
      "GULP1", "LSR", "PDK4", "ATP8B1", "PKP2", "PLS1", "RBM47", "EPB41L4B", "MANSC1", "ESRP1",
      "PPL", "PPP1R9A", "SYBU", "GPRC5C", "MYO5C", "RAB25", "MTUS1", "SQOR", "PLAAT4", "S100P",
      "CFB", "SCNN1A", "DEPTOR", "SLPI", "SORD", "SOX2", "DST", "SULT1A1", "TFAP2A", "NR2F2",
      "TNFAIP2", "TPD52L1", "C1orf116", "ALDH5A1", "FA2H", "C1orf115", "GRTP1", "ERMP1", "CAVIN2",
      "MAP7", "SLC16A7", "CD9", "ARHGAP29", "TJP2", "GDF15", "HS3ST1", "FGFBP1", "CDH1"
    )
  ),

  PID_PI3KCI_AKT = list(
    title  = "PID_PI3KCI_AKT_PATHWAY",
    source = "MSigDB C2 (PID)",
    genes  = c(
      "CHUK", "PDPK1", "FOXO3", "TBC1D4", "RAF1", "HSP90AA1", "SRC", "SLC2A4", "PRKACA",
      "YWHAQ", "AKT1", "AKT2", "YWHAB", "SFN", "CDKN1A", "MTOR", "CDKN1B", "GSK3A", "GSK3B",
      "KPNA1", "CASP9", "YWHAG", "YWHAE", "YWHAZ", "PRKDC", "FOXO4", "YWHAH", "BCL2L1",
      "FOXO1", "RICTOR", "BAD", "MAP3K5", "MAPKAP1", "MLST8", "AKT3"
    )
  )
)


# ---- 1b. Signatures shown as log2 fold-change bar plots (Figure 6n) ---------

# EMT signature: top 50 markers of cluster 17 (EMT) of the single-cell
# colorectal cancer reference dataset.
EMT_SIGNATURE <- c(
  "MZB1", "IGHG3", "IGKC", "IGHG2", "IGLC2", "IGHG1", "IGJ", "LGALS1", "IGLC7", "IGLC3",
  "IGHM", "CD79A", "LUM", "NEAT1", "IGHG4", "IGHA2", "IGLV6-57", "PRDM1", "ANKRD28", "SRGN",
  "JUND", "IGLL1", "HSPA6", "FOSB", "RGS1", "LY96", "COL1A2", "POSTN", "SELM", "CD27",
  "EGR3", "VIM", "HLA-DPB1", "MIR155HG", "TNFRSF17", "PCED1B-AS1", "CYSLTR1", "FBLN1",
  "C4orf26", "VCAM1", "KIAA0125", "NR4A1", "C1S", "DCN", "HSPA1B", "CPS1", "LYZ", "DUSP5",
  "SYNGR1", "AL928768.3"
)

# Intestinal markers signature: read from file (column 'Genes')
INTESTINAL_MARKERS <- read.csv(INTESTINAL_MARKERS_FILE)$Genes



################################################################################
# Section 2. Packages and helper functions
################################################################################

# ---- 2a. Load packages (stop with a clear message if any are missing) -------

required_packages <- c(
  "DESeq2",           # differential expression, PDOs
  "sva",              # ComBat-seq batch correction
  "edgeR",            # CPM normalisation, DGEList
  "limma",            # voom / treat, TCGA
  "clusterProfiler",  # KEGG GSEA, gene ID conversion
  "org.Hs.eg.db",     # human gene annotation
  "AnnotationDbi",    # mapIds()
  "fgsea",            # gene-set GSEA
  "EnhancedVolcano",  # volcano plots
  "ggplot2",          # plotting
  "pheatmap"          # heatmaps
)

is_installed <- vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)

if (any(!is_installed)) {
  stop("Please install these packages first: ",
       paste(required_packages[!is_installed], collapse = ", "))
}

suppressPackageStartupMessages({
  for (pkg in required_packages) {
    library(pkg, character.only = TRUE)
  }
})


# ---- 2b. Output helpers ------------------------------------------------------

# Build a path inside results/ and create the folder if it does not exist yet.
# Example: result_path("A_PDO", "table.csv") -> "results/A_PDO/table.csv"
result_path <- function(...) {
  path <- file.path(RESULTS_DIR, ...)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  path
}

# Save a ggplot (or any printable plot object) as a 300 dpi PNG.
save_png <- function(plot, file, width = 8, height = 6) {
  png(file, width = width, height = height, units = "in", res = 300)
  on.exit(dev.off())
  print(plot)
  invisible(file)
}

# Print a section heading to the console so progress is easy to follow.
announce <- function(text) {
  message("\n", strrep("=", 70), "\n", text, "\n", strrep("=", 70))
}


# ---- 2c. Gene identifier helpers ---------------------------------------------

# Ensembl gene IDs -> gene symbols (NA where no symbol exists)
ensembl_to_symbol <- function(ensembl_ids) {
  symbols <- mapIds(org.Hs.eg.db,
                    keys      = ensembl_ids,
                    keytype   = "ENSEMBL",
                    column    = "SYMBOL",
                    multiVals = "first")
  unname(symbols)
}

# Entrez gene IDs -> gene symbols (NA where no symbol exists)
entrez_to_symbol <- function(entrez_ids) {
  symbols <- mapIds(org.Hs.eg.db,
                    keys      = as.character(entrez_ids),
                    keytype   = "ENTREZID",
                    column    = "SYMBOL",
                    multiVals = "first")
  unname(symbols)
}


# ---- 2d. Ranked gene lists for GSEA (PDO data) -------------------------------

# For KEGG GSEA: all genes ranked by the DESeq2 Wald statistic, named by
# Entrez ID (KEGG uses Entrez IDs).
rank_by_wald_stat <- function(deseq_result) {

  id_map <- bitr(rownames(deseq_result),
                 fromType = "ENSEMBL",
                 toType   = "ENTREZID",
                 OrgDb    = org.Hs.eg.db)

  stats <- data.frame(ENSEMBL = rownames(deseq_result),
                      stat    = deseq_result$stat)

  stats <- merge(stats, id_map, by = "ENSEMBL")
  stats <- stats[!is.na(stats$stat), ]

  ranked <- setNames(stats$stat, stats$ENTREZID)
  sort(ranked, decreasing = TRUE)
}

# For fgsea: all genes ranked by log2 fold change, named by gene symbol.
# When several Ensembl IDs share one symbol, the one with the largest
# absolute log2 fold change is kept.
rank_by_log2fc <- function(deseq_result) {

  id_map <- bitr(rownames(deseq_result),
                 fromType = "ENSEMBL",
                 toType   = "SYMBOL",
                 OrgDb    = org.Hs.eg.db)

  lfc <- data.frame(ENSEMBL = rownames(deseq_result),
                    log2FC  = deseq_result$log2FoldChange)

  lfc <- merge(lfc, id_map, by = "ENSEMBL")
  lfc <- lfc[order(abs(lfc$log2FC), decreasing = TRUE), ]
  lfc <- lfc[!duplicated(lfc$SYMBOL), ]

  ranked <- setNames(lfc$log2FC, lfc$SYMBOL)
  sort(ranked, decreasing = TRUE)   # sort() also removes NA values
}


# ---- 2e. Plotting helpers ----------------------------------------------------

# PCA of log2(counts + 1) for a chosen set of samples.
#   counts   : count matrix (genes x samples)
#   samples  : sample names to include
#   pc_x/y   : principal components to plot
#   colour_by: metadata column used for point colour
plot_pca <- function(counts, metadata, samples, pc_x, pc_y, colour_by, title,
                     add_ellipse = FALSE) {

  log_counts <- log2(counts[, samples] + 1)
  pca        <- prcomp(t(log_counts))

  percent_variance <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)

  plot_data <- data.frame(
    PCx = pca$x[, pc_x],
    PCy = pca$x[, pc_y],
    metadata[samples, ]
  )

  p <- ggplot(plot_data, aes(x = PCx, y = PCy, colour = .data[[colour_by]])) +
    geom_point(size = 3) +
    labs(
      title = title,
      x     = sprintf("PC%d (%s%% variance)", pc_x, percent_variance[pc_x]),
      y     = sprintf("PC%d (%s%% variance)", pc_y, percent_variance[pc_y])
    ) +
    theme_bw()

  if (add_ellipse) {
    p <- p + stat_ellipse()
  }

  p
}

# Dot plot of significant KEGG pathways: NES on the x-axis, coloured by
# adjusted p-value, dot size = -log10(adjusted p-value). (Figure 3e, 7g)
plot_kegg_dotplot <- function(kegg_table, title, max_pathways = 30) {

  top <- kegg_table[order(kegg_table$p.adjust), ]
  top <- head(top, max_pathways)

  ggplot(top, aes(x      = NES,
                  y      = reorder(Description, NES),
                  colour = p.adjust,
                  size   = -log10(p.adjust))) +
    geom_point() +
    scale_colour_viridis_c(direction = -1) +
    labs(title = title,
         x     = "Normalized Enrichment Score",
         y     = NULL,
         size  = "-log10(p.adjust)") +
    theme_bw()
}

# Horizontal bar plot of log2 fold changes for the genes of one signature.
# (Figure 6n)
#   log2fc  : named numeric vector (names = gene symbols)
#   genes   : signature genes to show
#   colours : c(colour for down-regulated, colour for up-regulated)
plot_signature_bars <- function(log2fc, genes, title, colours) {

  bar_data <- data.frame(Gene = names(log2fc), log2FC = unname(log2fc))
  bar_data <- bar_data[bar_data$Gene %in% genes, ]
  bar_data <- bar_data[!is.na(bar_data$log2FC), ]
  bar_data <- bar_data[!duplicated(bar_data$Gene), ]

  message("  ", title, ": ", nrow(bar_data), " of ", length(unique(genes)),
          " signature genes found")

  ggplot(bar_data, aes(x    = reorder(Gene, log2FC),
                       y    = log2FC,
                       fill = log2FC > 0)) +
    geom_col(show.legend = FALSE) +
    coord_flip() +
    scale_fill_manual(values = c(`FALSE` = colours[1], `TRUE` = colours[2])) +
    labs(title = title, x = "Gene", y = "Log2 Fold Change") +
    theme_minimal()
}



################################################################################
#
#                    PART A - PATIENT-DERIVED ORGANOIDS
#
################################################################################



################################################################################
# Section 3. Load PDO counts and sample sheet
################################################################################

announce("PART A - Patient-derived organoids: loading data")

# ---- 3a. Count matrix --------------------------------------------------------

pdo_counts_table <- read.csv(PDO_COUNTS_FILE, check.names = FALSE)

pdo_counts <- as.matrix(pdo_counts_table[, -1])
rownames(pdo_counts) <- pdo_counts_table[[1]]   # Ensembl gene IDs
mode(pdo_counts) <- "numeric"

message("Count matrix: ", nrow(pdo_counts), " genes x ", ncol(pdo_counts), " samples")


# ---- 3b. Sample sheet --------------------------------------------------------

pdo_metadata <- read.csv(PDO_METADATA_FILE, check.names = FALSE)

# Required columns
required_columns <- c("Sample", "Group", "Batch")
missing_columns  <- setdiff(required_columns, colnames(pdo_metadata))
if (length(missing_columns) > 0) {
  stop("Sample sheet is missing columns: ", paste(missing_columns, collapse = ", "))
}

# Every sample in the count matrix must be in the sample sheet
samples_not_in_sheet <- setdiff(colnames(pdo_counts), pdo_metadata$Sample)
if (length(samples_not_in_sheet) > 0) {
  stop("These count-matrix samples are not in the sample sheet: ",
       paste(samples_not_in_sheet, collapse = ", "))
}

# Put the sample sheet in the same order as the count-matrix columns
pdo_metadata <- pdo_metadata[match(colnames(pdo_counts), pdo_metadata$Sample), ]
rownames(pdo_metadata) <- pdo_metadata$Sample

# All four groups (control and 5FU, MG and MGCOL) must be present
missing_groups <- setdiff(EXPECTED_GROUPS, pdo_metadata$Group)
if (length(missing_groups) > 0) {
  stop("These groups are missing from the sample sheet: ",
       paste(missing_groups, collapse = ", "))
}

pdo_metadata$Group <- factor(pdo_metadata$Group, levels = EXPECTED_GROUPS)
pdo_metadata$Batch <- factor(pdo_metadata$Batch)   # Batch = PDO line

message("Samples per group and PDO line:")
print(table(Group = pdo_metadata$Group, PDO_line = pdo_metadata$Batch))



################################################################################
# Section 4. Batch correction and quality-control PCA
#
# The four PDO lines are separate patients, so expression differs between
# lines. ComBat-seq removes this line (batch) effect from the counts while
# keeping differences between groups.
#
# The batch-corrected counts are used ONLY for PCA plots. Differential
# expression (Section 5) uses the raw counts with PDO line in the model.
################################################################################

announce("PART A - Batch correction (ComBat-seq)")

pdo_counts_corrected <- ComBat_seq(counts = pdo_counts,
                                   batch  = pdo_metadata$Batch,
                                   group  = pdo_metadata$Group)

write.csv(pdo_counts_corrected,
          result_path("A_PDO", "batch_corrected_counts.csv"))

# Batch-corrected, TMM-normalised log2 counts per million (for reference)
corrected_dge    <- calcNormFactors(DGEList(counts = pdo_counts_corrected))
corrected_log2cpm <- cpm(corrected_dge, log = TRUE, prior.count = 1)

write.csv(corrected_log2cpm,
          result_path("A_PDO", "batch_corrected_log2cpm.csv"))


# ---- Quality control: PCA of all samples before and after correction --------

all_samples <- colnames(pdo_counts)

pca_before <- plot_pca(pdo_counts, pdo_metadata, all_samples,
                       pc_x = 1, pc_y = 2, colour_by = "Batch",
                       title = "All samples - before batch correction")

pca_after <- plot_pca(pdo_counts_corrected, pdo_metadata, all_samples,
                      pc_x = 1, pc_y = 2, colour_by = "Batch",
                      title = "All samples - after batch correction")

save_png(pca_before, result_path("A_PDO", "QC_PCA_before_batch_correction.png"))
save_png(pca_after,  result_path("A_PDO", "QC_PCA_after_batch_correction.png"))



################################################################################
# Section 5. Differential expression with DESeq2
#
# Model: ~ Batch + Group
#   Batch (PDO line) is included so that each comparison is made within
#   PDO lines, i.e. a paired analysis across the four patients.
################################################################################

announce("PART A - Differential expression (DESeq2, ~ Batch + Group)")

dds <- DESeqDataSetFromMatrix(countData = pdo_counts,
                              colData   = pdo_metadata,
                              design    = ~ Batch + Group)

dds <- DESeq(dds)

saveRDS(dds, result_path("A_PDO", "DESeq2_object.rds"))



################################################################################
# Section 6. Per-contrast results, plots and gene set enrichment
#
# For each contrast in CONTRASTS:
#   6a  PCA of the two groups (batch-corrected)       Figure 3b / 7e
#   6b  DEG tables
#   6c  Volcano plot                                  Figure 3c / 7f
#   6d  KEGG GSEA (ranked by Wald statistic)          Figure 3e / 7g
#   6e  Gene-set GSEA with fgsea (ranked by log2FC)   Figure 3f, 3g
################################################################################

# Log2 fold changes kept for the Figure 6n bar plots (Section 7)
pdo_signature_log2fc <- NULL

for (contrast_name in names(CONTRASTS)) {

  contrast <- CONTRASTS[[contrast_name]]
  out_dir  <- file.path("A_PDO", contrast_name)

  announce(paste("PART A - Contrast:", contrast$label))


  # ---- 6a. PCA of the two groups in this contrast ---------------------------

  contrast_samples <- rownames(pdo_metadata)[
    pdo_metadata$Group %in% c(contrast$numerator, contrast$denominator)
  ]

  pca_contrast <- plot_pca(pdo_counts_corrected, pdo_metadata, contrast_samples,
                           pc_x = 2, pc_y = 3, colour_by = "Group",
                           title = contrast$label, add_ellipse = TRUE)

  save_png(pca_contrast, result_path(out_dir, "PCA_PC2_PC3.png"))


  # ---- 6b. Differential expression tables -----------------------------------

  deseq_result <- results(dds, contrast = c("Group",
                                            contrast$numerator,
                                            contrast$denominator))

  # Add gene symbols; genes without a symbol keep their Ensembl ID
  gene_symbols <- ensembl_to_symbol(rownames(deseq_result))
  gene_symbols <- ifelse(is.na(gene_symbols), rownames(deseq_result), gene_symbols)

  deg_table <- data.frame(ensembl_gene_id = rownames(deseq_result),
                          symbol          = gene_symbols,
                          as.data.frame(deseq_result),
                          row.names       = NULL)

  deg_table <- deg_table[order(deg_table$padj), ]

  significant_degs <- subset(deg_table,
                             !is.na(padj) &
                               padj < PADJ_CUTOFF &
                               abs(log2FoldChange) > LFC_CUTOFF)

  write.csv(deg_table,        result_path(out_dir, "DEG_all_genes.csv"),   row.names = FALSE)
  write.csv(significant_degs, result_path(out_dir, "DEG_significant.csv"), row.names = FALSE)

  message("Significant DEGs (padj < ", PADJ_CUTOFF, ", |log2FC| > ", LFC_CUTOFF, "): ",
          nrow(significant_degs),
          " (", sum(significant_degs$log2FoldChange > 0), " up, ",
          sum(significant_degs$log2FoldChange < 0), " down)")


  # ---- 6c. Volcano plot -----------------------------------------------------

  volcano <- EnhancedVolcano(deg_table,
                             lab            = deg_table$symbol,
                             x              = "log2FoldChange",
                             y              = "padj",
                             title          = contrast$label,
                             subtitle       = NULL,
                             xlim           = contrast$volcano$xlim,
                             ylim           = contrast$volcano$ylim,
                             pCutoff        = PADJ_CUTOFF,
                             FCcutoff       = contrast$volcano$fc_cutoff,
                             pointSize      = 1.0,
                             labSize        = 3.5,
                             col            = c("grey30", "forestgreen", "royalblue", "red2"),
                             colAlpha       = 0.9,
                             legendPosition = "right",
                             legendLabSize  = 12,
                             legendIconSize = 4.0)

  save_png(volcano, result_path(out_dir, "Volcano.png"), width = 10, height = 8)


  # ---- 6d. KEGG GSEA ----------------------------------------------------------

  kegg_ranked_genes <- rank_by_wald_stat(deseq_result)

  set.seed(SEED)
  kegg_gsea <- gseKEGG(geneList     = kegg_ranked_genes,
                       organism     = "hsa",
                       nPerm        = KEGG_PERMUTATIONS,
                       minGSSize    = 10,
                       pvalueCutoff = PADJ_CUTOFF,
                       seed         = TRUE,
                       verbose      = FALSE)

  kegg_table <- as.data.frame(kegg_gsea)

  if (nrow(kegg_table) == 0) {

    message("No significant KEGG pathways")

  } else {

    message("Significant KEGG pathways: ", nrow(kegg_table))

    # Convert the core (leading-edge) Entrez IDs to gene symbols
    kegg_table$core_enrichment_symbol <- vapply(
      strsplit(kegg_table$core_enrichment, "/"),
      function(entrez_ids) paste(entrez_to_symbol(entrez_ids), collapse = "/"),
      character(1)
    )

    kegg_columns <- c("ID", "Description", "setSize", "NES", "pvalue", "p.adjust",
                      "core_enrichment", "core_enrichment_symbol")

    write.csv(kegg_table[, kegg_columns],
              result_path(out_dir, "GSEA_KEGG.csv"),
              row.names = FALSE)

    kegg_plot <- plot_kegg_dotplot(kegg_table,
                                   title = paste0("GSEA top pathways\n", contrast$label))

    save_png(kegg_plot, result_path(out_dir, "GSEA_KEGG_dotplot.png"), width = 9, height = 9)
  }


  # ---- 6e. Gene-set GSEA with fgsea -------------------------------------------
  #
  # Each gene set is tested on its own against the log2FC-ranked gene list,
  # as in the published analysis.

  log2fc_ranked_genes <- rank_by_log2fc(deseq_result)

  fgsea_rows <- list()

  for (set_name in names(FGSEA_GENE_SETS)) {

    gene_set <- FGSEA_GENE_SETS[[set_name]]

    set.seed(SEED)
    fgsea_result <- fgseaSimple(pathways = setNames(list(gene_set$genes), set_name),
                                stats    = log2fc_ranked_genes,
                                nperm    = FGSEA_PERMUTATIONS)

    if (nrow(fgsea_result) == 0) {
      message("  ", set_name, ": too few genes found - skipped")
      next
    }

    message(sprintf("  %-32s NES = %5.2f   adj. p = %s",
                    set_name, fgsea_result$NES, signif(fgsea_result$padj, 2)))

    # Running enrichment score plot
    enrichment_plot <- plotEnrichment(gene_set$genes, log2fc_ranked_genes) +
      ggtitle(sprintf("%s\n(NES=%.1f; adj. p=%s)",
                      gene_set$title,
                      fgsea_result$NES,
                      signif(fgsea_result$padj, 2)))

    save_png(enrichment_plot, result_path(out_dir, paste0("fgsea_", set_name, ".png")))

    # One summary row per gene set
    fgsea_rows[[set_name]] <- data.frame(
      gene_set     = set_name,
      title        = gene_set$title,
      source       = gene_set$source,
      genes_tested = fgsea_result$size,
      NES          = fgsea_result$NES,
      pval         = fgsea_result$pval,
      padj         = fgsea_result$padj,
      leading_edge = paste(fgsea_result$leadingEdge[[1]], collapse = "/")
    )
  }

  write.csv(do.call(rbind, fgsea_rows),
            result_path(out_dir, "fgsea_gene_sets.csv"),
            row.names = FALSE)


  # Keep fold changes for the Figure 6n bar plots
  if (contrast_name == SIGNATURE_CONTRAST) {
    pdo_signature_log2fc <- log2fc_ranked_genes
  }
}



################################################################################
# Section 7. Signature bar plots - PDOs (Figure 6n, right-hand side)
################################################################################

announce("PART A - Signature bar plots (Figure 6n, PDOs)")

pdo_emt_bars <- plot_signature_bars(pdo_signature_log2fc,
                                    genes   = EMT_SIGNATURE,
                                    title   = "EMT signature - MGCOL vs MG CRC PDOs",
                                    colours = c("lightblue", "darkblue"))

pdo_intestinal_bars <- plot_signature_bars(pdo_signature_log2fc,
                                           genes   = INTESTINAL_MARKERS,
                                           title   = "Intestinal markers - MGCOL vs MG CRC PDOs",
                                           colours = c("lightblue", "darkblue"))

save_png(pdo_emt_bars,
         result_path("Figure6n", "PDO_EMT_signature.png"), width = 7, height = 9)
save_png(pdo_intestinal_bars,
         result_path("Figure6n", "PDO_intestinal_markers.png"), width = 7, height = 9)



################################################################################
#
#                    PART B - TCGA COLON ADENOCARCINOMA
#
################################################################################



################################################################################
# Section 8. Load TCGA-COAD counts and convert Entrez IDs to gene symbols
################################################################################

announce("PART B - TCGA-COAD: loading data")

tcga_table <- read.csv(TCGA_COUNTS_FILE, check.names = FALSE)

if (!"Entrez_Gene_Id" %in% colnames(tcga_table)) {
  stop("The TCGA counts file needs a column called 'Entrez_Gene_Id'")
}

# Every column except the gene identifier columns is a patient sample
patient_columns <- setdiff(colnames(tcga_table), c("Entrez_Gene_Id", "Hugo_Symbol"))

# Map each row's own Entrez ID to its gene symbol
tcga_symbols <- entrez_to_symbol(tcga_table$Entrez_Gene_Id)

# Keep genes with a symbol; if a symbol occurs twice, keep the first row
keep_row <- !is.na(tcga_symbols) & !duplicated(tcga_symbols)

tcga_counts <- as.matrix(tcga_table[keep_row, patient_columns])
rownames(tcga_counts) <- tcga_symbols[keep_row]
mode(tcga_counts) <- "numeric"

message("Rows in file: ", nrow(tcga_table),
        " | genes kept (unique symbol): ", nrow(tcga_counts),
        " | patients: ", ncol(tcga_counts))


# Library-size normalised log2 counts per million
tcga_log2cpm <- cpm(tcga_counts, log = TRUE)

write.csv(tcga_log2cpm, result_path("B_TCGA", "TCGA_COAD_log2cpm.csv"))



################################################################################
# Section 9. Group patients by COL1A1 / COL1A2 expression (Figure 6a)
#
#   1. z-score COL1A1 and COL1A2 across patients
#   2. cap z-scores at -2 and +2
#   3. hierarchical clustering of patients (pheatmap defaults: Euclidean
#      distance, complete linkage) and cut into 3 clusters
#   4. name the clusters Low / Intermediate / High by their mean COL1A1 level
#
# Note: pheatmap(scale = "row") re-scales the capped z-scores. This step is
# kept because it was used to define the published groups.
################################################################################

announce("PART B - COL1A1/COL1A2 patient groups")

collagen_genes <- c("COL1A1", "COL1A2")

missing_collagen <- setdiff(collagen_genes, rownames(tcga_log2cpm))
if (length(missing_collagen) > 0) {
  stop("Not found in TCGA data: ", paste(missing_collagen, collapse = ", "))
}

# Steps 1 and 2: z-scores, capped at +/- 2
collagen_z <- t(scale(t(tcga_log2cpm[collagen_genes, ])))
collagen_z <- pmax(pmin(collagen_z, 2), -2)

# Step 3: cluster patients into 3 groups
collagen_heatmap  <- pheatmap(collagen_z, scale = "row", cutree_cols = 3, silent = TRUE)
patient_cluster   <- cutree(collagen_heatmap$tree_col, k = 3)

# Step 4: name clusters by mean COL1A1 (lowest = Low, highest = High)
cluster_mean_col1a1 <- tapply(collagen_z["COL1A1", names(patient_cluster)],
                              patient_cluster, mean)

clusters_low_to_high <- names(sort(cluster_mean_col1a1))
cluster_names <- setNames(c("Low", "Intermediate", "High"), clusters_low_to_high)

col1_group <- factor(cluster_names[as.character(patient_cluster)],
                     levels = c("Low", "Intermediate", "High"))
names(col1_group) <- names(patient_cluster)

# Same order as the count-matrix columns
col1_group <- col1_group[colnames(tcga_counts)]

message("Patients per COL1 group:")
print(table(col1_group))

write.csv(data.frame(patient = names(col1_group), COL1_group = col1_group),
          result_path("B_TCGA", "COL1_patient_groups.csv"),
          row.names = FALSE)


# ---- Heatmap with group annotation (Figure 6a) --------------------------------

group_labels <- c(Low = "Low Col", Intermediate = "Intermediate", High = "High Col")

heatmap_annotation <- data.frame(
  cluster   = factor(group_labels[as.character(col1_group)]),
  row.names = names(col1_group)
)

heatmap_colours <- list(
  cluster = c(Intermediate = "white", `Low Col` = "lightblue", `High Col` = "plum")
)

pheatmap(collagen_z,
         scale             = "row",
         cutree_cols       = 3,
         show_rownames     = TRUE,
         show_colnames     = TRUE,
         fontsize_row      = 10,
         annotation_col    = heatmap_annotation,
         annotation_colors = heatmap_colours,
         filename          = result_path("Figure6a", "TCGA_COL1_heatmap.png"),
         width             = 12,
         height            = 4)



################################################################################
# Section 10. Differential expression: High vs Low COL1A1 (limma-voom)
#
#   - voom models the mean-variance relationship of the counts
#   - the linear model includes all three groups; the contrast is High - Low
#   - treat() tests for |log2 fold change| > TREAT_LFC rather than > 0
################################################################################

announce("PART B - Differential expression (limma-voom, High vs Low)")

design <- model.matrix(~ 0 + col1_group)
colnames(design) <- levels(col1_group)

contrast_matrix <- makeContrasts(High_vs_Low = High - Low, levels = design)

voom_data <- voom(DGEList(counts = tcga_counts), design)

linear_fit   <- lmFit(voom_data, design)
contrast_fit <- contrasts.fit(linear_fit, contrasts = contrast_matrix)
bayes_fit    <- eBayes(contrast_fit)
treat_fit    <- treat(bayes_fit, lfc = TREAT_LFC)

tcga_deg_table <- topTreat(treat_fit, coef = 1, n = Inf)

write.csv(data.frame(Gene = rownames(tcga_deg_table), tcga_deg_table, row.names = NULL),
          result_path("B_TCGA", "DEG_High_vs_Low_COL1A1.csv"),
          row.names = FALSE)

message("Genes with |log2FC| significantly > ", TREAT_LFC, ":")
print(summary(decideTests(treat_fit)))



################################################################################
# Section 11. Signature bar plots - TCGA (Figure 6n, left-hand side)
################################################################################

announce("PART B - Signature bar plots (Figure 6n, TCGA)")

tcga_log2fc <- setNames(tcga_deg_table$logFC, rownames(tcga_deg_table))

tcga_emt_bars <- plot_signature_bars(tcga_log2fc,
                                     genes   = EMT_SIGNATURE,
                                     title   = "EMT signature - High vs low COL1A1 TCGA-COAD",
                                     colours = c("darkblue", "firebrick"))

tcga_intestinal_bars <- plot_signature_bars(tcga_log2fc,
                                            genes   = INTESTINAL_MARKERS,
                                            title   = "Intestinal markers - High vs low COL1A1 TCGA-COAD",
                                            colours = c("darkblue", "firebrick"))

save_png(tcga_emt_bars,
         result_path("Figure6n", "TCGA_EMT_signature.png"), width = 7, height = 9)
save_png(tcga_intestinal_bars,
         result_path("Figure6n", "TCGA_intestinal_markers.png"), width = 7, height = 9)



################################################################################
# Section 12. Session information (package versions used for this run)
################################################################################

writeLines(capture.output(sessionInfo()), result_path("sessionInfo.txt"))

announce(paste("Finished. Results are in", normalizePath(RESULTS_DIR)))
