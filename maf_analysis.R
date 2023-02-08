
# R 4 
# DO


# # Generating MAF files
# 
# -   For VCF files or simple tabular files, easy option is to use [vcf2maf](https://github.com/mskcc/vcf2maf) utility which will annotate VCFs, prioritize transcripts, and generates an MAF. Recent updates to gatk has also enabled 
# 
# -   If you're using [ANNOVAR](http://annovar.openbioinformatics.org/en/latest/) for variant annotations, maftools has a handy function `annovarToMaf` for converting tabular annovar outputs to MAF.
# 
# # MAF field requirements
# 
# MAF files contain many fields ranging from chromosome names to cosmic annotations. However most of the analysis in maftools uses following fields.
# 
# -   Mandatory fields: **Hugo_Symbol, Chromosome, Start_Position, End_Position, Reference_Allele, Tumor_Seq_Allele2, Variant_Classification, Variant_Type and Tumor_Sample_Barcode**.
# 
# # Installation


if (!require("BiocManager"))
  install.packages("BiocManager")
BiocManager::install("maftools")

# Reading and summarizing maf files

## Required input files

# -   an MAF file - can be gz compressed. Required.
# -   an optional but recommended clinical data associated with each sample/Tumor_Sample_Barcode in MAF.
# -   an optional copy number data if available. Can be GISTIC output or a custom table containing sample names, gene names and 
# copy-number status (`Amp` or `Del`).
# 
# ## Reading MAF files.
# 
# `read.maf` function reads MAF files, summarizes it in various ways and stores it as an MAF object. Even though MAF file is alone enough, 
# it is recommended to provide annotations associated with samples in MAF. One can also integrate copy number data if available.


library(maftools)

#subset maf (MODEL9)
tsb_high = c('TCGA-IB-7891-01A-11D-2201-08', 'TCGA-HZ-7925-01A-11D-2154-08',
             'TCGA-HZ-8002-01A-11D-2201-08', 'TCGA-HZ-7924-01A-11D-2154-08',
             'TCGA-HZ-7918-01A-11D-2154-08', 'TCGA-HZ-8003-01A-21D-2201-08')

tsb_low = c('TCGA-IB-7893-01A-11D-2201-08', 'TCGA-IB-7890-01A-12D-2201-08',
            'TCGA-IB-7645-01A-22D-2201-08', 'TCGA-IB-7887-01A-11D-2154-08',
            'TCGA-IB-7649-01A-11D-2154-08', 'TCGA-F2-7276-01A-11D-2154-08',
            'TCGA-IB-7886-01A-11D-2154-08', 'TCGA-HZ-8001-01A-11D-2201-08',
            'TCGA-HZ-8005-01A-11D-2201-08', 'TCGA-F2-7273-01A-11D-2154-08',
            'TCGA-IB-7897-01A-21D-2201-08')



paad = read.maf(maf = '../data/paad_all.maf', clinicalData = '../data/group_model9.csv')
maf_data_view <- paad@data
maf_gene.summary <- paad@gene.summary

# GOOD PRONO
paad_high <- subsetMaf(
  maf = paad,
  tsb = tsb_high)

# BAD PRONO
paad_low <- subsetMaf(
  maf = paad,
  tsb = tsb_low)


# Collect MAF data from TCGA 
#https://portal.gdc.cancer.gov/

# library(TCGAbiolinks)
# GDCprojects <- getGDCprojects()
# summary <- TCGAbiolinks:::getProjectSummary("TCGA-CHOL")
# 
# query.maf.hg19 <- GDCquery(project = "TCGA-CHOL", 
#                            data.category = "Simple nucleotide variation", 
#                            data.type = "Simple somatic mutation",
#                            access = "open", 
#                            file.type = "bcgsc.ca_CHOL.IlluminaHiSeq_DNASeq.1.somatic.maf",
#                            legacy = TRUE)
# 
# GDCdownload(query.maf.hg19)
# maf_TCGA_data <- GDCprepare(query.maf.hg19)
# 
# maf_TCGA <- read.maf(maf_TCGA_data,isTCGA = TRUE)

## MAF object
# 
# Summarized MAF file is stored as an MAF object. MAF object contains main maf file, summarized data and any associated sample annotations.
# 
# There are accessor methods to access the useful slots from MAF object.


#Typing PAAD shows basic summary of MAF file.
paad_low


#extract sample summary from maf.
getSampleSummary(x=paad)
#extract gene summary.
getGeneSummary(paad)
#extract annotation from maf OBJECT
getClinicalData(paad)
#extract available fields (data columns name) from MAF object
getFields(paad)
#Writes maf summary to an output file with basename laml.
write.mafSummary(maf = paad, basename = 'data_summaries/paad')



# Visualization

## Plotting MAF summary.

# We can use `plotmafSummary` to plot the summary of the maf file, which displays number of variants in each sample as a 
# stacked barplot and variant types as a boxplot summarized by Variant_Classification.


plotmafSummary(maf = paad_high, showBarcodes = T, rmOutlier = TRUE, addStat = 'median', titvRaw = T)
plotmafSummary(maf = paad_low, showBarcodes = T, rmOutlier = TRUE, addStat = 'median', titvRaw = T)

vc_cols = RColorBrewer::brewer.pal(n = 8, name = 'Set2')
names(vc_cols) = c(
  'Frame_Shift_Del',
  'Missense_Mutation',
  'Nonsense_Mutation',  
  'Multi_Hit',
  'Frame_Shift_Ins',
  'In_Frame_Ins',
  'Splice_Site',
  'In_Frame_Del'
)
print(vc_cols)

mafbarplot(maf = paad_high, color = vc_cols, borderCol = "#000000")
mafbarplot(maf = paad_low, color = vc_cols, borderCol = "#000000")


# Use `mafbarplot` for a minimal barplot of mutated genes.
# 
# ## Oncoplots
# 
# ### Drawing oncoplots
# 
# Better representation of maf file can be shown as oncoplots, also known as waterfall plots.



#oncoplot for top ten mutated genes.
oncoplot(maf = paad,top = 50, clinicalFeatures = 'pred_dicho', sortByAnnotation = TRUE, colors = vc_cols, draw_titv = TRUE, legend_height = 4)

# NOTE: Variants annotated as `Multi_Hit` are those genes which are mutated more than once in the same sample.
# 
# For more details on customisation see the [Customizing oncoplots](http://bioconductor.org/packages/devel/bioc/vignettes/maftools/inst/doc/oncoplots.html) vignette.
# 
# ## Transition and Transversions.
# 
# `titv` function classifies SNPs into [Transitions and Transversions](http://www.mun.ca/biology/scarr/Transitions_vs_Transversions.html) 
# and returns a list of summarized tables in various ways. Summarized data can also be visualized as a boxplot showing overall 
# distribution of six different conversions and as a stacked barplot showing fraction of conversions in each sample.


# list with name, logical vector (passed or not), numeric vector of transition transversion ratio in coding and non coding regions 
paad.titv.high = titv(maf = paad_high, plot = T, useSyn = TRUE)
paad.titv.low = titv(maf = paad_low, plot = T, useSyn = TRUE)

#plot titv summary
#plotTiTv(res = laml.titv)


## Lollipop plots for amino acid changes

# `lollipopPlot` function requires us to have amino acid changes information in the maf file.
# However MAF files have no clear guidelines on naming the field for amino acid changes, 
# with different studies having different field (or column) names for amino acid changes.
# By default, `lollipopPlot` looks for column `AAChange`, and if its not found in the MAF file,
# it prints all available fields with a warning message. For below example, MAF file contains 
# amino acid changes under a field/column name 'Protein_Change'. We will manually specify this using argument `AACol`.
# 
# By default lollipopPlot uses the longest isoform of the gene.

#lollipop plot for DNMT3A, which is one of the most frequent mutated gene in Leukemia.
lollipopPlot(
  maf = paad_high,
  gene = 'TTN',
  AACol = 'amino_acid_change',
  showMutationRate = TRUE
  #labelPos = 200 # sur quelle aa 
)
lollipopPlot(
  maf = paad_low,
  gene = 'TTN',
  AACol = 'amino_acid_change',
  showMutationRate = TRUE
  #labelPos = 223 # sur quelle aa 
)

#General protein domains can be drawn with the function `plotProtein`
plotProtein(gene = "TTN", refSeqID = "NM_133379")


## Rainfall plots

#Cancer genomes, especially solid tumors are characterized by genomic loci with localized hyper-mutations [5](#references).
# Such hyper mutated genomic regions can be visualized by plotting inter variant distance on a linear genomic scale.
# These plots generally called rainfall plots and we can draw such plots using `rainfallPlot`.
# If `detectChangePoints` is set to TRUE, `rainfall` plot also highlights regions where potential changes in inter-event distances are located.


# DRAW PLOT FOR MOST MUTATED SAMPLE
#  identify Kataegis loci
# Kataegis = six consecutive mutations with an average distance of 1000bd
brca <- system.file("extdata", "brca.maf.gz", package = "maftools")
brca <- read.maf(maf = brca, verbose = FALSE)

rainfallPlot(maf = brca, detectChangePoints = TRUE, pointSize = 0.7)
rainfallPlot(maf = paad_high, detectChangePoints = TRUE, pointSize = 0.7)
rainfallPlot(maf = paad_low, detectChangePoints = TRUE, pointSize = 0.7)



# "Kataegis" are defined as those genomic segments containing six or more consecutive mutations with an average inter-mutation distance of less than or equal to 1,00 bp [5](#references).
#   
#   ## Compare mutation load against TCGA cohorts
#   
#   `tcgaCompare` uses mutation load from TCGA [MC3](https://gdc.cancer.gov/about-data/publications/mc3-2017) for 
#   comparing muttaion burden against 33 TCGA cohorts. Plot generated is [similar](http://www.nature.com/nature/journal/v500/n7463/fig_tab/nature12477_F1.html) to the one described in Alexandrov et al [5](#references).
#     
#     
#     name list of cohort = https://www.cancer.gov/research/key-initiatives/ras/ras-central/blog/2017/kras.pdf

high.mutload = tcgaCompare(maf = paad_high, cohortName = 'HIGH_exemple', logscale = TRUE, capture_size = 38.5, col = c("black", "green"))
low.mutload = tcgaCompare(maf = paad_low, cohortName = 'LOW_exemple', logscale = TRUE, capture_size = 38.5, col = c("black", "red"))



















































