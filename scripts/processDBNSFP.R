
library(data.table)

disease_variants <- fread(input=paste0(getwd(), "/data/disease_variants.txt"), header=TRUE, na.strings=".")

# Drop variants with multiple gene annotation
disease_variants <- disease_variants[!(grepl(";", disease_variants$genename)), ]

# Subset to ClinVar pathogenic variants with a clinical significance of 5
disease_variants <- disease_variants[grepl(pattern="^5$|^(5\\|)+5$", x=disease_variants$clinvar_clnsig), ]

# Drop pathogenic variants that are common in large sequencing cohorts
disease_variants <- disease_variants[(disease_variants$`1000Gp3_AF` < 0.01 | is.na(disease_variants$`1000Gp3_AF`)) & (disease_variants$TWINSUK_AF < 0.01 | is.na(disease_variants$TWINSUK_AF)) & (disease_variants$ALSPAC_AF < 0.01 | is.na(disease_variants$ALSPAC_AF)) & (disease_variants$ESP6500_AA_AF < 0.01 | is.na(disease_variants$ESP6500_AA_AF)) & (disease_variants$ESP6500_EA_AF < 0.01 | is.na(disease_variants$ESP6500_EA_AF)) & (disease_variants$ExAC_AF < 0.01 | is.na(disease_variants$ExAC_AF)), ]

disease_variants <- as.data.frame(disease_variants, stringsAsFactors=FALSE)

# Sequentially number variants within the same gene as this will be required for sampling
disease_variants$variant <- as.numeric(ave(disease_variants$genename, disease_variants$genename, FUN=seq_along))

colnames(disease_variants)[colnames(disease_variants)== "genename"] <- "Gene"

colnames(disease_variants)[colnames(disease_variants) == "#chr"] <- "Chr"

disease_variants$GeneVariant <- paste(disease_variants$Gene, disease_variants$variant, sep="_")

disease_variants$Pathogenic <- 1

# Drop any genes in the benign variant set that are in the disease variant set. This will only be relevant if using later versions of Phenolyzer and dbNSFP as there will be no overlap when using stated versions.
benign_variants <- benign_variants[!(benign_variants$Gene %in% disease_variants$Gene), ]

patient_variants <- read.delim(file="data/patient_variants.txt", stringsAsFactors=FALSE, na.strings=".", check.names=FALSE)

patient_variants$variant_id <- 1:nrow(patient_variants)

colnames(patient_variants)[colnames(patient_variants) == "#chr"] <- "Chr"
