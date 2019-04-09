# VARiant Prioritisation by Phenotype (VARPP)

## Contents

- [Overview](#overview)
- [Repo Contents](#repo-contents)
- [System Requirements](#system-requirements)
- [Instructions for Use](#instructions-for-use)
- [License](./LICENSE)
- [Issues](https://github.com/deniando/VARPP/issues)
- [Citation](#citation)

## Overview

Whole genome and exome sequencing is a standard tool for the diagnosis of patients suffering from rare and other genetic disorders. The clinical interpretation of the tens of thousands of genetic variants returned from such tests remains a major challenge.  A common approach is to apply a series of filtering steps to discard variants unlikely to cause the disease.  Here we focus on the problem of prioritising variants in respect to the observed disease phenotype(s) after these filtering steps have been applied.  We do this by linking patterns of gene expression across multiple tissues to the phenotypes which allows us to further narrow down candidates and thereby aid in discovering disease causing variants.  Our method is called VARiant Prioritisation by Phenotype (VARPP), and it is used to prioritise variants in a personalised manner by making use of the patient's phenotype(s).

## Repo Contents

- [data](./data): Data required to run VARPP
- [scripts](./scripts): Scripts to pre-process data and run VARPP

## System Requirements

### Hardware Requirements

VARPP can be run on a standard computer with the following specifications:

RAM: 32+ GB  
CPU: 4+ cores, 1.8+ GHz/core

### Software Requirements {#software}

Before using the functions, users should have R version 3.4.3 or higher and RStudio version 1.1.456 or higher. Install the following packages:

```
install.packages(c("data.table","dplyr","ranger","precrec"))
```

The package versions used were:

```
data.table: 1.12.0
dplyr: 0.8.0.1
ranger: 0.11.2
precrec: 0.10
```

To clone the VARPP GitHub repository, firstly install Git large file storage (https://git-lfs.github.com/) . Then, within Git Bash, navigate to your desired working directory and:

```
git lfs install
git clone https://github.com/deniando/VARPP.git
```

You can then set up a new RStudio project in this directory and proceed with analysis.

## Instructions for Use

Load the required files into R:

```
source(file="scripts/loadData.R")
```

Run Phenolyzer<sup>[1]</sup> for phenotypes of interest.  Phenolyzer returns candidate gene lists associated with human diseases or phenotypes. We used the command line version available from https://github.com/WGLab/phenolyzer.  The below example shows how to use Phenolyzer for a patient with abnormality of cardiovascular system morphology (HP:0030680) and short stature (HP:0004322). The output file of interest is `out.seed_gene_list` and a copy of this file is in [data](./data).
```
perl disease_annotation.pl "HP:0030680;HP:0004322" -p -ph -logistic -addon DB_DISGENET_GENE_DISEASE_SCORE,DB_GAD_GENE_DISEASE_SCORE -addon_weight 0.25
```

Read this file into R and output a list of genes for annotation with dbNSFP<sup>[2]</sup>:
```
genes <- read.table("data/out.seed_gene_list", header=TRUE, stringsAsFactors=FALSE)
write.table(genes$Gene, file="data/genes_for_annotation.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)
```

Annotate variants within the disease gene list using dbNSFP v3.4a (hg38) available from https://sites.google.com/site/jpopgen/dbNSFP. The `disease_variants.txt` output file for this query is in [data](./data).

```
java search_dbNSFP34a -i data/genes_for_annotation.txt -o disease_variants.txt -w 1-12,61,78,121,133,135,137,139,141,188-191
```

Annotate the patient's VCF file using dbNSFP 2.9.3. (hg19). This version of dbNSFP was used because the patient VCF file is based on hg19. We recommend pre-filtering the VCF file to exclude common variants. An example VCF file `patient.vcf.gz`, where variants with a minor allele frequency greater than 0.01 have been filtered out, is available under [data](./data). The `patients_variants.txt` output file for this query is in [data](./data).

```
java search_dbNSFP293 -i data/patient.vcf.gz -o patient_variants.txt -w 1-6,11,49,77,100,102,104,106,108-110,112,114,116,132-134
```

Read in and process the two dbNSFP output files:

```
source(file="scripts/processDBNSFP.R")
```

Source the VARPP script and run the method. The arguments are:

* `ntree`: The number of trees to grow (2000 trees takes approximately 1 hour and 54 minutes)
* `expression`: The gene expression matrix to use in VARPP. We provide expression and specificity data for GTEx and FANTOM5 under [data](./data). These were loaded into the R session at the start of this section.
* `disease_variants`: dbNSFP annotated variants within the disease gene list returned by Phenolyzer
* `patient_variants`: dbNSFP annotated variants from a patient VCF file

This will return a list with the following data frames:

* `accuracy`: Out of bag predicted probabilities of pathogenicity from VARPP. These are used to assess the performance of VARPP for the queried phenotype. There are two classifiers; `CADD_expression` CADD scores in combination with expression data and `MetaSVM_expression` MetaSVM scores in combination with expression data. In addition, `CADD_raw_rankscore` and `MetaSVM_rankscore` are the CADD and MetaSVM scores, whose performance can be compared to VARPP.
* `varimp`: The variable importance measures from VARPP for the two aforementioned classifiers
* `patient_predictions`: Predicted probabilities of pathogenicity from VARPP for the patient variants. These predictions are in columns `CADD_expression` and `MetaSVM_expression`. The remaining columns are the dbNSFP annotations merged in from the `patient_variants` file. Some variants will not have predictions due to missing data in either CADD/MetaSVM scores or gene expression.

```
source(file="scripts/VARPP.R")
VARPP_out <- VARPP(ntree=2000, expression=GTEx_specificity, disease_variants=disease_variants, patient_variants=patient_variants)
```

Calculate performance (area under the ROC curve and area under the precision-recall curve) of VARPP and compare this to the performance of CADD or MetaSVM alone. Also calculate the proportion of ClinVar pathogenic variants in the top 100 predictions of pathogenicity and compare this to CADD or MetaSVM alone. The following script calculates these performance statistics using the `VARPP_out` object:

```
source(file="scripts/performanceStats.R")
```

VARPP can also be run at the command line. Clone the VARPP GitHub repository as described under [Software Requirements](#Software) then:

## Citation

[1] Yang, H, Robinson, PN, Wang, K (2015). Phenolyzer: phenotype-based prioritization of candidate genes for human diseases. Nat. Methods, 12, 9:841-3.

[2] Liu, X, Wu, C, Li, C, Boerwinkle, E (2016). dbNSFP v3.0: A One-Stop Database of Functional Predictions and Annotations for Human Nonsynonymous and Splice-Site SNVs. Hum. Mutat., 37, 3:235-41.
