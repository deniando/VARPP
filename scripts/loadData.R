
benign_variants <- read.table("data/benign_variants.txt", header=TRUE, stringsAsFactors=FALSE)

benign_variants$Pathogenic <- 0

FANTOM5_expression <- read.csv("data/FANTOM5_expression.csv", stringsAsFactors=FALSE)

FANTOM5_specificity <- read.csv("data/FANTOM5_specificity.csv", stringsAsFactors=FALSE)

GTEx_expression <- read.csv("data/GTEx_expression.csv", stringsAsFactors=FALSE)

GTEx_specificity <- read.csv("data/GTEx_specificity.csv", stringsAsFactors=FALSE)
