
VARPP <- function(ntree, expression, disease_variants, patient_variants) {
  
  require(dplyr)
  
  require(ranger)
  
  exclude_cols <- c("ensembl_gene_id","hgnc_symbol","chromosome_name","start_position","end_position","gene_biotype","Pathogenic","hgnc_id","locus_group","locus_type","location","chromosome")
  
  # These will be progressively filled in loop for random forest votes
  rf_trees <- list()
  
  rf_trees$CADD_expression <- data.frame(
    bind_rows(disease_variants[ , c("GeneVariant","Pathogenic")], benign_variants[ , c("GeneVariant","Pathogenic")]),
    predictNum=vector(length=nrow(disease_variants) + nrow(benign_variants), mode="numeric"),
    predictDenom=vector(length=nrow(disease_variants) + nrow(benign_variants), mode="numeric"))
  
  rf_trees$CADD_expression_varimp <- data.frame(
    Variable=c("CADD_raw_rankscore", colnames(expression)[!colnames(expression) %in% exclude_cols]),
    sumVarimp=vector(length=length(colnames(expression)[!colnames(expression) %in% exclude_cols]) + 1, mode="numeric"),
    ntree=vector(length=length(colnames(expression)[!colnames(expression) %in% exclude_cols]) + 1, mode="numeric"),
    stringsAsFactors=FALSE)
  
  rf_trees$MetaSVM_expression <- rf_trees$CADD_expression
  
  rf_trees$MetaSVM_expression_varimp <- data.frame(
    Variable=c("MetaSVM_rankscore", colnames(expression)[!colnames(expression) %in% exclude_cols]),
    sumVarimp=vector(length=length(colnames(expression)[!colnames(expression) %in% exclude_cols]) + 1, mode="numeric"),
    ntree=vector(length=length(colnames(expression)[!colnames(expression) %in% exclude_cols]) + 1, mode="numeric"),
    stringsAsFactors=FALSE)
  
  rf_trees$CADD_expression_patient <- data.frame(
    variant_id=patient_variants$variant_id,
    predictNum=vector(length=nrow(patient_variants), mode="numeric"),
    predictDenom=vector(length=nrow(patient_variants), mode="numeric"))
  
  rf_trees$MetaSVM_expression_patient <- rf_trees$CADD_expression_patient
  
  for (i in 1:ntree) {
    
    #### Pathogenic variants ####
    
    # sample the genes with replacement
    cls <- sample(unique(disease_variants$Gene), replace=TRUE)
    
    # subset on the sampled clustering factors
    sub <- lapply(cls, FUN=function(x) disease_variants[which(disease_variants$Gene == x), "GeneVariant"])
    
    # sample a single variant within each gene
    sub <- lapply(sub, FUN=function(x) sample(x, size=1))
    
    # join and return samples
    sub <- dplyr::combine(sub)
    
    # Where the sample clusters have been selected, ensure same observation is selected
    # i.e. sample with replacement for clusters but without replacement within clusters
    sub <- sub[match(gsub("_.*$", "", sub), gsub("_.*$", "", sub))]
    
    # Add pathogenic indicator and weight for ranger. Samples with a case weight of zero (OOB sample) are used for variable importance and prediction error
    sub <- bind_cols(GeneVariant=sub,
                     disease_variants[match(sub, disease_variants$GeneVariant), colnames(disease_variants) %in% c("Pathogenic","MetaSVM_rankscore","CADD_raw_rankscore")],
                     weight=rep(1, length(sub)))
    
    # Add expression data
    sub <- bind_cols(sub, expression[match(gsub("_.*$", "", sub$GeneVariant), expression$hgnc_symbol), !(colnames(expression) %in% exclude_cols)])
    
    #### Benign variants ####
    
    # sample the genes with replacement
    cls_benign <- sample(unique(benign_variants$Gene), replace=TRUE)
    
    # Subset benign_variants for faster sampling
    benign_variants_sub <- benign_variants[benign_variants$Gene %in% unique(cls_benign), c("Gene","GeneVariant","Pathogenic","MetaSVM_rankscore","CADD_raw_rankscore")]
    
    # subset on the sampled clustering factors. This is the slowest step in the function.
    sub_benign <- lapply(cls_benign, FUN=function(x) benign_variants_sub[which(benign_variants_sub$Gene == x), "GeneVariant"])

    # sample a single variant within each gene
    sub_benign <- lapply(sub_benign, FUN=function(x) sample(x, size=1))
    
    # join and return samples
    sub_benign <- dplyr::combine(sub_benign)
    
    # Where the sample clusters have been selected, ensure same observation is selected
    # i.e. sample with replacement for clusters but without replacement within clusters
    sub_benign <- sub_benign[match(gsub("_.*$", "", sub_benign), gsub("_.*$", "", sub_benign))]
    
    # Add pathogenic indicator and weight for ranger. Samples with a case weight of zero (OOB sample) are used for variable importance and prediction error
    sub_benign <- bind_cols(GeneVariant=sub_benign,
                            benign_variants_sub[match(sub_benign, benign_variants_sub$GeneVariant), c("Pathogenic","MetaSVM_rankscore","CADD_raw_rankscore")],
                            weight=rep(1, length(sub_benign)))
    
    # Add expression data
    sub_benign <- bind_cols(sub_benign, expression[match(gsub("_.*$", "", sub_benign$GeneVariant), expression$hgnc_symbol), !(colnames(expression) %in% exclude_cols)])
    
    #### OOB samples are those genes not selected in bootstrap sample. Not necessary to restrict to one variant per gene in test set.
    sub_test <- disease_variants[!(disease_variants$Gene %in% cls), c("GeneVariant","Pathogenic","MetaSVM_rankscore","CADD_raw_rankscore")]
    
    # Add expression data
    sub_test <- bind_cols(sub_test,
                          weight=rep(0, nrow(sub_test)),
                          expression[match(gsub("_.*$", "", sub_test$GeneVariant), expression$hgnc_symbol), !(colnames(expression) %in% exclude_cols)])
    
    sub_benign_test <- benign_variants[!(benign_variants$Gene %in% cls_benign), c("GeneVariant","Pathogenic","MetaSVM_rankscore","CADD_raw_rankscore")]
    
    # Add expression data
    sub_benign_test <- bind_cols(sub_benign_test,
                                 weight=rep(0, nrow(sub_benign_test)),
                                 expression[match(gsub("_.*$", "", sub_benign_test$GeneVariant), expression$hgnc_symbol), !(colnames(expression) %in% exclude_cols)])
    
    #### OOB samples are those patient genes not selected in bootstrap sample. Not necessary to restrict to one variant per gene in test set.
    sub_test_patient <- patient_variants[!(patient_variants$genename %in% cls) & !(patient_variants$genename %in% cls_benign), c("variant_id","genename","MetaSVM_rankscore","CADD_raw_rankscore")]
    
    # Add expression data
    sub_test_patient <- data.frame(bind_cols(sub_test_patient, 
                                             expression[match(sub_test_patient$genename, expression$hgnc_symbol), !(colnames(expression) %in% exclude_cols)]),
                                   check.names=FALSE)
    
    sub_test_patient <- sub_test_patient[ , colnames(sub_test_patient) != "genename"]
    
    dat_boot <- list(training=data.frame(bind_rows(sub, sub_test, sub_benign, sub_benign_test), check.names=FALSE),
                     patient=sub_test_patient)
    
    rm(sub, sub_benign, benign_variants_sub, sub_test, cls, cls_benign, sub_benign_test, sub_test_patient)
    
    dat_boot$training$Pathogenic <- factor(dat_boot$training$Pathogenic)
    
    # Random forest does not handle missing data
    dat_boot$training_CADD <- na.omit(dat_boot$training[ , colnames(dat_boot$training) != "MetaSVM_rankscore"])
    
    dat_boot$training_MetaSVM <- na.omit(dat_boot$training[ , colnames(dat_boot$training) != "CADD_raw_rankscore"])
    
    dat_boot$training <- NULL
    
    # Random forest does not handle missing data
    dat_boot$patient_CADD <- na.omit(dat_boot$patient[ , colnames(dat_boot$patient) != "MetaSVM_rankscore"])
    
    dat_boot$patient_MetaSVM <- na.omit(dat_boot$patient[ , colnames(dat_boot$patient) != "CADD_raw_rankscore"])
    
    dat_boot$patient <- NULL
    
    CADD_expression_predict <-
      ranger(
        data = dat_boot$training_CADD[ , !colnames(dat_boot$training_CADD) %in% c("GeneVariant", "weight")],
        dependent.variable.name = "Pathogenic",
        num.trees = 1,
        replace = FALSE,
        sample.fraction = 0.9999,
        importance = "permutation",
        case.weights = dat_boot$training_CADD$weight,
        scale.permutation.importance = FALSE,
        holdout = TRUE
      )
    
    rf_trees$CADD_expression[match(dat_boot$training_CADD$GeneVariant[dat_boot$training_CADD$weight == 0], rf_trees$CADD_expression$GeneVariant), "predictNum"] <- rf_trees$CADD_expression[match(dat_boot$training_CADD$GeneVariant[dat_boot$training_CADD$weight == 0], rf_trees$CADD_expression$GeneVariant), "predictNum"] + as.numeric(as.character(CADD_expression_predict$predictions[dat_boot$training_CADD$weight == 0]))
    
    rf_trees$CADD_expression[match(dat_boot$training_CADD$GeneVariant[dat_boot$training_CADD$weight == 0], rf_trees$CADD_expression$GeneVariant), "predictDenom"] <- rf_trees$CADD_expression[match(dat_boot$training_CADD$GeneVariant[dat_boot$training_CADD$weight == 0], rf_trees$CADD_expression$GeneVariant), "predictDenom"] + 1
    
    rf_trees$CADD_expression_varimp[match(names(CADD_expression_predict$variable.importance), rf_trees$CADD_expression_varimp$Variable), "sumVarimp"] <- rf_trees$CADD_expression_varimp$sumVarimp + CADD_expression_predict$variable.importance
    
    rf_trees$CADD_expression_varimp[match(names(CADD_expression_predict$variable.importance), rf_trees$CADD_expression_varimp$Variable), "ntree"] <- rf_trees$CADD_expression_varimp$ntree + ifelse(CADD_expression_predict$variable.importance == 0, 0, 1)
    
    MetaSVM_expression_predict <-
      ranger(
        data = dat_boot$training_MetaSVM[ , !colnames(dat_boot$training_MetaSVM) %in% c("GeneVariant", "weight")],
        dependent.variable.name = "Pathogenic",
        num.trees = 1,
        replace = FALSE,
        sample.fraction = 0.9999,
        importance = "permutation",
        case.weights = dat_boot$training_MetaSVM$weight,
        scale.permutation.importance = FALSE,
        holdout = TRUE
      )
    
    rf_trees$MetaSVM_expression[match(dat_boot$training_MetaSVM$GeneVariant[dat_boot$training_MetaSVM$weight == 0], rf_trees$MetaSVM_expression$GeneVariant), "predictNum"] <- rf_trees$MetaSVM_expression[match(dat_boot$training_MetaSVM$GeneVariant[dat_boot$training_MetaSVM$weight == 0], rf_trees$MetaSVM_expression$GeneVariant), "predictNum"] + as.numeric(as.character(MetaSVM_expression_predict$predictions[dat_boot$training_MetaSVM$weight == 0]))
    
    rf_trees$MetaSVM_expression[match(dat_boot$training_MetaSVM$GeneVariant[dat_boot$training_MetaSVM$weight == 0], rf_trees$MetaSVM_expression$GeneVariant), "predictDenom"] <- rf_trees$MetaSVM_expression[match(dat_boot$training_MetaSVM$GeneVariant[dat_boot$training_MetaSVM$weight == 0], rf_trees$MetaSVM_expression$GeneVariant), "predictDenom"] + 1
    
    rf_trees$MetaSVM_expression_varimp[match(names(MetaSVM_expression_predict$variable.importance), rf_trees$MetaSVM_expression_varimp$Variable), "sumVarimp"] <- rf_trees$MetaSVM_expression_varimp$sumVarimp + MetaSVM_expression_predict$variable.importance
    
    rf_trees$MetaSVM_expression_varimp[match(names(MetaSVM_expression_predict$variable.importance), rf_trees$MetaSVM_expression_varimp$Variable), "ntree"] <- rf_trees$MetaSVM_expression_varimp$ntree + ifelse(MetaSVM_expression_predict$variable.importance == 0, 0, 1)
    
    CADD_expression_predict_patient <- predict(CADD_expression_predict, data=dat_boot$patient_CADD[ , colnames(dat_boot$patient_CADD) != "variant_id"])
    
    rf_trees$CADD_expression_patient[match(dat_boot$patient_CADD$variant_id, rf_trees$CADD_expression_patient$variant_id), "predictNum"] <- rf_trees$CADD_expression_patient[match(dat_boot$patient_CADD$variant_id, rf_trees$CADD_expression_patient$variant_id), "predictNum"] + as.numeric(as.character(CADD_expression_predict_patient$predictions))
    
    rf_trees$CADD_expression_patient[match(dat_boot$patient_CADD$variant_id, rf_trees$CADD_expression_patient$variant_id), "predictDenom"] <- rf_trees$CADD_expression_patient[match(dat_boot$patient_CADD$variant_id, rf_trees$CADD_expression_patient$variant_id), "predictDenom"] + 1
    
    MetaSVM_expression_predict_patient <- predict(MetaSVM_expression_predict, data=dat_boot$patient_MetaSVM[ , colnames(dat_boot$patient_MetaSVM) != "variant_id"])
    
    rf_trees$MetaSVM_expression_patient[match(dat_boot$patient_MetaSVM$variant_id, rf_trees$MetaSVM_expression_patient$variant_id), "predictNum"] <- rf_trees$MetaSVM_expression_patient[match(dat_boot$patient_MetaSVM$variant_id, rf_trees$MetaSVM_expression_patient$variant_id), "predictNum"] + as.numeric(as.character(MetaSVM_expression_predict_patient$predictions))
    
    rf_trees$MetaSVM_expression_patient[match(dat_boot$patient_MetaSVM$variant_id, rf_trees$MetaSVM_expression_patient$variant_id), "predictDenom"] <- rf_trees$MetaSVM_expression_patient[match(dat_boot$patient_MetaSVM$variant_id, rf_trees$MetaSVM_expression_patient$variant_id), "predictDenom"] + 1
    
    cat(paste("Tree", i), sep="\n")
    
  }
  
  # Drop variants not seen by the classifier
  rf_trees$CADD_expression <- rf_trees$CADD_expression[rf_trees$CADD_expression$predictDenom != 0, ]
  
  rf_trees$MetaSVM_expression <- rf_trees$MetaSVM_expression[rf_trees$MetaSVM_expression$predictDenom != 0, ]

  accuracy_CADD_expression <- data.frame(GeneVariant=rf_trees$CADD_expression$GeneVariant,
                                         Pathogenic=rf_trees$CADD_expression$Pathogenic,
                                         CADD_expression=rf_trees$CADD_expression$predictNum/rf_trees$CADD_expression$predictDenom)
  
  accuracy_MetaSVM_expression <- data.frame(GeneVariant=rf_trees$MetaSVM_expression$GeneVariant,
                                            MetaSVM_expression=rf_trees$MetaSVM_expression$predictNum/rf_trees$MetaSVM_expression$predictDenom)
  
  accuracy <- merge(accuracy_CADD_expression, accuracy_MetaSVM_expression[ , c("GeneVariant","MetaSVM_expression")], by="GeneVariant", all=TRUE)
  
  # Add CADD and MetaSVM scores for comparison to VARPP
  accuracy <- merge(accuracy, bind_rows(disease_variants[ , c("GeneVariant", "CADD_raw_rankscore", "MetaSVM_rankscore")],
                                        benign_variants[ , c("GeneVariant", "CADD_raw_rankscore", "MetaSVM_rankscore")]),
                    by.x="GeneVariant", by.y="GeneVariant", all.x=TRUE, all.y=FALSE)
  
  varimp <- data.frame(Variable=rf_trees$CADD_expression_varimp$Variable,
                       CADD_expression=rf_trees$CADD_expression_varimp$sumVarimp/rf_trees$CADD_expression_varimp$ntree,
                       MetaSVM_expression=rf_trees$MetaSVM_expression_varimp$sumVarimp/rf_trees$MetaSVM_expression_varimp$ntree)
                       
  
  patient_predictions <- data.frame(variant_id=rf_trees$CADD_expression_patient$variant_id,
                                    CADD_expression=rf_trees$CADD_expression_patient$predictNum/rf_trees$CADD_expression_patient$predictDenom,
                                    MetaSVM_expression=rf_trees$MetaSVM_expression_patient$predictNum/rf_trees$MetaSVM_expression_patient$predictDenom)
  
  patient_predictions$CADD_expression[is.nan(patient_predictions$CADD_expression)] <- NA
  
  patient_predictions$MetaSVM_expression[is.nan(patient_predictions$MetaSVM_expression)] <- NA
  
  patient_predictions <- merge(patient_predictions, patient_variants, by.x="variant_id", by.y="variant_id", all=TRUE)
  
  list(accuracy=accuracy, varimp=varimp, patient_predictions=patient_predictions)
  
}
