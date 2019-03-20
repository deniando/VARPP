
library(precrec)

# Have to manually drop NAs because the precrec package does not have options to handle them properly when using argument for na_worst
print(
  evalmod(
    mmdata(scores=join_scores(VARPP_out$accuracy$CADD_expression[!is.na(VARPP_out$accuracy$CADD_expression)],
                              VARPP_out$accuracy$CADD_raw_rankscore[!is.na(VARPP_out$accuracy$CADD_raw_rankscore)],
                              VARPP_out$accuracy$MetaSVM_expression[!is.na(VARPP_out$accuracy$MetaSVM_expression)],
                              VARPP_out$accuracy$MetaSVM_rankscore[!is.na(VARPP_out$accuracy$MetaSVM_rankscore)],
                              chklen=FALSE),
           join_labels(VARPP_out$accuracy$Pathogenic[!is.na(VARPP_out$accuracy$CADD_expression)],
                       VARPP_out$accuracy$Pathogenic[!is.na(VARPP_out$accuracy$CADD_raw_rankscore)],
                       VARPP_out$accuracy$Pathogenic[!is.na(VARPP_out$accuracy$MetaSVM_expression)],
                       VARPP_out$accuracy$Pathogenic[!is.na(VARPP_out$accuracy$MetaSVM_rankscore)],
                       chklen=FALSE),
           modnames=c("CADD_expression", "CADD_raw_rankscore", "MetaSVM_expression", "MetaSVM_rankscore"),
           dsids=1:4)
  )
)

cat("Proportion of pathogenic variants in the top 100 predictions of pathogenicity", sep="\n")
cat("", sep="\n")

print(
  data.frame(Model=c("CADD_expression", "CADD_raw_rankscore", "MetaSVM_expression", "MetaSVM_rankscore"),
             PP100=c(sum(VARPP_out$accuracy[order(-(VARPP_out$accuracy$CADD_expression)), ][1:100, "Pathogenic"])/100,
                     sum(VARPP_out$accuracy[order(-(VARPP_out$accuracy$CADD_raw_rankscore)), ][1:100, "Pathogenic"])/100,
                     sum(VARPP_out$accuracy[order(-(VARPP_out$accuracy$MetaSVM_expression)), ][1:100, "Pathogenic"])/100,
                     sum(VARPP_out$accuracy[order(-(VARPP_out$accuracy$MetaSVM_rankscore)), ][1:100, "Pathogenic"])/100))
)
