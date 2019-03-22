#!/usr/bin/env Rscript
libraries <- c("data.table","dplyr","ranger","precrec")

toinst = c();

for(library in libraries)
{
    f = is.element(library, installed.packages()[,1])
    print(paste("Library",library, "is installed?", f))
    if(!f)
    {
        toinst <- c(toinst, library)
    }
}

if(length(toinst))
{
    message("")
    message("ERROR: one or more required R packages are not installed:")
    message(paste(as.character(toinst),collapse=" ",sep=" "))
    message("Install using:\n R -e 'install.packages(c(",paste("\"",as.character(toinst),"\"",collapse=", ",sep=""),"), repos=\"http://cran.us.r-project.org\")'")
    message("")
    message("")
    quit(status=1)
}

quit(status=0)
