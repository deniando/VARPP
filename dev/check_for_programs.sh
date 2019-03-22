#!/usr/bin/env bash
programs=(Rscript java)

printf "Running Sanity checks:\n";

for item in ${programs[*]}
do
    if which $item >/dev/null; then
        printf "%15s found.\n"  $item;
    else
        printf "\nERROR: %s not found!\n\n" $item;
        exit 1;
    fi
done
