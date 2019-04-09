AUTOMAKE_OPTIONS = serial-tests
AM_COLOR_TESTS=always
EXTRA_DIST = check_for_R_libs.sh check_for_programs.sh

TESTS = check_for_programs.sh check_for_R_libs.sh


testrun-local:
		../scripts/runVAARP.sh  -b ../data/benign_variants.txt -p ../data/disease_variants.txt  -t ../data/patient_variants.txt -n 1  -e ../data/GTEx_specificity.csv
