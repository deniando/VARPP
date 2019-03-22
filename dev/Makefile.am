AUTOMAKE_OPTIONS = serial-tests
AM_COLOR_TESTS=always
EXTRA_DIST = check_for_R_libs.sh check_for_programs.sh

TESTS_ENVIRONMENT = \
      myprog_abspath='$(top_srcdir)'/dev/; \
         export myprog_abspath;\
	      testdatafiledir='$(top_srcdir)'/dev; \
         export testdatafiledir;
	$(VALGRIND)

TESTS = check_for_programs.sh check_for_R_libs.sh 





