#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([VARPP], [m4_esyscmd_s([git describe --always])])
AC_CONFIG_AUX_DIR([.])
AM_INIT_AUTOMAKE
AC_PROG_CC
AC_PROG_INSTALL
AC_PROG_RANLIB
AM_SILENT_RULES([no])


AC_CANONICAL_HOST
AC_SUBST([DATE]) DATE="`date`"

AC_SUBST([USER]) USER="`whoami`@`hostname`"

AC_SUBST([UNAME_INFO]) UNAME_INFO=`uname -n`
# Checks for programs.

AC_CHECK_PROG(RSCRIPT, Rscript, yes)
if test x"$RSCRIPT" != x"yes" ; then
AC_MSG_ERROR([Rscript not found - install R!])
fi

AC_CHECK_PROG(JAVA, java, yes)
if test x"$JAVA" != x"yes" ; then
AC_MSG_ERROR([java not found!])
fi


# Checks for libraries.

# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.

AC_CONFIG_FILES([Makefile])
AC_CONFIG_FILES([dev/Makefile])
AC_OUTPUT
AC_MSG_NOTICE([---------------------------------------------------------])
AC_MSG_NOTICE([Configuration complete - $PACKAGE_NAME-$PACKAGE_VERSION])
AC_MSG_NOTICE([])
AC_MSG_NOTICE([General Information:])
AC_MSG_NOTICE([               Version: $VERSION])
AC_MSG_NOTICE([         Configured on: $DATE])
AC_MSG_NOTICE([         Configured by: $USER])
AC_MSG_NOTICE([           Host system: $host_cpu-$host_vendor-$host_os])
AC_MSG_NOTICE([     Uname information: $UNAME_INFO])
AC_MSG_NOTICE([    Installation point: $prefix])
AC_MSG_NOTICE([---------------------------------------------------------])

