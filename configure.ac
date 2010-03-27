AC_INIT([lmfit], [3:2:0], [http://www.messen-und-deuten.de/lmfit])
AC_CONFIG_AUX_DIR([build-aux])
AC_CONFIG_SRCDIR([lib/lmmin.c], [demo/curve1.c]) # assure we have source dir's
AC_CONFIG_HEADERS([config.h])
AM_INIT_AUTOMAKE([-Wall -Werror foreign])
AC_CONFIG_MACRO_DIR([m4])

LT_INIT
AC_PROG_CC
AC_PROG_LIBTOOL
AC_CONFIG_FILES([Makefile lib/Makefile demo/Makefile doc/Makefile])
AC_OUTPUT
