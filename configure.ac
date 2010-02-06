AC_INIT([lmfit], [<version>], [http://www.messen-und-deuten.de/lmfit])
AC_MSG_NOTICE([Now configuring lmfit])
AC_CONFIG_AUX_DIR([build-aux])
AC_CONFIG_SRCDIR([lib/lmmin.c]) # assure we have the right src dir
AC_CONFIG_HEADERS([config.h])
AM_INIT_AUTOMAKE([-Wall -Werror foreign])
AC_CONFIG_MACRO_DIR([m4])

LT_INIT
AC_PROG_CC
AC_PROG_LIBTOOL
AC_CONFIG_FILES([Makefile lib/Makefile demo/Makefile])
AC_OUTPUT
