AC_INIT([lmfit], [4:0], [http://joachimwuttke.de/lmfit])

# standard configuration for shared libraries
AC_CONFIG_AUX_DIR([build-aux])
AC_CONFIG_MACRO_DIR([m4])

AM_INIT_AUTOMAKE([foreign]) # don't insert GNU standard text files
LT_INIT([disable-static])

AC_PROG_CC
AC_SUBST(AM_CFLAGS,"-pedantic -Wall -Werror")
AC_CONFIG_HEADERS([config.h]) # to avoid endless -D options
AC_HEADER_STDC

# consistency check: is source code present?
AC_CONFIG_SRCDIR([lib/lmmin.c], [demo/curve1.c])

# make these Makefiles
AC_CONFIG_FILES([Makefile lib/Makefile demo/Makefile man/Makefile])
AC_OUTPUT
