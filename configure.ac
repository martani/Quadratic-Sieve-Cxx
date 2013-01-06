AC_INIT([Quadratic Sieve Cxx], [1.0], [martani.net@gmail.com])
AM_INIT_AUTOMAKE([foreign -Wall -Werror])
AC_PROG_CC
AC_PROG_CXX

AC_CHECK_LIB(stdc++, main)

#GMP checks
AC_CHECK_LIB(gmp, __gmpz_legendre, ,
       [AC_MSG_ERROR(
       [GNU MP not found, or not 4.1 or up, see http://gmplib.org/])])
AC_CHECK_LIB(gmpxx, __gmpz_init)

AC_CHECK_HEADER([gmp.h], , [AC_MSG_ERROR([GNU MP not found, see http://gmplib.org/])])

AC_LANG_PUSH(C++)
AC_CHECK_HEADER([gmpxx.h], , [AC_MSG_ERROR([GMP with --enable-cxx is required])])
AC_LANG_POP

AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([Makefile src/Makefile tests/Makefile])
AC_OUTPUT