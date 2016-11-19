#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.68])
AC_INIT([treepl], [1.0], [eebsmith@umich.edu])
AC_CONFIG_SRCDIR([tnc.h])
AC_CONFIG_HEADERS([config.h])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC

# Checks for libraries.
have_adolc=no
AC_SEARCH_LIBS([main], [adolc], [have_adolc=yes])
if test "x${have_adolc}" = xyes; then
   AC_CHECK_HEADERS([adolc/adouble.h],[],[have_adolc=no])
fi

if test "x${have_adolc}" = xno; then
   AC_MSG_WARN([
--------------------------------------------------------
The adolc library is helpful for parallel excution of the
cross validation procedure. Failure to have this library
may cause problems for convergence.
In the deps folder, unpack the adolc package, go into the
folder and run
autoreconf -fi
./configure --with-openmp-flag=-fopenmp --prefix=/usr
make
sudo make install
--------------------------------------------------------])
fi

AC_CHECK_LIB([m], [fabs])

have_nlopt=no
AC_SEARCH_LIBS([nlopt_urand], [nlopt_cxx], [have_nlopt=yes])

if test "x${have_nlopt}" = xno; then
   AC_MSG_ERROR([
--------------------------------------------------------
The nlopt library is required for optimization. In the 
deps folder, you can unpack nlopt and run
./configure
make
sudo make install
--------------------------------------------------------])
fi

# Checks for header files.
AC_CHECK_HEADERS([fenv.h stdlib.h string.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_SIZE_T

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_STRTOD
AC_CHECK_FUNCS([floor pow sqrt])

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
