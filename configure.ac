AC_INIT([EGAnet], [1.3.0]) # Replace with your package name and version
AC_PREREQ([2.69]) # Replace with your Autoconf version

AC_CHECK_LIB(pthread, pthread_create, 
             [AC_SUBST([PTHREAD_LIBS], [-lpthread])],
             [AC_CHECK_LIB(winpthread, pthread_create,
                           [AC_SUBST([PTHREAD_LIBS], [-lwinpthread])],
                           [AC_MSG_ERROR([libpthread or libwinpthread required but not found])])])

AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT
