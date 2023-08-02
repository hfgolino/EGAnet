AC_INIT([EGAnet], [1.3.0])
AC_PREREQ([2.69])

# Test for the pthread library in a platform-independent way
saved_LIBS="$LIBS"
LIBS="-lpthread $LIBS"
AC_LINK_IFELSE([AC_LANG_PROGRAM([[#include <pthread.h>]], [[pthread_create(0, 0, 0, 0);]])],
               [AC_SUBST([PTHREAD_LIBS], [-lpthread]) pthread_found="yes"],
               [pthread_found="no"])
LIBS="$saved_LIBS"

# If pthread was not found, check for winpthread
if test "x$pthread_found" = "xno"; then
  AC_CHECK_LIB(winpthread, pthread_create,
               [AC_SUBST([PTHREAD_LIBS], [-lwinpthread])],
               [AC_MSG_ERROR([libpthread or libwinpthread required but not found])])
fi

AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT
