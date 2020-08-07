/*
  to be called from FORTRAN to get the system time 
  use both conventions (with and without underscore)
*/

#include <unistd.h>
#include <sys/times.h>

void myclock_(double *utime) {
  static struct tms t;
  static double     tickps;
  tickps= (double) sysconf(_SC_CLK_TCK);
  times(&t);
  *utime= (double)((t.tms_utime+t.tms_stime)*100.0/tickps);
}

void myclock(double *utime) {
  static struct tms t;
  static double     tickps;
  tickps= (double) sysconf(_SC_CLK_TCK);
  times(&t);
  *utime= (double)((t.tms_utime+t.tms_stime)*100.0/tickps);
}

#ifdef NO_TICK_DEF
#undef NO_TICK_DEF
#undef CLK_TCK 
#endif
