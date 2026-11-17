/* getarg_/iargc_ for CERNLIB kuip when linking with gfortran */
extern void _gfortran_getarg_i4(int *, char *, int);
extern int  _gfortran_iargc(void);
void getarg_(int *n, char *s, int len) { _gfortran_getarg_i4(n, s, len); }
int  iargc_(void) { return _gfortran_iargc(); }
