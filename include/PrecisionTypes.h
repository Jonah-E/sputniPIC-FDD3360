#ifndef PRECISIONTYPES_H
#define PRECISIONTYPES_H

#ifdef USE_DOUBLE
// precision for particle quanties (position, velocity, ...): choose double or double
typedef double FPpart;
// precision for field quantities (Ex, Bx, ...): choose double or double
typedef double FPfield;
// precision for interpolated quanties (Ex, Bx, ...): choose double or double
typedef double FPinterp;

#else

// precision for particle quanties (position, velocity, ...): choose float or double
typedef float FPpart;
// precision for field quantities (Ex, Bx, ...): choose float or double
typedef float FPfield;
// precision for interpolated quanties (Ex, Bx, ...): choose float or double
typedef float FPinterp;
#endif

#endif
