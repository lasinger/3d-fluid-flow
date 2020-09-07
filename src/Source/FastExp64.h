static const double log2_inv = 1.0 / log(2.0);

// Approximation for exponential function, code from Bernardt Duvenhage (https://bduvenhage.me/performance/machine_learning/2019/06/04/fast-exp.html)
// on VS 2017 performance improvement is not that significant
//! Approximate exp adapted from Schraudolph, 1999 - double precision floating point version.
inline
double fast_exp_64(const double x) {
	// Based on Schraudolph 1999, A Fast, Compact Approximation of the Exponential Function.
	// - Adapted to use 64-bit integer; reduces staircase effect.
	// - Valid for x in approx range (-700, 700).
	union { double d_; int64_t i_; } uid; //This could be moved to the thread scope.
	//BBBD(sizeof(uid)!=8)
	uid.i_ = int64_t(double((int64_t(1) << 52) * log2_inv) * x + double((int64_t(1) << 52) * 1023 - 0)); //c=0 for 1.0 at zero.
	return uid.d_;
}

