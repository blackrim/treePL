// RAD package (Reverse Automatic Differentiation) --
// a package specialized for function and gradient evaluations.
// Written in 2004 by David M. Gay at Sandia National Labs, Albuquerque, NM.

#include <cstddef>
#include <math.h>

 class ADvar;
 class ADvari;
 class ADvar1;
 class ADvar2;
 class Derp;

 struct
ADmemblock {	// We get memory in ADmemblock chunks
		// and never give it back, but reuse it
		// after each ADcontext::Gradcomp() call.
	ADmemblock *next;
	double memblk[1000];
	};

 class
ADcontext {	// A singleton class: one instance in radops.c
	ADmemblock *Busy, *Free;
	char *Mbase;
	size_t Mleft;
	ADmemblock First;
	void *new_ADmemblock(size_t);
 public:
	ADcontext();
	void *Memalloc(size_t len);
	static void Gradcomp();
	};

 class
Derp {		// one derivative-propagation operation
 public:
	static Derp *LastDerp;
	Derp *next;
	const double *a;
	ADvari *b;
	ADvari *c;
	Derp(){};
	Derp(ADvari *);
	Derp(const double *, ADvari *);
	Derp(const double *, ADvari *, ADvari *);
	/* c->aval += a * b->aval; */
	};

 class
ADvari {	// implementation of an ADvar
 public:
	double Val;	// result of this operation
	double aval;	// adjoint -- partial of final result w.r.t. this Val
	void *operator new(size_t len) { return ADvari::adc.Memalloc(len); }
	void operator delete(void*) {} /*Should never be called.*/
	ADvari(double t) { Val = t; aval = 0; }
	ADvari() { Val = 0; aval = 0; }
	static ADcontext adc;
	};

 class
ADvar {		// an "active" variable
	ADvari *cv;
 public:
	ADvar() { /*cv = 0;*/ }
	ADvar(double);
	ADvar(const ADvar &x) { cv = x.cv; }
	ADvar(ADvari *x) { cv = x; }
#ifdef RAD_NO_EQ_ALIAS
	ADvar& operator=(const ADvar &x);
#else	/* allow aliasing v and w after "v = w;" */
	ADvar& operator=(const ADvar &x) { cv = x.cv; return *this; }
#endif
	ADvar& operator=(const double);
	friend ADvar  operator+(const ADvar&);
	friend ADvar  operator-(const ADvar&);
	friend ADvar  operator+ (const ADvar&, const ADvar&);
	ADvar& operator+=(const ADvar&);
	friend ADvar  operator+ (const ADvar&, double);
	ADvar& operator+=(double);
	friend ADvar  operator+ (double L, const ADvar &R);
	friend ADvar  operator- (const ADvar&, const ADvar&);
	ADvar& operator-=(const ADvar&);
	friend ADvar  operator- (const ADvar&, double);
	ADvar& operator-=(double);
	friend ADvar  operator- (double L, const ADvar &R);
	friend ADvar  operator* (const ADvar&, const ADvar&);
	ADvar& operator*=(const ADvar&);
	friend ADvar  operator* (const ADvar&, double);
	ADvar& operator*=(double);
	friend ADvar  operator* (double L, const ADvar &R);
	friend ADvar  operator/ (const ADvar&, const ADvar&);
	ADvar& operator/=(const ADvar&);
	friend ADvar  operator/ (const ADvar&, double);
	ADvar& operator/=(double);
	friend ADvar atan(const ADvar&);
	friend ADvar atan2(const ADvar&, const ADvar&);
	friend ADvar atan2(double, const ADvar&);
	friend ADvar atan2(const ADvar&, double);
	friend ADvar cos (const ADvar&);
	friend ADvar exp (const ADvar&);
	friend ADvar log (const ADvar&);
	friend ADvar pow (const ADvar&, const ADvar&);
	friend ADvar pow (double, const ADvar&);
	friend ADvar pow (const ADvar&, double);
	friend ADvar sin (const ADvar&);
	friend ADvar sqrt(const ADvar&);
	friend ADvar tan (const ADvar&);
	friend int operator<(const ADvar&, const ADvar&);
	friend int operator<(const ADvar&, double);
	friend int operator<(double, const ADvar&);
	friend int operator<=(const ADvar&, const ADvar&);
	friend int operator<=(const ADvar&, double);
	friend int operator<=(double,const ADvar&);
	friend int operator==(const ADvar&, const ADvar&);
	friend int operator==(const ADvar&, double);
	friend int operator==(double, const ADvar&);
	friend int operator!=(const ADvar&, const ADvar&);
	friend int operator!=(const ADvar&, double);
	friend int operator!=(double, const ADvar&);
	friend int operator>=(const ADvar&, const ADvar&);
	friend int operator>=(const ADvar&, double);
	friend int operator>=(double, const ADvar&);
	friend int operator>(const ADvar&, const ADvar&);
	friend int operator>(const ADvar&, double);
	friend int operator>(double, const ADvar&);

	operator ADvari*() { return cv; }

	double val() { return cv->Val; }
	double adj() { return cv->aval; }
	};

 class
ADvar1: public ADvari {	// simplest unary ops
 public:
	Derp d;
	ADvar1(double val1): ADvari(val1) {}
	ADvar1(double val1, ADvari *c1): d(c1) { Val = val1; }
	ADvar1(double val1, const double *a1, ADvari *c1): d(a1,this,c1) { Val = val1; }
	};

 class
ADvar1s: public ADvar1 { // unary ops with partial "a"
 public:
	double a;
	ADvar1s(double val1, double a1, ADvari *c1): ADvar1(val1,&a,c1), a(a1) {}
	};

 class
ADvar2: public ADvari {	// basic binary ops
 public:
	Derp dL, dR;
	ADvar2(double val1): ADvari(val1) {}
	ADvar2(double val1, ADvari *Lcv, const double *Lc, ADvari *Rcv, const double *Rc):
			ADvari(val1) {
		dR.next = Derp::LastDerp;
		dL.next = &dR;
		Derp::LastDerp = &dL;
		dL.a = Lc;
		dL.c = Lcv;
		dR.a = Rc;
		dR.c = Rcv;
		dL.b = dR.b = this;
		}
	};

 class
ADvar2q: public ADvar2 { // binary ops with partials "a", "b"
 public:
	double a, b;
	ADvar2q(double val1, double Lp, double Rp, ADvari *Lcv, ADvari *Rcv):
			ADvar2(val1), a(Lp), b(Rp) {
		dR.next = Derp::LastDerp;
		dL.next = &dR;
		Derp::LastDerp = &dL;
		dL.a = &a;
		dL.c = Lcv;
		dR.a = &b;
		dR.c = Rcv;
		dL.b = dR.b = this;
		}
	};


ADvar operator+ (double L, const ADvar &R);
ADvar operator* (double L, const ADvar &R);
ADvar atan(const ADvar&);
ADvar atan2(const ADvar&, const ADvar&);
ADvar atan2(double, const ADvar&);
ADvar atan2(const ADvar&, double);
ADvar cos (const ADvar&);
ADvar exp (const ADvar&);
ADvar log (const ADvar&);
ADvar pow (const ADvar&, const ADvar&);
ADvar pow (double, const ADvar&);
ADvar pow (const ADvar&, double);
ADvar sin (const ADvar&);
ADvar sqrt(const ADvar&);
ADvar tan (const ADvar&);

inline ADvar operator+(const ADvar &T) { ADvar rv(T.cv); return rv; }

inline int operator<(const ADvar &L, const ADvar &R) { return L.cv->Val < R.cv->Val; }
inline int operator<(const ADvar &L, double R) { return L.cv->Val < R; }
inline int operator<(double L, const ADvar &R) { return L < R.cv->Val; }

inline int operator<=(const ADvar &L, const ADvar &R) { return L.cv->Val <= R.cv->Val; }
inline int operator<=(const ADvar &L, double R) { return L.cv->Val <= R; }
inline int operator<=(double L, const ADvar &R) { return L <= R.cv->Val; }

inline int operator==(const ADvar &L, const ADvar &R) { return L.cv->Val == R.cv->Val; }
inline int operator==(const ADvar &L, double R) { return L.cv->Val == R; }
inline int operator==(double L, const ADvar &R) { return L == R.cv->Val; }

inline int operator!=(const ADvar &L, const ADvar &R) { return L.cv->Val != R.cv->Val; }
inline int operator!=(const ADvar &L, double R) { return L.cv->Val != R; }
inline int operator!=(double L, const ADvar &R) { return L != R.cv->Val; }

inline int operator>=(const ADvar &L, const ADvar &R) { return L.cv->Val >= R.cv->Val; }
inline int operator>=(const ADvar &L, double R) { return L.cv->Val >= R; }
inline int operator>=(double L, const ADvar &R) { return L >= R.cv->Val; }

inline int operator>(const ADvar &L, const ADvar &R) { return L.cv->Val > R.cv->Val; }
inline int operator>(const ADvar &L, double R) { return L.cv->Val > R; }
inline int operator>(double L, const ADvar &R) { return L > R.cv->Val; }

inline void *ADcontext::Memalloc(size_t len) {
		if (Mleft >= len)
			return Mbase + (Mleft -= len);
		return new_ADmemblock(len);
		}

inline Derp::Derp(ADvari *c1): c(c1) {
		next = LastDerp;
		LastDerp = this;
		}

inline Derp::Derp(const double *a1, ADvari *c1): a(a1), c(c1) {
		next = LastDerp;
		LastDerp = this;
		}

inline Derp::Derp(const double *a1, ADvari *b1, ADvari *c1): a(a1), b(b1), c(c1) {
		next = LastDerp;
		LastDerp = this;
		}
