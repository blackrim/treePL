// Support routines for the RAD package (Reverse Automatic Differentiation) --
// a package specialized for function and gradient evaluations.
// Written in 2004 by David M. Gay at Sandia National Labs, Albuquerque, NM.

#include "myrad.h"

Derp *Derp::LastDerp = 0;

ADcontext ADvari::adc;

double One = 1., negOne = -1.;

ADcontext::ADcontext() {
    First.next = 0;
    Busy = &First;
    Free = 0;
    Mbase = (char*)First.memblk;
    Mleft = sizeof(First.memblk); 
}

void* ADcontext::new_ADmemblock(size_t len){
    ADmemblock *x;
    
    if (x = Free)
	Free = x->next;
    else
	x = new ADmemblock;
    x->next = Busy;
    Busy = x;
    return (Mbase = (char*)x->memblk) +
	(Mleft = sizeof(First.memblk) - len);
}

void ADcontext::Gradcomp(){
    ADmemblock *mb, *mb0, *mb1, *mbf;
    Derp *d = Derp::LastDerp;
    d->b->aval = 1;
    for(; d; d = d->next)
	d->c->aval += *d->a * d->b->aval;
    Derp::LastDerp = 0;
    mb0 = &ADvari::adc.First;
    mbf =  ADvari::adc.Free;
    for(mb = ADvari::adc.Busy; 
	mb != mb0; mb = mb1) {
	mb1 = mb->next;
	mb->next = mbf;
	mbf = mb;
    }
    ADvari::adc.Free = mbf;
    ADvari::adc.Busy = mb;
    ADvari::adc.Mbase = (char*)ADvari::adc.First.memblk;
    ADvari::adc.Mleft = sizeof(ADvari::adc.First.memblk);
}

ADvar::ADvar(double d){
    ADvari *x = new ADvari(d);
    cv = x;
}

#ifdef RAD_NO_EQ_ALIAS
 ADvar&
ADvar::operator=(const ADvar &x)
{ cv = new ADvar1(x.cv->Val, &One, x.cv); return *this; }
#endif

 ADvar&
ADvar::operator=(const double d)
{ cv = new ADvari(d); return *this; }

 ADvar
operator-(const ADvar &T) {
	return ADvar(new ADvar1(-T.cv->Val, &negOne, T.cv));
	}

 ADvar
operator+(const ADvar &L, const ADvar &R) {
	ADvari *Lcv = L.cv, *Rcv = R.cv;
	return ADvar(new ADvar2(Lcv->Val + Rcv->Val, Lcv, &One, Rcv, &One));
	}

 ADvar&
ADvar::operator+=(const ADvar &R) {
	ADvari *Lcv = cv, *Rcv = R.cv;
	cv = new ADvar2(Lcv->Val + Rcv->Val, Lcv, &One, Rcv, &One);
	return *this;
	}

 ADvar
operator+(const ADvar &L, double R) {
	ADvari *tcv = L.cv;
	return ADvar(new ADvar1(tcv->Val + R, &One, tcv));
	}

 ADvar&
ADvar::operator+=(double R) {
	ADvari *tcv = cv;
	cv = new ADvar1(tcv->Val + R, &One, tcv);
	return *this;
	}

 ADvar
operator+(double L, const ADvar &R) {
	ADvari *Rcv = R.cv;
	return ADvar(new ADvar1(L + Rcv->Val, &One, Rcv));
	}

 ADvar
operator-(const ADvar &L, const ADvar &R) {
	ADvari *Lcv = L.cv, *Rcv = R.cv;
	return ADvar(new ADvar2(Lcv->Val - Rcv->Val, Lcv, &One, Rcv, &negOne));
	}

 ADvar&
ADvar::operator-=(const ADvar &R) {
	ADvari *Lcv = cv, *Rcv = R.cv;
	cv = new ADvar2(Lcv->Val - Rcv->Val, Lcv, &One, Rcv, &negOne);
	}

 ADvar
operator-(const ADvar &L, double R) {
	ADvari *tcv = L.cv;
	return ADvar(new ADvar1(tcv->Val - R, &One, tcv));
	}

 ADvar&
ADvar::operator-=(double R) {
	ADvari *tcv = cv;
	cv = new ADvar1(tcv->Val - R, &One, tcv);
	return *this;
	}

 ADvar
operator-(double L, const ADvar &R) {
	ADvari *Rcv = R.cv;
	return ADvar(new ADvar1(L - Rcv->Val, &negOne, Rcv));
	}

 ADvar
operator*(const ADvar &L, const ADvar &R) {
	ADvari *Lcv = L.cv, *Rcv = R.cv;
	return ADvar(new ADvar2(Lcv->Val * Rcv->Val, Lcv, &Rcv->Val, Rcv, &Lcv->Val));
	}

 ADvar&
ADvar::operator*=(const ADvar &R) {
	ADvari *Lcv = cv, *Rcv = R.cv;
	cv = new ADvar2(Lcv->Val * Rcv->Val, Lcv, &Rcv->Val, Rcv, &Lcv->Val);
	}

 ADvar
operator*(const ADvar &L, double R) {
	ADvari *Lcv = L.cv;
	return ADvar(new ADvar1s(Lcv->Val * R, R, Lcv));
	}

 ADvar&
ADvar::operator*=(double R) {
	ADvari *Lcv = cv;
	cv = new ADvar1s(Lcv->Val * R, R, Lcv);
	return *this;
	}

 ADvar
operator*(double L, const ADvar &R) {
	ADvari *Rcv = R.cv;
	return ADvar(new ADvar1s(L * Rcv->Val, L, Rcv));
	}

 ADvar
operator/(const ADvar &L, const ADvar &R) {
	ADvari *Lcv = L.cv, *Rcv = R.cv;
	double Lv = Lcv->Val, Rv = Rcv->Val, pL = 1. / Rv, q = Lv/Rv;
	return ADvar(new ADvar2q(q, pL, -q*pL, Lcv, Rcv));
	}

 ADvar&
ADvar::operator/=(const ADvar &R) {
	ADvari *Lcv = cv, *Rcv = R.cv;
	double Lv = Lcv->Val, Rv = Rcv->Val, pL = 1. / Rv, q = Lv/Rv;
	cv = new ADvar2q(q, pL, -q*pL, Lcv, Rcv);
	}

 ADvar
operator/(const ADvar &L, double R) {
	ADvari *Lcv = L.cv;
	return ADvar(new ADvar1s(Lcv->Val / R, 1./R, Lcv));
	}

 ADvar&
ADvar::operator/=(double R) {
	ADvari *Lcv = cv;
	cv = new ADvar1s(Lcv->Val / R, 1./R, Lcv);
	return *this;
	}

 ADvar
atan(const ADvar &v) {
	ADvari *tcv = v.cv;
	double t = tcv->Val;
	return ADvar(new ADvar1s(atan(t), 1./(1. + t*t), tcv));
	}

 ADvar
atan2(const ADvar &L, const ADvar &R) {
	ADvari *Lcv = L.cv, *Rcv = R.cv;
	double x = Lcv->Val, y = Rcv->Val, t = x*x + y*y;
	return ADvar(new ADvar2q(atan2(x,y), y/t, -x/t, Lcv, Rcv));
	}

 ADvar
atan2(double x, const ADvar &R) {
	ADvari *Rcv = R.cv;
	double y = Rcv->Val, t = x*x + y*y;
	return ADvar(new ADvar1s(atan2(x,y), -x/t, Rcv));
	}

 ADvar
atan2(const ADvar &L, double y) {
	ADvari *Lcv = L.cv;
	double x = Lcv->Val, t = x*x + y*y;
	return ADvar(new ADvar1s(atan2(x,y), y/t, Lcv));
	}

 ADvar
cos(const ADvar &v) {
	ADvari *tcv = v.cv;
	return ADvar(new ADvar1s(cos(tcv->Val), -sin(tcv->Val), tcv));
	}

 ADvar
exp(const ADvar &v) {
	ADvar1* rcv;
	ADvari *tcv = v.cv;
	ADvar rv(rcv = new ADvar1(exp(tcv->Val), tcv));
	rcv->d.a = &rcv->Val;
	rcv->d.b = rcv;
	return rv;
	}

 ADvar
log(const ADvar &v) {
	ADvari *tcv = v.cv;
	double x = tcv->Val;
	return ADvar(new ADvar1s(log(x), 1. / x, tcv));
	}

 ADvar
pow(const ADvar &L, const ADvar &R) {
	ADvari *Lcv = L.cv, *Rcv = R.cv;
	double x = Lcv->Val, y = Rcv->Val, t = pow(x,y);
	return ADvar(new ADvar2q(t, y*t/x, t*log(x), Lcv, Rcv));
	}

 ADvar
pow(double x, const ADvar &R) {
	ADvari *Rcv = R.cv;
	double t = pow(x,Rcv->Val);
	return ADvar(new ADvar1s(t, t*log(x), Rcv));
	}

 ADvar
pow(const ADvar &L, double y) {
	ADvari *Lcv = L.cv;
	double x = Lcv->Val, t = pow(x,y);
	return ADvar(new ADvar1s(t, y*t/x, Lcv));
	}

 ADvar
sin(const ADvar &v) {
	ADvari *tcv = v.cv;
	return ADvar(new ADvar1s(sin(tcv->Val), cos(tcv->Val), tcv));
	}

 ADvar
sqrt(const ADvar &v) {
	ADvari *tcv = v.cv;
	double t = sqrt(tcv->Val);
	return ADvar(new ADvar1s(t, 0.5/t, tcv));
	}

 ADvar
tan(const ADvar &v) {
	ADvari *tcv = v.cv;
	double t = cos(tcv->Val);
	return ADvar(new ADvar1s(tan(tcv->Val), 1./(t*t), tcv));
	}
