#include "f2c.h"
#include <vector>
#include <string>
/* Table of constant values */

static int c__9 = 9;
static int c__1 = 1;
static int c__3 = 3;
static int c__5 = 5;

int main (void)
{
    /* System generated locals */
    int i__1, i__2;
    double d__1;
    
    std::vector<int> o__1;
    std::vector<long> cl__1;

    /* Builtin functions */
    int s_wsle(std::vector<int>* ilist), do_lio(int *, int *, std::string *, ftnlen), 
	    e_wsle(void), s_rsle(cilist *), e_rsle(void);
    /* Subroutine */ int s_stop(char *, ftnlen);
    integer f_open(olist *), f_clos(cllist *);
    double log(doublereal), exp(doublereal);

    /* Local variables */
    static int i__, j, n, p;
    static double m1;
    static int jj, kkk, nop[10];
    static double dist[1000010]	/* was [10][100001] */;
    extern int sstmm_(void);
    static int ncycle;
    extern int genrand_(double *);
    extern double faculty_(int *);

    /* Fortran I/O blocks */
    static cilist io___2 = { 0, 6, 0, 0, 0 };
    static cilist io___3 = { 0, 5, 0, 0, 0 };
    static cilist io___5 = { 0, 6, 0, 0, 0 };
    static cilist io___6 = { 0, 5, 0, 0, 0 };
    static cilist io___8 = { 0, 6, 0, 0, 0 };
    static cilist io___9 = { 0, 5, 0, 0, 0 };
    static cilist io___17 = { 0, 21, 0, 0, 0 };
    static cilist io___18 = { 0, 21, 0, 0, 0 };
    static cilist io___19 = { 0, 21, 0, 0, 0 };


/* ccccccccccccccccccccccccccccccccccccccccccccccccc */
/*     Divide N Particles Among P Compartments    C */
/* ccccccccccccccccccccccccccccccccccccccccccccccccc */
/* ccccccccccccccccccccccccccccccccccccccccccc */
/*     Initialize Random Number Generator   C */
/* ccccccccccccccccccccccccccccccccccccccccccc */
    m1 = (doublereal) ((sstmm_() * 10 + 10) % 1000) * .001;
    if (m1 < .001) {
	m1 = .001;
    }
    if (m1 > .999) {
	m1 = .999;
    }
    genrand_(&m1);
/* cccccccccccccccccccccccccccccccccccccc */
/*     Read Info From Standard Input   C */
/* cccccccccccccccccccccccccccccccccccccc */
    s_wsle(&io___2);
    do_lio(&c__9, &c__1, "Number Of Particles    ? ", (ftnlen)25);
    e_wsle();
    s_rsle(&io___3);
    do_lio(&c__3, &c__1, (char *)&n, (ftnlen)sizeof(integer));
    e_rsle();
    s_wsle(&io___5);
    do_lio(&c__9, &c__1, "Number Of Compartments ? ", (ftnlen)25);
    e_wsle();
    s_rsle(&io___6);
    do_lio(&c__3, &c__1, (char *)&p, (ftnlen)sizeof(integer));
    e_rsle();
    s_wsle(&io___8);
    do_lio(&c__9, &c__1, "Number Of Cycles       ? ", (ftnlen)25);
    e_wsle();
    s_rsle(&io___9);
    do_lio(&c__3, &c__1, (char *)&ncycle, (ftnlen)sizeof(integer));
    e_rsle();
    if (p < 2 || p > 10 || n < 2 || n > 100000) {
	s_stop("", (ftnlen)0);
    }
/* ccccccccccccccccccc */
/*     Initialize   C */
/* ccccccccccccccccccc */
    i__1 = n;
    for (jj = 0; jj <= i__1; ++jj) {
	i__2 = p;
	for (j = 1; j <= i__2; ++j) {
	    dist[j + jj * 10 - 1] = 0.;
	    nop[j - 1] = 0;
	}
    }
/* ccccccccccccccccccccccccccccc */
/*     Loop Over All Cycles   C */
/* ccccccccccccccccccccccccccccc */
    i__1 = ncycle;
    for (i__ = 1; i__ <= i__1; ++i__) {
	for (kkk = 1; kkk <= 1000; ++kkk) {
/* cccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/*     Distribute Particles                            C */
/*                                                     C */
/*     1. Loop Over All Particles                      C */
/*     2. Generate A Random Compartment (1...P)        C */
/*     3. Put This Particle In The Compartment         C */
/*                                                     C */
/*     Nop(J) = Number Of Particles In Compartment J   C */
/* cccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/*     Start Modification */
/*     End   Modification */
/* ccccccccccccccccccccccccccccc */
/*     Make Histrogram        C */
/* ccccccccccccccccccccccccccccc */
	    i__2 = p;
	    for (j = 1; j <= i__2; ++j) {
		dist[j + nop[j - 1] * 10 - 1] += 1.;
		nop[j - 1] = 0;
	    }
	}
    }
/* ccccccccccccccccccccccccccccc */
/*     Write Results          C */
/* ccccccccccccccccccccccccccccc */
    o__1.oerr = 0;
    o__1.ounit = 21;
    o__1.ofnmlen = 10;
    o__1.ofnm = "output.dat";
    o__1.orl = 0;
    o__1.osta = 0;
    o__1.oacc = 0;
    o__1.ofm = 0;
    o__1.oblnk = 0;
    f_open(&o__1);
    i__1 = n;
    for (jj = 0; jj <= i__1; ++jj) {
	i__2 = p;
	for (j = 1; j <= i__2; ++j) {
	    if (dist[j + jj * 10 - 1] > .5) {
		dist[j + jj * 10 - 1] /= (doublereal) (ncycle * 1000);
	    }
	}
	s_wsle(&io___17);
	do_lio(&c__3, &c__1, (char *)&jj, (ftnlen)sizeof(integer));
	i__2 = p;
	for (j = 1; j <= i__2; ++j) {
	    do_lio(&c__5, &c__1, (char *)&dist[j + jj * 10 - 1], (ftnlen)
		    sizeof(doublereal));
	}
	e_wsle();
    }
    cl__1.cerr = 0;
    cl__1.cunit = 21;
    cl__1.csta = 0;
    f_clos(&cl__1);
/* ccccccccccccccccccccccccccccccc */
/*     Write Analytical Dist.   C */
/*     For P=2                  C */
/* ccccccccccccccccccccccccccccccc */
    if (p == 2) {
	o__1.oerr = 0;
	o__1.ounit = 21;
	o__1.ofnmlen = 14;
	o__1.ofnm = "analytical.dat";
	o__1.orl = 0;
	o__1.osta = 0;
	o__1.oacc = 0;
	o__1.ofm = 0;
	o__1.oblnk = 0;
	f_open(&o__1);
	i__1 = n;
	for (jj = 0; jj <= i__1; ++jj) {
	    s_wsle(&io___18);
	    do_lio(&c__3, &c__1, (char *)&jj, (ftnlen)sizeof(integer));
	    i__2 = n - jj;
	    d__1 = exp(faculty_(&n) - faculty_(&jj) - faculty_(&i__2) - (
		    doublereal) n * log((doublereal) p));
	    do_lio(&c__5, &c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
	    e_wsle();
	}
	cl__1.cerr = 0;
	cl__1.cunit = 21;
	cl__1.csta = 0;
	f_clos(&cl__1);
    } else {
	o__1.oerr = 0;
	o__1.ounit = 21;
	o__1.ofnmlen = 14;
	o__1.ofnm = "analytical.dat";
	o__1.orl = 0;
	o__1.osta = 0;
	o__1.oacc = 0;
	o__1.ofm = 0;
	o__1.oblnk = 0;
	f_open(&o__1);
	s_wsle(&io___19);
	e_wsle();
	cl__1.cerr = 0;
	cl__1.cunit = 21;
	cl__1.csta = 0;
	f_clos(&cl__1);
    }
/* ccccccccccccccccccccccccccc */
/*     End Of The Program   C */
/* ccccccccccccccccccccccccccc */
    s_stop("", (ftnlen)0);
    return 0;
} /* MAIN__ */

/* Main program alias */ int distribution_ () { MAIN__ (); return 0; }
