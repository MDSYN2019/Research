
#include "f2c.h"

/* Table of constant values */

static integer c__9 = 9;
static integer c__1 = 1;
static integer c__3 = 3;
static integer c__5 = 5;

/* Main program */ int MAIN__(void)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;
    olist o__1;
    cllist cl__1;

    /* Builtin functions */
    integer s_wsle(cilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_wsle(void), s_rsle(cilist *), e_rsle(void);
    /* Subroutine */ int s_stop(char *, ftnlen);
    double exp(doublereal);
    integer f_open(olist *), s_wsfe(cilist *), do_fio(integer *, char *, 
	    ftnlen), e_wsfe(void), f_clos(cllist *);

    /* Local variables */
    static integer i__, n;
    static doublereal beta, temp, norm, distri[10001];

    /* Fortran I/O blocks */
    static cilist io___1 = { 0, 6, 0, 0, 0 };
    static cilist io___2 = { 0, 5, 0, 0, 0 };
    static cilist io___4 = { 0, 6, 0, 0, 0 };
    static cilist io___5 = { 0, 5, 0, 0, 0 };
    static cilist io___11 = { 0, 21, 0, "(I8,F20.10)", 0 };


/* ccccccccccccccccccccccccccccccccccccccccccccccccc */
/*     Calculate The Boltzmann Distribution       C */
/* ccccccccccccccccccccccccccccccccccccccccccccccccc */
/* ccccccccccccccccccccccccccccccccccccccccc */
/*     Read Info From Disk                C */
/*     Maxn = Maximum Number Over Levels  C */
/* ccccccccccccccccccccccccccccccccccccccccc */
    s_wsle(&io___1);
    do_lio(&c__9, &c__1, "Number Of Energy Levels (2-10000) ? ", (ftnlen)36);
    e_wsle();
    s_rsle(&io___2);
    do_lio(&c__3, &c__1, (char *)&n, (ftnlen)sizeof(integer));
    e_rsle();
    s_wsle(&io___4);
    do_lio(&c__9, &c__1, "Temperature ? ", (ftnlen)14);
    e_wsle();
    s_rsle(&io___5);
    do_lio(&c__5, &c__1, (char *)&temp, (ftnlen)sizeof(doublereal));
    e_rsle();
    if (n < 2 || n >= 10000 || temp < 1e-7 || temp > 1e7) {
	s_stop("", (ftnlen)0);
    }
    beta = 1. / temp;
    norm = 0.;
/* ccccccccccccccccccccccccccccc */
/*     Loop Over All Levels   C */
/* ccccccccccccccccccccccccccccc */
    i__1 = n - 1;
    for (i__ = 0; i__ <= i__1; ++i__) {
	temp = exp(-beta * (doublereal) i__);
	if (temp < 1e-70) {
	    temp = 0.;
	}
	norm += temp;
	distri[i__] = temp;
    }
    norm = 1. / norm;
/* cccccccccccccccccccccc */
/*     Write Results   C */
/* cccccccccccccccccccccc */
    o__1.oerr = 0;
    o__1.ounit = 21;
    o__1.ofnmlen = 10;
    o__1.ofnm = "result.dat";
    o__1.orl = 0;
    o__1.osta = 0;
    o__1.oacc = 0;
    o__1.ofm = 0;
    o__1.oblnk = 0;
    f_open(&o__1);
    i__1 = n - 1;
    for (i__ = 0; i__ <= i__1; ++i__) {
	s_wsfe(&io___11);
	do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	d__1 = distri[i__] * norm;
	do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
	e_wsfe();
    }
    cl__1.cerr = 0;
    cl__1.cunit = 21;
    cl__1.csta = 0;
    f_clos(&cl__1);
    s_stop("", (ftnlen)0);
    return 0;
} /* MAIN__ */

/* Main program alias */ int boltzmann_ () { MAIN__ (); return 0; }
