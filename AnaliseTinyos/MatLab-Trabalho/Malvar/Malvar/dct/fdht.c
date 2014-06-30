#include "/home/vinicius/image/prog/include/image.h"


/* ----------------------------------------------------------------------- *
 *     FDHT.C  -  Split-Radix Fast Hartley Transform.                      *
 *                                                                         *
 *     This is a decimation-in-frequency version, using the steps          *
 *     of the algorithm as described in Section 2.5.                       *
 *                                                                         *
 *     Author:  Henrique S. Malvar.                                        *
 *     Date:     October 8,  1991.                                         *
 *                                                                         *
 *     Usage:    fdht(x, logm);    --  for direct or inverse DHT           *
 *                                                                         *
 *     Arguments:  x (float)  -   input and output vector, length M;       *
 *                  logm (int) -  log (base 2) of vector length M,         *
 *                                e.g., for M = 256 -> logm = 8.           *
 * ----------------------------------------------------------------------- */

/* -------------------------------------------------------------------- *
 *     Recursive part of the FHT algorithm.   Not externally            *
 *     callable.                                                        *
 * -------------------------------------------------------------------- */

static void dhtrec(float *x, int logm)
{
    static     int        m, m2, m4, m8, nel, n;
    static     float      *x1, *x2, *x3, *x4;
    static     float      *sn, *cpsn, *cmsn, *s3n, *cps3n, *cms3n;
    static     float      tmp1, tmp2, tmp3, ang, c, s;
    static     float      *tab [MAXLOGM];

/*  Check  range of  logm */
if  ((logm < 0) || (logm > MAXLOGM)) {
    printf("Error : FDHT : logm = %d is out of bounds [%d, %d]\n",
       logm, 0, MAXLOGM);
    exit(1);
}
/* Compute trivial cases */
if (logm <= 2) {
    if (logm == 2) {     /* length m = 4 */
       x1 = x; x2 = x + 2;
       tmp1 = *x1 + *x2;
       *x2  = *x1 - *x2;
       *x1  = tmp1;
       x1++; x2++;
       tmp1 = *x1 + *x2;
       *x2  = *x1 - *x2;
       *x1  = tmp1;
       x1 = x;  x2 = x + 1;
       tmp1 = *x1 + *x2;
       *x2  = *x1 - *x2;
       *x1  = tmp1;
       x1++;   x2++;
       x1++;   x2++;
       tmp1 = *x1 + *x2;
       *x2  = *x1 - *x2;
       *x1  = tmp1;
       return;
   }   else   if (logm == 1) {  /* length m = 2 */
       x2 = x + 1;
       tmp1 = *x + *x2;
       *x2  = *x - *x2;
       *x   = tmp1;
       return;
             }
    else if (logm == 0) return;   /* length m == 1 */
}

/*  Compute a few constants */
m = 1 << logm; m2 = m / 2; m4 = m2 / 2; m8 = m4 / 2;

/* Build tables of butterfly coefficients, if necessary */
if  ((logm >= 4) && (tab[logm-4] == NULL)) {

    /* Allocate memory for tables */
    nel = m8 - 1;
    if ((tab[logm-4] = (float *) calloc(6 * nel, sizeof(float)))
        == NULL) {
        printf("Error : FDHT : not enough memory for cosine tables.\n");
        exit(1);
    }
    /* Initialize pointers */
    sn  = tab[logm-4]; cpsn = sn + nel;  cmsn = cpsn + nel;
    s3n = cmsn + nel;  cps3n = s3n + nel; cms3n = cps3n + nel;

    /* Compute tables */
    for (n = 1; n < m8; n++) {
        ang = n * TWOPI / m;
        c = cos(ang); s = sin(ang);
        *sn++ = s; *cpsn++ = - (c  + s); *cmsn++ = c - s;
        ang = 3 * n * TWOPI / m;
        c = cos(ang); s = sin(ang);
        *s3n++ = s; *cps3n++ = - (c + s); *cms3n++ = c - s;
    }
}

/*  step  1 */
x1 = x; x2 = x1 + m2;
for (n = 0; n < m2; n++) {
    tmp1 = *x1 + *x2;
    *x2  = *x1 - *x2;
    *x1  = tmp1;
    x1++; x2++;
}

/*  Step  2 */
x1 = x + m2; x2 = x1 + m4;
for (n = 0; n < m8; n++) {
    tmp1 = *x1 + *x2;
    *x2  = *x1 - *x2;
    *x1  = tmp1;
    x1++; x2--;
}

/*  Step  3 */
x1 = x + m2 + m4 + 1; x2 = x + m - 1;
for (n = 1; n < m8; n++) {
    tmp1 = *x2 + *x1;
    *x2  = *x2 - *x1;
    *x1  =  tmp1;
    x1++; x2--;
}

    /* Step 4 */
    x[m2 + m8]      *= SQ2;
    x[m2 + m4 + m8] *= SQ2;

    /* Step   5*/
    if (logm >= 4) {
       nel = m8 - 1;
       sn  = tab[logm-4]; cpsn  = sn + nel;  cmsn  = cpsn + nel;
       s3n = cmsn + nel;  cps3n = s3n + nel; cms3n = cps3n + nel;
    }
    x1 = x + m2; x2 = x1 + m4; x3 = x2; x4 = x + m;
    x1++; x2--; x3++; x4--;
    for (n = 1; n < m8; n++) {
       tmp2 = *sn++ * (*x1 + *x4);
       *x1  = *cmsn++ * *x1 + tmp2;
       tmp3 = *cpsn++ * *x4 + tmp2;
       tmp2 = *s3n++ * (*x2 + *x3);
       *x4  = *cps3n++ * *x3 + tmp2;
       *x3  = *cms3n++ * *x2 + tmp2;
       *x2  = tmp3;
       x1++; x2--; x3++; x4--;
    }

    /* Call dhtrec again with half length */
    dhtrec(x, logm-1);

    /* Call dhtrec again twice with one quarter length.
       Constants have to be recomputed, because they are static! */
    m = 1 << logm; m2 = m / 2;
    dhtrec(x + m2, logm-2);
    m = 1 << logm; m4 = 3 * (m / 4);
    dhtrec(x + m4, logm-2);
}



/* -------------------------------------------------------------------- *
 *     Externally callable module.                                      *
 * -------------------------------------------------------------------- */

void  fdht(float *x, int logm)
{
    int       i, m;
    float     fac, *xp;

    /* Call recursive routine */
    dhtrec(x, logm);

    /* output array unshuffling using bit-reversed indices */
    if (logm > 1) {
       BR_permute(x, logm);
    }
    /* Normalization */
    m = 1 << logm;
    fac = sqrt(1.0 / m);
    xp = x;
    for (i = 0; i < m; i++) {
        *xp++ *= fac;
    }
}
