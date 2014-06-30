#include "/home/vinicius/image/prog/include/image.h"


/* ------------------------------------------------------------------------- *
 *    FDCTIV2.C  -  Fast Discrete Cosine Transform, type IV                  *
 *                                                                           *
 *    Version 2: with scaling factor.                                        *
 *                                                                           *
 *    The computation is done via a half-length FFT, according               *
 *    to the algorithm described in Section 2.5.                             *
 *                                                                           *
 *    Author:  Henrique S. Malvar.                                           *
 *    Date:    October  8,  1991.                                            *
 *                                                                           *
 *    Usage:   fdctiv2(x, logm, fac);  --  direct or inverse transform       *
 *                                                                           *
 *    Arguments:  x (float)      -input and output vector, length M.         *
 *                logm  (int)   - log (base 2) of vector length M,           *
 *                                e.g., for M = 256 -> logm = 8.             *
 *                fac  (float)  - scaling factor.  If fac == 1, an            *
 *                                orthogonal transform is performed;         *
 *                                otherwise, outputs are multiplied          *
 *                                by fac; see note below.                    *
 *                                                                           *
 *    Note: the scaling factor argument fac is only used the first           *
 *    time the routine is called with a new value for logm.                  *
 *    Subsequent calls will use the previously defined value of fac.         *
 *                                                                           *
 *    This version of the fdctiv is called by the fdct module.               *
 * ------------------------------------------------------------------------- */


/* -------------------------------------------------------------------- *
 *      DCT-IV    transform                                             *
 * -------------------------------------------------------------------- */

void  fdctiv2(float *x, int logm, float fac)
{
    int        i, m, m2, m4, n, nel;
    static     float      *y, *xp1, *xp2, *yp1, *yp2;
    static     float      *c4n, *spc4n, *smc4n, *cn, *spcn, *smcn;
    static     float      tmp1, tmp2, ang, c, s;
    static     float      *tab[MAXLOGM];
    static     float      *yt[MAXLOGM];

    /* Check range of logm */
    if  ((logm < 0) || (logm > MAXLOGM)) {
        printf("Error : FDCTIV : logm = %d is out of bounds [%d, %d]\n",
           logm, 0, MAXLOGM);
        exit(1);
    }

    /*  Trivial  cases, m = 1 and m = 2 */
    if (logm < 1) {
        *x *= fac;
        return;
   }

    if (logm == 1) {
        if (tab[0] == NULL) {
           if ((tab[0] = (float *) calloc(2, sizeof(float))) == NULL) {
               printf("Error : FDCTIV2 : not enough memory for cosine"
                  " tables.\n");
               exit(1);
        }
           yp1 = tab[0];
           *yp1++ = fac * COSM2;
           *yp1   = fac * SINM2;
   }
        xp2 = x + 1;
        yp1 = tab[0];  yp2 = yp1 + 1;
        tmp1 = *yp1 * *x + *yp2 * *xp2;
        *xp2 = *yp2 * *x - *yp1 * *xp2;
        *x   = tmp1;
        return;
}

    /*  Compute m */
    m = 1 << logm;
    m2 = m / 2;
    m4 = m2 / 2;

 /* Build tables of butterfly coefficients, if necessary */
 if (tab[logm-1] == NULL) {

     /* Allocate memory for tables */
     nel = m2;
     if ((tab[logm-1] = (float *) calloc(6 * nel, sizeof(float)))
        == NULL) {
        printf("Error : FDCTIV2: not enough memory for"
            " cosine tables.\n");
        exit(1);
     }

     /* Initialize pointers */
     c4n = tab[logm-1]; spc4n = c4n + nel; smc4n = spc4n + nel;
     cn  = smc4n + nel; spcn  = cn  + nel; smcn  = spcn  + nel;

     /*  Compute  tables.  For  the  tables  in  the  second   modulation
        stage, all entries should be multiplied by sqrt(2/m), in
        order to generate an orthogonal transform */
     fac *= sqrt(2.0 / m);
     for (n = 0; n < m2; n++) {
        ang = (n + 0.25) * PI / m;
        c = cos(ang); s = sin(ang);
        *c4n++ = c; *spc4n++ = - (s  + c); *smc4n++ = s - c;
        if (n == m4) {
           c = fac * SQHALF;  *cn++ = c;
        } else {
           ang = n * PI / m;
           c = fac * cos(ang); s = fac * sin(ang);
           *cn++ = c; *spcn++ = - (s + c); *smcn++ = s - c;
         }
    }
}

/* Allocate space for working vector, if necessary */
if (yt[logm-1] == NULL) {
     if ((yt[logm-1] = (float *) calloc(m2, sizeof(float))) == NULL) {
        printf("Error : FDCTIV2 : not enough memory for"
           " working vector.\n");
        exit(1);
    }
}

/* Define table pointers */
nel = m2;
c4n = tab[logm-1]; spc4n = c4n + nel; smc4n = spc4n + nel;
cn  = smc4n + nel; spcn  = cn  + nel; smcn  = spcn  + nel;
y   = yt[logm-1];

/* Step 1: input data reordering */
xp1 = x; xp2 = x + 1;
yp1 = x; yp2 = y + m2 - 1;
for (i = 0; i < m2; i++) {
    *yp2-- = *xp2++;
    *yp1++ = *xp1++;
    xp1++; xp2++;
}
/*  Step  2:  first  modulation  stage */
yp1 = x; yp2 = y;
for (i = 0; i < m2; i++) {
    tmp2 = *c4n++ * (*yp1 + *yp2);
    tmp1 = *spc4n++ * *yp1 + tmp2;
    *yp1 = *smc4n++ * *yp2 + tmp2;
    *yp2 = tmp1;
    yp1++; yp2++;
}

/*  Step 3: call FFT of half length */
srfft(x, y, logm-1);

/* Step 4: second modulation stage */
yp1 = x; yp2 = y;
for (i = 0; i < m2; i++) {
    if (i == m4)  {
       tmp1 = *cn * (*yp1 + *yp2);
       *yp2 = *cn++ * (*yp2 - *yp1);
       *yp1 = tmp1;
    } else {
       tmp2 = *cn++ * (*yp1 + *yp2);
       tmp1 = *spcn++ * *yp1 + tmp2;
       *yp1 = *smcn++ * *yp2 + tmp2;
       *yp2 = tmp1;
     }
    yp1++; yp2++;
}

/* Step 5: output data reordering */
xp1 = x; xp2 = &x[m-1];
yp1 = y; yp2 = y + m2;

xp1 = x + m - 2;  xp2 = x + m -1;
yp1 = x + m2 - 1; yp2 = y;
for (i = 0; i < m2; i++) {
    *xp1-- =   *yp1--;
    *xp2-- = - *yp2++;
    xp1--; xp2--;
}
}
