
   /* -------------------------------------------------------------------- *
    *    RSFFT.C  -  Split-Radix Fast Fourier Transform                    *
    *                 for real-valued inputs.                              *
    *                                                                      *
    *    This is a decimation-in-frequency version, using the steps        *
    *    of the algorithm as described in Section 2.5.                     *
    *                                                                      *
    *    Author:  Henrique S. Malvar.                                      *
    *    Date:    October  8,   1991.                                      *
    *                                                                      *
    *    Usage:   rsfft(x, logm);         for direct DFT                   *
    *             rsifft(x, logm);        for inverse DFT                  *
    *                                                                      *
    *    Arguments:  x (float)  -   input and output vector, length M;     *
    *                               in the time domain, it contains the    *
    *                               M-point data vector; in the frequency  *
    *                               domain, it contains the real parts of  *
    *                               X(C) to X(M/2) followed by the         *
    *                               imaginary parts of X(M/2+1) to         *
    *                               X(M-1).                                *
    *                logm (int)     log (base 2) of vector length    M,    *
    *                               e.g., for M = 256 -> logm = 8.         *
    *--------------------------------------------------------------------  */

#include  <stdio.h>
#include  <stdlib.h>
#include  <alloc.h>
#include  <math.h>

#define   MAXLOGM      11     /* max FFT length      2^MAXLOGM */
#define   TWOPI        6.28318530717958647692
#define   SQHALF       0.707106781186547524401

/*  -------------------------------------------------------------------- *
 *    Prototypes  of   external   functions.                             *
 *  -------------------------------------------------------------------- */

void BR_permute(float *x, int logm);
void srrec(float *xr, float *xi, int logm);

/*  -------------------------------------------------------------------- *
 *    Error exit for program abortion.                                   *
 *  -------------------------------------------------------------------- */


 static    void  error_exit(void)
{
    exit (1);
}
/* -------------------------------------------------------------------- *
 *     Recursive part of the RSFFT algorithm.       Not externally      *
 *     callable.                                                        *
 * -------------------------------------------------------------------- */

 static void  rsrec(float *x, int logm)
{
     static   int       m, m2, m4, m8, nel, n;
     static   float     *xr1, *xr2, *xi1;
     static   float     *cn, *spcn, *smcn;
     static   float     tmp1, tmp2, ang, c, s;
     static   float     *tab[MAXLOGM];

     /* Check range   of logm */
     if ((logm < 0) || (logm > MAXLOGM)) {
        printf("Error : RSFFT : logm = %d is out of bounds [%d, %d]\n",
            logm, 0, MAXLOGM);
        error-exit();
}

     /* Compute trivial cases */
     if (logm < 2) {
          if (logm == 1) {    /* length m = 2 */
            xr2  = x + 1;
            tmp1 = *x + *xr2;
            *xr2 = *x - *xr2;
            *x   =  tmp1;
            return;
          }
        else if (logm == 0) return;      /* length m = 1 */
    }

     /* Compute a few constants */
     m = 1 << logm; m2 = m / 2; m4 = m2 / 2; m8 = m4 / 2;

     /* Build tables of butterfly coefficients, if necessary */
     if ((logm >= 4) && (tab[logm-4] == NULL)) {

         /* Allocate memory for tables */
         nel = m4 - 2;

    if ((tab[logm-4] = (float *) calloc(3 * nel, sizeof(float)))
       == NULL) {
       printf("Error : RSFFT : not enough memory for cosine tables.\n");
       error_exit();
    }

    /* Initialize pointers */
    cn  = tab[logm-4]; spcn = cn + nel;  smcn = spcn + nel;

    /* Compute tables */
    for (n = 1; n < m4; n++) {
       if (n == m8) continue;
       ang = n * TWOPI / m;
       c = cos(ang);  s = sin(ang);
       *cn++ = c;  *spcn++ = - (s + c); *smcn++ = s - c;
   }
}

/*  Step  1 */
xr1 = x;  xr2 = xr1 + m2;
for (n = 0; n < m2; n++) {
    tmp1 = *xr1 + *xr2;
    *xr2 = *xr1 - *xr2;
    *xr1 = tmp1;
    xr1++; xr2++;
}

/*  Step  2        */
xr1 = x + m2 + m4;
for (n = 0; n < m4; n++) {
    *xr1 = - *xr1;
    xr1++;
}

/*  Steps 3 &  4 */
xr1 = x + m2; xi1 = xr1 + m4;
if (logm >= 4) {
    nel = m4 - 2;
    cn  = tab[logm-4]; spcn = cn + nel;  smcn = spcn + nel;
}

xr1++; xi1++;
for (n = 1; n < m4; n++) {
    if (n == m8) {
       tmp1 =  SQHALF * (*xr1 + *xi1);
       *xi1  = SQHALF * (*xi1 - *xr1);
       *xr1  = tmp1;
    }  else {
       tmp2 = *cn++ * (*xr1 + *xi1);
       tmp1 = *spcn++ * *xr1 + tmp2;
       *xr1 = *smcn++ * *xi1 + tmp2;
       *xi1 = tmp1;
       }
    xr1++; xi1++;
}

    /*  Call rsrec again with half DFT length */
    rsrec(x, logm-1);

    /* Call complex DFT routine, with quarter DFT length.
        Constants have to be recomputed, because they are static! */
    m = 1 << logm; m2 = m / 2; m4 = 3 * (m / 4);
    srrec(x + m2, x + m4, logm-2);

    /* Step 5: sign change & data reordering */
    m = 1 << logm; m2 = m / 2; m4 = m2 / 2; m8 = m4 / 2;
    xr1 = x + m2 + m4;
    xr2 = x + m - 1;
    for (n = 0; n < m8; n++) {
       tmp1   = *xr1;
       *xr1++ = - *xr2;
       *xr2-- = - tmp1;
    }
    xr1 = x + m2 + 1;
    xr2 = x + m - 2;
    for (n = 0; n < m8; n++) {
        tmp1   =   *xr1;
        *xrl++ = - *xr2;
        *xr2-- =   tmp1;
        xr1++;
        xr2--;
    }
    if (logm == 2) x[3] = -x[3];
}
/* -------------------------------------------------------------------- *
 *      Direct transform for real inputs                                *
 * -------------------------------------------------------------------- */

void  rsfft(float *x, int logm)
{
    /* Call recursive routine */
    rsrec(x, logm);

    /* Output array unshuffling using bit-reversed indices */
    if (logm > 1) {
        BR_permute(x, logm);
    }
}


/* --------------------------------------------------------------------- *
 *  Recursive part of the inverse RSFFT algorithm.  Not externally       *
 *  callable.                                                            *
 *  -------------------------------------------------------------------- */

 static  void  rsirec(float  *x,  int   logm)
{
     static    int       m, m2, m4, m8, nel, n;
     static    float     *xr1, *xr2, *xi1;
     static    float     *cn, *spcn, *smcn;
     static    float     tmp1, tmp2, ang, c, s;
     static    float     *tab [MAXLOGM];

     /* Check  range  of logm */
     if ((logm < 0) || (logm > MAXLOGM)) {
        printf("Error : RSFFT : logm = %d is out of bounds [%d, %d]\n",
           logm, 0, MAXLOGM);
        error_exit();
     }

     /*  Compute  trivial  cases */
     if (logm < 2) {
        if (logm == 1) {     /* length m = 2 */
           xr2  = x + 1;
           tmp1 = *x + *xr2;
           *xr2 = *x - *xr2;
           *x   = tmp1;
           return;
        }
       else if (logm == 0) return;       /* length m = 1 */
    }

     /* Compute a few constants */
     m = 1 << logm; m2 = m / 2; m4 = m2 / 2; m8 = m4 / 2;

     /* Build tables of butterfly    coefficients, if necessary */
     if ((logm >= 4) && (tab[logm-4] == NULL)) {

       /* Allocate memory for tables */
       nel = m4 - 2;
       if ((tab[logm-4] = (float *) calloc(3 * nel, sizeof(float)))
           == NULL) {
           printf("Error : RSFFT : not enough memory for cosine tables.\n");
           error_exit();
       }

       /*  Initialize   pointers */
       cn  = tab[logm-4] ; spcn = cn + nel; smcn = spcn + nel;

     /*  Compute  tables */
     for (n = 1; n < m4; n++) {
         if (n == m8) continue;
         ang = n * TWOPI / m;
         c = cos(ang); s = sin(ang);
         *cn++ = c; *spcn++ = - (s + c); *smcn++ = s - c;
    }
}
 /* Reverse Step 5: sign change & data reordering */
 if (logm == 2) x[3] = -x[3];
 xr1 = x + m2 + 1;
 xr2 = x + m - 2;
 for (n = 0; n < m8; n++) {
     tmp1   =   *xr1;
     *xr1++ =   *xr2;
     *xr2-- = - tmp1;
     xr1++;
     xr2--;
 }
 xr1 = x + m2 + m4;
 xr2 = x + m - 1;
 for (n = 0; n < m8; n++) {
     tmp1   =   *xr1;
     *xr1++ = - *xr2;
     *xr2-- = - tmp1;
}
 /*  Call   rsirec again with half DFT length */
 rsirec(x, logm-1);

 /* Call complex DFT routine, with quarter DFT length.
     Constants have to be recomputed, because they are static! */
 m = 1 << logm; m2 = m / 2; m4 = 3 * (m / 4);
 srrec(x + m4, x + m2, logm-2);

 /* Reverse Steps 3 & 4 */
 m = 1 << logm; m2 = m / 2; m4 = m2 / 2; m8 = m4 / 2;
 xr1 = x + m2; xi1 = xr1 + m4;
 if (logm >= 4) {
     nel = m4 - 2;
     cn  = tab[logm-4]; spcn = cn + nel; smcn = spcn + nel;
 }
 xr1++; xi1++;
 for (n = 1; n < m4; n++) {
     if (n == m8) {
         tmp1 = SQHALF * (*xr1 - *xi1);
         *xi1 = SQHALF * (*xi1 + *xr1);
         *xr1 = tmp1;
     }    else {
            tmp2 = *cn++ * (*xr1 + *xi1);
            tmp1 = *smcn++ * *xr1 + tmp2;
            *xr1 = *spcn++ * *xi1 + tmp2;
            *xi1 = tmp1;
         }
        xr1++; xi1++;
 }

 /* Reverse Step 2 */
 xr1 = x + m2 + m4;
 for (n = 0; n < m4; n++) {
     *xr1 = - *xr1;
     xr1++;
 }

 /* Reverse  Step  1 */
 xr1 = x; xr2 = xr1 + m2;
 for (n = 0; n < m2; n++) {
     tmp1 = *xr1 + *xr2;
     *xr2 = *xr1 - *xr2;
     *xr1 = tmp1;
     xr1++; xr2++;
 }
}

  /* -------------------------------------------------------------------- *
   *   Inverse  transform  for  real  inputs                              *
   *--------------------------------------------------------------------  */

void  rsifft(float *x, int logm)
{
   int       i, m;
   float     fac, *xp;

   /* Output array unshuffling using bit-reversed indices */
   if (logm > 1) {
      BR_permute(x, logm);
   }
   x[0l *= 0.5;
   if (logm > 0) x[1] *= 0.5;

   /*  Call  recursive  routine */
   rsirec(x, logm);

   /* Normalization */
   m = 1 << logm;
   fac = 2.0 / m;
   xp = x;

   for (i = 0; i < m; i++) {
       *xp++ *= fac;
   }
}
