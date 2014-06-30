#include "/home/vinicius/image/prog/include/image.h"


/* -------------------------------------------------------------------- *
 *   FDCT.C  -  Fast Discrete Cosine Transform.                         *
 *                                                                      *
 *   The computation is done via a half-length DCT and a half-length    *
 *   DCT-IV, according to the algorithm described in Section 2.5.       *
 *                                                                      *
 *   In order to minimize the number of multiplications, the            *
 *   transform is not orthogonal. All outputs are multiplied by         *
 *   sqrt(M).                                                           *
 *                                                                      *
 *   Author:    Henrique S. Malvar.                                     *
 *   Date:      October 8, 1991.                                        *
 *                                                                      *
 *   Usage:     fdct(x, logm);    --  for direct transform              *
 *              fidct(x, logm);   --  for inverse transform             *
 *                                                                      *
 *   Arguments:   x (float)  -     input and output vector, length M.   *
 *                  logm (int) -   log (base 2) of vector length M,     *
 *                                 e.g., for M = 256 -> logm = 8.       *
 *--------------------------------------------------------------------- */

/* ------------------------------------------------------------------- *
 *   Pointers for working vectors.                                     *
 *-------------------------------------------------------------------- */

static float *yt[MAXLOGM];


/* ------------------------------------------------------------------- *
 *    Recursive part of the DCT algorithm.  Not externally             *
 *    callable                                                         *
 *-------------------------------------------------------------------- */

static void dctrec(float *x, int logm)
{
     static    float   tmp, *x1, *x2;
     static    int     n, m, m2, m4;

     /* Stop recursion when m = 2 */
     if (logm == 1) {
        x2 = x + 1;
        tmp = (*x + *x2);
        *x2 = (*x - *x2);
        *x  = tmp;
        return;
     }

     m = 1 << logm;
     m2 = m / 2;
     m4 = m2 / 2;

 /* +/- butterflies (see Fig. 2.18) */
     x1 = x; x2 = x1 + m - 1;
     for (n = 0; n < m2; n++) {
        tmp = *x1 + *x2;
        *x2 = *x1 - *x2;
        *x1 = tmp;
        x1++;  x2--;
     }
     /* Swap entries of bottom half vector */
     x1 = x + m2; x2 = x + m - 1;
     for (n = 0; n < m4; n++) {
        tmp = *x2;
        *x2 = *x1;
        *x1 = tmp;
        x1++;  x2--;
     }

     /*   DCT-IV on the second half of x */
     fdctiv2(x + m2, logm-1, sqrt(m2));


     /*  DCT  on the first half of x */
     dctrec(x, logm-1);
}


/* -------------------------------------------------------------------- *
 *    Direct transform                                                  *
 * -------------------------------------------------------------------- */

void  fdct(float *x, int logm)
{
    static  float     tmp, *xp, *y, *yp;
    static  int       gr, n, m, m2, dx;

    /* Check range of logm */
    if ((logm < 0) || (logm > MAXLOGM)) {
      printf("Error : FDCT : logm   %d is out of bounds [%d, %d]\n",
          logm, 0, MAXLOGM);
      exit(1);
    }
    /* Trivial cases, m = 1 and m = 2 */
    if (logm < 1) return;

    if (logm == 1) {
      xp  = x + 1;
      tmp = (*x + *xp);
      *xp = (*x - *xp);
      *x  = tmp;
      return;
    }
    m = 1 << logm;
    m2 = m / 2;

    /* Allocate space for working vector, if necessary */
    if (yt[logm-1] == NULL) {
      if ((yt[logm-1] = (float *) calloc(m, sizeof(float))) == NULL) {
          printf("Error : FDCT : not enough memory for working vector.\n");
          exit(1);
      }
    }
    y = yt[logm-1];

    /* Copy x into y */
    memcpy(y, x, sizeof(float) * m);

    /* Call recursive module */
    dctrec(y, logm);

    /* Copy y back into x, with  data unshuffling */
    dx = 2;

    for (gr = 0; gr < logm; gr++) {
       xp = x +  (1 << gr); yp = y + m2;
       for (n = 0; n < m2; n++) {
          *xp = *yp++;
          xp += dx;
       }
       m2 >>= 1;
       dx <<= 1;
    }
    x[0]   = y[0];
    x[m/2] = y[1];
}

