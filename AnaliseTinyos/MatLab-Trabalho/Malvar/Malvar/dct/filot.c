#include "/home/vinicius/image/prog/include/image.h"


/* ---------------------------------------------------------------------- *
 *    FLOT.C  -  Fast Lapped Orthogonal Transform                         *
 *                                                                        *
 *    This is the type-II LOT, as described in Section 4.4.               *
 *                                                                        *
 *    In order to minimize the number of multiplications, the             *
 *    transform is not orthogonal.     The basis functions have their     *
 *    energies equal to  1 / M.                                           *
 *                                                                        *
 *    Author:   Henrique S. Malvar.                                       *
 *    Date:     October 8, 1991.                                          *
 *                                                                        *
 *    Usage:    flot(x, logm);   --  direct transform                     *
 *              filot(x, logm);  --  inverse transform                    *
 *                                                                        *
 *    Ar uments:  x (float)     - input and output vector, length M.      *
 *                 logm (int)  -  log (base 2) of vector length M,        *
 *                                e.g., for M = 256 -> logm = 8.          *
 *                                                                        *
 *    Note: this  is a causal version.  Thus, there is a delay of         *
 *    half a block in flot or filot. Following a call to flot with        *
 *    a call to filot will recover the original signal with one           *
 *    block delay.                                                        *
 * ---------------------------------------------------------------------- */


/* ------------------------------------------------------------------- *
 *     Inverse LOT                                                     *
 * ------------------------------------------------------------------- */

void filot(float *x, int logm)
{
    static   int     n, m, m2, m4;
    static   float   tmp, fac, *xp1, *xp2, *yp1, *yp2, *y;
    static   float   *yt[MAXLOGM];

    /* Check range of logm */
    if ((logm < 2) || (logm > MAXLOGM)) {
       printf("Error : FLOT : logm = %d is out of bounds [%d, %d]\n",
          logm, 2, MAXLOGM);
       error_exit();
    }

    /* Define m */
    m  = 1 << logm;
    m2 = m / 2;
    m4 = m2 / 2;

    /* Allocate space for working vector, if necessary */
    if (yt[logm-2] == NULL) {
       if ((yt[logm-2] = (float *) calloc(3 * m2, sizeof(float))) ==   NULL) {
          printf("Error : FLOT : not enough memory for working vector.\n");
          error_exit();
       }
    }
    y = yt[logm-2];

    /* Even/odd re-indexing, output in  y */
    xp1 = x;
    yp1 = y;
    yp2 = y + m2;
    for (n = 0; n < m2; n++) {
       *(yp1++) = *(xp1++);
       *(yp2++) = *(xp1++);
    }
    /* Length-(m/2)  DST-IV */
    fac = sqrt(2.0 / m);
    yp1 = y + m2; yp2 = y + m - 1;

   for (n = 0; n < m4; n++) {
      tmp  = *yp2;
      *yp2 = fac * *yp1;
      *yp1 = fac * tmp;
      yp1++; yp2--;
   }
   fdctiv(y + m2, logm-1);
   yp1 = y + m2 + 1;
   for (n = 0; n < m4; n++) {
      *yp1 = - *yp1;
      yp1++; yp1++;
   }

   /* Length-(n/2)  DCT */
   fdct(y + m2, logm-1);

   /* First butterflies with  +1/-1  factors, with  1/2  factor,
      output in x */
   xp1 = x;
   xp2 = x + m2;
   yp1 = y;
   yp2 = y + m2;
   for (n = 0; n < m2; n++) {
      *(xp1++) = 0.5 * ( *yp1     + *yp2 );
      *(xp2++) = 0.5 * ( *(yp1++) - *(yp2++) );
   }


   /* This  piece of code correspond to the  z^(-l)  delays indicated in
      Section 4.4. The stored values are in the last  m/2  samples of y */
   memcpy(y, y + m, m2 * sizeof(float));
   memcpy(y + m, x, m2 * sizeof(float));

   /* Copy first  m/2  coefficients of  y  in  x[0] */
   memcpy(x, y, m2 * sizeof(float));

   /* Second stage of  +1/-1  butterflies, output in  y */
   xp1 = x;
   xp2 = x + m2;
   yp1 = y;
   yp2 = y + m2;
   for (n = 0; n < m2; n++) {
     *(yp1++) = *xp1 + *xp2;
     *(yp2++) = *xp1 - *xp2;
     xp1++;   xp2++;
   }

   /* Copy  y[0]  on  even-indexed  x's  and
   y[m/2]  on odd-indexed  x's */
   xp1 = x;
   yp1 = y;
   yp2 = y + m2;
   for (n = 0; n < m2; n++) {
      *(xp1++) = *(yp1++);
      *(xp1++) = *(yp2++);
   }

    /*  Compute  IDCT */
    fidct(x, logm);
}



