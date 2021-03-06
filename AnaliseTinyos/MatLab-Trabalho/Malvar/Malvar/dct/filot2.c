#include "/home/vinicius/image/prog/include/image.h"


/* -------------------------------------------------------------------- *
 *     FLOT2.C  -  Fast Lapped Orthogonal Transform                     *
 *                                                                      *
 *     This is the type-II LOT, as described in Section 4.4.            *
 *                                                                      *
 *     Version 2: with buffer memory as an argument.                    *
 *                                                                      *
 *     Author:  Henrique S. Malvar.                                     *
 *     Date:    October  8,  1991.                                      *
 *                                                                      *
 *     Usage:    flot2(x,  y,  logm);  --  direct transform             *
 *               filot2(x, y, logm);   --  inverse transform            *
 *                                                                      *
 *     Arguments:  x (float)      - input and output vector, length M.  *
 *                 y (float)      - pointer to internal buffers; it     *
 *                                  should point to a previously        *
 *                                  allocated buffer with room for      *
 *                                  3 * M / 2  values.                  *
 *                  logm  (int)   - log (base 2) of vector length M,    *
 *                                  e.g., for M = 256 -> logm = 8.      *
 * -------------------------------------------------------------------- */

/* -------------------------------------------------------------------- *
 *  Inverse    LOT                                                      *
 * -------------------------------------------------------------------- */

void filot2(float *x, float *y, int logm)
{

     static      int    n, m, m2, m4;
     static      float  tmp, fac, *xp1, *xp2, *yp1, *yp2;

     /* Check range of logm */
     if ((logm < 2)  || (logm > MAXLOGM)) {
       printf("Error : FLOT : logm = %d is out of bounds [%d, %d]\n",
            logm, 2, MAXLOGM);
       error_exit();
     }

     /* Define m */
     m  = 1 << logm;
     m2 = m / 2;
     m4 = m2 / 2;

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
      output in  x */
   xp1 = x;
   xp2 = x + m2;
   yp1 = y;
   yp2 = y + m2;
   for (n = 0; n < m2; n++) {
      *(xp1++) = 0.5 * ( *yp1     + *yp2);
      *(xp2++) = 0.5 * ( *(yp1++) - *(yp2++) );
   }

   /* This piece of code correspond to the  z^(-1)  delays indicated in
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
      xp1++; xp2++;
   }
   /* Copy  y[0]  on  even-indexed  x's  and
      y[m/2]    on odd-indexed  x's */
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
