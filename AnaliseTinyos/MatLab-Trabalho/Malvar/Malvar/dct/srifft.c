#include "/home/vinicius/image/prog/include/image.h"


 /* ------------------------------------------------------------------- *
  *      SRFFT.C  -  Split-Radix Fast Fourier Transform                 *
  *                                                                     *
  * This is a decimation-in-frequency version, using the steps          *
  * of the algorithm as described in Section 2.5.                       *
  *                                                                     *
  * Author:  Henrique S. Malvar.                                        *
  * Date:     October 8, 1991.                                          *
  *                                                                     *
  * Usage:    srfft(xr, xi, logm);   --  for direct DFT                 *
  *           srifft(xr, xi, logm);  --  for inverse DFT                *
  *                                                                     *
  * Arguments:  xr (float)  - input and       output vector, length M,  *
  *                           real    part.                             *
  *             xi (float)  - imaginary part.                           *
  *             logm (int)  - log (base   2) of vector length M,        *
  *                           e.g., for   M = 256 -> logm = 8.          *
  * ------------------------------------------------------------------- */


/* -------------------------------------------------------------------- *
 *   Inverse transform.  Uses Duhamel's trick (Ref: P. Duhamel          *
 *   et. al., "On computing the inverse DFT", IEEE Trans- ASSP,         *
 *   Feb. 1988, pp. 285-286).                                           *
 * -------------------------------------------------------------------- */

void srifft(float *xr, float *xi, int logm)
{
   int   i, m;
   float fac, *xrp, *xip;

   /* Call direct FFT, swapping real & imaginary addresses */
   srfft(xi, xr, logm);

   /* Normalization */
   m = 1 << logm;
   fac = 1.0 / m;
   xrp = xr; xip = xi;

     for (i = 0; i < m; i++) {
         *xrp++ *= fac;
         *xip++ *= fac;
    }
}
