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
 *    Direct transform                                                  *
 * -------------------------------------------------------------------- */

void srfft(float *xr, float *xi, int logm)
{
   /* Call recursive routine */
   srrec(xr, xi, logm);

   /* output array unshuffling using bit-reversed indices */
   if (logm > 1) {
     BR_permute(xr, logm);
     BR_permute(xi, logm);
   }
}
