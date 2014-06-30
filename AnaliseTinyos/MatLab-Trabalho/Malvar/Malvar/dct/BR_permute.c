#include "/home/vinicius/image/prog/include/image.h"


/* -------------------------------------------------------------------- *
 *    Data unshuffling according to bit-reversed indexing.              *
 *                                                                      *
 *                                                                      *
 *    Bit reversal is done using Evan's algorithm (Ref: D. M. W.        *
 *    Evans, "An improved digit-reversal permutation algorithm...",     *
 *    IEEE Trans.  ASSP, Aug. 1987, pp. 1120-1125).                     *
 * -------------------------------------------------------------------- */

static    int   brseed[256];     /* Evans' seed table */
static    int     brsflg;         /* flag for table building */


void BR_permute(float *x, int logm)
{
   int       i, j, imax, lg2, n;
   int       off, fj, gno, *brp;
   float     tmp, *xp, *xq;

   lg2 =  logm >> 1;
   n = 1  << lg2;
   if  (logm & 1) lg2++;

   /*  Create seed table if not yet built */
   if  (brsflg != logm) {
       brsflg = logm;
       brseed[0] = 0;
       brseed[1] = 1;
       for  (j=2; j <= lg2; j++) {
          imax = 1 << (j-1);
          for (i = 0; i < imax; i++) {
             brseed[i] <<= 1;
             brseed[i + imax] = brseed[i] + 1;
          }
       }
   }

   /*  Unshuffling   loop */
   for (off = 1; off < n; off++) {
       fj = n * brseed[off]; i = off; j = fj;
       tmp = x[i]; x[i] = x[j]; x[j] = tmp;
       xp = &x[i];
       brp = &brseed[1];
       for (gno = 1; gno < brseed[off]; gno++) {
          xp += n;
          j = fj + *brp++;
          xq = x + j;
          tmp = *xp; *xp = *xq; *xq = tmp;
       }
   }
}
