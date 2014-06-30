#include "/home/vinicius/image/prog/include/image.h"

/* -------------------------------------------------------------------- *
 *      Recursive part of the SRFFT algorithm.                          *
 * -------------------------------------------------------------------- */

void srrec(float *xr, float *xi, int logm)
{
   static    int        m, m2, m4, m8, nel, n;
   static    float      *xr1, *xr2, *xi1, *xi2;
   static    float      *cn, *spcn, *smcn, *c3n, *spc3n, *smc3n;
   static    float      tmp1, tmp2, ang, c, s;
   static    float      *tab[MAXLOGM];




/* Check range of logm */
if ((logm < 0) || (logm > MAXLOGM)) {
   printf("Error : SRFFT : logm = %d is out of bounds [%d, %d]\n",
       logm, 0, MAXLOGM);
   exit(1);
}

/*  Compute trivial cases */
if (logm < 3) {
      if (logm == 2) {  /* length m = 4 */
       xr2  = xr + 2;
       xi2  = xi + 2;
       tmp1 = *xr + *xr2;
       *xr2 = *xr - *xr2;
       *xr  = tmp1;
       tmp1 = *xi + *xi2;
       *xi2 = *xi - *xi2;
       *xi  = tmp1;
       xr1  = xr + 1;
       xi1  = xi + 1;
       xr2++;
       xi2++;
       tmp1 = *xr1 + *xr2;
       *xr2 = *xr1 - *xr2;
       *xr1 = tmp1;
       tmp1 = *xi1 + *xi2;
       *xi2 = *xi1 - *xi2;
       *xi1 = tmp1;
       xr2  = xr + 1;
       xi2  = xi + 1;
       tmp1 = *xr + *xr2;
       *xr2 = *xr - *xr2;
       *xr  = tmp1;
       tmp1 = *xi + *xi2;
       *xi2 = *xi - *xi2;
       *xi  = tmp1;
       xr1  = xr + 2;
       xi1  = xi + 2;
       xr2  = xr + 3;
       xi2  = xi + 3;
       tmp1 = *xr1 + *xi2;
       tmp2 = *xi1 + *xr2;
       *xi1 = *xi1 - *xr2;
       *xr2 = *xr1 - *xi2;
       *xr1 =tmp1;
       *xi2 =tmp2;
       return;
}

    else  if (logm == 1) { /* length m = 2 */
       xr2   = xr +  1;
       xi2   = xi +  1;
       tmp1  = *xr + *xr2;
       *xr2  = *xr - *xr2;
       *xr   = tmp1;
       tmp1  = *xi + *xi2;
       *xi2  = *xi - *xi2;
       *xi   = tmp1;
       return;
    }
    else if (logm == 0) return;     /* length m = 1*/
}

/* Compute a few constants */
m = 1 << logm; m2 = m / 2; m4 = m2 / 2; m8 = m4 / 2;

/* Build tables of butterfly coefficients, if necessary */
if ((logm >= 4) && (tab[logm-4] == NULL)) {

    /* Allocate memory for tables */
    nel = m4 - 2;
    if ((tab[logm-4] = (float *) calloc(6 * nel, sizeof(float)))
       == NULL) {
       exit(1);
     }
    /* Initialize pointers */

    cn  = tab[logm-4]; spcn = cn + nel;  smcn = spcn + nel;
    c3n = smcn + nel; spc3n = c3n + nel; smc3n = spc3n + nel;


    /* Compute tables */
    for (n = 1; n < m4; n++) {
       if (n == m8) continue;
       ang = n * TWOPI / m;
       c = cos(ang); s = sin(ang);
       *cn++ = c; *spcn++ = - (s + c); *smcn++ = s - c;
       ang = 3 * n * TWOPI / m;
       c = cos(ang); s = sin(ang);
       *c3n++ = c; *spc3n++ = - (s + c); *smc3n++ = s - c;
   }
}


/*  Step 1 */
xr1 = xr;  xr2 = xr1  +  m2;
xi1 = xi;  xi2 = xi1  +  m2;

for (n = 0; n < m2; n++) {
   tmp1 = *xr1 + *xr2;
   *xr2 = *xr1 - *xr2;
   *xr1 = tmp1;
   tmp2 = *xi1 + *xi2;
   *xi2 = *xi1 - *xi2;
   *xi1 = tmp2;
   xr1++;  xr2++;  xi1++;  xi2++;
}
/*   Step 2  */
xr1 = xr + m2; xr2 = xr1 + m4;
xi1 = xi + m2; xi2 = xi1 + m4;
for (n = 0; n < m4; n++) {
   tmp1 = *xr1 + *xi2;
   tmp2 = *xi1 + *xr2;
   *xi1 = *xi1 - *xr2;
   *xr2 = *xr1 - *xi2;
   *xr1 = tmp1;
   *xi2 = tmp2;
   xr1++;  xr2++;  xi1++;  xi2++;
}

/*   Steps  3 & 4 */
xr1 = xr + m2; xr2 = xr1 + m4;
xi1 = xi + m2; xi2 = xi1 + m4;
if (logm >= 4) {
   nel = m4 - 2;
   cn  = tab[logm-4]; spcn  = cn + nel;  smcn  = spcn + nel;
   c3n = smcn + nel;  spc3n = c3n + nel; smc3n = spc3n + nel;
}
xr1++; xr2++; xi1++; xi2++;
for (n = 1; n < m4; n++) {
   if (n == m8) {
       tmp1 = SQHALF * (*xr1 + *xi1);
       *xi1 = SQHALF * (*xi1 - *xr1);
       *xr1 = tmp1;
       tmp2 = SQHALF * (*xi2 - *xr2);
       *xi2 = -SQHALF * (*xr2 + *xi2);
       *xr2 = tmp2;
}     else {
       tmp2 = *cn++ * (*xr1 + *xi1);
       tmp1 = *spcn++ * *xr1 + tmp2;
       *xr1 = *smcn++ * *xi1 + tmp2;
       *xi1 = tmp1;
       tmp2 = *c3n++ * (*xr2 + *xi2);
       tmp1 = *spc3n++ * *xr2 + tmp2;
       *xr2 = *smc3n++ * *xi2 + tmp2;
       *xi2 = tmp1;
}
     xr1++; xr2++; xi1++; xi2++;
}
   /* Call ssrec again with half DFT length  */
   srrec(xr, xi, logm-1);

   /* Call ssrec again twice with one quarter DFT length.
     Constants have to be recomputed, because they are static!*/
   m = 1 << logm; m2 = m / 2;
   srrec(xr + m2, xi + m2, logm-2);
   m = 1 << logm; m4 = 3 * (m / 4);
   srrec(xr + m4, xi + m4, logm-2);
}
