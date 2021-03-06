 
#include "Timer.h"
#include "RadioDctPhoto.h"
#include <math.h>


#define MAXLOGM 12
#define TWOPI   6.28318530717958647692
#define PI      3.14159265358979323846
#define SQHALF  0.707106781186547524401
#define COSM2   0.92387953251129
#define SINM2   0.38268343236509

/**
 * Implementacao da aplicacao RadioDctPhoto. RadioDctPhoto amostra a tensao no fotoresistor
 * da placa mda100 com uma frequencia de 1Hz. O valor da tensao do fotoresistor e armazenado e depois
 * de 32 amostras a DCT e calculada e os coeficientes sao transmitidos em pacotes AM.
 *
 * @author Savio Neves
 * @date   22/05/14
 */

module RadioDctPhotoC @safe(){
  uses {
    interface Leds;
    interface Boot;
    interface Timer<TMilli> as MilliTimer;
    interface Timer<TMilli> as MilliTimer2;
    interface Timer<TMilli> as MilliTimer3;
    interface Read<uint16_t> as Photo;
    interface SplitControl as RadioControl;
    interface AMSend as PhotoSend;
    interface Packet as PhotoPacket;
    interface AMSend as CoefSend;
    interface Packet as CoefPacket; 
    interface GeneralIO as PortG1;       
  }
}
implementation {
  
  message_t packett;
  message_t packetc;
  bool lockedPhoto = FALSE;
  bool lockedCoef = FALSE;
  int counter=0;  
  int k;

  float photo[NREAD];
  uint16_t coef;
  uint16_t coef2;


  static float *yt[MAXLOGM];
  static    int   brseed[256];     /* Evans' seed table */
  static    int     brsflg;         /* flag for table building */

  struct elemento *inicio;
  uint16_t grupo=0;

  void BR_permute(float *x, int logm){
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

  void srrec(float *xr, float *xi, int logm){
    static    int        m, m2, m4, m8, nel, n;
    static    float      *xr1, *xr2, *xi1, *xi2;
    static    float      *cn, *spcn, *smcn, *c3n, *spc3n, *smc3n;
    static    float      tmp1, tmp2, ang, c, s;
    static    float      *tab[MAXLOGM];

    /* Check range of logm */
    if ((logm < 0) || (logm > MAXLOGM)) {
	  call Leds.led2On();
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
      if ((tab[logm-4] = (float *) calloc(6 * nel, sizeof(float)))== NULL) {
	call Leds.led2On();
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
      }else {
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

  void srfft(float *xr, float *xi, int logm){
    /* Call recursive routine */
    srrec(xr, xi, logm);

    /* output array unshuffling using bit-reversed indices */
    if (logm > 1) {
      BR_permute(xr, logm);
      BR_permute(xi, logm);
    }
  } 
  
  void  fdctiv2(float *x, int logm, float fac){
    int        i, m, m2, m4, n, nel;
    static     float      *y, *xp1, *xp2, *yp1, *yp2;
    static     float      *c4n, *spc4n, *smc4n, *cn, *spcn, *smcn;
    static     float      tmp1, tmp2, ang, c, s;
    static     float      *tab[MAXLOGM];

    /* Check range of logm */
    if  ((logm < 0) || (logm > MAXLOGM)) {
	  call Leds.led2On();
    }

    /*  Trivial  cases, m = 1 and m = 2 */
    if (logm < 1) {
      *x *= fac;
      return;
    }

    if (logm == 1) {
      if (tab[0] == NULL) {
        if ((tab[0] = (float *) calloc(2, sizeof(float))) == NULL) {
	  call Leds.led2On();
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
      if ((tab[logm-1] = (float *) calloc(6 * nel, sizeof(float)))== NULL) {
	  call Leds.led2On();
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
        }else{
          ang = n * PI / m;
          c = fac * cos(ang); s = fac * sin(ang);
          *cn++ = c; *spcn++ = - (s + c); *smcn++ = s - c;
        }
      }
    }

    /* Allocate space for working vector, if necessary */
    if (yt[logm-1] == NULL) {
      if ((yt[logm-1] = (float *) calloc(m2, sizeof(float))) == NULL) {
	  call Leds.led2On();
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

  static void dctrec(float *x, int logm){
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

  void  fdct(float *x, int logm){
    static  float     tmp, *xp, *y, *yp;
    static  int       gr, n, m, m2, dx;

    /* Check range of logm */
    if ((logm < 0) || (logm > MAXLOGM)) {
	  call Leds.led2On();
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
	  call Leds.led2On();
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

  void InicializaFila (struct elemento **fila) {
    *fila = NULL;
    return;
  }

  int FilaVazia (struct elemento *fila) {
    if (fila == NULL)
      return 1;
    else
      return 0;
  }

  void InsereFila(struct elemento **fila, int dadonovo) {
    struct elemento *f1, *f2;
    f1 = malloc (sizeof (struct elemento));
    f1->dado = dadonovo;
    f1->prox = NULL;
    if (*fila == NULL)
      *fila = f1;
    else {
      f2 = *fila;
      while (f2->prox != NULL)
        f2 = f2->prox;
      f2->prox = f1;
    }
    return;
  }

  int RetiraFila(struct elemento **fila) {
    struct elemento *f1;
    int dado;
    f1 = *fila;
    *fila = f1->prox;
    dado = f1->dado;
    free (f1);
    return dado;
  }

  void EnviaPacotes(){
    radio_coef_msg_t* rsm;
    rsm = (radio_coef_msg_t*)call CoefPacket.getPayload(&packetc, sizeof(radio_coef_msg_t));
    if (rsm == NULL) {
      return;
    }

    grupo=(grupo+1);
    rsm->coefic=coef;
    rsm->coefic2=coef2;
    rsm->grupo=grupo;

    if (call CoefSend.send(AM_BROADCAST_ADDR, &packetc, sizeof(radio_coef_msg_t)) == SUCCESS) {
      lockedCoef = TRUE;
    }
  }

  event void Boot.booted() {
    call PortG1.makeOutput();
    call PortG1.set();
    call MilliTimer2.startOneShot(1);
  }

  event void RadioControl.startDone(error_t err) {
    if (err == SUCCESS) {
      EnviaPacotes();
      call Leds.led0Toggle();
    }
  }

  event void RadioControl.stopDone(error_t err) {
    if (err == SUCCESS) {
      call Leds.led0Toggle();
    }
  }
  
  event void MilliTimer.fired() {
    call Photo.read();
  }

  event void MilliTimer2.fired() {
      call PortG1.makeOutput();
      call PortG1.clr();
      call MilliTimer3.startOneShot(1);     
  }

  event void MilliTimer3.fired() {
      call PortG1.makeOutput();
      call PortG1.set();
      call MilliTimer.startPeriodic(1000);
  }

  event void Photo.readDone(error_t result, uint16_t data) {
    if (lockedPhoto) {
      return;
    }
    else {
      photo[counter]=data;
      counter++;
      if(counter==NREAD){
          counter=0;
	  coef=0;
	  coef2=0;

	  fdct(photo,LOG);
          InicializaFila(&inicio);
	  for (k=0;k<NREAD;k++){
	      if((int) photo[k] != 0){
		 if(k<16)
	     	    coef|=1<<k;
		 if((k>=16)&&(k<32))
	     	    coef2|=1<<(k-16);
		 InsereFila(&inicio,(int)photo[k]);	
              }
	  }
          call RadioControl.start();
      }
    }
  }

  event void PhotoSend.sendDone(message_t* bufPtr, error_t error) {
    if (&packett == bufPtr) {
      if (lockedPhoto == TRUE) {
          lockedPhoto = FALSE;
	  if(!FilaVazia(inicio)){
	    radio_sense_msg_t* rsm;
            rsm = (radio_sense_msg_t*)call PhotoPacket.getPayload(&packett, sizeof(radio_sense_msg_t));
            if (rsm == NULL) {
	      return;
            }
	    rsm->dataint=RetiraFila(&inicio);
	    rsm->grupo=grupo;
	    if (call PhotoSend.send(AM_BROADCAST_ADDR, &packett, sizeof(radio_sense_msg_t)) == SUCCESS) {
	      lockedPhoto = TRUE;
            }
	  }
	  if(FilaVazia(inicio)){
	    call RadioControl.stop();
	  }
      }
    }
  }

  event void CoefSend.sendDone(message_t* bufPtr, error_t error) {
    if (&packetc == bufPtr) {
      if (lockedCoef == TRUE) {
         radio_sense_msg_t* rsm;
  
         lockedCoef = FALSE;


          rsm = (radio_sense_msg_t*)call PhotoPacket.getPayload(&packett, sizeof(radio_sense_msg_t));
          if (rsm == NULL) {
	      return;
          }
	  rsm->dataint= RetiraFila(&inicio);
	  rsm->grupo=grupo;
	  if (call PhotoSend.send(AM_BROADCAST_ADDR, &packett, sizeof(radio_sense_msg_t)) == SUCCESS) {
	      lockedPhoto = TRUE;
          }

      }
    }
  }
}
