 
#include "Timer.h"
#include "RadioPhoto.h"

/**
 * Implementacao da aplicacao RadioPhoto. RadioPhoto amostra a tensao no fotoresistor
 * da placa mda100 em uma frequencia de 1Hz. O valor da tensao no fotoresistor e 
 * transmitido em pacotes AM.
 *
 * @author Savio Neves
 * @date   22/05/14
 */

module RadioPhotoC @safe(){
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
    interface GeneralIO as PortG1;    
  }
}
implementation {
  
  message_t packett;
  bool lockedPhoto = FALSE;

  struct elemento *inicio;
  uint16_t photo[NREAD];
  int counter=0;  
  int k;
 
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
    radio_sense_msg_t* rsm;
    rsm = (radio_sense_msg_t*)call PhotoPacket.getPayload(&packett, sizeof(radio_sense_msg_t));
    if (rsm == NULL) {
      return;
    }
    rsm->data= RetiraFila(&inicio);
	
    if (call PhotoSend.send(AM_BROADCAST_ADDR, &packett, sizeof(radio_sense_msg_t)) == SUCCESS) {
      lockedPhoto = TRUE;
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
	InicializaFila(&inicio);
	
	for (k=0;k<NREAD;k++){
	  InsereFila(&inicio,photo[k]);	
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
	  rsm->data= RetiraFila(&inicio);
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

}
