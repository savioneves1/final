 
#include "RadioPhoto.h"

/**
 * Configuracao da aplicacao RadioPhoto. RadioPhoto amostra a tensao no fotoresistor
 * da placa mda100 em uma frequencia de 1Hz. O valor da tensao no fotoresistor e transmitido
 * em pacotes AM.
 *
 * @author Savio Neves
 * @date   22/05/14
 * @param AM_RADIO_SENSE_MSG Identificacao do pacote AM
 */

configuration RadioPhotoAppC {}
implementation {
  components MainC, RadioPhotoC as App, LedsC, new PhotoC();
  components ActiveMessageC;
  components new AMSenderC(AM_RADIO_SENSE_MSG) as PhotoSenderC;
  components new TimerMilliC() as Tempo1C;
  components new TimerMilliC() as Tempo2C;
  components new TimerMilliC() as Tempo3C;
  components HplAtm128GeneralIOC;

  App.Boot -> MainC.Boot;
  App.RadioControl -> ActiveMessageC;
  App.Leds -> LedsC;
  App.MilliTimer -> Tempo1C;
  App.MilliTimer2 -> Tempo2C;
  App.MilliTimer3 -> Tempo3C;
  App.Photo -> PhotoC;
  App.PhotoSend -> PhotoSenderC;
  App.PhotoPacket -> PhotoSenderC;  
  App.PortG1 -> HplAtm128GeneralIOC.PortG1;

}
