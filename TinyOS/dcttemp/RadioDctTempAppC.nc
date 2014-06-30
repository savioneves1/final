 
#include "RadioDctTemp.h"

/**
 * Configuracao para a aplicacao RadioDctTemp. RadioDctTemp amostra o sensor de temperatura
 * da placa mda100 com uma frequencia de 1Hz. O valor da temperatura e armazenado e depois 
 * de 32 amostras a DCT e calculada e os coeficientes sao transmistidos em pacotes AM.
 *
 * @author Savio Neves
 * @date   22/05/14
 * @param AM_RADIO_SENSE_MSG Identificacao do pacote AM 
 * @param AM_RADIO_COEF_MSG Identificacao do pacote AM
 */

configuration RadioDctTempAppC {}
implementation {
  components MainC, RadioDctTempC as App, LedsC, new TempC();
  components ActiveMessageC;
  components new AMSenderC(AM_RADIO_SENSE_MSG) as TempSenderC;
  components new AMSenderC(AM_RADIO_COEF_MSG) as CoefSenderC;
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
  App.Temp -> TempC;
  App.TempSend -> TempSenderC;
  App.TempPacket -> TempSenderC;
  App.CoefSend -> CoefSenderC;
  App.CoefPacket -> CoefSenderC;
  App.PortG1 -> HplAtm128GeneralIOC.PortG1;  
}
