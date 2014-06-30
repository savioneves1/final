
#ifndef RADIO_DCT_TEMP_H
#define RADIO_DCT_TEMP_H


enum {
  AM_RADIO_SENSE_MSG = 7,
  AM_RADIO_COEF_MSG = 9,
  NREAD=32,
  LOG=5,
};

typedef nx_struct radio_sense_msg {
  
  nx_int16_t dataint;
  nx_uint16_t grupo;

} radio_sense_msg_t;

typedef nx_struct radio_coef_msg {
  
  nx_uint16_t coefic;
  nx_uint16_t coefic2;
  nx_uint16_t grupo;

} radio_coef_msg_t;

struct elemento{
  int dado;
  struct elemento *prox;
};


#endif
