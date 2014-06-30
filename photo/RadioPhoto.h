
#ifndef RADIO_PHOTO_H
#define RADIO_PHOTO_H


enum {
  AM_RADIO_SENSE_MSG = 7,
  NREAD=32,
};

typedef nx_struct radio_sense_msg {
  
  nx_uint16_t data;
  
} radio_sense_msg_t;

struct elemento{
  int dado;
  struct elemento *prox;
};


#endif
