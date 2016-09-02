//
// Created by Zhengyu on 8/30/16.
//

#ifndef INC_2DFLUIDSTRUCTURE_INTEGER_LOG_H
#define INC_2DFLUIDSTRUCTURE_INTEGER_LOG_H



inline int Integer_Log(unsigned int v) // this works for any v, but it is slower
{int c=0;
    if(v&0xffff0000){v>>=16;c|=16;}
    if(v&0xff00){v>>=8;c|=8;}
    if(v&0xf0){v>>=4;c|=4;}
    if(v&0xc){v>>=2;c|=2;}
    if(v&2)c|=1;return c;}


#endif //INC_2DFLUIDSTRUCTURE_INTEGER_LOG_H
