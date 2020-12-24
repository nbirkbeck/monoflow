#ifndef UTIL_TIMER_H
#define UTIL_TIMER_H

#include<sys/time.h>

#warning utigl/timer is depricated, use common/misc/timer

class Timer{
 protected:
  struct timeval tm;
 public:
  Timer(){
    updateTime();
  }
  void updateTime(){
    gettimeofday(&tm,0);
  }
  Timer operator-(const Timer & t){
    Timer result;
    timersub(&tm,&(t.tm),&(result.tm));
    return result;
  }
  double getDoubleTime(){
    return ((double)tm.tv_sec)+((double)tm.tv_usec)*1e-6;
  }
  operator double (){
    return getDoubleTime();
  }
};

#endif
