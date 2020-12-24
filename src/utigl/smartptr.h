#ifndef SMART_PTR_H
#define SMART_PTR_H

#include <stdlib.h>

//new gcc needs this, it is actually declard in globals.h
char * combine(char * fmt,...);

template <class T>
class SmartPtr{
 private:
  T * ptr;
  int * ref;
 public:
  SmartPtr(T * ptr=0){
    this->ptr=ptr;
    ref=new int;
    (*ref)=0;
    incRefCnt();
  }
  SmartPtr(const SmartPtr & other){
    //twLog("copy\n");
    ref=0;
    ptr=0;
    operator=(other);
  }
  void operator=(const SmartPtr & other){
    if(this==&other)return;
    decRefCnt();
    this->ref=other.ref;
    ptr=other.ptr;
    incRefCnt();
  }
  ~SmartPtr(){
    decRefCnt();
  }
  void incRefCnt(){
    if(ref)
      (*ref)++;
  }
  void decRefCnt(){
    if(ref){
      (*ref)--;
      if((*ref)==0){
	if(ptr)delete ptr;
	delete ref;
	ref=0;
	ptr=0;
      }
      else if(*ref<0){
	//fprintf(stderr,combine("references less than 0:%d",*ref));
	exit(1);
      }	
    }
  }
  const T & operator*() const{
    return *ptr;
  }
  T & operator*(){
    return *ptr;
  }
  T * operator->(){
    return ptr;
  }
  const T * operator->() const {
    return ptr;
  }
  T * get() {
    return ptr;
  }
  const T * get() const {
    return ptr;
  }
  /*int operator!(){
    return empty();
    }*/
  int empty() const {
    if(!ptr || !ref || !(*ref))return 1;
    return 0;
  }
  int getRefCnt() const {
    if(ref)return *ref;
    return -1;
  }
};

#endif
