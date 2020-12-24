#include "simplesequence.h"


template<> double SimpleSequence<float>::s_subfrac = 1;
template<> int SimpleSequence<float>::s_subdiv = 0;


template<> double SimpleSequence<unsigned char>::s_subfrac = 1;
template<> int SimpleSequence<unsigned char>::s_subdiv = 0;
