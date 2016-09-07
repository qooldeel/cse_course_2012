#ifndef SOME_GENERALLY_USED_CONSTANTS_HH
#define SOME_GENERALLY_USED_CONSTANTS_HH

namespace lebnum{

  template<class T>
  class Constants{
  public:
    static const T TinyNumber;

  };

  template<class T> const T Constants<T>::TinyNumber = 1.e-12;

}

#endif
