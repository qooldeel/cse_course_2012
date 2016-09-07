#ifndef STUFF_FOR_NUMERICAL_TYPES_HH
#define STUFF_FOR_NUMERICAL_TYPES_HH

#include <complex>

namespace lebnum{

  /**
   * \brief Selects the right numerical type
   *  For example, the (2-)norm for both a double and 
   *  a complex-valued vector is a double value. In such cases
   *  this tiny class guarantees that for different vector types, the 
   *  meaningful base type is always chosen.
   */
  template<class T>
  class SelectProperType{
  public:
    typedef T Type;
    typedef T BaseType;
  };

  //! partial specialization
  template<class T>
  class SelectProperType<std::complex<T> >{
  public:
    typedef std::complex<T> Type;
    typedef T BaseType;
  };

} //end namespace

#endif 
