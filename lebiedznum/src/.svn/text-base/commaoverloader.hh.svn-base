#ifndef COMMA_OVERLOADING_IN_C_PLUS_PLUS_HH
#define COMMA_OVERLOADING_IN_C_PLUS_PLUS_HH

#include <stdio.h>
#include <stdlib.h>
#include <iostream>

namespace lebnum{

  /**
   * \brief Overload commas; code adapted from [1]
   * 
   * Example:
   * \code 
         Vector<double> v(4);
         v <<= 0.5, -2.25, 1.5, 0.75; 
   * \endcode
   *
   * Reference:
   *
   * [1] [T. Veldhuizen, "Techniques for Scientific C++", Tech. Report, Indiana University Computer Science, 2000]
   *
   */
  template<class T, class ITER>
  class CommaOverloading{
  public:
    typedef std::size_t SizeType;
    typedef T value_type;
    typedef ITER iter_type;
    
    //! \param iter iterator of random access container
    //! \param n length of container (for size check)
    //! \param count counter that counts commata (for size check)  
    CommaOverloading(ITER iter, const SizeType& n, SizeType count):it_(iter),length_(n),count_(count){}

    CommaOverloading operator,(const T& t){ //! overload comma here
      count_++;
      //!check if there are more entries than lenght_ are assigned (less are 
      //! not problematic, since following entries are zero)
#ifndef NDEBUG    
      if(count_ >= length_){ 
	std::cerr << "****Error: CommaOverloading<T,ITER>::operator, : Length of rac and number of assigned types disagree:\n   There are more elements than container can store." << std::endl;
	abort();
      }
#endif
      *it_ = t;
      return CommaOverloading(it_+1,length_,count_);
    }

  private:
    ITER it_;
    const SizeType& length_;   //! length of random access container
    SizeType count_;          //! count entries to the right of assignment 
                               //! operator (assume that count_ = 0)

    CommaOverloading(){}  //no default constructor
  };

} //end namespace 

#endif
