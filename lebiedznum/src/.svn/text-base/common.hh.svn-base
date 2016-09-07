#ifndef COMMON_FUNCTIONALITIES_HH
#define COMMON_FUNCTIONALITIES_HH

#include <cmath>
#include <cassert>
#include <complex>
#include <sstream>
#include <string>
#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <stdlib.h>

#include "type.hh"

namespace lebnum{

  //! if a value (may be complex) is in absolute value smaller than a 
  //! prescribed positive threshold then it is considered (approx.) zero
  template<class T>
  inline bool is_zero(const T& val, const typename SelectProperType<T>::BaseType& eps = 1.e-13){
    typedef typename SelectProperType<T>::BaseType BaseType;
    assert(eps > BaseType() && eps < 1.e-08); //!guarantee sufficiently small value
    return ( (std::abs(val) <= eps) ? true : false );
  }

  template<class T>
  inline T Sqr(const T& x){
    return x*x;
  }

  template<class T>
  inline T Cub(const T& x){
    return Sqr(x)*x;
  }

  template<class T>
  inline T Absval(const T& x){
    return ( (x>= T()) ? x : -x);
  }

  //! specialisation for complex numbers
  template<class T>
  inline T Absval(const std::complex<T>& z){
    return std::sqrt(Sqr(z.real()) + Sqr(z.imag()));
  }

  template<class T> 
  inline T NaturalPow(const T& x, int e){
    int p = Absval(e);
    T pow(1);

    for(int i = 0; i < p; ++i)
      pow *= x;

    return ((e >= 0) ? pow : T(1)/pow);
  }

  template<class T>  //!for non complex values, just return the value
  inline const T& Conj(const T& x){
    return x;
  }
  
  //! partial specialization for complex numbers
  template<class T>
  inline std::complex<T> Conj(const std::complex<T>& x){
    return std::conj(x);
  }
    
  template<class T> 
  inline std::string num2str(const T& number, int prec = 14){
    std::stringstream giveback;
    giveback << std::setprecision(prec)  << number; 
    return giveback.str();
  }

  //!not necessarily the best formula to solve \f[ax^2 + bx + c = 0 \f], since cancellation may occur for slightly 
  //!different input values \f$a, b, c \f$
  template<class T>
  inline std::pair<std::complex<T>,std::complex<T> > quadratic_solution(const T& a, const T& b, const T& c){
    typedef std::pair<std::complex<T>,std::complex<T> > PairType;
    PairType p;
    std::complex<T> discr = b*b - 4*a*c, // now you can call complex sqrt!!
      st = std::sqrt(discr);
    p.first = (-b + st)/(2*a);
    p.second = (-b - st)/(2*a);
    return p;
  }


  //!convert string into number (int, long int, double,...): 
  //!invoke it via str2num<double>(string). For non-rounded output use, 
  //! e.g. 'std::cout.precision(17)' etc..
  template<class T> 
  inline T str2num(const std::string& s){
    return T(atof(s.c_str()));
  } 

  template<class T>  //last two arguments are dummy arguments
  inline const T& write_complex(const T& d, char isymbol = 'i',int prec = 14){
    return d;
  }
  
  //! output when \f$z\f$ is complex: 'a + bi'
  template<class T>
  inline std::string write_complex(const std::complex<T>& z, char isymbol = 'i',int prec = 14){
    std::string oper;
    if(z.imag() >= T())
      oper = "+";
    std::string s = num2str(z.real())+oper+num2str(z.imag())+isymbol;
    return s;
  }

} //end namespace 

#endif
