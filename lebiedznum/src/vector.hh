#ifndef LEBIEDZ_NUMERICS_VECTOR_HH
#define LEBIEDZ_NUMERICS_VECTOR_HH

/**
 * @ file a simple vector class.
 * Deriving from the STL has the great advantage that we can simply use
 * std::vector's member functions, e.g.
 * \code 
    lebnum::Vector<double> v(5);
    std::cout << "length of v: "<< v.size() << std::endl;
    // copy construction
    double af[] = {0.5,-2.15,0.75,0.85,-1.5};
    lebnum::Vector<double> w(af,af+5);
    std::cout << "w = "<< w << std::endl;
    lebnum::Vector<double> x(w);
    std::cout << "y = "<< y << std::endl;
 * \endcode 
 * We only have to define members for 
 * e.g. vector-vector addition, scalar-vec multiplication and so forth,
 * thus expanding std::vector's versatility. 
 *
 * TODO: Feel free to improve the class to your liking
*/
#include <stdio.h>
#include <stdlib.h>

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cassert>
#include <iomanip>
#include <cmath>
#include <algorithm>

#include "common.hh"
#include "commaoverloader.hh"
#include "type.hh"

namespace lebnum{  

  template<class T>
  class Vector: public std::vector<T>{
  public:
    //!aliases for different types 
    typedef T value_type;
    typedef std::size_t size_type;
    typedef std::vector<T> BaseVecType;
    typedef std::string StringType;
    typedef typename BaseVecType::iterator iterator;
    typedef typename BaseVecType::const_iterator const_iterator;

    enum{outputprecision = 12};

    //! constructor -- serves also as default vec
    Vector(size_type n = 0, const T& val = T()):BaseVecType(n,val){}
    
    //! 2nd constructor -- use std::vector's constructors
    template<class ITER>
    Vector(ITER itstart, ITER itend):BaseVecType(itstart,itend){}
    
    //! quite cool MATLAB(r)-style value assignment
    CommaOverloading<T,iterator> operator=(const T& val){
      size_type count = 0;
      (*this)[0] = val;
      return CommaOverloading<T,iterator>((*this).begin()+1,(*this).size(),count);
    }
    
    //! you can use the paranthesis operator in conjunction to []
    T& operator()(size_type i){
      assert(i < (*this).size());
      return this->operator[](i);
    }

    //read only access
    const T& operator()(size_type i) const{
      assert(i < (*this).size());
      return (*this)[i];
    }

    //! the std::vector has no output functionality
    //! the output can then be easily reused in MATLAB(r) 
    friend std::ostream& operator<<(std::ostream& os, const Vector& v){
      os << " [  "; 
      for(size_type i = 0; i < v.size(); ++i){
	os << std::setprecision(outputprecision) << write_complex(v[i],'i',outputprecision) << "  "; 
      }
      os << "]'" << std::endl;
      return os;
    }

    std::ofstream& print(std::ofstream& of){
      of << std::setprecision(outputprecision);
      for(size_type i = 0; i < (*this).size(); ++i){
	of << write_complex((*this)[i],'i',outputprecision) << "  ";
      }
      of << std::endl;
      return of;
    }

    //!###################################################################
    //!##################  OPERATIONS ON VECTORS #########################
    //!###################################################################
    
    //! \f$ v = 2*v \f$
    Vector& operator*=(const T& val){
      for(size_t i = 0; i < (*this).size(); ++i) //same as this->size()
	this->operator[](i) *= val;
      return *this;
    }
    
    //! \f$ v = v/val \f$, note that vector is overwritten!
    Vector& operator/=(const T& val){
      assert(!is_zero(val));  //! prevent division by (nearly) zero values
      for(size_t i = 0; i < (*this).size(); ++i) //same as this->size()
	this->operator[](i) /= val;
      return *this;
    }

    //! \f$ v = v+w\f$
    Vector& operator+=(const Vector& w){
      assert((*this).size() == w.size()); //!have to be same length
      for(size_t i = 0; i < (*this).size(); ++i) //same as this->size()
	this->operator[](i) += w[i];
      return *this;
    }
    
    //! \f$ v = v-w\f$
    Vector& operator-=(const Vector& w){
      assert((*this).size() == w.size()); //!have to be same length
      for(size_t i = 0; i < (*this).size(); ++i) //same as this->size()
	this->operator[](i) -= w[i];
      return *this;
    }

    Vector operator+(const Vector& w) const {
      assert((*this).size() == w.size());
      Vector y(*this); //!copy construct temporary object
      y += w;          //!apply += from above
      return y;
    }

    Vector operator-(const Vector& w) const {
      assert((*this).size() == w.size());
      Vector y(*this); //!copy construct temporary object
      y -= w;          //!apply -= from above
      return y;
    }

   
    //! usage: -x, i.e. x *= -1 
    Vector& operator-(){
      return ( (*this) *= -1 ); 
    }

  }; //end of class

  //! some additional vector operations: right multiplication with a scalar
  template<class T>
  inline Vector<T> operator*(const T& val, const Vector<T>& v){
    Vector<T> w(v);  //! copy construction
    for(std::size_t i = 0; i < v.size(); ++i)
      w[i] *= val;
    return w;
  }

   //! some additional vector operations: left multiplication with a scalar
  template<class T>
  inline Vector<T> operator*(const Vector<T>& v, const T& val){
    return val*v;  //use the stuff from above
  }

  template<class T>
  inline Vector<T> operator/(const Vector<T>& v, const T& val){
    assert(Absval(val) > 1.e-13);
    Vector<T> w(v);  //! copy construction
    for(std::size_t i = 0; i < v.size(); ++i)
      w[i] /= val;
    return w;
    
  }

  //! give back norm of a vector (default: 2-norm)
  //! Other possibilities: 1-norm: '1', infty-norm: 'i' or 'I' 
  template<class T> 
  inline typename SelectProperType<T>::BaseType norm(const Vector<T>& v, char c = '2'){
    typedef typename SelectProperType<T>::BaseType BaseType;   
    BaseType nm = BaseType(); 
    if(c == '1'){
      for(std::size_t i = 0; i < v.size(); ++i)
	nm += std::abs(v[i]);
    }
    else if(c == '2'){
      for(std::size_t i = 0; i < v.size(); ++i){
	nm += Sqr(std::abs(v[i]));
      }
      nm = std::sqrt(nm);
    }
    else if(c == 'i' || c == 'I'){
      for(std::size_t i = 0; i < v.size(); ++i){
	nm = std::max(nm,std::abs(v[i]));
      }
    }
    else{
      std::cerr << "Nothing defined for c = " << c << std::endl;
      abort();
    }
    return nm;
  }

  //! the dot product, i.e. dt = v'*w, also for complex numbers!
  //! Note that Conj(·) just returns the value when not being a complex number!
  template<class T>
  inline T dot(const Vector<T>& v, const Vector<T>& w){
    assert(v.size() == w.size());
    T dt = T();
    for(std::size_t i = 0; i < v.size(); ++i)
      dt += v[i]*Conj(w[i]);
    return dt;
  }

  

} //end namespace

#endif
