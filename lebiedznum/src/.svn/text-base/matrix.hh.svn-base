#ifndef CONTAINER_FOR_A_DENSE_MATRIX_HH
#define CONTAINER_FOR_A_DENSE_MATRIX_HH

#include <iostream>
#include <iomanip>
#include <vector>
#include <cassert>
#include <algorithm>

#include "common.hh"
#include "commaoverloader.hh"
#include "consts.hh"
#include "vector.hh"

namespace lebnum{

  /**
   * \brief <BB>Dense</BB> matrix object. 
   *
   * Example how to create a 3 x 3 matrix:
   * \code 
     lebnum::Matrix<double> AA(3,3);
     AA = 1, 4, 7,          //just written in a matrix-look-alike-style    
          2, 5, 8,
          3, 6, 10;
   *\endcode
   */
  template<class T>
  class Matrix{
  public:
    typedef T value_type;
    typedef std::size_t size_type;
    typedef std::vector<T> VType;
    typedef typename VType::iterator DataTypeIterator;
    typedef typename VType::const_iterator ConstDataTypeIterator;
    typedef typename SelectProperType<T>::BaseType BaseType;

    typedef value_type* TpointerType;
    typedef const value_type* ConstTpointerType;

    enum{outputprecision = 12};

    //! 1st constructor; serves also as default constructor
    Matrix(size_type rows = 0, size_type cols = 0, const T& val = T()):rows_(rows),cols_(cols),mtx_(rows_*cols_,val){}

    Matrix(size_type rows, size_type cols, T** ptr):rows_(rows),cols_(cols),mtx_(rows*cols){
      mtx_.reserve(rows*cols);
      for(size_type i =  0; i < rows_; ++i) {      
	for (size_type j = 0; j < cols_; ++j) 
	  mtx_[offset(i,j)] = ptr[i][j];
      }
    }

    //! quite cool MATLAB(r)-style value assignment
    //! Usage:
    //! \code
    //!  Matrix<double> A(3,4);
    //! 
    //!  A = 0.5,  2.75, -1,   0.75,
    //!      1.25, -1.5, 3.25, 2.45,
    //!      -3,   -7.5, 0.65, -1.85;
    //!\endcode
    CommaOverloading<T,DataTypeIterator> operator=(const T& val){
      size_type count = 0;
      mtx_[0] = val;
      return CommaOverloading<T,DataTypeIterator>(mtx_.begin()+1,mtx_.size(),count);
    }

    //!copy constr.
    Matrix(const Matrix& m):rows_(m.rows_),cols_(m.cols_),mtx_(m.mtx_){}

    //!copy assignment
    Matrix& operator=(const Matrix& m){
      if(this != &m){ //! no self-assignment
	rows_ = m.rows_;
	cols_ = m.cols_;
	mtx_ = m.mtx_;
      }
      return *this;
    }

    //! read-only access
    const T& operator()(size_type i, size_type j) const{
      assert(i < rows_ && j < cols_);
      return mtx_[offset(i,j)];
    }

    //!rw access
    T& operator()(size_type i, size_type j){
      assert(i < rows_ && j < cols_);
      return mtx_[offset(i,j)];
    }

    //!operator[], you can use A[i][j] to access element \f$a_{ij}\f$
     ConstDataTypeIterator operator[](size_type i) const{
      assert(i < rows_);
      return mtx_.begin() + i*cols_;
    }

    //read write access
    DataTypeIterator operator[](size_type i){
      assert(i < rows_);
      return mtx_.begin() + i*cols_;
    }
    


    //! return pointer to matrix storage container mtx_
    TpointerType data(){ return &mtx_[0];}
    ConstTpointerType data() const {return &mtx_[0];}

    //! well not the nicest output, but you can just copy the console output
    //! to Matlab(r)
    friend std::ostream& operator<<(std::ostream& os, const Matrix& m){
      os  << " [  ";
      for(size_type i = 0; i < m.rows_; ++i){
	for(size_type j = 0; j < m.cols_; ++j){
	  os << std::setprecision(outputprecision) << write_complex(m(i,j),'i',outputprecision) << "  ";
	}
	if(i != m.rows()-1)
	  os << ";" << std::endl;
	else 
	  os << "]" << std::endl;
      }
      return os;
    }

    const size_type rows() const {return rows_;}
    const size_type cols() const {return cols_;}

    Matrix& add_unity_matrix(){
      assert(rows_ == cols_);
      for(size_type i = 0; i < rows_; ++i)
        this->operator()(i,i) += 1;
      return *this;
    }

    Matrix& operator*=(const T& val){
      for(size_type i = 0; i < rows_; ++i)
	for(size_type j = 0; j < cols_; ++j)
	  (*this)(i,j) *= val;
      return *this;
    }

    Matrix& operator-(){
      return ( (*this) *= -1 ); 
    }
    
     Matrix& operator/=(const T& val){
       assert(!is_zero(val));
       for(size_type i = 0; i < rows_; ++i)
	for(size_type j = 0; j < cols_; ++j)
	  (*this)(i,j) /= val;
      return *this;
    }

    Matrix& operator+= (const Matrix& B)
    {
      for(size_type i=0; i < rows_; ++i) 
        for(size_type j=0; j < cols_; ++j) 
	  (*this)(i,j) += B(i,j);
      return *this;
    }

    
    Matrix& operator-= (const Matrix& B)
    {
      for(size_type i=0; i < rows_; ++i) 
        for(size_type j=0; j < cols_; ++j) 
	  (*this)(i,j) -= B(i,j);
      return *this;
    }

    Matrix operator/(const T& val) const{
      assert(!is_zero(val));
      Matrix a(*this);
      return a/=val;
    }

    Matrix operator+(const Matrix& B) const {
      assert((*this).rows() == B.rows() && (*this).cols() == B.cols());
      Matrix C(*this); //!copy construct temporary object
      C += B;          //!apply += from above
      return C;
    }

    Matrix operator-(const Matrix& B) const {
      assert((*this).rows() == B.rows() && (*this).cols() == B.cols());
      Matrix C(*this); //!copy construct temporary object
      C -= B;          //!apply += from above
      return C;
    }

    // u = M*v
    Vector<T> operator*(const Vector<T>& v) const{
      assert((*this).cols() == v.size());
      Vector<T> u((*this).rows()); 
      for(size_type i=0; i < (*this).rows(); ++i){
	for(size_type j=0; j < (*this).cols(); ++j){
	  u[i] += (*this)(i,j)*v[j];
	}
      }
      return u;
    }

    //! Direct linear solver based on  PLU decomposition, i.e.
    //! Gaussian elimination with partial pivoting. We solve the lin. sys.
    //! \f[ A\cdot x = b, \qquad A \in R^{n\times n}\f]
    //! note: both *this and b will be overwritten: A = PLU, b holds the 
    //  solution and is returned
    //! Complete pivoting is more robust but due to increased computational
    //! costs it is hardly used in practice
    //! special feature: the determinant is calculated simultaneously!
    //! Note that the det is hardly used, and it just gives you some info
    //! about the regularity of the matrix under consideration
    //!\tparam V some STL-compliant random access container
    //! \param b rhs 
    //! \param caldet calculates determinant from PLU (default val: false)
    template<class V>
    V& solve(V& b, bool calcdet = false){
      if(rows_ != cols_){ 
	std::cerr << ("Matrix<T>::solve: Only square matrices allowed.") << std::endl;
	abort();
      }
      if(cols_ != b.size()){
	std::cerr << "Matrix<T>::solve: Vector size and coefficient matrix dimension do not match."<<std::endl;
	abort();
      }
	
      //!compute A = PLU, where pvt contains the column indices for 
      //!the permutation matrix P
      int nrows = static_cast<int>(rows_),
	ncols = static_cast<int>(cols_);
      
      int sign(1);   //signum for determinant
     
      Vector<int> pvt(nrows);       //vector containing pivot elements 
      for(int k = 0; k < nrows; ++k) 
	pvt[k] = k;
      
      for (int k = 0; k < nrows-1; ++k) {  // main loop
	
	// find the pivot in column k in rows pvt[k], 
	// pvt[k+1], ..., pvt[n-1]
	int pc  = k; 
	BaseType aet = Absval((*this)(pvt[k],k));
	for (int i = k+1; i < nrows; ++i) {
	  if (Absval((*this)(pvt[i],k)) > aet) {
	    aet = Absval((*this)(pvt[i],k)); 
	    pc = i;
	  }
	}
	if (!aet){ 
	  std::cerr << "Matrix<T>::solve: pivot is zero" << std::endl; 
	  abort();
	}
	if (pc != k){ 
	  std::swap(pvt[k], pvt[pc]);
	  sign *= -1;             //determinant: even # of permutations: +
	}                         //             odd # of permutations:  -
	int pvtk = pvt[k];                  // pivot row
	T pivot = (*this)(pvtk,k);            // pivot

	// now eliminate column entries logically 
	// below (*this)(pvt[k],k)
	for (int i = k + 1; i < nrows; ++i) {
	  int pvti = pvt[i];
	  if ((*this)(pvti,k) != T()) {
	    T mult = (*this)(pvti,k)/pivot;
	    (*this)(pvti,k) = mult;
	    for (int j = k + 1; j < ncols; ++j) 
	      (*this)(pvti,j) -= mult*(*this)(pvtk,j);
	  }
	}
      }

      //std::cout << "PLU = " << (*this) <<  std::endl; //!some output
      //std::cout << "piv = "<< pvt << std::endl;
      

      //!actual determinant calculation -- if <TT>calcdet = false</TT>,
      //!no determinant will be calculated to save time
      calc_determinant(pvt,sign,calcdet);
      
      //!Construct solution
      // forward substitution for L y = Pb.
      for (int i = 1; i < nrows; ++i)  
	for (int j = 0; j < i; ++j) 
	  b[pvt[i]] -= (*this)(pvt[i],j)*b[pvt[j]];
      
      // back substitution for Ux = y
      V x(nrows);  // x stores solution in correct order
      for (int i = nrows - 1; i >= 0; --i) {
	for (int j = i+1; j < ncols; ++j) 
	  b[pvt[i]] -= (*this)(pvt[i],j)*x[j];
	x[i] = b[pvt[i]] / (*this)(pvt[i],i);
      }
      
      b = x;             //assign solution to b
      
      return b;
    }

    //! compute A' resp. A* (complex conjugate matrix, i.e. the analogon to
    //! real matrix transposition
    Matrix transpose() const{
      Matrix Trans(cols_,rows_);
      for(size_type i=0; i < rows_; ++i){
	for(size_type j=0; j< cols_; ++j) Trans(j,i)= Conj((*this)(i,j)); 
      }
      return Trans;
    }

    //! note that an actual matrix transposition (and hence copying)
    //! can be avoided. Element access can be established by using
    //! <TT> A.transposed_offset(i,j) </TT> where the dimensions of 
    //! the fictitious transposed are given by the members below
    const size_type transposed_rows() const {return cols_;}
    const size_type transposed_cols() const {return rows_;}
    
    size_type transposed_offset(size_type i, size_type j) const{
      assert(i < (*this).transposed_rows() && j < (*this).transposed_cols());
      return i + j*rows_;
    }

    T trace() const{
      assert(rows_ == cols_);
      T tr = T();
      for(size_type i = 0; i < rows_; ++i)
	tr += (*this)(i,i);
      return tr;
    }


    //! power method, see [GOLUB,VAN LOAN, "Matrix computations", 3rd ed. 1996]
    //! pp. 330
    
    
  private:
    size_type rows_, cols_;
    VType mtx_;
   
    //! Assumption: rows are stored one after another (so-called) <I>row major</I> ordering (which is the usual C/C++-style, in constrast to the column-major FORTRAN/Matlab(r)-style).
    size_type offset(size_type i, size_type j) const{
      return i*cols_ + j;
    }

    //! don't allow the user to invoke this function -- gives only some info
    T calc_determinant(const Vector<int>& piv, const int sign, bool calcdet){
      assert(piv.size() == rows_);
      T det(1);
      if(calcdet){
	for(size_type l = 0; l < rows_; ++l)
	  det *= (*this)(piv[l],l);
	det *= sign;
	std::cout << "DETERMINANT = "<< det << std::endl;
      }
      return det;
    }
  }; //END OF CLASS MATRIX

  

  //!additional operators defined outside class Matrix
  template<class T>
  inline Matrix<T> operator*(const T& val, const Matrix<T>& mtx){
    Matrix<T> m(mtx);
    m *= val;
    return m;
  }

  template<class T>
  inline Matrix<T> operator*(const Matrix<T>& mtx, const T& val){
    return val*mtx;
  }

  // s = v'M
  template<class T>
  inline Vector<T> operator*(const Vector<T>& w, const Matrix<T>& M){
    assert(w.size() == M.rows());
    typedef typename Vector<T>::size_type size_type;
    Vector<T> s(M.cols());
    for(size_type i = 0; i < M.rows(); ++i){
      for(size_type j = 0; j < M.cols(); ++j){
	s[j]+= w[i]*M(i,j);
      }
    }
  return s;
  }

  //! matrix-matrix multiplication
  template<class T>
  inline Matrix<T> operator*(const Matrix<T>& A, const Matrix<T>& B){
    if(A.cols() != B.rows()){
      std::cerr << "Matrix dimensions do not match! " << std::endl;
      abort();
    }
    typedef typename Matrix<T>::size_type size_type;

    Matrix<T> C(A.rows(), B.cols());
    for(size_type i=0; i<A.rows(); ++i){
      for(size_type j=0; j<B.cols(); ++j){ 
	for(size_type k=0; k<A.cols(); ++k){
	  C(i,j) += A(i,k)*B(k,j); 
	}
      }
    }  
    return C;
  }

  template<class T>
  inline typename SelectProperType<T>::BaseType norm(const Matrix<T>& a, char c = '1'){
    typedef typename SelectProperType<T>::BaseType BaseType;   
    BaseType nm = BaseType();
    typedef typename Vector<T>::size_type size_type;

    if(c == '1'){
      for(size_type j = 0; j< a.cols(); ++j){
	BaseType tmp = BaseType();
	for(size_type i =0; i < a.rows(); ++i) tmp += Absval(a(i,j));
	nm = std::max(nm, tmp);
      }
    }
    else if (c == '2'){
      std::cerr << "Well, the matrix 2-norm needs the calculation of eigen values. Feel free to do so ;)"<< std::endl;
      abort();
    }
    else if(c == 'i' || c == 'I'){ //max norm
      for(size_type i = 0;  i< a.rows(); ++i){
	BaseType tmp = BaseType();
	for(size_type j =0; j < a.cols(); ++j) tmp += Absval(a(i,j));
	nm = std::max(nm, tmp);
      }
    }
    else if(c == 'f' || c == 'F'){ //Frobenius norm
      for(size_type i=0; i < a.rows(); ++i) {
	for(size_type j=0; j < a.cols(); ++j) {
	  nm += Sqr(Absval(a(i,j))); 
	}
      }
      nm = std::sqrt(nm);
    }
    else{
      std::cerr << "Nothing defined for c = " << c << std::endl;
      abort();
    }
    return nm;
  }


  template<class T>         //!dyadic product, i.e. D = v*w'
  inline Matrix<T> dyad(const Vector<T>& v, const Vector<T>& u){
    Matrix<T> D(v.size(), u.size());
    typedef typename Vector<T>::size_type size_type;
    for(size_type i = 0; i < v.size(); ++i){
      for(size_type j = 0; j < u.size(); ++j){
	D(i,j) += v[i]*u[j];   //update A = A + vu^T
      }
    }
    return D;
  }
  

}//end namespace 

#endif
