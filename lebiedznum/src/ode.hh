#ifndef NUMERICAL_METHODS_FOR_SOLVING_ODES_HH
#define NUMERICAL_METHODS_FOR_SOLVING_ODES_HH

#include <cmath>
#include <cassert>
#include <string>
#include <fstream>

namespace lebnum{

/**
 * @file Autonomous representation of ODEs because every time-dependent ODE can
 * be converted to an autonomous one by introducing a new state representing time.
*/

  template<class MODEL>
  class ClassicalRK{
  public:
    typedef typename MODEL::value_type value_type;
    typedef typename MODEL::size_type size_type;
    typedef typename MODEL::VType VType;

    ClassicalRK(MODEL& ode):ode_(ode){}

    template<class TSTEP, class V>
    VType& step(const TSTEP& deltaT, const V& x){
      k1_ = deltaT*ode_(x);
      k2_ = deltaT*ode_(x+0.5*k1_);
      k3_ = deltaT*ode_(x+0.5*k2_);
      k4_ = deltaT*ode_(x + k3_);

      next_ = x + 1./6.*(k1_ + 2.*k2_ + 2.*k3_ + k4_);
      return next_;
    }
    

  private:
    MODEL& ode_;  //reference
    VType next_, k1_,k2_, k3_,k4_;
  };


  //! should be the same as above
 template<class MODEL>
  class RungeKutta4{
  public:
    typedef typename MODEL::value_type value_type;
    typedef typename MODEL::size_type size_type;
    typedef typename MODEL::VType VType;

    RungeKutta4(MODEL& ode):ode_(ode){}

    template<class TSTEP, class V>
    VType& step(const TSTEP& deltaT, const V& x){
      k1_ = ode_(x);
      k2_ = ode_(x+0.5* deltaT*k1_);
      k3_ = ode_(x+0.5* deltaT*k2_);
      k4_ = ode_(x +  deltaT*k3_);

      next_ = x + 1./6.*deltaT*(k1_ + 2.*k2_ + 2.*k3_ + k4_);
      return next_;
    }
    

  private:
    MODEL& ode_;  //reference
    VType next_, k1_,k2_, k3_,k4_;
  };


  template<class MODEL>
  class ImplicitEuler{
  public:
    typedef typename MODEL::value_type value_type;
    typedef typename MODEL::size_type size_type;
    typedef typename MODEL::VType VType;
    typedef typename MODEL::MtxType MtxType;

    ImplicitEuler(MODEL& ode, const value_type& eps, const value_type& tol, const value_type& cfac = 1.):ode_(ode),eps_(eps),tol_(tol),Cfac_(cfac),h_(value_type()){}

    VType integrate(const value_type& t0, const value_type& tend,const VType& u0,const value_type& hestim, int maxit, const std::string& fname = "BackwardEuler.dat"){
      std::ofstream of(fname.c_str());

      h_ = hestim;
      value_type time(t0);
      Uprev_ = u0;
      U_ = u0;

      unsigned k(0);
      VType err;

      while(time < tend){
	std::cout << k << ".)   time: "<< time << "  dt = " << h_ << std::endl;
	of << time << "    ";
	U_.print(of);
	
	//Newton Iteration
	U_ = Uprev_;
	for(int l = 1; l <= maxit; ++l){
	  b_ = -(U_ - Uprev_ - h_*ode_(U_));
	  A_ = -h_*ode_.jacobian(U_);
	  A_.add_unity_matrix();
	  A_.solve(b_);
	  U_ += b_;
	  if(norm(b_) <= eps_)
	    break;
	  if(l == maxit)
	    std::cout << "NEWTON ITERATION DID NOT CONVERGE WITHIN "<< maxit << " ITERATIONS..." << std::endl;
	}
	err = U_ - Uprev_;  //STPZCTRL

	//stpszctrl
	value_type nmerr = norm(err);
	if(nmerr <= tol_)
	  h_ = tol_*h_/(Cfac_*nmerr);
	else
	  h_ *= 0.5;

	Uprev_ = U_;
	time += h_;
	k++;
      }
      of.close();
      return U_;
    }

  private:
    MODEL& ode_;  //reference
    value_type eps_, tol_, Cfac_, h_;  
    VType U_, Uprev_, b_;
    MtxType A_;
  };


   //! WORKS LIKE SHIT ON A CAROUSEL
  template<class MODEL>
  class Explicit4Stiff{
  public:
    typedef typename MODEL::value_type value_type;
    typedef typename MODEL::size_type size_type;
    typedef typename MODEL::VType VType;

    Explicit4Stiff(MODEL& ode, const value_type& TOL, const value_type& eps, const value_type& sigma, const value_type& c, const value_type& hmin):ode_(ode),TOL_(TOL),eps_(eps),sigma_(sigma),c_(c), hmin_(hmin),htilde_(value_type()),h_(value_type()){}

    VType& integrate(const value_type& t0, const value_type& tend,const VType& u0, const value_type& hestim, int maxit, const std::string& fname = "Explicit4Stiff.dat"){
      assert(maxit > 2);
      
      std::ofstream of(fname.c_str());

      Uprev_ = u0;
      
      value_type time(0.), L(0.), nm(0.), nmprev(0.),
	h_ = hestim;
      int k(0);
      bool flag(false);
      VType Ulp;

      while(time < tend){  //TIME LOOP
	std::cout << k << ".)   time: "<< time << "  h = " << h_ << std::endl;
	of << time << "    ";
	U_.print(of);
	
	if(k > 0){
	  R_ = (U_ - Uprev_)/h_ - ode_(Uprev_); //continuous residuum
	   htilde_ = TOL_/(sigma_*norm(R_));
	  // std::cout << "htilde = "<< htilde_ << std::endl;
	   h_ = (2.*htilde_*h_)/(htilde_+h_); //k_n <= 2*k_n-1
	   //h_ = 0.9*h_*(htilde_);
	}
	
	
	 // if(h_ < hmin_)  //beware of too small stepsizes
	 //   h_ = hmin_;

	
	//Functional iteration
	U_ = Uprev_;
       
	for(int l = 1; l <= maxit; ++l){
	  Ulp = U_;
	  U_  = Uprev_ + h_*ode_((U_ + Uprev_)/2.);
	  r_ = (U_ - Uprev_)/h_ - ode_((U_ + Uprev_)/2.);
	  nm = norm(r_);
	  
	  if(nm <= eps_){
	    flag = true;
	    //Uprev_ = U_;
	    break;
	  }
	  if(l == maxit-1){
	    nmprev = nm;
	  }

	  if(l == maxit){ //no convergence detected
	    flag = false;
	    std::cout << "FUNCTIONAL ITERATION DID NOT CONVERGE WITHIN "<< maxit << " ITERATIONS..." << std::endl;
	  }
	}
	Uprev_ = U_;

	if(!flag){
	  std::cout << "nm = "<< nm << "    nmprev = "<< nmprev << std::endl;
	  L = 2./h_*(nm/nmprev);
	  int m = (int)std::ceil(log(h_*L)/log(2));//(int)ceil(Absval(log2(h_*L)));
	  value_type H = c_/L;   //compute new stepsize
	  //value_type H = h_/2.;
	  std::cout << "L = "<< L << "    m = "<<m << "    H = "<< H << std::endl;
	  // if(H < hmin_)
	  //   H = hmin_;
	  
	  // value_type H = 1.e-08;
	  // int m = 5;
	  for(int s = 1; s <= m; ++s){  //advance m steps with new stepsize H
	    U_ += H*ode_(U_);
	  }
	  time += m*H;
	  k += m;
	}
	else{ //has converged
	  time += h_;
	  k++;
	}
	if(time >= tend)
	  break;

      } //end while
      of.close();
      return U_;
    }
    

  private:
    MODEL& ode_;  //reference
    value_type TOL_, eps_, sigma_,c_, hmin_,
      htilde_, h_;
    VType U_, Uprev_, R_, r_;
  };



}//end namespace

#endif
