#ifndef MODELS_FOR_ODES_HH
#define MODELS_FOR_ODES_HH

#include <iostream>
#include <cassert>


namespace lebnum{

  template<class T>
  class DavisSkodje{
  public:
    typedef T value_type;
    typedef std::size_t size_type;
    typedef Vector<T> VType;

    //!default gamma is taken to be 60, but feel free to use other vals 
    DavisSkodje(const T& gamma = 60):rhs_(2),gamma_(gamma){}

    size_type domain_dim() const {return 2;}
    size_type range_dim() const {return 2;}

    VType& operator()(const VType& y){
      rhs_[0] = -y[0];
      rhs_[1] = -gamma_*y[1] + ((gamma_-1)*y[0] + gamma_*y[0]*y[0])/((1+y[0])*(1+y[0]));
      return rhs_;
    }
    

  private:
    VType rhs_;
    T gamma_;
  };


  //!examples: mechanical clock or vibrating string
  template<class T>
  class VanDerPolOscillator{
  public:
    typedef T value_type;
    typedef std::size_t size_type;
    typedef Vector<T> VType;

    VanDerPolOscillator(const T& mu):rhs_(2),mu_(mu){}

    size_type domain_dim() const {return 2;}
    size_type range_dim() const {return 2;}

    VType& operator()(const VType& y){
      rhs_[0] = y[1];
      rhs_[1] = mu_*(1 - y[0]*y[0])*y[1] - y[0];
      
      return rhs_;
    }
      

  private:
    VType rhs_;
    T mu_;
  };


  //! predetor-prey model: a simple but not very realistic example of 
  //! two rivaling species
  template<class T>
  class LotkaVolterra{
  public:
    typedef T value_type;
    typedef std::size_t size_type;
    typedef Vector<T> VType;

    LotkaVolterra(const T& a, const T& b, const T& c, const T& d):rhs_(2),alpha_(a),beta_(b),gamma_(c),delta_(d){}

    size_type domain_dim() const {return 2;}
    size_type range_dim() const {return 2;}


    VType& operator()(const VType& y){
      rhs_[0] = y[0]*(alpha_ -beta_*y[1]);     //!rate of number of rabbits
      rhs_[1] = -y[1]*(gamma_ - delta_*y[0]);  //!rate of number of foxes
      
      return rhs_;
    }

  private:
    VType rhs_;
    T alpha_,beta_,gamma_,delta_;
  };

  
  /**
   * \brief Brusselator: an easy example of an oscillating chemical reaction
   *
   * This model has been taken from
   *
   *  [S.H.STROGAT, "Nonlinear Dynamics and Chaos", Addison-Wesely,1994,p. 290]
   */
  template<class T>
  class Brusselator{
  public:
    typedef T value_type;
    typedef std::size_t size_type;
    typedef Vector<T> VType;

    Brusselator(const T& a, const T& b):rhs_(2),a_(a),b_(b){}

    size_type domain_dim() const {return 2;}
    size_type range_dim() const {return 2;}
    
    VType& operator()(const VType& w){      
      rhs_[0] = 1 - (b_+1)*w[0] + a_*w[0]*w[0]*w[1];
      rhs_[1] = b_*w[0] - a_*w[0]*w[0]*w[1];

      return rhs_;
    }
    
    void change_a(const T& a){
      a_ = a;
    }

    void change_b(const T& b){
      b_ = b;
    }

    T nullclineX(const T& x) const{
      return ((b_+1)*x -1)/(a_*x*x);
    }

    T nullclineY(const T& x) const{
      return b_/(a_*x);
    }
    
  private:
    VType rhs_;
    T a_,b_;  
  };

  /**
   * \brief The famous Lorenz attractor which displays some chaotic features
   */
  template<class T>
  class LorenzEquations{
  public:
    typedef T value_type;
    typedef std::size_t size_type;
    typedef Vector<T> VType;

    //! input in order of occurence in Lorenz system
    //! if \f$ \rho < 1\f$ there is exactly 1 equilibrium point which is at the 
    //! origin
    //! if \f$ \rho = 1\f$ we have a saddle node bifurcation 
    //! if \f$ \rho > 1\f$ 2 additional critical points appear
    //! for \f$ \sigma = 10, \rho = 28, \beta = \frac{8}{3}\f$, chaotic 
    //! solutions emerge.
    LorenzEquations(const T& s, const T& r, const T& b):rhs_(3),sigma_(s),rho_(r),beta_(b){
      assert(s > T() && r > T() && b > T());
    }

    size_type domain_dim() const {return 3;}
    size_type range_dim() const {return 3;}


    VType& operator()(const VType& w){
      rhs_[0] = sigma_*(w[1]-w[0]);
      rhs_[1] = w[0]*(rho_-w[2])-w[1];
      rhs_[2] = w[0]*w[1] - beta_*w[2];
      
      return rhs_;
    }

  private:
    VType rhs_;
    T sigma_,rho_,beta_;
  };


  
  /**
   * \brief Example from plant physiology, taken from [1]
   *
   *  [1] [ERIKSSON, JOHNSON and LOGG, "Explicit Time-Stepping for Stiff ODEs"]
   */
 template<class T>
 class HIRESProblem{
  public:
   typedef T value_type;
   typedef std::size_t size_type;
   typedef Vector<T> VType;
   typedef Matrix<T> MtxType;
    
   HIRESProblem():rhs_(8),jac_(8,8){}

    size_type domain_dim() const {return 8;}
    size_type range_dim() const {return 8;}

   VType& starting_point(){
     rhs_[0] = 1.;
     rhs_[1] = rhs_[2] = rhs_[3] = rhs_[4] = rhs_[5] = rhs_[6] = 0.;
     rhs_[7] = 0.0057;
   }

   MtxType& jacobian(const VType& u){
     jac_ = -1.71,  0.43, 8.32, 0, 0, 0, 0, 0,
       1.71, -8.75, 0,    0, 0, 0, 0, 0,
       0,0, -10.03, 0.43, 0.035, 0, 0, 0,
       0, 8.32, 1.71, -1.12, 0, 0, 0,0,
       0,0,0,0,-1.745, 0.43, 0.43, 0,
       0,0,0,0.69, 1.71, -280.*u[7], 0.69, -280.*u[5],
       0,0,0,0,0,280.*u[7],-1.81, 280*u[5],
       0,0,0,0,0,-280.*u[7],1.81, -280.*u[5];

     return jac_;
   }


    VType& operator()(const VType& u){
      rhs_[0] = -1.71*u[0] + 0.43*u[1] + 8.32*u[2] + 0.0007;
      rhs_[1] = 1.71*u[0] - 8.75*u[1];
      rhs_[2] = -10.03*u[2] + 0.43*u[3] + 0.035*u[4];
      rhs_[3] = 8.32*u[1] + 1.71*u[2] - 1.12*u[3];
      rhs_[4] = -1.745*u[4] + 0.43*u[5] + 0.43*u[6];
      rhs_[5] = -280.*u[5]*u[7] + 0.69*u[3] + 1.71*u[4] - 0.43*u[5] + 0.69*u[6];
      rhs_[6] = 280.*u[5]*u[7] - 1.81*u[6];
      rhs_[7] = -280.*u[5]*u[7] + 1.81*u[6];
      
      return rhs_;
    }

  private:
   VType rhs_;
   MtxType jac_;
  };



  /**
   * Robertson's simple chemical reaction, see. [1] 
   *  [1, HAIRER,NORSETT,WANNER, Vol. II, p.3]
   *
   * Starting point: [1,0,0]
   */
  template<class T>
  class ChemicalReactionRobertson{
  public:
    typedef T value_type;
    typedef std::size_t size_type;
    typedef Vector<T> VType;
    typedef Matrix<T> MtxType;

    ChemicalReactionRobertson():rhs_(3),jac_(3,3){}

    size_type domain_dim() const {return 3;}
    size_type range_dim() const {return 3;}

    VType& starting_point() {
      rhs_[0] = 1.;
      rhs_[1] = 0.;
      rhs_[2] = 0.;
      return rhs_;
    }

    VType& operator()(const VType& y){
      rhs_[0] = -0.04*y[0] + 1.e+04*y[1]*y[2];
      rhs_[1] =  0.04*y[0] - 1.e+04*y[1]*y[2] - 3.e+07*y[1]*y[1];
      rhs_[2] = 3.e+07*y[1]*y[1];
      return rhs_;
    }
    

    MtxType& jacobian(const VType& y){
      jac_[0][0] = -0.04; 
      jac_[0][1] = 1.e+04*y[2]; 
      jac_[0][2] = 1.e+04*y[1];
      
      jac_[1][0] = 0.04; 
      jac_[1][1] = -1.e+04*y[2]-6.e+07*y[1]; 
      jac_[1][2] = -1.e+04*y[1];
      
      jac_[2][0] = 0.; 
      jac_[2][1] = 6.e+07*y[1]; 
      jac_[2][2] = 0.;
      
      return jac_;
    }

  private:
    VType rhs_;
    MtxType jac_;
  };



} //end namespace

#endif
