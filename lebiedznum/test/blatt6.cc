#include <fstream>
#include <string>
#include "lebnum.hh"

#include "model.hh"


//! exact solution of the Davis Skodje model in dependence of c1, c2 and gamma
template<class T>
class ExactSolutionOfDavisSkodje{
public:
  typedef lebnum::Vector<T> VType; 

  ExactSolutionOfDavisSkodje(const T& c1 = T(), const T& c2 = T(), const T& gamma = T()):sol_(2),c1_(c1),c2_(c2),gamma_(gamma){}

  VType& operator()(const T& t){
    sol_[0] = c1_*std::exp(-t);
    sol_[1] = c2_*std::exp(-gamma_*t) + c1_/(c1_ + std::exp(t));
    return sol_;
  }

  const T& operator[](std::size_t i) const {return sol_[i];}

private:
  VType sol_;
  T c1_,c2_,gamma_; 
};


//! Analytical Jacobian of Davis Skodje. For the first entry in
//! the second row, the quotient rule for derivatives has been used. 
//! Check it out!
//!
//!  --                   --
//! | -1                  0 |
//! |                       | 
//! | (1+g)y1 + g -1        |      where g abbreviates gamma 
//! | ---------------    -g |
//! |   (1 + y1)³           |
//!  --                   --  
template<class T>
class JacobianOfDavisSkodje{
public:
  typedef lebnum::Matrix<T> MtxType;

  JacobianOfDavisSkodje(const T& gamma):jac_(2,2),gamma_(gamma){}

  template<class V>
  MtxType& operator()(const V& x){
    T df2dx1 = ((1+gamma_)*x[0] + gamma_ -1)/lebnum::Cub(1+x[0]);
    jac_(0,0) = -1.;     jac_(0,1) = 0.;
    jac_(1,0) = df2dx1;  jac_(1,1) = -gamma_;
    
    return jac_;
  }

private:
  MtxType jac_;
  T gamma_;
  
};

//! If you don't want to  write here 'std::cout' or 'lebnum::Vector<double>' 
//! all the time, just tell the compiler to use the following namespaces 
using namespace std;
using namespace lebnum;

// Here you go...
//MAIN FILE · MAIN FILE · MAIN FILE · MAIN FILE · MAIN FILE · MAIN FILE ·
int main(int argc, char** argv){ 

  //! the executable will get an argument which corresponds to the different 
  //! parts of the exercise sheet
  if(argc != 2){
    cerr << "***** Error: "<<argv[0] << " takes ONE input argument" << endl;
    abort();
  }

  int example = atoi(argv[1]);

  if(example == 0){
    cout << "Exercise 6.2 a) + b)"<<endl;
   
    const double gamma = 60;
    const double Tol = 1.e-08;


    Vector<double> U(2);
    U = 3., 1.5;

    //! c1, and c2 can be determined from the exact solution at t = 0 and 
    //! the initial values above
    double c1 = 3,  //since c1·1 = 3
      c2 = 0.75;    //since c2·1 + 3/(3+1) = 1.5
    ExactSolutionOfDavisSkodje<double> exact(c1,c2,gamma);

    double time = 0, 
      tend = 1.;

    double step = 3.e-02; //1.e-05; //0.0476;

    DavisSkodje<double> F(gamma);
    JacobianOfDavisSkodje<double> DF(gamma);
    
    string fexpl = "DS_explicit.dat",
      fimpl = "DS_implicit.dat";

    ofstream ofile(fexpl.c_str());
    ofile << "# Zeit         U[0]          U[1]       u[0]     u[1]"<<endl;
  
    ofstream ofileIM(fimpl.c_str());
    ofileIM << "# Zeit         Uim[0]          Uim[1]      u[0]     u[1]"<<endl;
  
    size_t count = 0;
    //! note that the the use of classes yields another advantage:
    //! we can nearly write 'mathematical' code (remember the usage of Vector
    //! and the functor DavisSkodje!
    // Vector<double> Utilde(2), f(2);

    size_t maxNewt = 5;
    Vector<double> Uimprev(U), Uim(Uimprev), b(2);
    Matrix<double> IterationMatrix(2,2);

    while(time < tend){
      exact(time);
      ofile << setprecision(14) << time << "  " << U[0] << "  " << U[1] << "  "<< exact[0] << "  "<< exact[1] << endl;
      ofileIM << setprecision(14) << time << "  " << Uim[0] << "  " << Uim[1] << "  "<< exact[0] << "  "<< exact[1] << endl;
      cout << count++ << ".)   dt = "<< step << "   time = " << time << endl; 
      //! f = F(U);
      //! Utilde = U + step*f;  //explicit Euler
      //! U += step/2.*(f + F(Utilde)); //explicit trapezoidal
    
      //! EXPLICIT EULER
      U += step*F(U);
   
  
      //!IMPLICIT EULER
      Uim = Uimprev;
      for(size_t l = 0; l < maxNewt; ++l){ //Newton iteration
	b = -(Uim - Uimprev - step*F(Uim));
	IterationMatrix = -step*DF(Uim);  
	IterationMatrix.add_unity_matrix();
	IterationMatrix.solve(b);
	Uim += b;    // b contains solution now
	if(norm(b) <= Tol)  //o.k. convergence within tolerance achieved
	  break;

	if(l == maxNewt-1){
	  cout << "+++++ WARNING: Newton iteration did not converge within "<< maxNewt << " iterations" << endl;
	  break;  //leave anyway
	}
      } // end of Newton iteration
      Uimprev = Uim;  //store current approximation

      time += step; // increment time
    }
  
    ofile.close();
    ofileIM.close();
  
    int sys = system(("./davisskodje.sh "+fexpl+" Explicit_Euler; ./davisskodje.sh "+fimpl+" Implicit_Euler").c_str());
    cout << "Status: "<< sys << endl;
      
    
  }

  else if(example == 1){
    /**
     *  fixed point: (0,0)
     *  classification:
     *   0 < mu < 2: unstable focus
     *  -2 < mu < 0: stable focus
     *   mu > 2: unstable node
     *   mu < -2: stable node
     *   mu = 2: unstable star
     *   mu = -2: stable star
     *   mu = 0: elliptic
     */
    cout << endl << "Exercise 6.1 a) VAN DER POL OSCILLATOR" <<endl;
    typedef VanDerPolOscillator<double> ODEModel;

  
    for(double mu = -2.5; mu <= 2.5; mu += 0.5){
      //const double mu = 1; //0.25;  //! damping factor
      ODEModel VDP(mu);

      ClassicalRK<ODEModel> RK4(VDP);
  
      double st = 0.1;          //!decrep.: Strogatz, p. 181
  
      double tend = 40, //30
	time = 0;
      Vector<double> U(2);
      U = 1, 1; //1., 0.01;
  
      size_t count = 0;
      string fname = "VanDerPol.dat";

      ofstream ofVan(fname.c_str());
      while(time < tend){
	cout << count++ << ".)     dt = " << st << "    time = "  << time << endl;  
	ofVan << setprecision(14) << time << "  " << U[0] << "  " << U[1] << endl;
	U = RK4.step(st,U);
	time += st;
      }
      ofVan.close();

      int syst = system(("./vanderpol.sh "+fname+" "+num2str(mu)+"; ./timevolution.sh "+fname).c_str());
      cout << "system status: "<< syst << endl;
    } //end mu-loop
  }

  else if(example == 2){
    cout << "Exercise 6.1 b)  LOTKA-VOLTERRA" <<endl;
    //Parameter for LV model
    double aval = 1.5, //1, //1.5, 
      bval = 1., //0.2,    //1,
      cval = 3.,  //0.04,   //3,
      dval = 0.9; // 0.5;    //1;

  
    typedef LotkaVolterra<double> ODE1;
    ODE1 LV(aval,bval,cval,dval);

    ClassicalRK<ODE1> RUKU(LV);
    //time horizon
    double tend = 30, //35;
      time = 0;
  
    //step size, if it gets nasty, try 0.01 etc
    double st = 0.1;

  
    //initial condition, first entry corresponds to prey, second to predator
    Vector<double> U(2);
    U = 4.5, 1.4; 
    //U = 5, 2;
    cout << "U = "<< U << endl;
  
    string filename = "LotkaVolterra.dat"; 

    ofstream ofLV(filename.c_str());  //convert to C-string, needed here
    Vector<double> dirvec(2);  // tangent, i.e. x' = f(x)
    double scal = 0.2;

    size_t count(0); //same as count = 0;
    ofLV << "# t       y1(t)     y2(t)    d/dty1(t)     d/dty2(t)" << endl;
    while(time < tend){
      cout << count++ << ".)     dt = " << st << "    time = "  << time << endl;  
      dirvec = LV(U);  //tangent, given by x' = f(x)
      ofLV << setprecision(14) << time << "  " << U[0] << "  " << U[1] << "   "<< scal*dirvec[0] << "   "<< scal*dirvec[1] << endl;
      //<< scal*(Absval(dirvec[0]-U[0]) << "   "<< scal*Absval(dirvec[1]-U[1]) << endl;
      U = RUKU.step(st,U);  //approximate solution
   
      time += st;  //iterate time 
    }
    ofLV.close();
 
    //write commands to shell
    int isys = system(("./lotkavolt.sh "+filename+"; ./direction.sh "+filename+"; ./timevolution.sh "+filename).c_str());  //again: convert to C-string

    cout << endl<< "Status: "<< isys << endl;
  }
  else{
    cerr << "****ERROR: executable not defined yet"<<endl;
    abort();
  }
}
