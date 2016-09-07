#include <fstream>
#include <string>
#include "lebnum.hh"

#include "model.hh"

template<class T>
class HopfEasy{
public:
  typedef T value_type;
  typedef std::size_t size_type;
  typedef lebnum::Vector<T> VType;

  HopfEasy(const T& mu):rhs_(2),mu_(mu){}

  size_type domain_dim() const {return 2;}
  size_type range_dim() const {return 2;}

  VType& operator()(const VType& u){
    rhs_[0] = -u[1] + mu_*u[0] + u[0]*u[1]*u[1];
    rhs_[1] = u[0] + mu_*u[1] - u[0]*u[0];
    
    return rhs_;
  }

  void change_mu(const T& muNew){mu_ = muNew;}

  private:
    VType rhs_;
    T mu_;
};

using namespace std;
using namespace lebnum;
int main(int argc, char** argv){
  
  //! the executable will get an argument which corresponds to the different 
  //! parts of the exercise sheet
  if(argc != 2){
    cerr << "***** Error: "<<argv[0] << " takes ONE input argument" << endl;
    abort();
  }

  int example = atoi(argv[1]);
  
  if(example == 0){
    ofstream ld("lambda_mu.dat");
    ld << "# mu      imag   imag(conj)"<<endl;
    for(double mu = -1.5; mu <= 1.5; mu += 0.05){
      ld << mu << "  "<< +1 << "   " << -1 << endl;
    }
    ld.close();

    int isys = system(("./lamdaplot.sh lambda_mu.dat"));
    cout << "System status: "<< isys <<endl;


    double t0 = 0.,
      tend = 100.,
      dt = 0.1;
    double Mu = -0.2;
    typedef HopfEasy<double> ModelType;
    ModelType F(Mu);

    ClassicalRK<ModelType> RK(F);
    Vector<double> init(2);
    init = 1., 1;
    Vector<double> U(init);

    double time(t0);
    ofstream hf("hopfeasy.dat");
    size_t count(0);
    while(time < tend){
      cout << count++ << ".)     dt = " << dt << "    time = "  << time << endl;  
      hf << setprecision(14) << time << "  " << U[0] << "  " << U[1] << endl;
      U = RK.step(dt,U);
      time += dt;
    }
    hf.close();
    string addstr = num2str(Mu); //+"  IV: " + num2str(init[0])+", "+num2str(init[1]);
    int syst = system(("./hopf.sh hopfeasy.dat " +addstr).c_str());
    cout << "system status: "<< syst << endl;
  }

  else if(example == 1){
    int syst = system("./brussnullclines.sh");
    cout << "system status: "<< syst << endl;

    double a = 1, b = 2.5;
    typedef Brusselator<double> ModelType;
    ModelType Bruss(a,b);

    Vector<double> eval(2), between1(2), onydot1(2), under(2), left(2);
    eval[0] = 0.9;
    eval[1] = 5;

    between1[0] = 1.3;
    between1[1] = 2.4;
 
    onydot1[0] = 1.3183;
    onydot1[1] = 2.27567;

    under = 1,1;
    left = 0.4, 4.5;
    cout << "Bruss(above) = "<< Bruss(eval) << endl;
    cout << "Bruss(between1) = "<< Bruss(between1) << endl;
    cout << "Bruss(onydot1) = "<< Bruss(onydot1) << endl;
    cout << "Bruss(under) = "<< Bruss(under) << endl;
    cout << "Bruss(left) = "<< Bruss(left) << endl;
 

    cout << endl<< "BRUSSELATOR for fixed a and changing b"<<endl;
    double t0 = 0.,
      tend = 100,
      dt = 0.01;

    
    ClassicalRK<ModelType> RK(Bruss);
   
    Vector<double> init(2);
    init = 1., 1;
    Vector<double> U(init);

    double time(t0);
    ofstream hf("bruss.dat");
    hf << "# time       x      y     nullclineX     nullclineY"<<endl;
    size_t count(0);
    while(time < tend){
      cout << count++ << ".)     dt = " << dt << "    time = "  << time << endl;  
      hf << setprecision(14) << time << "  " << U[0] << "  " << U[1] << "  " <<
	Bruss.nullclineX(U[0]) << "  " << Bruss.nullclineY(U[0]) << endl;
      U = RK.step(dt,U);
      time += dt;
    }
    hf.close();

    syst = system(("./oszillator_bruss.sh bruss.dat " + num2str(a)+" "+num2str(b)).c_str());
    cout << "system status: "<< syst << endl;
  }

  //! default case
  else{
    cerr << "****ERROR: executable not defined yet"<<endl;
    abort();
  }

} //end main
