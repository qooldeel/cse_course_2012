#include <fstream>
#include <string>
#include "lebnum.hh"

/**
 * Solve y'' + 4y' + 3y = 0, y(0) = 2, y'(0) = -4
*/
template<class T>
class SecondOrderODE{
public:
  typedef T value_type;
  typedef std::size_t size_type;
  typedef lebnum::Vector<T> VType;

  SecondOrderODE():rhs_(2){}

  size_type domain_dim() const {return 2;}
  size_type range_dim() const {return 2;}

  VType& operator()(const VType& u){
    rhs_[0] = u[1];
    rhs_[1] = -3*u[0] - 4*u[1];
    
    return rhs_;
  }


  private:
    VType rhs_;
};

using namespace std;
using namespace lebnum;
int main() {
  
  typedef SecondOrderODE<double> ModelType;
  ModelType dotx;
  //ClassicalRK<ModelType> RK4(dotx);
  RungeKutta4<ModelType> RK4(dotx);
  Vector<double> x(2);
  x = 2.,-4.;

  double time(0), tend(2), dt = 0.01;
  ofstream hf("Test.dat"); 
  size_t count(0);
  while(time < tend){
    cout << count++ << ".)     dt = " << dt << "    time = "  << time << endl;  
    hf << setprecision(14) << time << "  " << x[0] << "  " << x[1] << "  "<< endl;
    x = RK4.step(dt,x);
    time += dt;
  }
  hf.close();
  
  cout << "Classical RK"<<endl;
  string add;
  int syst = system(("./plottest.sh Test.dat "+add).c_str());
  cout << "system status: "<< syst << endl;


  return 0;
}
