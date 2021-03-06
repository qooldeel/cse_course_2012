#include <fstream>
#include <string>
#include "lebnum.hh"

#include "model.hh"

using namespace std;
using namespace lebnum;
int main(){

  double sigma = 10.,
    beta = 8./3.,
    rho = 28.;

  typedef LorenzEquations<double> ModelType;
  
  ModelType Lorenz(sigma,rho,beta);
  double t0 = 0.,
      tend = 100.,
      dt = 0.01;
  ClassicalRK<ModelType> RK4(Lorenz);

  //!starting point
  Vector<double> U(Lorenz.range_dim());
  //U = 1.,1.,1.;   //Jochen
  // U = 0., 1., 0.;//wiki
  U = 1.5,2.,3.; //Peter Bastian
  

  double time(t0);
  ofstream hf("LorenzAttractor.dat"); 
  size_t count(0);
  while(time < tend){
    cout << count++ << ".)     dt = " << dt << "    time = "  << time << endl;  
    hf << setprecision(14) << time << "  " << U[0] << "  " << U[1] << "  "<< U[2]<< endl;
    U = RK4.step(dt,U);
    time += dt;
  }
  hf.close();
  
  cout << "integrated LORENZ equations"<<endl;
  string add;
  int syst = system(("./lorenz.sh LorenzAttractor.dat "+add).c_str());
  cout << "system status: "<< syst << endl;

  return 0;
}
