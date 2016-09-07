#include <fstream>
#include <string>
#include "lebnum.hh"

#include "model.hh"

//============= TODO =========================================================
#define ODE_INSTANCE 2  //select ode instance at compile time
//============================================================================

using namespace std;
using namespace lebnum;
int main() {

//============= TODO =========================================================
  const int choose = 2;   //1 = impl. Euler, 2 = expl4stiff, 3 = both
//============================================================================

  int dim;
  Vector<double> u0;
  double t0,tend, hestim;
  string fname, fname2, add, add2;

#if ODE_INSTANCE == 1
  typedef ChemicalReactionRobertson<double> OdeType;
  OdeType reaction;
  dim = 3;
  u0.resize(dim);
  u0 = 1.,0.,0.;

  t0 = 0.;
  tend = 0.3;
  hestim = 1.25e-03;
  fname = "ImplicitEulerRobertson.dat";
  add = "Robertson";
  fname2 = "E4StiffRobertson.dat";
  add2 = "Robertson-reaction system solved with explicit method";
#endif

#if ODE_INSTANCE == 2
  typedef HIRESProblem<double> OdeType;
  OdeType reaction;
  dim = 8;
  u0.resize(dim);
  u0 = 1.,0.,0.,0.,0.,0.,0.,0.0057;

  t0 = 0.;
  tend = 321.8122;
  hestim = 0.0045;
  fname = "IE_HIRES.dat";
  add = "HIRES";
  fname2 = "E4SHIRES.dat";
  add2 = "HIRES explicit";
#endif



  if(choose == 1 || choose == 3){
    double stpsztol = 1.e-05;
    ImplicitEuler<OdeType> IE(reaction,1.e-06,stpsztol);
    IE.integrate(t0,tend,u0,hestim,5,fname);
    int sys = system(("./versus.sh "+fname+" "+ num2str(dim)+" "+add).c_str());
    cout << "Status: "<< sys << endl;
  }
 
  if(choose == 2 || choose == 3){
    double hmin = 1.e-06;
    int maxit = 20;
    Explicit4Stiff<OdeType> E4S(reaction,1.e-04,1.e-04,0.65,0.9995,hmin);
    E4S.integrate(t0,tend,u0,1.e-05,maxit,fname2);
    int sys = system(("./versus.sh "+fname2+" "+ num2str(dim)+" "+add2).c_str());
    cout << "Status: "<< sys << endl;
  }
  
  cout << "Method used for the integration the stiff system:"<<endl;
  switch(choose){
  case 1: cout << " IMPLICIT EULER" <<endl; break;
  case 2: cout << " EXPLICIT4STIFF" <<endl; break;
  case 3: cout << " both methods" <<endl; break;
  default: cout << " method not defined" <<endl;    
  }
  
  

 return 0;
}
