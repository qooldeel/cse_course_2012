#include <iostream>
#include <fstream>
#include <complex>
#include <cmath>

#include "lebnum.hh"


template<class T>
inline T discr(const T& mu){ return mu*mu - 4; }


using namespace std;
using namespace lebnum;

int main(){
  
  typedef pair<complex<double>,complex<double> > PairType;
  PairType lambda;

  string fname = "classific.dat";
  ofstream ofs(fname.c_str());
  ofs << "## mu      mu² - 4     Re(lambda1)     Im(lambda1)     Re(lambda2)    Im(lambda2)" << endl; 
  for(double mu = -2.5; mu <= 2.5; mu += 0.05){
    cout << "mu = "<< mu << endl;
    lambda = quadratic_solution(1.,-mu,1.);  //solve characteristic equation
    ofs << mu << "   " << discr(mu) << "   "  <<  "   "<< lambda.first.real() << "   "<< lambda.first.imag() <<"   "<< lambda.second.real() << "   "<< lambda.second.imag() << endl; 
  }
  ofs.close();

  int isys = system(("./spezi.sh " + fname).c_str());
  cout << "system's status: "<< isys << endl;

  //! complex solution
  // lambda = quadratic_solution(1.,-0.5,3.);
  // cout<< "lambda1 = "<< lambda.first << "     lambda2 = "<< lambda.second << endl;

}
