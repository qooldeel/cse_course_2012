#include "lebnum.hh"


/**
 * @ file compilation and execution via use of a 'Makefile' which automates
 * compiling and linking commands
 *  $ make clean;make; ./test
*/

//! faciliates readability by allowing to skip std::, lebnum::, etc.
using namespace std;
using namespace lebnum;

int main(){

  Vector<double> v(5);
  cout << "v.size() = " << v.size() << endl;
  // copy construction
  double af[] = {0.5,-2.15,0.75,0.85,-1.5};
  Vector<double> w(af,af+5);
  cout << "w = "<< w << endl;
  Vector<double> x(w);
  cout << "x = "<< x << endl;

  cout << "This is awsome:"<<endl;
  Vector<double> y(4), def;
  y = -3.5, 2., 0.5, 4.;
  def = y;
  cout << "def = "<< def << endl;

  Vector<double> one(1);
  one = 1;
  cout << "one = "<< one << endl;

  y *= 2;
  cout << "y *= scal: " << y << endl;

  Vector<double> V(4), W(4), X(4), Y(4);
  V = 1, 2,3,4;
  W = 0.5, -1, -6, 1.5;
  X = 0.5*V;  
  Y = V*0.75;

  cout << "X = "<< X << "   Y = "<< Y << endl;
  V += W;  //overwrites V
  cout << "V = V + W: "<< V << endl;
  V -= X; //overwrites V again
  cout << " V = V - Z" << V << endl;
  Vector<double> Z =  X + Y, A = X - Y;
  cout << "Z = X + Y: "<<  Z   << endl;
  cout << "A = X - Y: "<<  A   << endl;
  
  cout << "X +2.*Y - 3.75*V: " << X +2.*Y - 3.75*V << endl;

  -X;
  cout << "neg(X) = "<< X << endl;
   
  
  cout << "Test various norms:"<<endl;
  double d1 = 0.5, d2 = -0.75, d3 = 1.5;
  std::complex<double> z1(3.,-4), z2(-0.5,-0.75), z3(1.5,2.);
  
  Vector<double> vn(3);
  vn[0] = d1; vn[1] = d2; vn[2] = d3;
  Vector<complex<double> > vz(3);
  vz[0] = z1; vz[1] = z2; vz[2] = z3;

  cout << "1-norm(vn) = "<< norm(vn,'1') << "   1-norm(vz) = "<< norm(vz,'1') << endl;
  cout << "2-norm(vn) = "<< norm(vn) << "   2-norm(vz) = "<< norm(vz) << endl;
  cout << "inf-norm(vn) = "<< norm(vn,'i') << "   inf-norm(vz) = "<< norm(vz,'i') << endl;
  

  
  cout << endl << "Test some matrix operations:"<<endl;
  Matrix<double> M1(2,3), M2(2,3);
  M1 = 1, 2, 3,
       4, 5, 6; 

  M2 = -0.5, -1, 4,
    1.75, 0.25, -1.5;

  Matrix<double> M3 = M1 + M2/4.;
  cout << "M1 = " << M1 << endl;
  cout << "M3 = " << M3 << endl;
  

  Matrix<double> M1prime;
  M1prime = M1.transpose();
  cout << "M1prime.rows() = "<< M1prime.rows() << "   M1prime.cols() = " << M1prime.cols() << endl;
  cout << "M1prime = "<< M1prime << endl;


  Matrix<double> M4(4,3), M5(3,2);
  M4 = 0.5,  0.75,  -0.24,  3.15,
       0.01, -2.15, -3.05,  0.8,
    4.34, -6.76,  2.096,   -1.67;

  M5 = 0.8, -1.9,
    1.,  0,
    4.5, -7.8;

  cout << "M4 = "<< M4 << endl;
  cout << "M5 = "<< M5 << endl;

  Matrix<double> Mult(M4*M5);
  cout << "Mult.rows() = "<< Mult.rows() << "  Mult.cols() = "<<Mult.cols() <<endl; 
  cout << "M4*M5 = "<< Mult << endl;
       

  Vector<double> bv(3), cv(4);
  bv = -3.65, -1.5, 2.3;
  cv = 0.75, -0.98, 1.45, 0.55;

  cout << "M4.cols() = "<< M4.cols() << "  bv.size() = "<< bv.size() << endl;   
  cout << "M4*bv = "<< M4*bv <<endl;
  Vector<double> vM = cv*M4;
  cout << "cv*M4 = "<< vM << endl;
  

  cout << "1-norm of M4: "<< norm(M4,'1') << endl;
  cout << "inf-norm of M4: "<< norm(M4,'i') << endl;
  cout << "frob-norm of M4: "<< norm(M4,'F') << endl;

  complex<double> zz1(3,-0.5), zz2(0.5,0.75), zz3(-1,-2),
    zz4(8,-4.65), zz5(2,-3), zz6(-1.25,3);

  Matrix<complex<double> > AZ(3,2);
  AZ = zz1, zz2, zz3,
    zz4, zz5, zz6;

  cout << "AZ = "<< AZ << endl;

  cout << "AZ' = "<< AZ.transpose() << endl;

  cout << "Test operator [i][j] on matrices" << endl;
  for(size_t i = 0; i < 4; ++i){
    for(size_t j = 0; j < 3; ++j){
      if(i==j)
	M4[i][j] = 1000;
      cout << M4[i][j] << " ";
    }
    cout << endl;
  }

  Matrix<double> AA(3,3);
  AA = 1, 4, 7,     //see Golub and Van Loan, Ex. 3.2.2, p. 99
    2,  5, 8,
    3, 6, 10;

  cout << "AA = "<< AA << endl;

  Vector<double> bb(3);

  bb = -1.5, 0.5, -4.85;

  cout << "solution = "<< AA.solve(bb) << endl;
 
  cout << endl << "Complex matrices:"<< endl;
  Matrix<complex<double> > MCpx(3,3);
  
  MCpx = z1, z2, z3,
    zz1, zz2, zz3,
    zz4, zz5, zz6;

  cout << "MCpx = "<< MCpx <<endl;
  
  complex<double> c1(0.5,-1),c2(-1.45,3),c3 = 5;
  Vector<complex<double> > bz(3);
  bz = c1,c2,c3;
  Vector<complex<double> > xz = MCpx.solve(bz);
  cout << "complex solution = "<< xz << endl;


  cout << endl << "Test dot product:"<<endl;
  Vector<double> WD1(3),WD2(3);

  WD1 = 0.5,-1.45,4.5;
  WD2 = -0.75,3.25,-0.25;

  cout << "WD1'*WD2 = "<< dot(WD1,WD2) << endl;
  
  Vector<complex<double> > ZV1(3), ZV2(3);
  ZV1 = z1, z2, z3;
  ZV2 = zz4, zz1, zz3;

  cout << "ZV1 = "<< ZV1 << endl;
  cout << "ZV2 = "<< ZV2 << endl;
  

  cout << "ZV1'*ZV2 = "<< dot(ZV1,ZV2) << endl; 
  

  Matrix<double> MUE3(3,3);
  Vector<double> VUE3(3);

  MUE3 = 6, -4, 7,
    -12, 5, -12,
    18, 0, 22;

  VUE3 = 41./12., -22./3., 29./2.;

  cout << "solution = "<< MUE3.solve(VUE3) << endl;
 
  pair<complex<double>,complex<double> > pay = quadratic_solution(1.,2.,1.);
  cout << "Quadrat. solution:  x1 = "<< pay.first << "    x2 = "<< pay.second << endl;
 // cout << endl << "Test me" <<endl;
  
  // double do1 = -1.24;

  // cout << NaturalPow(do1,3) << endl;
  // cout << NaturalPow(do1,0) << endl;
  // cout << NaturalPow(do1,-3) << endl;

  // cout << "complex" <<endl;
  // complex<double> zo1(3,-4);
  // cout << NaturalPow(zo1,3) << endl;
  // cout << NaturalPow(zo1,0) << endl;
  // cout << NaturalPow(zo1,-3) << endl;
}
