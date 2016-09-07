#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>

/**
 * Konsole:      1.) Oeffne Terminal 
 * Compilation:  2.) c++ -o trapez trapez.cc
 * Ausfuehren:   3.) ./trapez
 * Graphik:      4.) gnuplot
 *               5.) plot 'FreierFall.dat' using 1:2 with lines lw 2    
*/
template<class T>
inline T fun(const T& v, const T& m , const T& g, const T& k){
  return (g - k/m*v*v);
}

template<class T>
inline T exact(const T& t, const T& m , const T& g, const T& k){
  return ( std::tanh(t*std::sqrt(k*g/m)+1)*std::sqrt(g)/std::sqrt(k/m) );
} 

using namespace std;
int main(){
  
  double m = 120,       //!Masse des Felix Baumgartner mit Raumanzug [kg]
    g = 9.81,           //!Erdbeschleunigung  [m/s²]
    k = 0.73,           //!Luftwiderstand     [kg/m] 
    t0 = 0.,            //!Anfangszeit (Sprung zur Erde)
    tend = 300.;        //!Endzeit (tend-t0 = so lange dauerte der freie Fall)
    
  //**** TODO: aendere Zeitschritt auf 0.1, 0.01,...
  double step = 1;     //!Zeitschritt
  
  double time(t0),      //!Zeitpunkt
    v(0.),              //!Diskrete Loesung, d.h. \f$ v_{n+1}\f$ 
                        //! Start: \f$ v_0 = 0 \f$
    v_old(0.),          //!Loesung zuvor, d.h. \f$ v_{n}\f$ 
    f(0.),             //! Speichere Auswertung der Rechten Seite
    err(0.),           //! Fehler bzgl. diskreter und exakter Loesung
    ex(0.);            //! exacte Lösung

  ofstream ofile;
  ofile.open("FreierFall.dat");
  ofile << "# Zeit       num. v    exakt. v   Fehler" << endl;
    
  while(time < tend){
    v_old = v;
    f = fun(v_old, m,g,k);
    v = v_old + step/2*(f + fun(v_old + step*f, m,g,k));
    
    //! das ist quasi nur fuer die Ausgabe
    ex = exact(time,m,g,k);
    err = std::abs(ex-v);
    cout<< setprecision(15) << "t = "<< time << "  v = " <<  v << " |v(t_i) - v_i| = " <<  err << endl;
    ofile << time << "  " << v << "  " << ex << "  " << err << endl; //! schreibe in File
    time += step;  //! inkrementiere Zeit um festen Zeitschritt
  }

  ofile.close(); 
}
