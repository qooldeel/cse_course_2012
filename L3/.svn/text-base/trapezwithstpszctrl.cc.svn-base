#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <limits>

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


using namespace std;
int main(){
  
  double m = 120,       //!Masse des Felix Baumgartner mit Raumanzug [kg]
    g = 9.81,           //!Erdbeschleunigung  [m/s²]
    k = 0.73,           //!Luftwiderstand     [kg/m] 
    t0 = 0.,            //!Anfangszeit (Sprung zur Erde)
    tend = 300.;        //!Endzeit (tend-t0 = so lange dauerte der freie Fall)
    
  //**** TODO: aendere Zeitschritt auf 0.1, 0.01,... *************************
  double step = 2.,     //!Zeitschritt
    tol = 1.e-02;
  //**************************************************************************

  double time(t0),      //!Zeitpunkt
    v(0.),              //!Diskrete Loesung, d.h. \f$ v_{n+1}\f$ 
                        //! Start: \f$ v_0 = 0 \f$
    v_old(0.),          //!Loesung zuvor, d.h. \f$ v_{n}\f$ 
    f(0.),             //! Speichere Auswertung der Rechten Seite
    vlow(0.),          //! Methode geringerer Ordnung (hier: 1. Ordnung)
    discrerr(0.);      //! Fehler zw. diskr. Loesungen verschiedener Ordnungen

  ofstream ofile;
  ofile.open("Schrittweitensteuerung.dat");
  ofile << "# i     Zeit    dt     v_i" << endl;
    
  int count = 0;  //!Zaehle Schritte
  while(time < tend){
    v_old = v;
    f = fun(v_old, m,g,k);
    vlow = v_old + step*f;                     //!Verfahren 1. Ordnung
    v = v_old + step/2*(f + fun(vlow, m,g,k)); //!Verfahren 2. Ordnung
    
    //! Schrittweitensteuerung
    discrerr = std::abs(v - vlow);
    // step = step*std::pow((std::numeric_limits<double>::epsilon()*std::abs(vlow)/step)/discrerr,0.5);
    step = step*(std::pow(tol/discrerr,0.5));

    //! das ist nur fuer die Ausgabe
    cout<< setprecision(15) << count << ".)   t = "<< time << "   dt = "<< step << "  v = " <<  v << endl;
    ofile << setprecision(15) << count << "  "<< time << "  "<< step << "  " <<  v << endl;  //schreib's auch in 'ne File
    count++; //!Inkrementiere Schrittzaehler

    time += step;  //! inkrementiere Zeit um festen Zeitschritt
    
    //! o.k. bastle etwas herum, so dass wir auch den exakten Endpunkt erwischen
    // if (time >= tend)
    //   break;
    // else if (time+step > tend)
    //   step = tend - time;

  }

  ofile.close(); 
}
