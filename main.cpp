#include <iostream>
#include <matplot/matplot.h>
#include <cmath>
#include <iomanip>
#include <complex>
using namespace std;
namespace mp = matplot;
vector<complex<double>> dft(vector<double> sygnal) {
    vector<double> Re,Im;
    vector<complex<double>> complex_sygnal;
    for(double wartosc : sygnal) {
        complex_sygnal.emplace_back(wartosc, 0.0);
    }
    vector<complex<double>> output;
    complex<double> suma;
    int N = sygnal.size();
    
    for(int f = 0; f < N; f++) {
        suma = complex<double>(0.0, 0.0);
        for(int k = 0; k < N; k++) {
            double re = cos(((2*M_PI)/N)*k*f);
            double im = sin(((2*M_PI)/N)*k*f);
            complex<double> zeta (re, -im);
            suma += complex_sygnal[k] * zeta;
        }
        output.push_back(suma);
    }
    return output;
}
vector<double> cosinus(double amplituda, double czestotliwosc, double t_min, double t_max, double samples) {
    vector<double> wynik;
    // wynik.resize(samples);
    vector<double> os_x = mp::linspace(t_min, t_max, samples);
    for(int i = 0; i < samples; i++) {
        cout << setw(5) << os_x[i] << endl;
    }
    for(int i = 0; i < samples; i++) {
        wynik.push_back(amplituda * cos(2 * M_PI * os_x[i] * czestotliwosc));
    }
    return wynik;
}
vector<double> pila(double amplitude, double czestotliwosc, double t_min, double t_max, double samples, double faza) {
    vector<double> os_x = mp::linspace(t_min, t_max, samples);
    vector<double> wynik;
    for(int i = 0 ; i < samples; i++) {
        wynik.push_back(2 * (((os_x[i] + faza)*czestotliwosc) - floor(0.5 + ((os_x[i] + faza)*czestotliwosc))));
    }
    if(faza == 0) {
        mp::plot(os_x, wynik);
        mp::show();

    }
    return wynik;
}
vector <double> prostokatny(double amplitude, double czestotliwosc, double t_min, double t_max, double samples) {
    vector<double> os_x = mp::linspace(t_min, t_max, samples);
    vector<double> wynik;
    for(int i = 0 ; i < samples; i++) {
        if(sin(2 * M_PI * czestotliwosc * os_x[i]) >= 0) wynik.push_back(1);
        else wynik.push_back(0);
    }
//   double okres = 1/czestotliwosc;
//   double polowa_okresu = okres/2;
//   for (int i = 0; i<samples; i++){
//   double obecnyCzas = fmod(os_x[i] ,okres);
//   if(obecnyCzas<polowa_okresu) wynik.push_back(amplitude);
//   else wynik.push_back(0);
//  }
    return wynik;
}
vector<double> sinus(double amplituda, double czestotliwosc, double t_min, double t_max, double samples) {
    vector<double> wynik;
    // wynik.resize(samples);
    vector<double> os_x = mp::linspace(t_min, t_max, samples);
    for(int i = 0; i < samples; i++) {
        cout << setw(5) << os_x[i] << endl;
    }
    for(int i = 0; i < samples; i++) {
        wynik.push_back(amplituda * sin(2 * M_PI * os_x[i] * czestotliwosc));
    }
    for(int i = 0; i < samples; i++) {
        // cout<< setw(5)<< os_x[i];
        // cout<<wynik[i]<<endl;
    }
    return wynik;
}
int main() {
    vector<complex<double>> poDft = dft(sinus(1,20,0,100,100));
    vector<double> rePoDft, imPoDft, magPoDft;
    for (int i = 0; i<poDft.size();i++){
        rePoDft.push_back(poDft[i].real());
        imPoDft.push_back(poDft[i].imag());
        magPoDft.push_back(sqrt(imPoDft[i]*imPoDft[i] + rePoDft[i] *rePoDft[i]));
    }
    vector<double> os;
    for (int i = 0; i<sinus(1,20,0,100,100).size();i++){
      os.push_back(i);
    }
    mp::plot(os,magPoDft);
    mp::show();
    return 0;
    

}
