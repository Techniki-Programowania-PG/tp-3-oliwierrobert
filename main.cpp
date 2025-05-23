#include <iostream>
#include <matplot/matplot.h>
#include <cmath>
#include <iomanip>
#include <complex>
using namespace std;
namespace mp = matplot;
vector<complex<double>> dft(vector<double> sygnal, int samples) {
    vector<double> Re, Im;
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
            double re = cos(((2 * M_PI) / N) * k * f);
            double im = sin(((2 * M_PI) / N) * k * f);
            complex<double> zeta(re, im);
            suma += complex_sygnal[k] * zeta;
        }
        output.push_back(suma);
    }
    vector<double> rePoDft, imPoDft, magPoDft;
    for(int i = 0; i < output.size(); i++) {
        rePoDft.push_back(output[i].real());
        imPoDft.push_back(output[i].imag());
        magPoDft.push_back(sqrt(imPoDft[i]*imPoDft[i] + rePoDft[i] *rePoDft[i]));
    }
    vector<double> os_x = mp::linspace(0, sygnal.size(), samples);
    mp::plot(os_x, magPoDft);
    mp::show();
    return output;
}
vector<double> idft(vector<complex<double>> sygnal, int samples) {
    vector<double> Re, Im;
    
    vector<complex<double>> output;
    complex<double> suma;
    double N = sygnal.size();

    for(int f = 0; f < N; f++) {
        suma = complex<double>(0.0, 0.0);
        for(int k = 0; k < N; k++) {
            double re = cos(((2 * M_PI) / N) * k * f);
            double im = sin(((2 * M_PI) / N) * k * f);
            complex<double> zeta(re, -im);
            suma +=sygnal[k] * zeta;
        }
        output.push_back((1 / N)*suma);
    }
    vector<double> rePoDft, imPoDft, magPoDft;
    for(int i = 0; i < output.size(); i++) {
        rePoDft.push_back(output[i].real());
        imPoDft.push_back(output[i].imag());
        magPoDft.push_back(sqrt(imPoDft[i]*imPoDft[i] + rePoDft[i] *rePoDft[i]));
    }
    double koniec = rePoDft.size();
    vector<double> os_x = mp::linspace(0, koniec, samples);
    mp::plot(os_x, rePoDft);
    mp::show();
    return rePoDft;
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
        if(sin(2 * M_PI * czestotliwosc * os_x[i]) >= 0) wynik.push_back(amplitude);
        else wynik.push_back(0);
    }
//   double okres = 1/czestotliwosc;
//   double polowa_okresu = okres/2;
//   for (int i = 0; i<samples; i++){
//   double obecnyCzas = fmod(os_x[i] ,okres);
//   if(obecnyCzas<polowa_okresu) wynik.push_back(amplitude);
//   else wynik.push_back(0);
//  }
    mp::plot(os_x, wynik);
    mp::show();
    return wynik;
}

vector<double> sinnnnnnn(double amplituda, double okres, double czestotliwosc, double sampleCzestotliwosc) {
    vector<double> os_x = mp::linspace(0, okres, okres * sampleCzestotliwosc);
    vector<double> wynik;
    for(int i = 0; i < os_x.size(); ++i) {
        wynik.push_back(sin(2 * M_PI * czestotliwosc * os_x[i]));
    }
    mp::plot(os_x, wynik);
    mp::show();
    return wynik;
}
vector<double> cosssss(double amplituda, double okres, double czestotliwosc, double sampleCzestotliwosc){
    vector<double> os_x = mp::linspace(0, okres, okres * sampleCzestotliwosc);
    vector<double> wynik;
    for(int i = 0; i < os_x.size(); ++i) {
        wynik.push_back(cos(2 * M_PI * czestotliwosc * os_x[i]));
    }
    mp::plot(os_x, wynik);
    mp::show();
    return wynik;
}
int main() {
//    idft(dft(sinus(1,20,0,500,256),256),256);
//  prostokatny(5,20,0,100,250);
    vector<double> sygnal = cosssss(1, 1, 5, 1000);
    idft(dft(sygnal, 256), 256);
    return 0;


}
