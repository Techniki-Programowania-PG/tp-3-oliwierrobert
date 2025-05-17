#include <iostream>
#include <matplot/matplot.h>
#include <cmath>
#include <iomanip>
using namespace std;
namespace mp = matplot;
void cosinus(double amplituda, double czestotliwosc, double t_min, double t_max, double samples) {
    vector<double> wynik;
    // wynik.resize(samples);
    vector<double> os_x = mp::linspace(t_min, t_max, samples);
    for(int i = 0; i < samples; i++) {
        cout << setw(5) << os_x[i] << endl;
    }
    for(int i = 0; i < samples; i++) {
        wynik.push_back(amplituda * cos(2 * M_PI * os_x[i] * czestotliwosc));
    }
    mp::plot(os_x, wynik);
    mp::show();
}
void pila(double amplitude, double czestotliwosc, double t_min, double t_max, double samples){
    vector<double> os_x = mp::linspace(t_min, t_max, samples);
    vector<double> wynik;
    for( int i=0 ; i <samples; i++){
        wynik.push_back(2*((os_x[i]*czestotliwosc )- floor(0.5 + (os_x[i]*czestotliwosc))));
    }
    mp::plot(os_x,wynik);
    mp::show();
}
void sinus(double amplituda, double czestotliwosc, double t_min, double t_max, double samples) {
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
    mp::plot(os_x, wynik);
    mp::show();
}
int main() {
    pila(5,10,0,5000,5000);
    return 0;
}
