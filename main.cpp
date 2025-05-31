// Projekt.cpp: definiuje punkt wejścia dla aplikacji.
//
#define _USE_MATH_DEFINES
#include <matplot/matplot.h>
#include <vector>
#include <complex>
#include <cmath>
#include <iomanip>
#include <pybind11/pybind11.h>
using namespace std;
namespace mp = matplot;
vector<complex<double>> dft(vector<double> sygnal, int samples) {
	vector<double> Re, Im;
	vector<complex<double>> complex_sygnal;
	for (double wartosc : sygnal) {
		complex_sygnal.emplace_back(wartosc, 0.0);
	}
	vector<complex<double>> output;
	complex<double> suma;
	int N = sygnal.size();

	for (int f = 0; f < N; f++) {
		suma = complex<double>(0.0, 0.0);
		for (int k = 0; k < N; k++) {
			double re = cos(((2 * M_PI) / N) * k * f);
			double im = sin(((2 * M_PI) / N) * k * f);
			complex<double> zeta(re, im);
			suma += complex_sygnal[k] * zeta;
		}
		output.push_back(suma);
	}
	vector<double> rePoDft, imPoDft, magPoDft;
	for (int i = 0; i < output.size(); i++) {
		rePoDft.push_back(output[i].real());
		imPoDft.push_back(output[i].imag());
		magPoDft.push_back(sqrt(imPoDft[i] * imPoDft[i] + rePoDft[i] * rePoDft[i]));
	}
	vector<double> os_x = mp::linspace(0, sygnal.size(), samples);
	mp::plot(os_x, magPoDft);
	return output;
}
vector<double> idft(vector<complex<double>> sygnal, int samples, double okres) {
	vector<double> Re, Im;

	vector<complex<double>> output;
	complex<double> suma;
	double N = sygnal.size();

	for (int f = 0; f < N; f++) {
		suma = complex<double>(0.0, 0.0);
		for (int k = 0; k < N; k++) {
			double re = cos(((2 * M_PI) / N) * k * f);
			double im = sin(((2 * M_PI) / N) * k * f);
			complex<double> zeta(re, -im);
			suma += sygnal[k] * zeta;
		}
		output.push_back((1 / N) * suma);
	}
	vector<double> rePoDft, imPoDft, magPoDft;
	for (int i = 0; i < output.size(); i++) {
		rePoDft.push_back(output[i].real());
		imPoDft.push_back(output[i].imag());
		magPoDft.push_back(sqrt(imPoDft[i] * imPoDft[i] + rePoDft[i] * rePoDft[i]));
	}
	double koniec = rePoDft.size();
	vector<double> os_x = mp::linspace(0, okres, samples);
	mp::plot(os_x, rePoDft);
	mp::show();
	return rePoDft;
}

vector<double> pila(double amplitude, double czestotliwosc, double t_min, double t_max, double samples, double faza) {
	vector<double> os_x = mp::linspace(t_min, t_max, samples);
	vector<double> wynik;
	for (int i = 0; i < samples; i++) {
		wynik.push_back(amplitude * (2 * (((os_x[i] + faza) * czestotliwosc) - floor(0.5 + ((os_x[i] + faza) * czestotliwosc)))));
	}
	if (faza == 0) {
		mp::xlabel("Czas [s]");
		mp::ylabel("y(t)");
		mp::plot(os_x, wynik);
		mp::show();

	}
	return wynik;
}
vector <double> prostokatny(double amplitude, double czestotliwosc, double t_min, double t_max, double samples) {
	vector<double> os_x = mp::linspace(t_min, t_max, samples);
	vector<double> wynik;
	for (int i = 0; i < samples; i++) {
		if (sin(2 * M_PI * czestotliwosc * os_x[i]) >= 0) wynik.push_back(amplitude);
		else wynik.push_back(0);
	}
	//   double okres = 1/czestotliwosc;
	//   double polowa_okresu = okres/2;
	//   for (int i = 0; i<samples; i++){
	//   double obecnyCzas = fmod(os_x[i] ,okres);
	//   if(obecnyCzas<polowa_okresu) wynik.push_back(amplitude);
	//   else wynik.push_back(0);
	//  }
	mp::xlabel("Czas [s]");
	mp::ylabel("y(t)");
	mp::plot(os_x, wynik);
	mp::show();
	return wynik;
}

vector<double> sinnnnnnn(double amplituda, double okres, double czestotliwosc, double sampleCzestotliwosc) {
	vector<double> os_x = mp::linspace(0, okres, okres * sampleCzestotliwosc);
	vector<double> wynik;
	for (int i = 0; i < os_x.size(); ++i) {
		wynik.push_back(amplituda * sin(2 * M_PI * czestotliwosc * os_x[i]));
	}
	mp::xlabel("Czas [s]");
	mp::ylabel("y(t)");
	mp::plot(os_x, wynik);
	//mp::show();
	return wynik;
}
vector<double> cosssss(double amplituda, double okres, double czestotliwosc, double sampleCzestotliwosc) {
	vector<double> os_x = mp::linspace(0, okres, okres * sampleCzestotliwosc);
	vector<double> wynik;
	for (int i = 0; i < os_x.size(); ++i) {
		wynik.push_back(amplituda * cos(2 * M_PI * czestotliwosc * os_x[i]));
	}
	mp::xlabel("Czas [s]");
	mp::ylabel("y(t)");
	mp::plot(os_x, wynik);
	//mp::show();
	return wynik;
}
vector<double> czestotliwoscUsun(vector<double> sygnal, double sample, double cutoff, double okres) {
	vector<complex<double>> wynik = dft(sygnal, sample);
	vector<double> mag;
	for (int i = 0; i < wynik.size(); i++) {
		mag.push_back(sqrt(wynik[i].real() * wynik[i].real() + wynik[i].imag() * wynik[i].imag()));

	}
	for (int j = 0; j < mag.size(); j++) {
		if (j > cutoff) {
			wynik[j] = (0.0, 0.0);
		}
	}
	vector<double> doWykresu = idft(wynik, sample, okres);
	vector<double> os_x = mp::linspace(0,doWykresu.size(), 1000);
	mp::plot(os_x, doWykresu);
	mp::show();
	return doWykresu;

}

vector<double> filtracja1D(vector<double> sygnal, vector<double> h)
{
	size_t N = sygnal.size(), M = h.size();
	vector<double> sygnalPo(N + M - 1, 0.0);

	for (size_t n = 0; n < sygnalPo.size(); n++)
	{
		double sum = 0;

		for (size_t k = 0; k < M; k++)
		{
			if (n >= k && (n - k) < N)
				sum += h[k] * sygnal[n - k];
		}
		sygnalPo[n] = sum;
	}
	return sygnalPo;
}

PYBIND11_MODULE(sygnaly, m) {
	m.doc() = "test";
	m.def("plotSin", &sinnnnnnn, "sinus",
		pybind11::arg("amplituda"), pybind11::arg("okres"), pybind11::arg("czestotliwosc"), pybind11::arg("sampleCzestotliwosc"));
	m.def("plotCos", &cosssss, "cosinus",
		pybind11::arg("amplituda"), pybind11::arg("okres"), pybind11::arg("czestotliwosc"), pybind11::arg("sampleCzestotliwosc"));
	m.def("plotSquare", &prostokatny, "prostokatny",
		pybind11::arg("amplitude"), pybind11::arg("czestotliwosc"), pybind11::arg("t_min"), pybind11::arg("t_max"), pybind11::arg("samples"));
	m.def("plotSaw", &pila, "pila",
		pybind11::arg("amplitude"), pybind11::arg("czestotliwosc"), pybind11::arg("t_min"), pybind11::arg("t_max"), pybind11::arg("samples"), pybind11::arg("faza"));
	m.def("DFT", &dft, "dyskretna transformata fouriera",
		pybind11::arg("sygnal"), pybind11::arg("samples"));
	m.def("IDFT", &idft, "odwrotna dyskretna transformata fouriera",
		pybind11::arg("sygnal"), pybind11::arg("samples"), pybind11::arg("okres"));
	m.def("filter1D", &filtracja1D, "filtracja 1D",
		pybind11::arg("sygnal"), pybind11::arg("h"));
	m.def("czestotliwoscUsun", &czestotliwoscUsun, "czestotliwoscUsun",
		pybind11::arg("sygnal"), pybind11::arg("sample"), pybind11::arg("cutoff"), pybind11::arg("okres"));
}



int main() {
	//cosssss	(5, 3, 2, 1000);
	//prostokatny(5, 3, 0, 2, 1000);
	//pila(5, 3, 0, 2, 1000, 0);
	//idft(dft(sinnnnnnn(1, 20, 0, 500, 256), 256), 256);
	//prostokatny(5, 20, 0, 100, 250);
	









	//vector<double> sygnal = cosssss(1, 1, 5, 200);
	//vector<double> sygnal2 = cosssss(2, 3, 5, 200);
	//vector<double> outSygnal;
	//auto t = mp::linspace(0, 1, 200);

	//for (int i = 0; i < sygnal2.size(); i++) {

	//	outSygnal.push_back(sygnal[i] + sygnal2[i]);
	//}



	//auto fig = mp::figure(true);
	//fig->size(800, 1000);
	//mp::grid(true);

	//mp::tiledlayout(4, 1);

	//mp::nexttile();
	//mp::title("Sygnał 1");
	//mp::xlabel("czas(t)");
	//mp::ylabel("y(t)");
	//mp::plot(t, sygnal);

	//mp::nexttile();
	//mp::title("Sygnał 2");
	//mp::xlabel("czas(t)");
	//mp::ylabel("y(t)");
	//mp::plot(t, sygnal2);

	//mp::nexttile();
	//mp::title("Suma sygnalow");
	//mp::xlabel("czas(t)");
	//mp::ylabel("y(t)");
	//mp::plot(t, outSygnal);


	//mp::nexttile();
	//mp::title("Sygnal po usunięciu czestotliwosci");
	//mp::xlabel("czas(t)");
	//mp::ylabel("y(t)");
	//czestotliwoscUsun(sygnal, 256, 200, 1);
	//mp::show();







	const int Amplituda = 5;
	const int N = 1024;
	auto t = mp::linspace(0, 1, N);
	auto y1 = sinnnnnnn	(5, 3, 6, N);

	auto y2 = sinnnnnnn(5, 3, 2, N);

	vector <double>suma_sygnalow_y(N);
	for (int i = 0; i < N; i++)
	{
		suma_sygnalow_y[i] = y1[i] + y2[i];
	}

	auto fig = mp::figure(true);

	fig->size(800, 1000);
	mp::grid(true);
	
	mp::tiledlayout(5, 1);

	mp::nexttile();
	mp::title("Sygnał 1");
	mp::xlabel("czas(t)");
	mp::ylabel("y(t)");
	mp::plot(t, y1);

	mp::nexttile();
	mp::title("Sygnał 2");
	mp::xlabel("czas(s)");
	mp::ylabel("y(t)");
	mp::plot(t, y2);


	mp::nexttile();
	mp::title("Sygnał");
	mp::xlabel("czas(s)");
	mp::ylabel("y(t)");
	mp::plot(t, suma_sygnalow_y);

	mp::nexttile();
	mp::title("Sygnał po usunieciu");
	mp::xlabel("czas(s)");
	mp::ylabel("y(t)");
	czestotliwoscUsun(suma_sygnalow_y, N, 3, 1);

	
	//mp::nexttile();
	//mp::title("DFT");
	//mp::xlabel("czestotliwosc(Hz)");
	//dft(suma_sygnalow_y, N);

	//mp::nexttile();
	//mp::title("IDFT");
	//mp::xlabel("czas(s)");
	//mp::ylabel("y(t)");
	//idft(dft(suma_sygnalow_y, N), N, 1);
	










	//const size_t N = 500;
	//vector<double> os_x(N), sygnal(N);
	//for (size_t i = 0; i < N; i++)
	//{
	//	os_x[i] = i * 0.01;
	//	sygnal[i] = sin(2 * M_PI * 2 * os_x[i]) + 0.5 * ((rand() / (double)RAND_MAX) - 0.5);
	//}

	//const size_t M = 21;
	//vector<double> h(M, 1.0 / M);

	//vector<double> SygnalPO = filtracja1D(sygnal, h);

	//auto fig = mp::figure(true);
	//fig->size(800, 400);
	//mp::tiledlayout(2, 1);


	//mp::nexttile();
	//mp::plot(os_x, sygnal);
	//mp::title("Sygnał oryginalny");
	//mp::xlabel("Czas [s]");
	//mp::ylabel("Amplituda");


	////vector<double> y_cut(N), t_cut(N);
	////size_t offset = M / 2;
	////for (size_t i = 0; i < N; ++i) {
	////	y_cut[i] = SygnalPO[i + offset];
	////	t_cut[i] = os_x[i];
	////}
	//mp::nexttile();
	//mp::plot(os_x, SygnalPO);
	////mp::plot(os_x, SygnalPO);
	//mp::title("Sygnał przefiltrowany");
	//mp::xlabel("Czas [s]");
	//mp::ylabel("Amplituda");

	//mp::show();


	return 0;
}
























