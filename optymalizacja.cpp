// GNU GPL license v3, Janek Kozicki
#include "RandomizeArguments.hpp"
#include "FunctionToMinimize.hpp"
#include "FindMinimum.hpp"

#include <cmath>
#include <vector>
#include <iostream>
#include <iomanip>

// gdybym chciał wszystkie `double` zamienić poniżej na `float` albo `long double`, to wystarczy zmienić tylko tego typedefa
typedef long double Real;

// ta funkcja jest tylko po to, żeby działało wypisywanie na ekran std::vector<double>, dla N-wymiarowych funkcji
template<typename T>
std::ostream& operator<<(std::ostream& return_stream,const std::vector<T>& s) {
	for(const auto& a : s) return_stream << a << " ";
	return return_stream;
}

void printSpins(const std::vector<Real>& spins, const Real& E, int printEveryNth) {
	static const std::vector<std::string> arrows{"↑","↗","→","↘","↓","↙","←","↖"};
	static int callNum{0};
	if(++callNum % printEveryNth != 0) { return; }

	int N=std::sqrt(spins.size());
	if(N*N != int(spins.size())) throw std::length_error("const std::vector<Real>& spins should be a pseudo-square table.");

	std::ios_base::fmtflags previousFlags( std::cout.flags() );
	std::cout << std::setw(4) << std::setfill('0') << std::left << std::setprecision(2) << std::fixed << std::showpos;
	for(bool what : {true,false}) {
		int count=0;
		for(const Real& spin : spins) {
			int pos = (unsigned int)(std::floor(Real(8.0*fmod(spin,2*M_PI)/(2*M_PI))))%8;
			if(pos<0 or pos>=8) {
				std::cerr << "Error: " << spin << " " << pos << "\n";
				throw std::logic_error("Arrows calculated wrong.");
			}
			//std::cout << arrows[pos] << ":" << spin << ":" << pos << " ";
			if(what) {
				std::cout << spin << " ";
			} else {
				std::cout << arrows[pos] << " ";
			}
			if( (++count % N) ==0) { std::cout << "\n"; }
		}
		std::cout << "\n";
	}
	std::cout.flags( previousFlags );
	std::cout << std::setprecision(std::numeric_limits<Real>::digits10+1) << "Energy of the system " << E << "\n\n";
}

void metropolis_parabola1D(int steps, Real temperature) {
	RealWalk1D<Real>		real_walk(0.5); // będzie zmieniał losowo o -0.5 do +0.5
	Parabola1D<Real, Real>		parabola;
	//          ↓↓↓ będzie liczył Parabola1D z pliku FunctionToMinimize.hpp
	//                                   ↓↓↓ oraz losował według RealWalk1D z pliku Random.hpp
	Metropolis< Parabola1D<Real, Real> , RealWalk1D  <Real> > metropolis1D( parabola , real_walk , temperature /* tutaj pomija drukowanie (funkcja Metropolis::m_printer) */ );

	// szukamy minimum                      ↓↓ zacznij szukanie od wartości -10
	Real minimum1D = metropolis1D.findMinimum(-10, steps);
	// drukujemy wynik
	std::cout << "metropolis 1D minimum: " << minimum1D << " has energy = " << parabola(minimum1D) << "\n\n\n";
}

void metropolis_parabolaND(int steps, Real temperature) {
	RealWalkND<Real>		real_walk_ND(2.5);    // będzie zmieniał losowo o 2.5 w N-wymiarach
	ParabolaCosND<Real, Real>	parabola_ND;          // poszukajmy minimum takiej funkcji parabola wymieszana z cosinusem
	Metropolis<ParabolaCosND<Real, Real>,RealWalkND<Real>> metropolisND( parabola_ND , real_walk_ND , temperature ,
		// tutaj drukuje przy pomocy funkcji lambda
		[](const std::vector<Real>& arg, const Real& E){
			std::cout << arg << " has energy=" << E << "\n";
		}
	);

	// poszukaj minimum w paraboli 4D         ↓↓↓↓↓↓↓↓↓↓↓↓ zacznij w tym punkcie
	auto minimum4D = metropolisND.findMinimum({10,10,20,-50}, steps);
	// drukujemy wynik
	std::cout << "metropolis 4D minimum: " << minimum4D << " has energy = " << parabola_ND(minimum4D) << "\n\n\n";
}

void metropolis_ferromagnetyk(int steps, std::vector<Real> temps, int printEveryNth) {
	RealWalkND<Real>		real_walk_ND(0.2);    // będzie zmieniał losowo o -0.2 N-wymiarach
	FerromagnetND<Real, Real>	ferromagnet_ND;       // poszukujemy minimum takiej funkcji
	using namespace std::placeholders;                    // żeby móc pisać _1 zamiast std::placeholders::_1 jako pierwszy argument dla funkcji
	Metropolis<FerromagnetND<Real, Real>,RealWalkND<Real>> metropolisND( ferromagnet_ND , real_walk_ND , temps[0] , std::bind(printSpins, _1 , _2 , printEveryNth)); // używamy takiego algorytmu szukania minimum
	// trzeba zainicjalizować tablicę ze spinami jakoś losowo
	std::vector<Real> ferroSpins(100,0);
	Random_0_1<Real>		uniform_0_1; // z pliku Random.hpp
	for(Real& s : ferroSpins) { s = uniform_0_1()*2.0*M_PI; };
	// szukamy minimum
	for(Real T : temps ) {
		std::cout << "Using temperature: " << T << "\n";
		metropolisND.setTemperature(T);
		ferroSpins = metropolisND.findMinimum(ferroSpins, steps);
	}
	// drukujemy wynik
	std::cout << "Finished: metropolis ferromagnet minimum:\n";
	printSpins(ferroSpins , ferromagnet_ND(ferroSpins) , 1);
}

using namespace std;

template <typename T> vector <T> create_random_2PI_table(int n)
{
	/*
	creating a starting point of the system:
	n-dimentional table with random values
	values range: [0, 2*pi)
	*/
	vector <T> tab(n);
	for(int i = 0; i<n; i++)
	{
		tab[i] = fmod(rand(), 2*M_PI);
	}
	return tab;
}

void amebsa_parabola_ND(int n, int steps, int printEveryNth)
{
    ParabolaND<double, double> parabola_ND;
    Amebsa<ParabolaND<double, double>, double, double> a ( parabola_ND );
    vector<double> v = create_random_2PI_table<double>(n);
    a.NMAX = steps;
    a.temperature = {0};
    a.iter_period = printEveryNth;
    a.findMinimum(v);
}

void amebsa_parabola_cos_ND(int n, int steps, int printEveryNth)
{
    ParabolaCosND<double, double> parabola_cos_ND;
    Amebsa<ParabolaCosND<double, double>, double, double> a ( parabola_cos_ND );
    vector<double> v = create_random_2PI_table<double>(n);
    a.NMAX = steps;
    a.temperature = {0};
    a.iter_period = printEveryNth;
    a.findMinimum(v);
}

void amebsa_ferromagnetyk(int n, int steps, int printEveryNth, vector<Real> temps)
{
    FerromagnetND<Real, Real> ferromagnet_ND;
    Amebsa<FerromagnetND<Real, Real>, Real, Real> a ( ferromagnet_ND );
    vector<Real> v = create_random_2PI_table<Real>(n);
    a.NMAX = steps;
    a.temperature = temps;
    a.iter_period = printEveryNth;
    a.arrows_table_printing = true;
    a.findMinimum(v);
}

int main()
{
    int n = 400;
	int steps=10000000;
	int printEveryNth=1000;
    //amebsa_parabola_ND(n, steps, printEveryNth);
    //amebsa_parabola_cos_ND(n, steps, printEveryNth);
	amebsa_ferromagnetyk(n, steps, printEveryNth, {1, 1e-10, 0});
    return 0;
}

