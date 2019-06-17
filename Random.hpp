// GNU GPL license v3, Janek Kozicki
#ifndef RANDOM_HEAD
#define RANDOM_HEAD

// w tym pliku są dwie przydatne klasy do generowania liczb losowych, w zakresie [0,1] i [-1,1]
// przykład: deklarujemy zmienną:
//	Random_1_1<double> r1;
// losujemy liczbę:
//      double losowa_liczba = r1();

#include <boost/random/uniform_real.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/random_device.hpp>

template<typename RealType>
class RandomMinMax {
	protected:
		boost::random_device		m_rd;
		boost::uniform_real<RealType>	m_uni_dist;
		boost::variate_generator<boost::random_device&, boost::uniform_real<RealType> > m_uni_min_max;
	public:
		RandomMinMax(RealType min, RealType max) : m_uni_dist(min,max) , m_uni_min_max(m_rd,m_uni_dist) {};
		RealType operator()() { return m_uni_min_max(); };
};

// ta klasa zawsze da liczbę losową w przedziale [0,1]
template<typename RealType>
struct Random_0_1 : public RandomMinMax<RealType> {
	Random_0_1() : RandomMinMax<RealType>(0,1) {};
};

// ta klasa zawsze da liczbę losową w przedziale [-1,1]
// przykład: deklarujemy zmienną:
//	Random_1_1<double> r1;
// losujemy liczbę:
//      double losowa_liczba = r1();
template<typename RealType>
struct Random_1_1 : public RandomMinMax<RealType> {
	Random_1_1() : RandomMinMax<RealType>(-1,1) {};
};

// klasa do losowania liczb całkowitych z zakresu [min,max], taki rzut kostką, np. od 1 do 6:
// przykład:
//	RandomMinMax<int> r2(1,6);
// losujemy liczbę:
//      int losowa_liczba = r2();
template<>
class RandomMinMax<int> {
	protected:
		boost::random_device		m_rd;
		boost::uniform_int<int>		m_uni_dist;
		boost::variate_generator<boost::random_device&, boost::uniform_int<int> > m_uni_min_max;
	public:
		RandomMinMax(int min, int max) : m_uni_dist(min,max) , m_uni_min_max(m_rd,m_uni_dist) {};
		int operator()() { return m_uni_min_max(); };
};

#endif

