// GNU GPL license v3, Janek Kozicki
#ifndef RANDOMIZE_ARGUMENTS_HEAD
#define RANDOMIZE_ARGUMENTS_HEAD

// W tym pliku wszystkie klasy używane przez metropolis, każda z nich daje funkcję randomize(…) która w jakiś sposób losowo zmienia argument.
#include "Random.hpp"

#include <vector>

// to jest ogólna postać tej klasy, tylko określa interfejs, tzn, że wszystkie one muszą mieć funkcję randomize(…)
// oprócz tego jest nieprzydatna.
template<typename ArgumentType>
class RandomizeArguments {
	public:
		ArgumentType randomize(const ArgumentType& prev);
};


template<typename ArgumentType>
//            ↓↓ jest 1D, tzn do problemów jednowymiarowych.
//        Walk - tzn. że losowo chodzi, modyfikuje poprzednią wartość.    
class RealWalk1D : public RandomizeArguments<ArgumentType> {
	private:
		mutable Random_1_1<ArgumentType>	m_uni_real;
		ArgumentType				m_scale;
	public:
		RealWalk1D(ArgumentType scale) : m_scale(scale) {};
		// losuje liczbę w przedziale [-1,1]*m_scale i dodaje do poprzednie wartości
		ArgumentType randomize(const ArgumentType& prev) {
			return prev + m_uni_real()*m_scale;
		}
};

template<typename ArgumentType>
//            ↓↓ jest ND, tzn do problemów N-wymiarowych.
//        Walk - tzn. że losowo chodzi, modyfikuje poprzednią wartość.    
class RealWalkND : public RandomizeArguments<ArgumentType> {
	private:
		mutable Random_1_1<ArgumentType>	m_uni_real;
		ArgumentType				m_scale;
	public:
		RealWalkND(ArgumentType scale) : m_scale(scale) {};
		// losuje N liczb w przedziale [-1,1]*m_scale i dodaje do poprzednich wartości
		std::vector<ArgumentType> randomize(std::vector<ArgumentType> ret) {
			for(auto& a : ret) a += m_uni_real()*m_scale;
			return ret;
		}
};


#endif

