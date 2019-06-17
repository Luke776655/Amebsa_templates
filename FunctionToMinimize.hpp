// GNU GPL license v3, Janek Kozicki
#ifndef FUNCTION_TO_MINIMIZE_HEAD
#define FUNCTION_TO_MINIMIZE_HEAD

#include <cmath>
#include <stdexcept>

using namespace std;

// w tym pliku są wzory funkcji które chcemy zminimalizować

// to jest ogólna postać tej klasy, tylko określa interfejs, tzn, że wszystkie one muszą mieć operator()(), tzn. żeby można było ją wywołać pisząc po nazwie zmiennej tylko nawiasy.
// oprócz tego jest nieprzydatna.
template<typename ReturnType, typename ArgumentType>
struct FunctionToMinimize {
	typedef ReturnType	ret_type;
	typedef ArgumentType	arg_type;
	ReturnType operator()(const ArgumentType& arg) const;
};



// parabola jednowymiarowa
template<typename ReturnType, typename ArgumentType>
struct Parabola1D : public FunctionToMinimize<ReturnType,ArgumentType> {
	ReturnType operator()(const ArgumentType& x) const {
		return x*x;
	};
};

// parabola N-wymiarowa
template<typename ReturnType, typename ArgumentType>
struct ParabolaND : public FunctionToMinimize<ReturnType,vector<ArgumentType>> {
	ReturnType operator()(const vector<ArgumentType>& xx) const {
		ReturnType ret=0;
		for(const auto& x : xx) ret+=x*x;
		return ret;
	};
};

// parabola N-wymiarowa zamieszana z cosinuesem.
template<typename ReturnType, typename ArgumentType>
struct ParabolaCosND : public FunctionToMinimize<ReturnType,vector<ArgumentType>> {
	ReturnType operator()(const vector<ArgumentType>& xx) const {
		ReturnType ret=0;
		for(const auto& x : xx) ret+=x*x - 10*cos(x);
		return ret;
	};
};

// ferromagnetyk N-wymiarowy
template<typename ReturnType, typename ArgumentType>
struct FerromagnetND : public FunctionToMinimize<ReturnType,vector<ArgumentType>> {
	ReturnType operator()(const vector<ArgumentType>& A) const {
	ReturnType s = 0;
	int N = A.size();
    int k = sqrt(N);
	for (int i = 0; i<N; i++)
	{
		if(i>=k) //top boundary
			s = s-cos(abs(A[i]-A[i-k]))+1;
		if(i<N-k) //bottom boundary
			s = s-cos(abs(A[i]-A[i+k]))+1;
		if(i%k!=0) //left boundary
			s = s-cos(abs(A[i]-A[i-1]))+1;
		if(i%k!=(k-1)) //right boundary
			s = s-cos(abs(A[i]-A[i+1]))+1;
	}
	return s;
    };
};

#endif

