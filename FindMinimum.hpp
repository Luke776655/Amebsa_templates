// GNU GPL license v3, Janek Kozicki
#ifndef FIND_MINIMUM_HEAD
#define FIND_MINIMUM_HEAD

#include <cstddef>
#include <functional>
#include <cstdio>
#include <iostream>
#include <cstdlib>
#include <vector>
#include <cmath>
#include <limits>

// w tym pliku jest rodzina klas wyszukujących minimum
// to jest ogóla postać tej klasy, tylko określa interfejs, tzn, że wszystkie one muszą mieć funkcję findMinimum()
// oprócz tego jest nieprzydatna.
template<class FunctionToMinimizeType>
class FindMinimum {
	protected:
		const FunctionToMinimizeType& m_func;
	public:
		FindMinimum(const FunctionToMinimizeType& func) // konstruktorowi trzeba podać JAKĄ funkcję chcemy zminimalizować, i to wszystko.
			: m_func(func)
		{};
		typename FunctionToMinimizeType::arg_type findMinimum();
};



// zaimplementowany algorytm metropolis do wyszukiwania minimum funkcji
template<class FunctionToMinimizeType, class RandomizeArgumentsType, typename TypeReal = typename FunctionToMinimizeType::ret_type>
class Metropolis : public FindMinimum<FunctionToMinimizeType> {
	private:
		typedef					typename FunctionToMinimizeType::ret_type	ret_type;
		typedef					typename FunctionToMinimizeType::arg_type	arg_type;
		typedef					TypeReal					real_type;
		RandomizeArgumentsType&			m_randomize_arguments;
		real_type				m_temperature;
		Random_0_1<real_type>			m_uniform_0_1;
		std::function<void(arg_type,ret_type)>	m_printer;
	public:
		Metropolis(const FunctionToMinimizeType& function_to_minimize, RandomizeArgumentsType& r, TypeReal t, std::function<void(arg_type,ret_type)> printer=nullptr)
		// konstruktorowi trzeba podać
			: FindMinimum<FunctionToMinimizeType>(function_to_minimize) // JAKĄ funkcję chcemy zminimalizować, jest ona przekazana do klasy macierzystej FindMinimum
			, m_randomize_arguments(r) // trzeba podać specjalny parametr - w jaki sposób losować nowe parametry, bo algorytm metropolis opiera się na losowaniu
			, m_temperature(t)         // trzeba podać temperaturę wyżarzania.
			, m_printer(printer)       // opcjonalnie można podać funkcję do drukowania na ekranie.
		{};
		void setTemperature(TypeReal t) { m_temperature=t; };
		// findMinimum() według alg. metropolis działa tak, że losuje nowe parametry i sprawdza czy są lepsze.
		// jeśli są lepsze (lub minimalnie gorsze - zależy od m_temperature) to się na nie przełącza
		typename FunctionToMinimizeType::arg_type findMinimum(arg_type the_arg, int steps) {
			ret_type the_energy = FindMinimum<FunctionToMinimizeType>::m_func(the_arg);
			while(--steps>0) {
				arg_type new_arg	= m_randomize_arguments.randomize(the_arg);
				ret_type new_energy	= FindMinimum<FunctionToMinimizeType>::m_func(new_arg);
				if(std::log(m_uniform_0_1()) < (-1.0*(new_energy - the_energy)/m_temperature) ) {
					the_energy	= std::move(new_energy);
					the_arg		= std::move(new_arg   );
				}
				if(m_printer) { m_printer(the_arg, the_energy); };
			}
			return the_arg;
		}
};

using namespace std;

template<class FunctionToMinimizeType, class arg_type, class ret_type>
class Amebsa : public FindMinimum<FunctionToMinimizeType>
{
	public:
		Amebsa(const FunctionToMinimizeType& function_to_minimize)
		:FindMinimum<FunctionToMinimizeType>(function_to_minimize)
		{	};
		ret_type findMinimum(vector <arg_type> &pt);
		vector <double> temperature = {1, 1e-10, 0};
		ret_type fmin = numeric_limits<ret_type>::max();
		vector <arg_type> pmin;
		int NMAX = 10000000; //maximum allowed number of function evaluations
		int delta = 1; //displacement
		double ftol = 1e-16; //fractional convergence tolerance
		double TINY = 1e-7; //tiny value preventing from dividing by 0
		bool arrows_table_printing = false;
		bool show_iter_output = true;
		int iter_period = 1000;
	private:
		int N; //dimention of side of pseudo-square table
		void print_array(vector <arg_type> tab);
		void print_table(vector <arg_type> tab);
		void print_arrow(arg_type s);
		void print_arrow_table(vector <arg_type> tab);
		void print_result();
		vector <arg_type> get_psum();
		ret_type amotsa(ret_type fac);
		int amebsa_alg();
		vector<arg_type> point;
		int ndim;
		vector <vector <arg_type> > p;
		vector <ret_type> y;
		int ilo;
		int ihi;
		int inhi;
		ret_type ynhi;
		ret_type yhi;
		ret_type ylo;
		ret_type ytry;
		ret_type rtol;
		int mpts;
		int iter;
		int nfunc;
		ret_type tt;
		vector <arg_type> psum;
};

template<class FunctionToMinimizeType, class arg_type, class ret_type> void Amebsa<FunctionToMinimizeType, arg_type, ret_type>::print_array(vector <arg_type> tab)
{
	/*
	printing values of the array
	*/
	for(int i = 0; i<N; i++)
	{
		cout<<tab[i]<<"  ";	
	}
	printf("\n");
}

template<class FunctionToMinimizeType, class arg_type, class ret_type> void Amebsa<FunctionToMinimizeType, arg_type, ret_type>::print_table(vector <arg_type> tab)
{
	/*
	printing values of the table
	*/
	for(int i = 0; i<N; i++)
	{
		if(i%(int)sqrt(N)==0)
		{
			printf("\n");
		}
		cout<<tab[i]<<"  ";	
	}
	printf("\n");
}

template<class FunctionToMinimizeType, class arg_type, class ret_type> void Amebsa<FunctionToMinimizeType, arg_type, ret_type>::print_arrow(arg_type s)
{
	/*
	printing value as an arrow with proper slope
	*/
	if(s<(M_PI/8+0*M_PI/4) || s>=(M_PI/8+7*M_PI/4))
		printf("↑ ");
	else if(s>=(M_PI/8+0*M_PI/4) && s<(M_PI/8+1*M_PI/4))
		printf("↗ ");
	else if(s>=(M_PI/8+1*M_PI/4) && s<(M_PI/8+2*M_PI/4))
		printf("→ ");
	else if(s>=(M_PI/8+2*M_PI/4) && s<(M_PI/8+3*M_PI/4))
		printf("↘ ");
	else if(s>=(M_PI/8+3*M_PI/4) && s<(M_PI/8+4*M_PI/4))
		printf("↓ ");
	else if(s>=(M_PI/8+4*M_PI/4) && s<(M_PI/8+5*M_PI/4))
		printf("↙ ");
	else if(s>=(M_PI/8+5*M_PI/4) && s<(M_PI/8+6*M_PI/4))
		printf("← ");
	else if(s>=(M_PI/8+6*M_PI/4) && s<(M_PI/8+7*M_PI/4))
		printf("↖ ");
	else
		printf("x ");
}

template<class FunctionToMinimizeType, class arg_type, class ret_type> void Amebsa<FunctionToMinimizeType, arg_type, ret_type>::print_arrow_table(vector <arg_type> tab)
{
	/*
	printing values of the table as arrows with proper slope
	*/
	for(int i=0; i<N; i++)
	{
		if(i%(int)sqrt(N)==0)
			printf("\n");
		print_arrow(fmod(fmod(tab[i], 2*M_PI)+2*M_PI, 2*M_PI));
	}
	printf("\n\n");
}

template<class FunctionToMinimizeType, class arg_type, class ret_type> void Amebsa<FunctionToMinimizeType, arg_type, ret_type>::print_result()
{
	/*
	printing result of optimisation
	*/
	arg_type z = 0;
	vector <ret_type> k = y;
	vector <vector <arg_type> > q = p;

	z = k[0];
	k[0] = k[ilo];
	k[ilo] = z;

	ret_type g;
	for(int i = 0; i<ndim; i++)
	{
		g = q[0][i];
		q[0][i] = q[ilo][i];
		q[ilo][i] = g;
		pmin[i] = q[0][i];
	}
	fmin=k[0];
	printf("\nIteration %d\tTemperature %g\n", iter, -tt);
	if(arrows_table_printing == true)
	{
		print_table(pmin);
		print_arrow_table(pmin);
	}
	else
		print_array(pmin);
	printf("Energy of the system %g\n\n", fmin);
}

template<class FunctionToMinimizeType, class arg_type, class ret_type> vector <arg_type> Amebsa<FunctionToMinimizeType, arg_type, ret_type>::get_psum()
{
	/*
	counting partial sum
	*/
	vector <arg_type> ps(ndim);
	for (int j=0; j<ndim; j++)
	{
		arg_type sum=0;
		for(int i=0; i<mpts; i++)
		{
			sum += p[i][j];
		}
		ps[j]=sum;
	}
	return ps;
}

template<class FunctionToMinimizeType, class arg_type, class ret_type> ret_type Amebsa<FunctionToMinimizeType, arg_type, ret_type>::amotsa(ret_type fac)
{
	/*
	extrapolation by a factor fac through the face of the simplex across from the high point
	replacing the high point if the new point is better
	*/
	vector <arg_type> ptry(ndim);
	ret_type fac1=(1.0-fac)/ndim;
	ret_type fac2=fac1-fac;
	for(int j=0; j<ndim; j++)
	{
		ptry[j]=psum[j]*fac1-p[ihi][j]*fac2;
	}
	ytry=FindMinimum<FunctionToMinimizeType>::m_func(ptry);
	ret_type yflu=ytry-tt*log(0.5*(ret_type)rand()/RAND_MAX);
	if(yflu < yhi)
	{
		y[ihi]=ytry;
		yhi=yflu;
		for(int j=0; j<ndim; j++)
		{
			psum[j] += ptry[j]-p[ihi][j];
			p[ihi][j]=ptry[j];
		}
	}
	return yflu;
}

template<class FunctionToMinimizeType, class arg_type, class ret_type> ret_type Amebsa<FunctionToMinimizeType, arg_type, ret_type>::findMinimum(vector <arg_type> &pt)
{
	srand (time(NULL));
	point = pt;
	N = pt.size();
	ndim = N;
	pmin.resize(ndim);

	//creating delta values table
	vector <int> delta_tab;
	for(int i=0; i<ndim; i++)
	{
		delta_tab.push_back(delta);
	}

	//adding delta values to the point table with extended dimention as p simplex

	for(int i=0; i<ndim+1; i++)
	{
		vector<arg_type> k;
		for(int j=0; j<ndim; j++)
		{
			k.push_back(point[j]);
		}
		p.push_back(k);
		if(i!=0)
		{
			p[i][i-1]+=delta_tab[i-1];
		}
	}
	//getting y table of solutions
	y.resize(ndim+1);
	yhi = 0;
	mpts = p.size(); //number of rows
	for(int i = 0; i<mpts; i++)
	{
		vector <arg_type> x;
		for(int j=0; j<ndim; j++)
		{
			x.push_back(p[i][j]);
		}
		y[i] = FindMinimum<FunctionToMinimizeType>::m_func(x);
	}
	//parameters for iterating
	nfunc = 0;
	psum = get_psum();
	for(int unsigned i=0; i<temperature.size(); i++)
	{
		int w = -1;
		iter = 0;
		nfunc = 0;
		tt = -temperature[i];
		while(w<0)
		{
			w = amebsa_alg();
		}
	}
	return fmin;
}

template<class FunctionToMinimizeType, class arg_type, class ret_type> int Amebsa<FunctionToMinimizeType, arg_type, ret_type>::amebsa_alg()
{
	ilo=0;
	ylo=y[ilo]+tt*log(0.5*(ret_type)rand()/RAND_MAX);
	ihi=1;
	yhi=y[ihi]+tt*log(0.5*(ret_type)rand()/RAND_MAX);

	if(ylo>yhi)
	{
		inhi = 1;
		ynhi = yhi;
		ihi = 0;
		yhi = ylo;
	}
	else
	{
		inhi = 0;
		ynhi = ylo;
	}

	for (int i=0; i<mpts; i++)
	{
		ret_type yi = y[i]+tt*log(0.5*(ret_type)rand()/RAND_MAX);
		if(yi<=ylo)
		{
			ilo = i;
			ylo = yi;
		}
		if(yi>yhi)
		{
			inhi = ihi;
			ynhi = yhi;
			ihi = i;
			yhi = yi;
		}
		else if(yi > ynhi && i != ihi)
		{
			inhi = i;
			ynhi = yi;
		}
	}

	rtol=2.0*abs(y[ihi]-y[ilo])/(abs(y[ihi])+abs(y[ilo])+TINY);
	if (rtol < ftol)
	{
		print_result();
		printf("Optimisation succeeded\n");
		return 0;
	}

	if (nfunc >= NMAX)
	{
		print_result();
		printf("Maximal number of iterations exceeded\n");
		return 1;
	}
	nfunc += 2;

	ytry = amotsa(-1.0); //simplex reflection
	if(ytry <= ylo)
	{
		ytry = amotsa(2.0); //simplex extrapolation
	}
	else if (ytry >= ynhi)
	{
		ret_type ysave = yhi;
		ytry = amotsa(0.5); //simplex contraction
		if(ytry >= ysave)
		{
			for (int i=0;i<mpts;i++) {
				if (i != ilo) {
					for (int j=0;j<ndim;j++)
					{
						psum[j]=0.5*(p[i][j]+p[ilo][j]);
						p[i][j]=psum[j];
					}
					y[i]=FindMinimum<FunctionToMinimizeType>::m_func(psum);
				}
			}
			nfunc += ndim;
			psum = get_psum();
		}
	} else nfunc = nfunc-1;
	if(show_iter_output == true)
	{
		if(iter%iter_period == 0)
		{
			print_result();
		}
	}
	iter++;
	return -1;
}



#endif
