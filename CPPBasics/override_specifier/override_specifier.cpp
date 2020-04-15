#include <iostream>
#include <iostream>
#include <string>
#include <iterator>
#include <functional>
#include <algorithm>
#include <iostream>
#include <vector>
#include <numeric>
#include <string>
#include <functional>

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

#pragma pack(1)

struct Empty {};
struct EmptyVirt { virtual ~EmptyVirt() {} };
struct NotEmpty { int m_i; };
struct NotEmptyVirt
{
    virtual ~NotEmptyVirt() {}
    int m_i;
};
struct NotEmptyNonVirt
{
    void foo() const {}
    int m_i;
};


class Parent
{
public:
	// Этот метод getThis() возвращает указатель на класс Parent
	virtual Parent* getThis() { std::cout << "called Parent::getThis()\n"; return this; }
	virtual void printType() { std::cout << "returned a Parent\n"; }
};

class Child : public Parent
{
public:
	// Обычно, типы возврата переопределений и виртуальных функций родительского класса должны совпадать
	// Однако, поскольку Child наследует класс Parent, то следующий метод может возвращать Child* вместо Parent*
	virtual Child* getThis() { std::cout << "called Child::getThis()\n";  return this; }
	void printType() { std::cout << "returned a Child\n"; }
};

class A {
public:
    A() { std::cout << "A()" << std::endl; }
    ~A() { std::cout << "~A()" << std::endl; }
};

class B : public A {
public:
    B() { std::cout << "B()" << std::endl; }
    ~B() { std::cout << "~B()" << std::endl; }
};


struct Size {
    int width, height;
};


class Spline;
class SplineImpl;


/***** Spline factory functions *****/

/* Linear interpolation */
Spline LinearSpline(const std::vector<double>& X, const std::vector<double>& Y);
Spline LinearSpline(int N, const double* X, const double* Y);

/* Shifted linear interpolation (see Blu, Thevenaz, and Unser, 2004) */
Spline ShiftedLinearSpline(const std::vector<double>& X, const std::vector<double>& Y, double tau = 0.2);
Spline ShiftedLinearSpline(int N, const double* X, const double* Y, double tau = 0.2);

/* Natural cubic spline */
Spline CubicSpline(const std::vector<double>& X, const std::vector<double>& Y);
Spline CubicSpline(int N, const double* X, const double* Y);


/***************************************************************
 * Spline
 *
 * Generic spline wrapper class.
 ***************************************************************/
class Spline {
public:
    /* Default to cubic spline */
    Spline(const std::vector<double>& X, const std::vector<double>& Y);
    Spline(int N, const double* X, const double* Y);

    Spline();
    Spline(SplineImpl* impl, bool clone = false);
    Spline(const Spline& F);
    ~Spline();

    Spline& operator=(const Spline& F);

    double Evaluate(double x) const;
    double EvaluateDerivative(double x) const;

    double operator()(double x) const { return Evaluate(x); }

    /* Find a local maximum or minimum of the interpolated function */
//    double FindMaximum(double xguess, double* ymax = 0);
//    double FindMinimum(double xguess, double& ymin = 0);

protected:
    SplineImpl* impl;   // internal spline implementation
};


/**************************************************
 * SplineImpl
 *
 * Base class for internal spline implementations.
 **************************************************/
struct SplineImpl {
    SplineImpl() {}
    virtual ~SplineImpl() {}

    virtual double y(double x) const = 0;
    virtual double dydx(double x) const = 0;
    virtual SplineImpl* clone() const = 0;

    double xmin, xmax;    // domain
};


int intrpolation(void)
{
    int i;
    double xi, yi, x[10], y[10];

    printf("#m=0,S=17\n");

    for (i = 0; i < 10; i++)
    {
        x[i] = i + 0.5 * sin(i);
        y[i] = i + cos(i * i);
        printf("%f %f\n", x[i], y[i]);
    }

    printf("#m=1,S=0\n");

    {
        gsl_interp_accel* acc
            = gsl_interp_accel_alloc();
        gsl_spline* spline
            = gsl_spline_alloc(gsl_interp_linear, 10);

        gsl_spline_init(spline, x, y, 10);

        for (xi = x[0]; xi < x[9]; xi += 0.01)
        {
            yi = gsl_spline_eval(spline, xi, acc);
            printf("%f %f\n", xi, yi);
        }
        gsl_spline_free(spline);
        gsl_interp_accel_free(acc);
    }
    return 0;
}

int main()
{
    B b;
	Child ch;
	Parent* p = &ch;
	ch.getThis()->printType(); // вызывается Child::getThis(), возвращается Child*, вызывается Child::printType
	p->getThis()->printType(); // вызывается Child::getThis(), возвращается Parent*, вызывается Parent::printType
	auto apthis = p->getThis();
	apthis->printType();
	((Child *)apthis)->printType();
	p->printType();

    std::cout << sizeof(Empty) << std::endl;
    std::cout << sizeof(EmptyVirt) << std::endl;
    std::cout << sizeof(NotEmpty) << std::endl;
    std::cout << sizeof(NotEmptyVirt) << std::endl;
    std::cout << sizeof(NotEmptyNonVirt) << std::endl;
    A* pA = new B;
    delete pA;

    std::string A = "Privet! Kak, dela?.,,A?", B;
    std::copy_if(A.begin(), A.end(), std::inserter(B, B.begin()), std::not_fn(std::ref(ispunct)));
    std::cout << "Result: " << B << std::endl;

    std::string s1 = "Privet! Kak, dela?.,,A?", s2;
    std::remove_copy_if(s1.begin(), s1.end(), std::back_inserter(s2), ispunct);
    std::cout << "Result: " << s2 << std::endl;

    copy_if(A.begin(), A.end(), back_inserter(B), [](char a) {return (!ispunct(a)); });

    std::vector<int> v{ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 };

    int sum = std::accumulate(v.begin(), v.end(), 0);

    int product = std::accumulate(v.begin(), v.end(), 1, std::multiplies<int>());

    auto dash_fold = [](std::string a, int b) {
        return std::move(a) + '-' + std::to_string(b);
    };

    std::string s = std::accumulate(std::next(v.begin()), v.end(),
        std::to_string(v[0]), // start with first element
        dash_fold);

    // Right fold using reverse iterators
    std::string rs = std::accumulate(std::next(v.rbegin()), v.rend(),
        std::to_string(v.back()), // start with last element
        dash_fold);

    std::cout << "sum: " << sum << '\n'
        << "product: " << product << '\n'
        << "dash-separated string: " << s << '\n'
        << "dash-separated string (right-folded): " << rs << '\n';


    std::vector<std::pair<int, double>> v2{ std::make_pair(1, 1.0), std::make_pair(2,2.0), std::make_pair(3, 3.0), std::make_pair(4,4.0), std::make_pair(5, 5.0), std::make_pair(6, 6.0), std::make_pair(7,7.0), std::make_pair(8, 8.0), std::make_pair(9, 9.0), std::make_pair(10,10.0) };
   // std::vector<std::pair<int, double>> v2{ std::make_pair(1, 1.0), std::make_pair(1,1.0), std::make_pair(2, 2.3), std::make_pair(1,1.0), std::make_pair(1, 1.0), std::make_pair(1, 1.0), std::make_pair(1,1.0), std::make_pair(1, 1.0), std::make_pair(1, 1.0), std::make_pair(1,1.0) };
    std::vector<std::tuple<int, double>> v3{ std::make_tuple(1, 1.0), std::make_tuple(2,2.0), std::make_tuple(3, 3.0), std::make_tuple(4,4.0), std::make_tuple(5, 5.0), std::make_tuple(6, 6.0), std::make_tuple(7,7.0), std::make_tuple(8, 8.0), std::make_tuple(9, 9.0), std::make_tuple(10,10.0) };

    std::vector<double> vd2;
    std::cout << "vd2 size: " << vd2.size() << '\n';
    std::transform(v2.begin(), v2.end(), std::back_inserter(vd2), [](const std::pair<int, double>& p) -> double { return p.second; });


    std::vector<double> vd3;
    std::cout << "vd3 size: " << vd3.size() << '\n';
    std::transform(v3.begin(), v3.end(), std::back_inserter(vd3), [](const std::tuple<int, double>& p) -> double { return std::get<1>(p); });
    std::vector<std::tuple<int, double>> usedPoints;
    std::copy(v3.begin() + 3, v3.end(), std::back_inserter(usedPoints));

    auto sum_d = std::accumulate(vd2.begin(), vd2.end(), 0.0);
    if (vd2.size()) sum_d /= vd2.size();
    std::cout << "vd2 size2: " << vd2.size() << '\n';
    std::cout << "sum: " << sum_d << '\n';
    auto minel = std::min_element(vd2.begin(), vd2.end());
    std::cout << "min: " << *minel << '\n';
    std::vector<double> maxd;
    std::copy_if(vd2.begin(), vd2.end(), back_inserter(maxd), [sum_d](const double a) {return a < 2 * sum_d; });
    auto pointsMinLastsTwo = *std::min_element(vd2.end()-2, vd2.end());
    std::vector<double> pointsFirstTwo = { *vd2.begin(),  *(vd2.begin() + 1) };

    /*
    std::vector<int> v;
    std::transform(test_vector.begin(), test_vector.end(), std::back_inserter(v),
        [&v](const std::pair<int, int>& p)
        { v.push_back(p.first);
    return p.second; });
*/
 //   std::cout << "Result: " << s2 << std::endl;
 //   auto sum2 = std::accumulate(v2.begin(), v2.end(), std::make_pair(0, 0.0));

 //   auto product2 = std::accumulate(v.begin(), v.end(), std::make_pair(1,1.0), std::multiplies<std::pair<int,double>>());
/*
    auto dash_fold2 = [](std::string a, int b) {
        return std::move(a) + '-' + std::to_string(b);
    };

    std::string su2 = std::accumulate(std::next(v.begin()), v.end(),
        std::to_string(v[0]), // start with first element
        dash_fold2);

    // Right fold using reverse iterators
    std::string r2s = std::accumulate(std::next(v.rbegin()), v.rend(),
        std::to_string(v.back()), // start with last element
        dash_fold2);
*/
 //   std::cout << "sum: " << product2.first << '\n';
    /*
    std::cout << "sum: " << sum2 << '\n'
        << "product: " << product2 << '\n'
        << "dash-separated string: " << su2 << '\n'
        << "dash-separated string (right-folded): " << r2s << '\n';
        */



    std::vector<Size> sizes = { {4, 1}, {2, 3}, {1, 2} };

    decltype(sizes)::iterator minEl, maxEl;
    std::tie(minEl, maxEl) = std::minmax_element(begin(sizes), end(sizes),
        [](Size const& s1, Size const& s2)
        {
            return s1.width < s2.width;
        });

    std::cout << "Minimum (based on width): "
        << minEl->width << "," << minEl->height << std::endl;

    std::cout << "Maximum (based on width): "
        << maxEl->width << "," << maxEl->height << std::endl;


    intrpolation();
    return EXIT_SUCCESS;
}


