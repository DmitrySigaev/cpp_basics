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
#include <assert.h>
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


/* Preconditions: X[0] <= x <= X[N-1], X[i] < X[i+1] for 0 <= i <= N-2
 * Postcondition: 0 <= i <= N-2 */
static int LookupIndex(double x, int N, const double* X) {
    assert(!isinf(x) && !isnan(x) && X[0] <= x && x <= X[N - 1]);

    /* Potentially pathological case */
    if (x == X[N - 1])
        return N - 2;

    /* First calculate what the index would be for uniformly spaced points */
    double I = (N - 1) * (x - X[0]) / (X[N - 1] - X[0]);
    int i = (int)I;

    /* Check this index and its immediate neighbors, in case the points are
     * actually uniform (or close to it) */
    if (X[i] <= x && x <= X[i + 1])
        return i;
    else if (X[i - 1] <= x && x < X[i])
        return i - 1;
    else if (X[i + 1] < x && x <= X[i + 2])
        return i + 1;
    else {
        /* Finally resort to a binary search starting from the initial guess */
        int ilo = 0;
        int ihi = N - 2;
        while (true) {
            if (x < X[i]) {
                ihi = i - 1;
                i = (ilo + ihi) / 2;
            }
            else if (x > X[i + 1]) {
                ilo = i + 1;
                i = (ilo + ihi) / 2;
            }
            else
                break;
        }
        return i;
    }
}

static double* new_array(int n) {
    return (double*)malloc(n * sizeof(double));
}

static double* new_array(int n, const double* src) {
    double* dest = (double*)malloc(n * sizeof(double));
    memcpy(dest, src, n * sizeof(double));
    return dest;
}


/*****************
 * Linear spline *
 *****************/

struct LinearSplineImpl : public SplineImpl {
    int N;
    double* X;
    double* Y;

    LinearSplineImpl() {
    }

    LinearSplineImpl(int nn, const double* xx, const double* yy) {
        N = nn;
        X = new_array(N, xx);
        Y = new_array(N, yy);

        xmin = X[0];
        xmax = X[N - 1];
    }

    virtual ~LinearSplineImpl() {
        free(X);
        free(Y);
    }

    double y(double x) const {
        int i = LookupIndex(x, N, X);
        assert(0 <= i && i <= N - 2);
        double x1 = X[i], x2 = X[i + 1], y1 = Y[i], y2 = Y[i + 1];
        return y1 + (y2 - y1) * (x - x1) / (x2 - x1);
    }

    double dydx(double x) const {
        int i = LookupIndex(x, N, X);
        assert(0 <= i && i <= N - 2);
        double x1 = X[i], x2 = X[i + 1], y1 = Y[i], y2 = Y[i + 1];
        return (y2 - y1) / (x2 - x1);
    }

    double max(double& xmax) const {
        double ymax = -1e100;
        for (int i = 0; i < N; i++) {
            if (Y[i] > ymax) {
                xmax = X[i];
                ymax = Y[i];
            }
        }
        return ymax;
    }

    double min(double& xmin) const {
        double ymin = 1e100;
        for (int i = 0; i < N; i++) {
            if (Y[i] < ymin) {
                xmin = X[i];
                ymin = Y[i];
            }
        }
        return ymin;
    }

    LinearSplineImpl* clone() const {
        return new LinearSplineImpl(N, X, Y);
    }
};

Spline LinearSpline(const std::vector<double>& X, const std::vector<double>& Y) {
    assert(X.size() == Y.size());
    return Spline(new LinearSplineImpl(X.size(), &X[0], &Y[0]));
}

Spline LinearSpline(int N, const double* X, const double* Y) {
    return Spline(new LinearSplineImpl(N, X, Y));
}


/**********************************************
 * Shifted linear spline
 *
 * Uses constant extrapolation at boundaries.
 *********************************************/

struct ShiftedLinearSplineImpl : public LinearSplineImpl {
    double tau;

    ShiftedLinearSplineImpl(int nn, const double* xx, const double* yy, double tau) {
        xmin = xx[0];
        xmax = xx[nn - 1];
        double T = xx[1] - xx[0];

#ifdef DEBUG
        /* (Shifted linear interpolation only works for uniformly spaced points.) */
        for (int i = 0; i < nn - 1; i++)
            assert(fabs(xx[i + 1] - xx[i] - T) < 1e-10);
#endif

        N = nn + 1;
        X = new_array(N);
        Y = new_array(N);
        X[0] = xx[0] - (1 - tau) * T;
        Y[0] = yy[0];
        for (int i = 0; i < nn; i++) {
            X[i + 1] = xx[i] + tau * T;
            Y[i + 1] = -tau / (1 - tau) * Y[i] + 1 / (1 - tau) * yy[i];
        }
    }
};

Spline ShiftedLinearSpline(const std::vector<double>& X, const std::vector<double>& Y, double tau) {
    assert(X.size() == Y.size());
    return Spline(new ShiftedLinearSplineImpl(X.size(), &X[0], &Y[0], tau));
}

Spline ShiftedLinearSpline(int N, const double* X, const double* Y, double tau) {
    return Spline(new ShiftedLinearSplineImpl(N, X, Y, tau));
}


/************************************************************
 * Cubic spline
 *
 * Implementation copied from netlib's fmm/spline.f program.
 ************************************************************/

struct CubicSplineImpl : public SplineImpl {
    /* y(x) = y[i] + b (x-x[i]) + c (x-x[i])^2 + d (x-x[i])^3  for  x[i] <= x < x[i+1] */
    int N;
    double* X;
    double* Y;
    double* b;
    double* c;
    double* d;

    CubicSplineImpl(int nn, const double* xx, const double* yy) {
        N = nn;
        assert(N >= 2);
        X = new_array(N, xx);
        Y = new_array(N, yy);
        b = new_array(N);
        c = new_array(N);
        d = new_array(N);

        xmin = X[0];
        xmax = X[N - 1];

        if (N == 2) {
            b[0] = b[1] = (Y[1] - Y[0]) / (X[1] - X[0]);
            c[0] = c[1] = 0;
            d[0] = d[1] = 0;
        }
        else {
            /* Set up tridiagonal system: b = diagonal, d = offdiagonal, c = right hand side */
            d[0] = X[1] - X[0];
            c[1] = (Y[1] - Y[0]) / d[0];
            for (int i = 1; i < N - 1; i++) {
                d[i] = X[i + 1] - X[i];
                b[i] = 2 * (d[i - 1] + d[i]);
                c[i + 1] = (Y[i + 1] - Y[i]) / d[i];
                c[i] = c[i + 1] - c[i];
            }

            /* End conditions: third derivatives at X[0] and X[N-1] obtained from divide differences */
            b[0] = -d[0];
            b[N - 1] = -d[N - 2];
            c[0] = 0;
            c[N - 1] = 0;
            if (N > 3) {
                c[0] = c[2] / (X[3] - X[1]) - c[1] / (X[2] - X[0]);
                c[N - 1] = c[N - 2] / (X[N - 1] - X[N - 3]) - c[N - 3] / (X[N - 2] - X[N - 4]);
                c[0] = c[0] * d[0] * d[0] / (X[3] - X[0]);
                c[N - 1] = -c[N - 1] * d[N - 1] * d[N - 1] / (X[N - 1] - X[N - 4]);
            }

            /* Forward elimination */
            for (int i = 1; i < N; i++) {
                double t = d[i - 1] / b[i - 1];
                b[i] -= t * d[i - 1];
                c[i] -= t * c[i - 1];
            }

            /* Back substitution */
            c[N - 1] = c[N - 1] / b[N - 1];
            for (int j = 0; j < N - 1; j++) {
                int i = N - j - 2;
                c[i] = (c[i] - d[i] * c[i + 1]) / b[i];
            }

            /* Compute polynomial coefficients */
            b[N - 1] = (Y[N - 1] - Y[N - 2]) / d[N - 2] + d[N - 2] * (c[N - 2] + 2 * c[N - 1]);
            for (int i = 0; i < N - 1; i++) {
                b[i] = (Y[i + 1] - Y[i]) / d[i] - d[i] * (c[i + 1] + 2 * c[i]);
                d[i] = (c[i + 1] - c[i]) / d[i];
                c[i] = 3 * c[i];
            }
            c[N - 1] = 3 * c[N - 1];
            d[N - 1] = d[N - 2];
        }
    }

    virtual ~CubicSplineImpl() {
        free(X);
        free(Y);
        free(b);
        free(c);
        free(d);
    }

    double y(double x) const {
        int i = LookupIndex(x, N, X);
        assert(0 <= i && i <= N - 2);
        double u = x - X[i];
        return Y[i] + b[i] * u + c[i] * u * u + d[i] * u * u * u;
    }

    double dydx(double x) const {
        int i = LookupIndex(x, N, X);
        assert(0 <= i && i <= N - 2);
        double u = x - X[i];
        return b[i] + 2 * c[i] * u + 3 * d[i] * u * u;
    }

    //    double max(int i, double& xmax) const {
    //        return 0;
    //    }

    CubicSplineImpl* clone() const {
        return new CubicSplineImpl(N, X, Y);
    }
};

Spline CubicSpline(const std::vector<double>& X, const std::vector<double>& Y) {
    assert(X.size() == Y.size());
    return Spline(new CubicSplineImpl(X.size(), &X[0], &Y[0]));
}

Spline CubicSpline(int N, const double* X, const double* Y) {
    return Spline(new CubicSplineImpl(N, X, Y));
}


/**********
 * Spline *
 **********/

Spline::Spline(const std::vector<double>& X, const std::vector<double>& Y) {
    assert(X.size() == Y.size());
    impl = new LinearSplineImpl(X.size(), &X[0], &Y[0]);
}

Spline::Spline(int n, const double* X, const double* Y) {
    impl = new LinearSplineImpl(n, X, Y);
}

Spline::Spline() {
    impl = NULL;
}

Spline::Spline(SplineImpl* otherimpl, bool clone) {
    if (clone && otherimpl != NULL)
        impl = otherimpl->clone();
    else
        impl = otherimpl;
}

Spline::Spline(const Spline& S) {
    if (S.impl)
        impl = S.impl->clone();
    else
        impl = NULL;
}

Spline::~Spline() {
    delete impl;
}

Spline& Spline::operator=(const Spline& S) {
    delete impl;
    impl = (S.impl) ? S.impl->clone() : NULL;
    return *this;
}

double Spline::Evaluate(double x) const {
    if (!impl || x < impl->xmin || x > impl->xmax)
        return 0;
    else
        return impl->y(x);
}

double Spline::EvaluateDerivative(double x) const {
    if (!impl || x < impl->xmin || x > impl->xmax)
        return 0;
    else
        return impl->dydx(x);
}

int intrpolation(void)
{
    int i;
    double xi, yi, x[10], y[10];
    std::vector<double> xPoints;
    std::vector<double> yPoints;

    printf("#m=0,S=17\n");

    for (i = 0; i < 10; i++)
    {

        x[i] = i + 0.5 * sin(i);
        xPoints.push_back(x[i]);
        y[i] = i + cos(i * i);
        yPoints.push_back(y[i]);
        printf("%f %f\n", x[i], y[i]);
    }

    printf("#m=1,S=0\n");

    {
        gsl_interp_accel* acc
            = gsl_interp_accel_alloc();
        gsl_spline* spline
            = gsl_spline_alloc(gsl_interp_linear, 10);

        gsl_spline_init(spline, x, y, 10);

        Spline sl = LinearSpline(xPoints, yPoints);

        for (xi = x[0]; xi < x[9]; xi += 0.01)
        {
            yi = gsl_spline_eval(spline, xi, acc);
            
            printf("%f %f %f\n", xi, yi, sl(xi));
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


    std::unordered_map<int, double> map1;
    std::unordered_map<int, double> *p_map = &map1;
 //   std::shared_ptr< std::unordered_map<int, double>> sp_map = std::make_shared(p_map);
    std::transform(sizes.begin(), sizes.end(), std::inserter(map1, map1.end()), [](Size const& s) -> std::unordered_map<int, double>::value_type { return std::unordered_map<int, double>::value_type( s.width, s.height); });
  //  From < https://stackoverflow.com/questions/4375180/how-to-insert-into-stdmap> 

    std::vector<std::pair<int, double>> points;
    auto factResPressures = map1;
  //  if (factResPressures)
     std::copy(factResPressures.begin(), factResPressures.end(), std::back_inserter(points));

     std::shared_ptr<std::unordered_map<int, double>> dmap(new std::unordered_map<int, double>);
     std::transform(sizes.begin(), sizes.end(), std::inserter(*dmap, dmap->end()), [](Size const& s) -> std::unordered_map<int, double>::value_type { return std::unordered_map<int, double>::value_type(s.width, s.height); });

     std::vector<std::pair<int, double>> points2;
     const auto factResPressures2 = dmap;
     if (factResPressures2)
     std::copy(factResPressures2->begin(), factResPressures2->end(), std::back_inserter(points2));
    return EXIT_SUCCESS;
}


