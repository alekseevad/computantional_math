#include <iostream>
#include <cstdlib>
#include <vector>

using std::pow;
using std::endl;
using std::vector;
using std::copy;
using std::cout;

double function(double x);
vector<double>& getNodes(double k, uint16_t n, double step);
double lagrange_interpolation(double x, const vector<double>& x_i, const vector<double>& f_xi, int n, int m);
void hermite_cubic_spline_value(int nn, double xn[], double fn[],
    double dn[], int n, double x[], double f[], double d[], double s[],
    double t[]);
void r8vec_bracket3(int n, double t[], double tval, int* left);
void hermite_cubic_value(double x1, double f1, double d1, double x2,
    double f2, double d2, int n, double x[], double f[], double d[],
    double s[], double t[]);
double derivative_value(double x);

void main()
{
    vector<double> nodes = getNodes( 0.0, 21, 1.0 );
    vector<double> functionValues = { };
    double* arr_dn = new double[nodes.size()];

    cout << "Values:" << endl;

    for (uint16_t i = 0; i < 21; ++i) {
        functionValues.push_back(function(nodes[i]));
        cout << "(" << nodes[i] << ";" << functionValues[i] << ")" << endl;
        arr_dn[i] = derivative_value(nodes[i]);
    }

    cout << endl << endl;

    cout << "Lagranzh: " << endl;

    for (uint16_t i = 0; i < 21; ++i) {
        cout << "(" << nodes[i] << ";" << lagrange_interpolation(nodes[i], nodes, functionValues, i, 21) << ")" << endl;
    }

    cout << endl << endl << "Derivative values: " << endl;

    for (uint16_t i = 0; i < 21; ++i) {
        cout << "(" << nodes[i] << ";" << arr_dn[i] << ")" << endl;
    }

    double* arr_xn = new double[nodes.size()];
    double* arr_fn = new double[functionValues.size()];
    double* arr_f = new double[21];
    double* arr_d = new double[21];
    double* arr_sd = new double[21];
    double* arr_td = new double[21];

    copy(nodes.begin(), nodes.end(), arr_xn);
    copy(functionValues.begin(), functionValues.end(), arr_fn);

    cout << endl << endl;
    
    cout << "Hermite cubic spline: " << endl;
    hermite_cubic_spline_value(21, arr_xn, arr_fn, arr_dn, 21, arr_xn, arr_f, arr_d, arr_sd, arr_td);
    for (uint16_t i = 0; i < 21; ++i) {
        cout << "(" << nodes[i] << ";" 
            <<  arr_f[i]
            << ")" << endl;
    }

    delete[] arr_xn;
    delete[] arr_dn;
    delete[] arr_fn;
}

double function(double x) {
    return 1 / (1 + 25 * pow(x, 2));
}

double derivative_value(double x) {
    return -50.0 * x / pow(25.0 * x * x + 1.0, 2.0);
}

vector<double>& getNodes(double k, uint16_t n, double step) {
    vector<double>* nodes = new vector<double>();
    while (n > 0) {
        nodes->push_back(-1 + 0.1 * k);
        k += step;
        --n;
    }
    return *nodes;
}

double lagrange_interpolation(double x, const vector<double>& x_i, const vector<double>& f_xi, int n, int m) {

    int h = 0;
    for (int i = 0; i < n; i++) {
        if ((x_i[i] > x) && (i >= m)) {
            h = i + 1 - m;
            break;
        }
    }

    double q = 0.0;
    for (int i = h; i < h + m; i++) {
        double p = 1.0;
        for (int j = h; j < h + m; j++) {
            if (j != i) {
                p *= x - x_i[j];
                p /= x_i[i] - x_i[j];
            }
        }
        p *= f_xi[i];
        q += p;
    }

    return q;
}

void hermite_cubic_spline_value(int nn, double xn[], double fn[],
    double dn[], int n, double x[], double f[], double d[], double s[],
    double t[]) {
        {
            int i;
            int left;

            left = n / 2;

            for (i = 0; i < n; i++)
            {
                r8vec_bracket3(nn, xn, x[i], &left);

                hermite_cubic_value(xn[left], fn[left], dn[left], xn[left + 1],
                    fn[left + 1], dn[left + 1], 1, x + i, f + i, d + i, s + i, t + i);
            }
            return;
        }
}

void r8vec_bracket3(int n, double t[], double tval, int* left) {
    int high;
    int low;
    int mid;
    //
    //  Check the input data.
    //
    if (n < 2)
    {
        std::cerr << "\n";
        std::cerr << "R8VEC_BRACKET3 - Fatal error!\n";
        std::cerr << "  N must be at least 2.\n";
        std::exit(1);
    }
    //
    //  If *LEFT is not between 0 and N-2, set it to the middle value.
    //
    if (*left < 0 || n - 2 < *left)
    {
        *left = (n - 1) / 2;
    }
    //
    //  CASE 1: TVAL < T[*LEFT]:
    //  Search for TVAL in (T[I],T[I+1]), for I = 0 to *LEFT-1.
    //
    if (tval < t[*left])
    {
        if (*left == 0)
        {
            return;
        }
        else if (*left == 1)
        {
            *left = 0;
            return;
        }
        else if (t[*left - 1] <= tval)
        {
            *left = *left - 1;
            return;
        }
        else if (tval <= t[1])
        {
            *left = 0;
            return;
        }
        //
        //  ...Binary search for TVAL in (T[I],T[I+1]), for I = 1 to *LEFT-2.
        //
        low = 1;
        high = *left - 2;

        for (; ; )
        {
            if (low == high)
            {
                *left = low;
                return;
            }

            mid = (low + high + 1) / 2;

            if (t[mid] <= tval)
            {
                low = mid;
            }
            else
            {
                high = mid - 1;
            }
        }
    }
    //
    //  CASE 2: T[*LEFT+1] < TVAL:
    //  Search for TVAL in (T[I],T[I+1]) for intervals I = *LEFT+1 to N-2.
    //
    else if (t[*left + 1] < tval)
    {
        if (*left == n - 2)
        {
            return;
        }
        else if (*left == n - 3)
        {
            *left = *left + 1;
            return;
        }
        else if (tval <= t[*left + 2])
        {
            *left = *left + 1;
            return;
        }
        else if (t[n - 2] <= tval)
        {
            *left = n - 2;
            return;
        }
        //
        //  ...Binary search for TVAL in (T[I],T[I+1]) for intervals I = *LEFT+2 to N-3.
        //
        low = *left + 2;
        high = n - 3;

        for (; ; )
        {

            if (low == high)
            {
                *left = low;
                return;
            }

            mid = (low + high + 1) / 2;

            if (t[mid] <= tval)
            {
                low = mid;
            }
            else
            {
                high = mid - 1;
            }
        }
    }
    //
    //  CASE 3: T[*LEFT] <= TVAL <= T[*LEFT+1]:
    //  T is just where the user said it might be.
    //
    else
    {
    }

    return;
}

void hermite_cubic_value(double x1, double f1, double d1, double x2,
    double f2, double d2, int n, double x[], double f[], double d[],
    double s[], double t[])

    //****************************************************************************80
    //
    //  Purpose:
    //
    //    HERMITE_CUBIC_VALUE evaluates a Hermite cubic polynomial.
    //
    //  Discussion:
    //
    //    The input arguments can be vectors.
    //
    //  Licensing:
    //
    //    This code is distributed under the GNU LGPL license.
    //
    //  Modified:
    //
    //    13 February 2011
    //
    //  Author:
    //
    //    John Burkardt
    //
    //  Reference:
    //
    //    Fred Fritsch, Ralph Carlson,
    //    Monotone Piecewise Cubic Interpolation,
    //    SIAM Journal on Numerical Analysis,
    //    Volume 17, Number 2, April 1980, pages 238-246.
    //
    //  Parameters:
    //
    //    Input, double X1, F1, D1, the left endpoint, function value
    //    and derivative.
    //
    //    Input, double X2, F2, D2, the right endpoint, function value
    //    and derivative.
    //
    //    Input, int N, the number of evaluation points.
    //
    //    Input, double X[N], the points at which the Hermite cubic
    //    is to be evaluated.
    //
    //    Output, double F[N], D[N], S[N], T[N], the value and first
    //    three derivatives of the Hermite cubic at X.
    //
{
    double c2;
    double c3;
    double df;
    double h;
    int i;

    h = x2 - x1;
    df = (f2 - f1) / h;

    c2 = -(2.0 * d1 - 3.0 * df + d2) / h;
    c3 = (d1 - 2.0 * df + d2) / h / h;

    for (i = 0; i < n; i++)
    {
        f[i] = f1 + (x[i] - x1) * (d1
            + (x[i] - x1) * (c2
                + (x[i] - x1) * c3));
        d[i] = d1 + (x[i] - x1) * (2.0 * c2
            + (x[i] - x1) * 3.0 * c3);
        s[i] = 2.0 * c2 + (x[i] - x1) * 6.0 * c3;
        t[i] = 6.0 * c3;
    }
    return;
}
