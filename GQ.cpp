#include "GQ.h"
#include <cmath>        // std::fabs
#include <iostream>

Gauss_integration::Gauss_integration(int N) {
    Num_of_points = N;
    p = new double*[N+1];
    dp = new double*[N+1];
    xi = new double*[N+1];
    wi = new double*[N+1];
    point_range = new double*[N+1];
}

Gauss_integration::~Gauss_integration() {
    for (int n = 0; n < Num_of_points; n++) {
        delete [] p[n];
        delete [] dp[n];
        delete [] xi[n];
        delete [] wi[n];
        delete [] point_range[n];
    }
    delete [] p;
    delete [] dp;
    delete [] xi;
    delete [] wi;
    delete [] point_range;

    p = nullptr;
    dp = nullptr;
    xi = nullptr;
    wi = nullptr;
    point_range = nullptr;
}

double Gauss_integration::polyn(const int n, double x) {
    if (n == 0)
        return 1.0;
    else if (n == 1)
        return x;
    else
        return (double)(((2*n-1)*x*polyn(n-1, x)-(n-1)*polyn(n-2, x))/n);
}

double Gauss_integration::dpolyn(const int n, double x) {
    if (n == 0)
        return 0.0;
    else if (n == 1)
        return 1.0;
    else
        return (double)((x*polyn(n, x)-polyn(n-1, x))*n/(x*x-1.0));
}

void Gauss_integration::data_base(const int N, double x_axis[]) {
    for (int n = 0; n < Num_of_points+1; n++) {
        p[n] = new double[N];
        if ((Num_of_points) >= 0)
            dp[n] = new double[N-2];
        else
            dp[n] = NULL;

        if (n != 0) {
            xi[n] = new double[n];
            wi[n] = new double[n]; 
        } else {
            xi[n] = NULL;
            wi[n] = NULL;
        }

        for (int i = 0; i < N; i++) {
            p[n][i] = polyn(n, x_axis[i]);
            if ((1 <= i) && (i <= (N-2)))
                dp[n][i-1] = dpolyn(n, x_axis[i]);
        } 
    }
}

void Gauss_integration::points_range(const int N, double x_axis[]) {
    for (int n = 0; n < (Num_of_points+1); n++) {
        if (n != 0) {
           if (n == 1) {
               point_range[n] = new double[n];
           }  
            else {
                point_range[n] = new double[n*2];
                int xr = 0;
                for (int i = 0; i < (N-1); i++) {
                    if (p[n][i]*p[n][i+1] < 0) {
                        point_range[n][xr] = x_axis[i];
                        point_range[n][xr+1] = x_axis[i+1];
                        xr += 2;
                    }    
                }
            } 
        }                      
    }
}

void Gauss_integration::sort(double arr[], const int n) {
    for (int i = 0; i < n; i++) {
        for (int j = i; j < n; j++) {
            if (arr[i] > arr[j]) {
                int temp = arr[i];
                arr[i] = arr[j];
                arr[j] = temp;
            }
        }
    }
}

void Gauss_integration::newton_method() {
    double TOL = 1E-8;
    for (int n = 1; n < Num_of_points; n++) {
        int m = 0;
        int j = 0;

        for (int i = 0; i < (n+1); i++) {
            double xo = point_range[n+1][m];

            while (true) {
                xo -= polyn(n+1, xo)/dpolyn(n+1, xo);
                if (std::fabs(polyn(n+1, xo)) <= TOL) {
                    xi[n+1][j] = xo;
                    j += 1;
                    m += 2;
                    break;
                }   
            }
        sort(xi[n], n);
        }   
    }
}

double Gauss_integration::w(const int n, const double xi) {
    double a = (1.0-xi*xi)*(dpolyn(n, xi))*(dpolyn(n, xi));
    return (double)2.0/a;
}

void Gauss_integration::calculate_weight() {
    for (int n = 1; n < (Num_of_points+1); n++) {
        for (int i = 0; i < n; i++) {
            wi[n][i] = this->w(n, xi[n][i]);
        }
    }  
}

double Gauss_integration::calculate_Func(double (*f)(double), const int a, const int b) {
    double I = 0.0;
    for (int i = this->Num_of_points; i < (this->Num_of_points+1); i++) {
        for (int j = 0; j < this->Num_of_points; j++) {
            double x = (b-a)*xi[i][j]/2.0+(b+a)/2.0;
            I += wi[i][j]*f(x);
        }
    }    
    return (double)(b-a)*I/2.0;
}
