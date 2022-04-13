#include <iostream>
#include <cmath>
#include "GQ.cpp"

double f(double x);
void linspace(const int N, double axis[]);

int main() {
    int Num_of_points = 5;
    int a = 0;
    int b = 4;
    const int N = 101; // 1001

    double x_axis[N] = {0};
    linspace(N, x_axis);

    Gauss_integration gi(Num_of_points);
    gi.data_base(N, x_axis);
    gi.points_range(N, x_axis);
    gi.newton_method();
    gi.calculate_weight();

    double I = gi.calculate_Func(f, a, b);
    std::cout << "f(x) = " << I << std::endl;
    return 0;
}

double f(double x) {
    return x*std::exp(x);
}

void linspace(const int N, double axis[]) {
    int a = -1, b = 1;
    double d = (double)(b-a)/(N-1);
    for (int i = 0; i < N; i++) {
        axis[i] = (double)(a+i*d);
    }
}     