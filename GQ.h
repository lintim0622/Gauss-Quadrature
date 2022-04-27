#ifndef GQ_H
#define GQ_H

class Gauss_integration {
private:
    int Num_of_points;
    double **p = nullptr;
    double **dp = nullptr;
    double **xi = nullptr;
    double **wi = nullptr;
    double **point_range = nullptr;

public:
    Gauss_integration(const int N);
    ~Gauss_integration();

    double polyn(const int n, double x);
    double dpolyn(const int n, double x);
    double w(const int n, const double xi);
    double calculate_Func(double (*f)(double x), const int a, const int b);

    void points_range(const int N, double x_axis[]);
    void data_base(const int N, double x_axis[]);
    void newton_method();
    void calculate_weight();
    void sort(double arr[], const int n);
};

#endif
