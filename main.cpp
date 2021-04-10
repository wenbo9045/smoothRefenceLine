#include <iostream>
#include <vector>
#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;

using namespace std;

double basis_function(vector<std::size_t> knots, std::size_t j, std::size_t d, double u) {
    // Compute B[j,d](u) (all basis functions)
    if (d == 1) {
        //  When d = 1:
        // 	B[j,1](u) = (t[j] <= u < t[j+1]) ? 1 : 0
        if (knots.at(j) <= u && u < knots.at(j + 1)) return 1.0;
        else return 0.0;
    } else // if (d > 1)
    {
        // 	B[j,d](u) = ((u - t[j])/(t[j+d-1] - t[j])) * B[j,d-1](u) + ((t[j+d] - u)/(t[j+d] - t[j+1])) * B[j+1,d-1](u)
        double div1 = 1.0;
        double num1 = (u - knots[j]);
        double denum1 = (double) (knots[j + d - 1] - knots[j]);

        if (std::fabs(denum1) < 1.0e-8) div1 = 0.0;        // define (u/0) as 0 (to avoid devision by zero)
        else div1 = num1 / denum1;

        double div2 = 1.0;
        double num2 = (knots[j + d] - u);
        double denum2 = (double) (knots[j + d] - knots[j + 1]);

        if (std::fabs(denum2) < 1.0e-8) div2 = 0.0;
        else div2 = num2 / denum2;

        // definition of B-Spline basis function
        return (div1 * basis_function(knots, j, d - 1, u) + div2 * basis_function(knots, j + 1, d - 1, u));
    }
}

std::vector<std::pair<double, double>> Spline(const std::vector<std::pair<double, double>>& xy_lists_in, const std::size_t degree, const double intervel_u) {

    auto point_size = xy_lists_in.size();
    auto n = point_size - 1;    // number of control points - 1
    auto d = degree + 1;        // degree of the polynomial + 1

    auto knots_size = n + d + 1;                       // knots_size
    vector<std::size_t> knots;
    // knots between 'start' and 'm' will be updated, others require no updating
    for (std::size_t j = 0; j < knots_size; j++) {
        std::size_t value = 0;
        if (j < d) value = 0;        // if (j < d) 	t[j] = 0;
        else if (d <= j && j <= n) value = j - d + 1;    // if (d <= j && j <= n)	t[j] = j-d+1;
        else if (j > n) value = n - d + 2;    // if (j > n)				t[j] = n-d+2;
        knots.emplace_back(value);
    }

    // variable 'u' must be between t[d-1] and t[n+1])
    std::vector<std::pair<double, double>> xy_lists_out;
    for (double u = (double) knots[d - 1]; u <= (double) knots[n + 1]; u += intervel_u) {
        double x = 0.0; double y = 0.0;
        // definition of the B-Spline curve
        for (std::size_t j = 0; j < point_size; j++) {
            double b = basis_function(knots, j, d, u);
            x += xy_lists_in[j].first * b;
            y += xy_lists_in[j].second * b;
        }
        xy_lists_out.emplace_back(std::make_pair(x, y)); // add new segment to curve that will be drawn
    }

    return xy_lists_out;
}

int main()
{
    vector<double> x_set = {-1.0, 1.0, 3.0, 3.0, 3.0, 3.0, 5.0, 6.0, 8.0, 9.0};
    vector<double> y_set = {0.0, 0.0, 0.0, 1.0, 2.0, 3.0, 3.0, 3.0, 3.5, 5.0};
    std::vector<std::pair<double, double>> xy_lists_in;
    for(auto i=0; i<x_set.size(); i++){
        xy_lists_in.emplace_back(std::make_pair(x_set[i], y_set[i]));
    }
    plt::plot(x_set, y_set, "ob");

    x_set.clear();
    y_set.clear();
    auto xy_lists = Spline(xy_lists_in, 3,0.05);
    for (const auto& xy_point : xy_lists) {
        x_set.emplace_back(xy_point.first);
        y_set.emplace_back(xy_point.second);
    }
    plt::plot(x_set, y_set, "r");
    plt::axis("equal");
    plt::show();

}