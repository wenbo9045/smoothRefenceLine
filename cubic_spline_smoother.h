#pragma once

#include <cmath>
#include <array>
#include <vector>

class CubicSplineSmoother{

    typedef struct {
        double a0;
        double a1;
        double a2;
        double a3;
        double b0;
        double b1;
        double b2;
        double b3;
        double length;
    } CubicSplineParams;

public:
    CubicSplineSmoother(const std::vector<ReferencePoint_T> &ref_points);
    ~CubicSplineSmoother() = default;
    void smooth();
    std::vector<ReferencePoint_T> send_smooth_points();

private:
    std::vector<ReferencePoint_T> g2_continuity_connect_by_poly(const double x0,
                                                              const double y0,
                                                              const double theta0,
                                                              const double kappa0,
                                                              const double x1,
                                                              const double y1,
                                                              const double theta1,
                                                              const double kappa1,
                                                              const double resolution);

    double spline_length_calc(const double a0, const double a1, const double a2, const double a3,
                              const double b0, const double b1, const double b2, const double b3,
                              const double lower, const double upper);

private:
    std::vector<ReferencePoint_T> ref_points_;
    std::size_t num_of_points_;
    std::vector<ReferencePoint_T> smooth_points_;

    std::array<double, 10> gauss_points_;
    std::array<double, 10> gauss_point_weights_;

};
