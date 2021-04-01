
#include "cubic_spline_smoother.h"
#include "Eigen/Dense"

#include <iostream>

using namespace xpilot::reference_line::smoother;

CubicSplineSmoother::CubicSplineSmoother(const std::vector<ReferencePoint_T> &ref_points)
{
    ref_points_ = ref_points;
    num_of_points_ = ref_points_.size();

    gauss_points_[0] = 1.48874338981631210881e-01;
    gauss_points_[1] = -1.48874338981631210881e-01;
    gauss_points_[2] = 4.33395394129247190794e-01;
    gauss_points_[3] = -4.33395394129247190794e-01;
    gauss_points_[4] = 6.79409568299024406207e-01;
    gauss_points_[5] = -6.79409568299024406207e-01;
    gauss_points_[6] = 8.65063366688984510759e-01;
    gauss_points_[7] = -8.65063366688984510759e-01;
    gauss_points_[8] = 9.73906528517171720066e-01;
    gauss_points_[9] = -9.73906528517171720066e-01;

    gauss_point_weights_[0] = 2.95524224714752870187e-01;
    gauss_point_weights_[1] = 2.95524224714752870187e-01;
    gauss_point_weights_[2] = 2.69266719309996355105e-01;
    gauss_point_weights_[3] = 2.69266719309996355105e-01;
    gauss_point_weights_[4] = 2.19086362515982044000e-01;
    gauss_point_weights_[5] = 2.19086362515982044000e-01;
    gauss_point_weights_[6] = 1.49451349150580593150e-01;
    gauss_point_weights_[7] = 1.49451349150580593150e-01;
    gauss_point_weights_[8] = 6.66713443086881375920e-02;
    gauss_point_weights_[9] = 6.66713443086881375920e-02;

}

void CubicSplineSmoother::smooth()
{
    this->smooth_points_.clear();
    //insert a virtual points at end
    auto pt1 = this->ref_points_.back();
    auto pt2 = this->ref_points_[num_of_points_ - 2];
    ReferencePoint_T virtual_pt = {2 * pt1.x - pt2.x, 2 * pt1.y - pt2.y, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    this->ref_points_.emplace_back(virtual_pt);
    //fitting cubic B spline
    std::vector<CubicSplineParams> spline_params;
    for (std::size_t i = 0; i < num_of_points_ - 3; i++) {
        CubicSplineParams tmp;
        tmp.a3 =
                -(ref_points_[i].x - 3 * ref_points_[i + 1].x + 3 * ref_points_[i + 2].x -
                  ref_points_[i + 3].x)
                / 6;
        tmp.a2 = (ref_points_[i].x - 2 * ref_points_[i + 1].x + ref_points_[i + 2].x) / 2;
        tmp.a1 = -(ref_points_[i].x - ref_points_[i + 2].x) / 2;
        tmp.a0 = (ref_points_[i].x + 4 * ref_points_[i + 1].x + ref_points_[i + 2].x) / 6;
        tmp.b3 =
                -(ref_points_[i].y - 3 * ref_points_[i + 1].y + 3 * ref_points_[i + 2].y -
                  ref_points_[i + 3].y)
                / 6;
        tmp.b2 = (ref_points_[i].y - 2 * ref_points_[i + 1].y + ref_points_[i + 2].y) / 2;
        tmp.b1 = -(ref_points_[i].y - ref_points_[i + 2].y) / 2;
        tmp.b0 = (ref_points_[i].y + 4 * ref_points_[i + 1].y + ref_points_[i + 2].y) / 6;
        tmp.length = spline_length_calc(tmp.a0, tmp.a1, tmp.a2, tmp.a3, tmp.b0, tmp.b1, tmp.b2, tmp.b3, 0.0, 1.0);
        spline_params.emplace_back(tmp);
    }
    this->ref_points_.pop_back();
    std::size_t points_num = 1;//(std::size_t) (param.length * 20);
    //computing raw reference points
    for (std::size_t i = 0; i < spline_params.size(); i++) {
        auto param = spline_params[i];

        double delta_s = param.length / static_cast<double>(points_num);
        for (std::size_t j = 0; j <= points_num; j++) {
            //if current is not the last param, don't compute last point
            if ((i != spline_params.size() - 1) && j == points_num) {
                continue;
            }
            //std::cout<<spline_params.size()<< " " <<i<<" "<<points_num<< " " <<j<<std::endl;
            double t = static_cast<double>(j) / static_cast<double>(points_num);
            double t2 = t * t;
            double t3 = t2 * t;
            double x = param.a0 + param.a1 * t + param.a2 * t2 + param.a3 * t3;
            double y = param.b0 + param.b1 * t + param.b2 * t2 + param.b3 * t3;
            double z = t * (ref_points_[i + 2].z - ref_points_[i + 1].z) + ref_points_[i + 1].z;
            double dx = param.a1 + 2 * param.a2 * t + 3 * param.a3 * t2;
            double d2x = 2 * param.a2 + 6 * param.a3 * t;
            double d3x = 6 * param.a3;
            double dy = param.b1 + 2 * param.b2 * t + 3 * param.b3 * t2;
            double d2y = 2 * param.b2 + 6 * param.b3 * t;
            double d3y = 6 * param.b3;
            double theta = std::atan2(dy, dx);
            double tmp = dx * dx + dy * dy;
            double kappa = (dx * d2y - d2x * dy) / tmp / std::sqrt(tmp);
            double a = dx * d2y - dy * d2x;
            double b = dx * d3y - dy * d3x;
            double c = dx * d2x + dy * d2y;
            double d = dx * dx + dy * dy;
            double dkappa = (b * d - 3.0 * a * c) / (d * d * d);
            double ddkappa = 0.0;
            double s = 0;
            if (!this->smooth_points_.empty()) {
                s = this->smooth_points_.back().s + delta_s;
            }
            ReferencePoint_T tmp_pt = {x, y, z, theta, kappa, dkappa, ddkappa, s};
            this->smooth_points_.emplace_back(tmp_pt);
        }
    }
    double resolution = 0.1 / static_cast<double>(points_num);
    //connect the point
    auto connected_points = g2_continuity_connect_by_poly(ref_points_.front().x, ref_points_.front().y,
                                                          ref_points_.front().theta, ref_points_.front().kappa,
                                                          smooth_points_.front().x, smooth_points_.front().y,
                                                          smooth_points_.front().theta, smooth_points_.front().kappa, resolution);
    for (std::size_t i = 0; i < connected_points.size(); ++i) {
        double t = static_cast<double>(i) / static_cast<double>(connected_points.size());
        connected_points[i].z = t * (smooth_points_.front().z - ref_points_.front().z) + ref_points_.front().z;
        connected_points[i].s += ref_points_.front().s;
    }
    //insert connecting points
    for (auto &pt: this->smooth_points_) {
        pt.s += connected_points.back().s;
    }
    connected_points.pop_back();
    this->smooth_points_.insert(smooth_points_.begin(), connected_points.begin(), connected_points.end());
}


std::vector<ReferencePoint_T> CubicSplineSmoother::g2_continuity_connect_by_poly(const double x0,
                                                                             const double y0,
                                                                             const double theta0,
                                                                             const double kappa0,
                                                                             const double x1,
                                                                             const double y1,
                                                                             const double theta1,
                                                                             const double kappa1,
                                                                             const double resolution) {
    //first convert frame to the start point and compute new points
    double sin_theta = std::sin(theta0);
    double cos_theta = std::cos(theta0);
    double xf = cos_theta * (x1 - x0) + sin_theta * (y1 - y0);
    double yf = -sin_theta * (x1 - x0) + cos_theta * (y1 - y0);
    double theta1_ = theta1 - theta0;
    //compute 5th poly params
    double d2y0 = kappa0;
    //end x, y, dy, d2y
    double xf2 = xf * xf;
    double xf3 = xf2 * xf;
    double xf4 = xf3 * xf;
    double xf5 = xf4 * xf;
    double dyf = std::tan(theta1_);
    double tmp = 1.0 + dyf * dyf;
    double d2yf = kappa1 * tmp * std::sqrt(tmp);
    //compute path parameters
    double a0 = 0;
    double a1 = 0;
    double a2 = 0.5 * d2y0;
    Eigen::Matrix3d A;
    A << xf3, xf4, xf5,
            3.0 * xf2, 4.0 * xf3, 5.0 * xf4,
            6.0 * xf, 12.0 * xf2, 20.0 * xf3;
    Eigen::Vector3d b;
    b << yf - 0.5 * d2y0 * xf2, dyf - d2y0 * xf, d2yf - d2y0;
    auto ans = A.inverse() * b;
    double a3 = ans(0);
    double a4 = ans(1);
    double a5 = ans(2);
    //get points
    int points_num = std::max(1, (int) (xf / resolution));
    double point_distance = std::sqrt(xf * xf + yf * yf) / static_cast<double>(points_num);
    std::vector<ReferencePoint_T> result;
    for (int i = 0; i <= points_num; i++) {
        double t = xf * static_cast<double>(i) / static_cast<double>(points_num);
        double t2 = t * t;
        double t3 = t2 * t;
        double t4 = t3 * t;
        double t5 = t4 * t;
        double x_ = t;
        double y_ = a0 + a1 * t + a2 * t2 + a3 * t3 + a4 * t4 + a5 * t5;
        double dy_ = a1 + 2 * a2 * t + 3 * a3 * t2 + 4 * a4 * t3 + 5 * a5 * t4;
        double d2y_ = 2 * a2 + 6 * a3 * t + 12 * a4 * t2 + 20 * a5 * t3;
        double d3y_ = 6 * a3 + 24 * a4 * t + 60 * a5 * t2;
        double tmp_ = 1 + dy_ + dy_;
        double kappa_ = d2y_ / tmp / std::sqrt(tmp);
        double dkappa_ = (d3y_ * tmp_ - 3.0 * d2y_ * dy_ * d2y_) / (tmp_ * tmp_ * tmp_);
        double theta_ = std::atan(dy_);
        double s = 0;
        if (!result.empty()) {
            s = result.back().s + point_distance;
        }
        ReferencePoint_T tmp_pt = {x_, y_, 0.0, theta_, kappa_, dkappa_, 0.0, s};
        result.emplace_back(tmp_pt);
    }
    //transfer to previous frame
    for (auto &pt : result) {
        double new_x = cos_theta * pt.x - sin_theta * pt.y + x0;
        double new_y = sin_theta * pt.x + cos_theta * pt.y + y0;
        double new_theta = pt.theta + theta0;
        pt.x = new_x;
        pt.y = new_y;
        pt.theta = new_theta;
    }
    return result;
}

double CubicSplineSmoother::spline_length_calc(const double a0, const double a1, const double a2,
                                             const double a3, const double b0, const double b1,
                                             const double b2, const double b3, const double lower,
                                             const double upper) {

    double total_length = 0;
    for (int i = 0; i < 10; ++i) {
        double t = 0.5 * (upper - lower) * gauss_points_[i] + 0.5 * (upper + lower);
        double dx = a1 + 2 * a2 * t + 3 * a3 * t * t;
        double dy = b1 + 2 * b2 * t + 3 * b3 * t * t;
        total_length += 0.5 * (upper - lower) * gauss_point_weights_[i] * std::hypot(dx, dy);
    }
    return total_length;
}

std::vector<ReferencePoint_T> CubicSplineSmoother::send_smooth_points() {
    return this->smooth_points_;
}