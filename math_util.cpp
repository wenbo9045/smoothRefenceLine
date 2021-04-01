#include "math_util.h"
#include <cmath>
#include <algorithm>

double correctAngle(double last, double cur) {
    double resAngle = cur;
    double delta_theta = cur - last;
    while (std::fabs(delta_theta) > M_PI) {
      if (delta_theta > 0) {
        resAngle -= 2 * M_PI;
      } else {
        resAngle += 2 * M_PI;
      }
      delta_theta = resAngle - last;
    }

    return resAngle;
}

double normalizeAngle(double angle) {
    double radian = angle / 180.0 * M_PI;
    double a = std::fmod(radian + M_PI, 2.0 * M_PI);
    if (a < 0.0) {
        a += (2.0 * M_PI);
    }
    return a - M_PI;
}

double normalizeRadian(double radian) {
    double a = std::fmod(radian + M_PI, 2.0 * M_PI);
    if (a < 0.0) {
        a += (2.0 * M_PI);
    }
    return a - M_PI;
}


std::vector<ReferencePoint_T> ComputePathProfileWithXYCoordinate(const std::vector<Eigen::Vector2d> &xy_points, const std::vector<ReferencePoint_T> &refpoints) {
    std::size_t points_size = xy_points.size();
    std::vector<double> thetas;
    // Get finite difference approximated dx and dy for heading
    for (std::size_t i = 0; i < points_size; ++i) {
        double x_delta = 0.0;
        double y_delta = 0.0;
        if (i == 0) {
            x_delta = (xy_points[i + 1].x() - xy_points[i].x());
            y_delta = (xy_points[i + 1].y() - xy_points[i].y());
        } else if (i == points_size - 1) {
            x_delta = (xy_points[i].x() - xy_points[i - 1].x());
            y_delta = (xy_points[i].y() - xy_points[i - 1].y());
        } else {
            x_delta = 0.5 * (xy_points[i + 1].x() - xy_points[i - 1].x());
            y_delta = 0.5 * (xy_points[i + 1].y() - xy_points[i - 1].y());
        }
        thetas.emplace_back(std::atan2(y_delta, x_delta));
    }

    // Get linear interpolated s
    std::vector<double> accumulated_s;
    double distance = 0.0;
    for (std::size_t i = 0; i < points_size; ++i) {
        if (i == 0) {
            accumulated_s.emplace_back(distance);
        } else {
            auto p0 = xy_points[i - 1];
            auto p1 = xy_points[i];
            double ds = std::hypot(p1.y() - p0.y(), p1.x() - p0.x());
            distance += ds;
            accumulated_s.emplace_back(distance);
        }
    }

    // Get finite difference approximated first derivative of y and x respective to s
    std::vector<double> y_over_s_first_derivatives;
    std::vector<double> x_over_s_first_derivatives;
    for (std::size_t i = 0; i < points_size; ++i) {
        double xds = 0.0;
        double yds = 0.0;
        if (i == 0) {
            xds = (xy_points[i + 1].x() - xy_points[i].x()) / (accumulated_s[i + 1] - accumulated_s[i]);
            yds = (xy_points[i + 1].y() - xy_points[i].y()) / (accumulated_s[i + 1] - accumulated_s[i]);
        } else if (i == points_size - 1) {
            xds = (xy_points[i].x() - xy_points[i - 1].x()) / (accumulated_s[i] - accumulated_s[i - 1]);
            yds = (xy_points[i].y() - xy_points[i - 1].y()) / (accumulated_s[i] - accumulated_s[i - 1]);
        } else {
            xds = (xy_points[i + 1].x() - xy_points[i - 1].x()) / (accumulated_s[i + 1] - accumulated_s[i - 1]);
            yds = (xy_points[i + 1].y() - xy_points[i - 1].y()) / (accumulated_s[i + 1] - accumulated_s[i - 1]);
        }
        x_over_s_first_derivatives.push_back(xds);
        y_over_s_first_derivatives.push_back(yds);
    }

    // Get finite difference approximated second derivative of y and x respective to s
    std::vector<double> y_over_s_second_derivatives;
    std::vector<double> x_over_s_second_derivatives;
    for (std::size_t i = 0; i < points_size; ++i) {
        double xdds = 0.0;
        double ydds = 0.0;
        if (i == 0) {
            xdds = (x_over_s_first_derivatives[i + 1] - x_over_s_first_derivatives[i])
                   / (accumulated_s.at(i + 1) - accumulated_s.at(i));
            ydds = (y_over_s_first_derivatives[i + 1] - y_over_s_first_derivatives[i])
                   / (accumulated_s.at(i + 1) - accumulated_s.at(i));
        } else if (i == points_size - 1) {
            xdds = (x_over_s_first_derivatives[i] - x_over_s_first_derivatives[i - 1])
                   / (accumulated_s.at(i) - accumulated_s.at(i - 1));
            ydds = (y_over_s_first_derivatives[i] - y_over_s_first_derivatives[i - 1])
                   / (accumulated_s.at(i) - accumulated_s.at(i - 1));
        } else {
            xdds = (x_over_s_first_derivatives[i + 1] - x_over_s_first_derivatives[i - 1])
                   / (accumulated_s.at(i + 1) - accumulated_s.at(i - 1));
            ydds = (y_over_s_first_derivatives[i + 1] - y_over_s_first_derivatives[i - 1])
                   / (accumulated_s.at(i + 1) - accumulated_s.at(i - 1));
        }
        x_over_s_second_derivatives.push_back(xdds);
        y_over_s_second_derivatives.push_back(ydds);
    }

    std::vector<double> kappas;
    for (std::size_t i = 0; i < points_size; ++i) {
        double xds = x_over_s_first_derivatives[i];
        double yds = y_over_s_first_derivatives[i];
        double xdds = x_over_s_second_derivatives[i];
        double ydds = y_over_s_second_derivatives[i];
        double kappa = (xds * ydds - yds * xdds) / (std::sqrt(xds * xds + yds * yds) * (xds * xds + yds * yds) + 1e-6);
        kappas.push_back(kappa);
    }

    // Dkappa calculation
    std::vector<double> dkappas;
    for (std::size_t i = 0; i < points_size; ++i) {
        double dkappa = 0.0;
        if (i == 0) {
            dkappa = (kappas.at(i + 1) - kappas.at(i)) / (accumulated_s.at(i + 1) - accumulated_s.at(i));
        } else if (i == points_size - 1) {
            dkappa = (kappas.at(i) - kappas.at(i - 1)) / (accumulated_s.at(i) - accumulated_s.at(i - 1));
        } else {
            dkappa = (kappas.at(i + 1) - kappas.at(i - 1)) / (accumulated_s.at(i + 1) - accumulated_s.at(i - 1));
        }
        dkappas.push_back(dkappa);
    }

    std::vector<ReferencePoint_T> smooth_points;
    for (std::size_t i = 0; i < points_size; ++i) {
        ReferencePoint_T point;
        auto ref_point = refpoints[i];
        point.x = xy_points[i].x();
        point.y = xy_points[i].y();
        point.z = ref_point.z;
        point.theta = thetas[i];
        point.kappa = kappas[i];
        point.dkappa = dkappas[i];
        point.s = accumulated_s[i];

        point.link_id = ref_point.link_id;
        point.left_boundary_exist = false;
        point.right_boundary_exist = false;
        point.left_foot_point = ref_point.left_foot_point;
        point.right_foot_point = ref_point.right_foot_point;
        point.left_boundary = ref_point.left_boundary;
        point.right_boundary = ref_point.right_boundary;

        smooth_points.push_back(point);
    }
    return smooth_points;
}

void handleReferencePoints(std::vector<ReferencePoint_T>& refpoints)
{
    if(refpoints.size() < 2)
        return;

    //1. handle reference points data
    auto num_of_points = static_cast<int>(refpoints.size());
    //theta and s 处理
    double cur_s = 0.0;
    for (int i = 1; i < num_of_points; i++) {
        auto &cur_point = refpoints[i];
        auto last_point = refpoints[i - 1];
        cur_point.theta = std::atan2(cur_point.y-last_point.y, cur_point.x-last_point.x);
        cur_s += std::hypot(cur_point.x-last_point.x, cur_point.y-last_point.y);
        cur_point.s = cur_s;
    }
    refpoints.front().theta = refpoints[1].theta;//front
    refpoints.front().s = 0.0;

    for (int i = 1; i < num_of_points; i++) {
        auto last_point = refpoints[i - 1];
        auto &cur_point = refpoints[i];
        cur_point.theta = correctAngle(last_point.theta, cur_point.theta);
    }

    //kappa处理
    double kappa_bound = 0.2;
    for (int i = 1; i < num_of_points; i++) {
        auto &cur_point = refpoints[i];
        auto last_point = refpoints[i - 1];
        double dtheta = cur_point.theta - last_point.theta;
        double ds = cur_point.s - last_point.s;
        //kappa = dtheta / ds;
        double kappa = dtheta / ds;
        kappa = fabs(kappa) > kappa_bound ? (kappa / fabs(kappa)) * kappa_bound : kappa;
        cur_point.kappa = kappa;
    }
    refpoints.front().kappa = refpoints[1].kappa;//front
}

std::vector<ReferencePoint_T> uniformReferencePoints(std::vector<ReferencePoint_T>& refpoints, double delta_s)
{
    //1. correct angle
    if(refpoints.size() < 2)
        return refpoints;

    handleReferencePoints(refpoints);
    //2. Uniform reference points data
    std::vector<double> vector_s;
    for (const auto& p : refpoints) {
        vector_s.emplace_back(p.s);
    }
    std::size_t N = std::floor(vector_s.back() / delta_s);

    std::vector<ReferencePoint_T> traj_out;
    double s = 0;
    for ( std::size_t i=0; i<=N; i++ )
    {
        s = delta_s * i;
        auto index = findIndex(vector_s, s);
        auto p0 = refpoints[index.first];
        auto p1 = refpoints[index.second];
        auto refpoint = InterpolatePointUsingLinearApproximation(p0, p1, s);
        traj_out.emplace_back((refpoint));
    }
    traj_out.emplace_back(refpoints.back());//不要最后一个点

    return traj_out;
}
