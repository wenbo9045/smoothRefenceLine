#include "math_util.h"
#include <cmath>
#include <algorithm>

namespace xpilot {
namespace reference_line {
namespace math {

Eigen::Vector2d RotateVector2d(Eigen::Vector2d v_in, double theta) {
    double cos_theta = std::cos(theta);
    double sin_theta = std::sin(theta);

    auto x = cos_theta * v_in.x() - sin_theta * v_in.y();
    auto y = sin_theta * v_in.x() + cos_theta * v_in.y();

    return Eigen::Vector2d(x, y);
}

Eigen::Vector2d orthogonalOperation(Eigen::Vector2d v0, Eigen::Vector2d v1) {
    Eigen::Vector2d v;
    // multiply v1 by the dot product of this and v1 then divide it by v1's length
    v = v0 - v1 * InnerProd(v0, v1) / (v1.norm() * v1.norm());
    return v;
}

double CrossProd(Eigen::Vector2d v0, Eigen::Vector2d v1) {
    return v0.x() * v1.y() - v1.x() * v0.y();
}

double InnerProd(Eigen::Vector2d v0, Eigen::Vector2d v1) {
    return v0.x() * v1.x() + v0.y() * v1.y();
}


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

double lerp(double x0, double s0, double x1, double s1, double s) {
    double r = (s - s0) / (s1 - s0);
    double x = x0 + r * (x1 - x0);
    return x;
}

double slerp(double x0, double s0, double x1, double s1, double s) {
    double a0_n = normalizeRadian(x0);
    double a1_n = normalizeRadian(x1);
    double d = a1_n - a0_n;
    if (d > M_PI) {
        d = d - 2 * M_PI;
    } else if (d < -M_PI) {
        d = d + 2 * M_PI;
    }

    double r = (s - s0) / (s1 - s0);
    double x = a0_n + d * r;
    return normalizeRadian(x);
}


std::vector<double> ExtendVectorWithLinearInterpolation(const std::vector<double> &x, int extend_size) {
    // interplation example:
    // x: [x0, x1, x2], extend_size: 3
    // output: [y0(x0), y1, y2, y3(x1), y4, y5, y6(x2)]
    size_t origin_last = x.size() - 1;
    std::vector<double> res(origin_last * extend_size + 1, 0.0);

    for (size_t i = 0; i < origin_last * extend_size; ++i) {
        size_t idx0 = i / extend_size;
        size_t idx1 = idx0 + 1;
        double w = static_cast<double>(i % extend_size) / static_cast<double>(extend_size);
        res[i] = x[idx0] * (1.0 - w) + x[idx1] * w;
    }

    res.back() = x.back();
    return res;
}

Waypoints_T
InterpolatePointUsingLinearApproximation(const Waypoints_T &p0, const Waypoints_T &p1, double s) {
    Waypoints_T waypoint;
    waypoint.Location[0] = static_cast<float>(lerp(p0.Location[0], p0.TempS, p1.Location[0], p1.TempS, s));
    waypoint.Location[1] = static_cast<float>(lerp(p0.Location[1], p0.TempS, p1.Location[1], p1.TempS, s));
    waypoint.Location[2] = static_cast<float>(lerp(p0.Location[2], p0.TempS, p1.Location[2], p1.TempS, s));
    waypoint.Heading = static_cast<float>(slerp(p0.Heading, p0.TempS, p1.Heading, p1.TempS, s));
    waypoint.TempS = static_cast<float>(s);
    waypoint.Kappa = p0.Kappa;
    return waypoint;
}

ReferencePoint_T
InterpolatePointUsingLinearApproximation(const ReferencePoint_T &p0, const ReferencePoint_T &p1, double s) {
    ReferencePoint_T path_point;
    path_point.x = lerp(p0.x, p0.s, p1.x, p1.s, s);
    path_point.y = lerp(p0.y, p0.s, p1.y, p1.s, s);
    path_point.z = lerp(p0.z, p0.s, p1.z, p1.s, s);
    path_point.theta = slerp(p0.theta, p0.s, p1.theta, p1.s, s);
    path_point.kappa = lerp(p0.kappa, p0.s, p1.kappa, p1.s, s);
    path_point.dkappa = lerp(p0.dkappa, p0.s, p1.dkappa, p1.s, s);
    path_point.s = s;
    path_point.link_id = p0.link_id;
    path_point.left_boundary_exist = false;
    path_point.right_boundary_exist = false;
    path_point.left_foot_point.first = lerp(p0.left_foot_point.first, p0.s, p1.left_foot_point.first, p1.s, s);
    path_point.left_foot_point.second = lerp(p0.left_foot_point.second, p0.s, p1.left_foot_point.second, p1.s, s);
    path_point.right_foot_point.first = lerp(p0.right_foot_point.first, p0.s, p1.right_foot_point.first, p1.s, s);
    path_point.right_foot_point.second = lerp(p0.right_foot_point.second, p0.s, p1.right_foot_point.second, p1.s, s);
    path_point.left_boundary = lerp(p0.left_boundary, p0.s, p1.left_boundary, p1.s, s);
    path_point.right_boundary = lerp(p0.right_boundary, p0.s, p1.right_boundary, p1.s, s);
    return path_point;
}

std::pair<std::size_t, std::size_t> findIndex(std::vector<double> vector_s, double s) {
    auto upper_bound = std::upper_bound(vector_s.begin(), vector_s.end(), s);
    std::size_t dis_low = static_cast<std::size_t>(std::distance(vector_s.begin(), upper_bound)) - 1;
    std::size_t dis_up = static_cast<std::size_t>(std::distance(vector_s.begin(), upper_bound));
    if (dis_low < 0)
        dis_low = 0;
    if (dis_up > vector_s.size() - 1)
        dis_up = vector_s.size() - 1;
    if (dis_low == dis_up)
        dis_low = dis_low - 1;

    return std::make_pair(dis_low, dis_up);
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

void handleWayPoints(std::vector<Waypoints_T>& waypoints)
{
    if(waypoints.size() < 2)
        return;

    //1. handle reference points data
    float cur_s = 0.0;
    auto num_of_points = static_cast<int>(waypoints.size());
    for (int i = 1; i < num_of_points; i++) {
        auto &cur_point = waypoints[i];
        auto last_point = waypoints[i - 1];
        cur_point.Heading = std::atan2(cur_point.Location[1]-last_point.Location[1], cur_point.Location[0]-last_point.Location[0]);
        cur_s += std::hypot(cur_point.Location[0]-last_point.Location[0], cur_point.Location[1]-last_point.Location[1]);
        cur_point.TempS = cur_s;
    }
    waypoints.front().Heading = waypoints[1].Heading;//front
    waypoints.front().TempS = 0.0;

    // handle angle
    for (int i = 1; i < num_of_points; i++) {
        auto last_point = waypoints[i - 1];
        auto &cur_point = waypoints[i];
        cur_point.Heading = correctAngle(last_point.Heading, cur_point.Heading);
    }
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


std::vector<Waypoints_T> uniformWayPoints(std::vector<Waypoints_T>& waypoints, double delta_s)
{
    //1. correct angle
    if(waypoints.size() < 2)
        return waypoints;

    //1. correct angle
    handleWayPoints(waypoints);

    //2. Uniform reference points data
    std::vector<double> vector_s;
    for (const auto& p : waypoints) {
        vector_s.emplace_back(p.TempS);
    }
    std::size_t N = std::floor(vector_s.back() / delta_s);

    std::vector<Waypoints_T> traj_out;
    double s = 0;
    for ( std::size_t i=0; i<=N; i++ )
    {
        s = delta_s * i;
        auto index = findIndex(vector_s, s);
        auto p0 = waypoints[index.first];
        auto p1 = waypoints[index.second];
        auto refpoint = InterpolatePointUsingLinearApproximation(p0, p1, s);
        traj_out.emplace_back((refpoint));
    }
    traj_out.emplace_back(waypoints.back());//不要最后一个点

    return traj_out;
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

void InterpolateReferencePointsBoundary(std::vector<ReferencePoint_T>& refpoints, const std::vector<ReferencePoint_T>& original_points)
{
    //1. Uniform reference points data
    std::vector<double> vector_s;
    for (const auto& p : original_points) {
        vector_s.emplace_back(p.s);
    }

    for (auto& point : refpoints)
    {
        double s = point.s;
        auto index = findIndex(vector_s, s);
        auto p0 = original_points[index.first];
        auto p1 = original_points[index.second];

        point.link_id = p0.link_id;
        point.left_boundary_exist = false;
        point.right_boundary_exist = false;
        point.left_foot_point.first = lerp(p0.left_foot_point.first, p0.s, p1.left_foot_point.first, p1.s, s);
        point.left_foot_point.second = lerp(p0.left_foot_point.second, p0.s, p1.left_foot_point.second, p1.s, s);
        point.right_foot_point.first = lerp(p0.right_foot_point.first, p0.s, p1.right_foot_point.first, p1.s, s);
        point.right_foot_point.second = lerp(p0.right_foot_point.second, p0.s, p1.right_foot_point.second, p1.s, s);
        point.left_boundary = lerp(p0.left_boundary, p0.s, p1.left_boundary, p1.s, s);
        point.right_boundary = lerp(p0.right_boundary, p0.s, p1.right_boundary, p1.s, s);
    }
}


void
LineFitWithLeastSquares(const std::vector<double> &data_x, const std::vector<double> &data_y, double *a, double *b, double *c) {

    double A = 0.0, B = 0.0, C = 0.0;

    int xjPcount = data_x.size();
    double meanX = 0.0, meanY = 0.0;
    double sumXX = 0.0, sumXY = 0.0, sumYY = 0.0;
    for (int i = 0; i < xjPcount; i++)
    {
        meanX += data_x[i];
        meanY += data_y[i];

        sumXX += data_x[i] * data_x[i];
        sumXY += data_x[i] * data_y[i];
        sumYY += data_y[i] * data_y[i];
    }
    meanX /= xjPcount;
    meanY /= xjPcount;

    sumXX -= xjPcount * meanX*meanX;
    sumXY -= xjPcount * meanX*meanY;
    sumYY -= xjPcount * meanY*meanY;
    if (sumXX < MATH_FLAG_EPS) {
        A = 1.0;
        B = 0.0;
    }
    else
    {
        double ev = (sumXX + sumYY + sqrt((sumXX - sumYY)*(sumXX - sumYY) + 4 * sumXY*sumXY)) / 2.0;
        A = -sumXY;
        B = ev - sumYY;
        double norm = sqrt(A*A + B * B);
        A /= norm;
        B /= norm;
    }
    C = -(A * meanX + B * meanY);

    *a = A;
    *b = B;
    *c = C;
}

void CircleFitWithLeastSquares(const std::vector<double> &data_x, const std::vector<double> &data_y, double *radius) {

    double sum_x = 0.0, sum_y = 0.0;
    double sum_x2 = 0.0, sum_y2 = 0.0;
    double sum_x3 = 0.0, sum_y3 = 0.0;
    double sum_xy = 0.0, sum_x1y2 = 0.0, sum_x2y1 = 0.0;

    std::size_t data_n = data_x.size();
    for (std::size_t i = 0; i < data_n; i++) {
        double x = data_x[i];
        double y = data_y[i];
        double x2 = x * x;
        double y2 = y * y;
        sum_x += x;
        sum_y += y;
        sum_x2 += x2;
        sum_y2 += y2;
        sum_x3 += x2 * x;
        sum_y3 += y2 * y;
        sum_xy += x * y;
        sum_x1y2 += x * y2;
        sum_x2y1 += x2 * y;
    }

    double C, D, E, G, H;
    double a, b, c;

    C = data_n * sum_x2 - sum_x * sum_x;
    D = data_n * sum_xy - sum_x * sum_y;
    E = data_n * sum_x3 + data_n * sum_x1y2 - (sum_x2 + sum_y2) * sum_x;
    G = data_n * sum_y2 - sum_y * sum_y;
    H = data_n * sum_x2y1 + data_n * sum_y3 - (sum_x2 + sum_y2) * sum_y;
    a = (H * D - E * G) / (C * G - D * D + MATH_FLAG_EPS);
    b = (H * C - E * D) / (D * D - G * C + MATH_FLAG_EPS);
    c = -(a * sum_x + b * sum_y + sum_x2 + sum_y2) / data_n;

    *radius = sqrt(a * a + b * b - 4 * c) / 2;
    double center_x = a / (-2);
    double center_y = b / (-2);
}

bool IntersectionOfTwoCircles(const Eigen::Vector2d& c1, const double r1, const Eigen::Vector2d& c2, const double r2, Eigen::Vector2d* P1, Eigen::Vector2d* P2)
{
    double x1 = c1.x();
    double y1 = c1.y();

    double x2 = c2.x();
    double y2 = c2.y();

    double distance = std::hypot(x2-x1, y2-y1);
    if (r1 + r2 < distance)
        return false;

    double r1r1 = r1*r1;
    double x1x1 = x1*x1;
    double y1y1 = y1*y1;

    double x2x2 = x2*x2;
    double y2y2 = y2*y2;
    double r2r2 = r2*r2;

    double subs1 = x1x1 - 2 * x1*x2 + x2x2 + y1y1 - 2 * y1*y2 + y2y2;
    double subs2 = -r1r1 * x1 + r1r1 * x2 + r2r2 * x1 - r2r2 * x2 + x1x1*x1 - x1x1 * x2 - x1*x2x2 + x1*y1y1 - 2 * x1*y1*y2 + x1*y2y2 + x2x2*x2 + x2*y1y1 - 2 * x2*y1*y2 + x2*y2y2;
    double subs3 = -r1r1 * y1 + r1r1 * y2 + r2r2 * y1 - r2r2 * y2 + x1x1*y1 + x1x1 * y2 - 2 * x1*x2*y1 - 2 * x1*x2*y2 + x2x2 * y1 + x2x2 * y2 + y1y1*y1 - y1y1 * y2 - y1*y2y2 + y2y2*y2;
    double sigma = std::sqrt((r1r1 + 2 * r1*r2 + r2r2 - x1x1 + 2 * x1*x2 - x2x2 - y1y1 + 2 * y1*y2 - y2y2)*(-r1r1 + 2 * r1*r2 - r2r2 + subs1));
    P1->x() = (subs2 - sigma*y1 + sigma*y2) / (2 * subs1 + MATH_FLAG_EPS);
    P1->y() = (subs3 + sigma*x1 - sigma*x2) / (2 * subs1 + MATH_FLAG_EPS);

    P2->x() = (subs2 + sigma*y1 - sigma*y2) / (2 * subs1 + MATH_FLAG_EPS);
    P2->y() = (subs3 - sigma*x1 + sigma*x2) / (2 * subs1 + MATH_FLAG_EPS);

    return true;
}

std::pair<Eigen::Vector2d, double> CalcMinDistance(const std::pair<Eigen::Vector2d, Eigen::Vector2d> segment_line,
                                                   const Eigen::Vector2d pt) {
    Eigen::Vector2d p1 = segment_line.first;
    Eigen::Vector2d p2 = segment_line.second;
    // a = y2 - y1; b = x1 - x2; c = x2y1 - x1y2;
    // ax + by + c = 0; -bx + ay + d = 0
    double coeff_a = p2.y() - p1.y();
    double coeff_b = p1.x() - p2.x();
    double coeff_c = p2.x() * p1.y() - p1.x() * p2.y();

    double sqrt_temp = std::sqrt(pow(coeff_a, 2) + pow(coeff_b, 2));
    double dis_temp = (coeff_a * pt.x() + coeff_b * pt.y() + coeff_c);

    double min_dis = std::abs(dis_temp) / sqrt_temp;

    double foot_point_x = pt.x() - coeff_a * dis_temp / pow(sqrt_temp, 2);
    double foot_point_y = pt.y() - coeff_b * dis_temp / pow(sqrt_temp, 2);

    return std::make_pair(Eigen::Vector2d(foot_point_x, foot_point_y), min_dis);
}


bool FootPointIsInLine(const std::pair<Eigen::Vector2d, Eigen::Vector2d> segment_line,
                       const Eigen::Vector2d foot_point) {
    Eigen::Vector2d p1 = segment_line.first;
    Eigen::Vector2d p2 = segment_line.second;
    // PA.PB < 0
    double value_temp = (p1.x() - foot_point.x()) * (p2.x() - foot_point.x())
                        + (p1.y() - foot_point.y()) * (p2.y() - foot_point.y());

    if (value_temp <= 0) {
        return true;
    } else {
        return false;
    }
}


uint8_t CalcPointDirection(const std::pair<Eigen::Vector2d, Eigen::Vector2d> direction_segment,
                           const Eigen::Vector2d pt) {

    Eigen::Vector2d p_s = direction_segment.first;  //start p
    Eigen::Vector2d p_e = direction_segment.second; //end p
    //AC X BC >0
    double value_temp = (pt.x() - p_s.x()) * (p_e.y() - p_s.y()) - (p_e.x() - p_s.x()) * (pt.y() - p_s.y());

    if (value_temp >= 0) {
        return 1;   // on the right side
    } else {
        return 2;   // on the left side
    }
}

double CalcRelativeAngle(const std::pair<Eigen::Vector2d, Eigen::Vector2d> vector1,
                         const std::pair<Eigen::Vector2d, Eigen::Vector2d> vector2) {

    Eigen::Vector2d v1 = Eigen::Vector2d(vector1.second.x() - vector1.first.x(), vector1.second.y() - vector1.first.y());
    Eigen::Vector2d v2 = Eigen::Vector2d(vector2.second.x() - vector2.first.x(), vector2.second.y() - vector2.first.y());

    double a = v1.x() * v2.x() + v1.y() * v2.y();
    double b = std::sqrt(pow(v1.x(), 2.0) + pow(v1.y(),2.0)) * std::sqrt(pow(v2.x(), 2.0) + pow(v2.y(),2.0));

    return acos(a/b);
}

int RectClipSegment(Rect_T Rect,LineSeg_T Segm)
{
  int rt = 0;	//   输出初始化
  Eigen::Vector2d RectCentPt;
  Eigen::Vector2d e0,e1,tempVec1,tempVec2;
  double halfLen_e0, halfLen_e1;
  double x0,y0,x1,y1,t0,t1,xmin,ymin,xmax,ymax,dx,dy,x,y;
  double p1,p2,p3,p4,q1,q2,q3,q4;
  int rt1,rt2,rt3,rt4;

  RectCentPt.x() = (Rect.pt_fl.x() + Rect.pt_rr.x())*0.5;
  RectCentPt.y() = (Rect.pt_fl.y() + Rect.pt_rr.y())*0.5;
  e0.x() = Rect.pt_fr.x() - Rect.pt_fl.x();
  e0.y() = Rect.pt_fr.y() - Rect.pt_fl.y();
  e1.x() = Rect.pt_fr.x() - Rect.pt_rr.x();
  e1.y() = Rect.pt_fr.y() - Rect.pt_rr.y();
  halfLen_e0 = Get_Norm_Vec2(e0)*0.5;
  halfLen_e1 = Get_Norm_Vec2(e1)*0.5;
  Normalize_Vec2(e0);
  Normalize_Vec2(e1);

  //  将线段的两个端点转换到矩形中心轴对称坐标系下

  tempVec1.x() = Segm.pt1.x() - RectCentPt.x();
  tempVec1.y() = Segm.pt1.y() - RectCentPt.y();
  tempVec2.x() = Segm.pt2.x() - RectCentPt.x();
  tempVec2.y() = Segm.pt2.y() - RectCentPt.y();

  x0 = Get_Dot_Vec2(tempVec1,e0); //   点积
  y0 = Get_Dot_Vec2(tempVec1,e1);
  x1 = Get_Dot_Vec2(tempVec2,e0);
  y1 = Get_Dot_Vec2(tempVec2,e1);

  // 计算矩形四条边对线段的裁剪
  t0 = 0; t1 = 1; //  线段端点参数初始化为无效值
  xmin = -halfLen_e0; xmax = halfLen_e0;  //   矩形边界
  ymin = -halfLen_e1; ymax = halfLen_e1;
  dx = x1 - x0;
  dy = y1 - y0;
  p1 = -dx; q1 = x0 - xmin;
  p2 = dx; q2 = xmax - x0;
  p3 = -dy; q3 = y0 - ymin;
  p4 = dy; q4 = ymax - y0;

  rt1 = SegClip( p1, q1, &t0, &t1);
  rt2 = SegClip( p3, q3, &t0, &t1);
  if(rt1 == 1 && rt2 == 1)
  {
    rt3 = SegClip( p2, q2, &t0, &t1);
    rt4 = SegClip( p4, q4, &t0, &t1);
    if(rt3 == 1 && rt4 == 1)
    {
      rt = 1;
    }
  }
  return rt;
}

double Get_Norm_Vec2(Eigen::Vector2d vec)
{
  return std::sqrt(pow(vec.x(),2)+pow(vec.y(),2));
}

double Get_Dot_Vec2(Eigen::Vector2d vec1, Eigen::Vector2d vec2)
{
  return (vec1.x()*vec2.x()+vec1.y()*vec2.y());
}

int SegClip(double p,double q,double* t0, double* t1)
{
  double r;
  int rst = 1;  //  初始化
  if (p < 0 && q < 0)
  {
    r = q / p;
    if (r > *t1)
    {
      rst = 0;
    }
    else if(r > *t0 )
    {
      *t0 = r; //   更新
    }
  }
  else if (p > 0 && q < p)
  {
    r = q/p;
    if (r < *t0)
    {
      rst = 0;
    }
    else if (r < *t1)
    {
      *t1 = r;
    }
  }
  else if (p == 0 && q < 0)
  {
    rst = 0;
  }
  return rst;
}


void Normalize_Vec2(Eigen::Vector2d &vec)
{
  double VecLen = Get_Norm_Vec2(vec);
  if (VecLen < 1e-5f) {
    vec.x() = 0;
    vec.y() = 0;
  }
  else {
    vec.x() /= VecLen;
    vec.y() /= VecLen;
  }
}


void ThreeOrderSpline(double start_state[], double end_state[], double c, std::array<double,8> &coeff) {
  coeff[0] = 2.0 * start_state[0] - 2.0 * end_state[0] + c*cos(start_state[2]) + c*cos(end_state[2]);
  coeff[1] = -3.0 * start_state[0] + 3.0 * end_state[0] - 2.0*c*cos(start_state[2]) - c*cos(end_state[2]);
  coeff[2] = c*cos(start_state[2]);
  coeff[3] = start_state[0];
  coeff[4] = 2.0 * start_state[1] - 2.0 * end_state[1] + c*sin(start_state[2]) + c*sin(end_state[2]);
  coeff[5] = -3.0 * start_state[1] + 3.0 * end_state[1] - 2.0*c*sin(start_state[2]) - c*sin(end_state[2]);
  coeff[6] = c*sin(start_state[2]);
  coeff[7] = start_state[1];
/*
  coeff(1) = 2 * start_point(1) - 2 * final_point(1) + c*cos(start_point(3)) + c*cos(final_point(3));
  coeff(2) = -3 * start_point(1) + 3 * final_point(1) - 2*c*cos(start_point(3)) - c*cos(final_point(3));
  coeff(3) = c*cos(start_point(3));
  coeff(4) = start_point(1);
  coeff(5) = 2 * start_point(2) - 2 * final_point(2) + c*sin(start_point(3)) + c*sin(final_point(3));
  coeff(6) = -3 * start_point(2) + 3 * final_point(2) - 2*c*sin(start_point(3)) - c*sin(final_point(3));
  coeff(7) = c*sin(start_point(3));
  coeff(8) = start_point(2);
  */
}

void RLSPreProcess(const std::vector<double> &value_in, std::vector<double> &value_out) {

  double ry = 0.855;
  double sy = std::sqrt(ry);
  double k = 0.0f;
  double value_last;
  double value_new;

  value_out.clear();
  value_out.shrink_to_fit();

  value_last = value_in[0];


  for (int i = 0; i < value_in.size(); ++i) {
    k = sy * sy / (sy * sy + ry);
    value_new = value_last + k * (value_in[i] - value_last);
    sy = (sy - k * sy / (1 + sqrtf(ry / (sy * sy + ry)))) / sqrtf(ry);


    value_out.push_back(value_new);

    value_last = value_new;
  }

}

std::string ResetString(std::string str) {
  std::string str_out = str;
  std::string::iterator it;
  std::size_t str_size = str_out.size();
  std::size_t j = 0;
  for (it = str_out.begin(); it < str_out.end(); it++)
  {
    j++;
    if (j == str_size - 3)
    {
      str_out.erase(it, str_out.end());
      break;
    }
  }
  return str_out;
}

}
}
}
