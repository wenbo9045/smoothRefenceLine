#pragma once
#include <array>
#include <vector>
#include "Eigen/Dense"
#include "../../common/Rte_AP_REFERENCE_POINT_Types.h"
#include "../../common/Rte_AP_PKMAP_Types.h"

namespace xpilot {
namespace reference_line {
namespace math {

typedef struct _Rect_T
{
  Eigen::Vector2d pt_fr;
  Eigen::Vector2d pt_fl;
  Eigen::Vector2d pt_rl;
  Eigen::Vector2d pt_rr;
}Rect_T;

typedef struct _LineSeg_T
{
  Eigen::Vector2d pt1;
  Eigen::Vector2d pt2;
}LineSeg_T;

static const double MATH_FLAG_EPS = 1e-8;

constexpr double kMathEpsilon = 1e-10;

Eigen::Vector2d RotateVector2d(Eigen::Vector2d v_in, double theta);

Eigen::Vector2d orthogonalOperation(Eigen::Vector2d v0, Eigen::Vector2d v1);

double CrossProd(Eigen::Vector2d v0, Eigen::Vector2d v1);

double InnerProd(Eigen::Vector2d v0, Eigen::Vector2d v1);

double correctAngle(double last, double cur);

double normalizeAngle(double angle);

double normalizeRadian(double radian);

double lerp(double x0, double s0, double x1, double s1, double s);

double slerp(double x0, double s0, double x1, double s1, double s);

std::pair<std::size_t, std::size_t> findIndex(std::vector<double> vector_s, double s);

std::vector<double> ExtendVectorWithLinearInterpolation(const std::vector<double>& x, int extend_size);

Waypoints_T InterpolatePointUsingLinearApproximation(const Waypoints_T& p0, const Waypoints_T& p1, double s);

ReferencePoint_T InterpolatePointUsingLinearApproximation(const ReferencePoint_T& p0, const ReferencePoint_T& p1, double s);

void handleWayPoints(std::vector<Waypoints_T> &waypoints);

void handleReferencePoints(std::vector<ReferencePoint_T> &refpoints);

std::vector<Waypoints_T> uniformWayPoints(std::vector<Waypoints_T>& waypoints, double delta_s);

std::vector<ReferencePoint_T> uniformReferencePoints(std::vector<ReferencePoint_T>& refpoints, double delta_s);

void InterpolateReferencePointsBoundary(std::vector<ReferencePoint_T>& refpoints, const std::vector<ReferencePoint_T>& original_points);

std::vector<ReferencePoint_T> ComputePathProfileWithXYCoordinate(const std::vector<Eigen::Vector2d> &xy_points, const std::vector<ReferencePoint_T> &refpoints);

void LineFitWithLeastSquares(const std::vector<double> &data_x, const std::vector<double> &data_y, double *a, double *b, double *c);

void CircleFitWithLeastSquares(const std::vector<double> &data_x, const std::vector<double> &data_y, double *radius);

bool IntersectionOfTwoCircles(const Eigen::Vector2d& c1, const double r1, const Eigen::Vector2d& c2, const double r2, Eigen::Vector2d* P1, Eigen::Vector2d* P2);

std::pair<Eigen::Vector2d, double> CalcMinDistance(const std::pair<Eigen::Vector2d, Eigen::Vector2d> segment_line, const Eigen::Vector2d pt);

bool FootPointIsInLine(const std::pair<Eigen::Vector2d, Eigen::Vector2d> segment_line, const Eigen::Vector2d foot_point);

uint8_t CalcPointDirection(const std::pair<Eigen::Vector2d, Eigen::Vector2d> direction_segment, const Eigen::Vector2d pt);

double CalcRelativeAngle(const std::pair<Eigen::Vector2d, Eigen::Vector2d> vector1, const std::pair<Eigen::Vector2d, Eigen::Vector2d> vector2);

int RectClipSegment(Rect_T Rect,LineSeg_T Segm);

void Normalize_Vec2(Eigen::Vector2d &vec);

int SegClip(double p,double q,double* t0, double* t1);

double Get_Norm_Vec2(Eigen::Vector2d vec);

double Get_Dot_Vec2(Eigen::Vector2d vec1, Eigen::Vector2d vec2);

void ThreeOrderSpline(double start_state[], double end_state[], double c, std::array<double,8> &coeff);
void RLSPreProcess(const std::vector<double> &value_in, std::vector<double> &value_out);
std::string ResetString(std::string str);

/**
 * @brief Compute squared value.
 * @param value The target value to get its squared value.
 * @return Squared value of the input value.
 */
template <typename T>
inline T Square(const T value) {
    return value * value;
}

/**
 * @brief Clamp a value between two bounds.
 *        If the value goes beyond the bounds, return one of the bounds,
 *        otherwise, return the original value.
 * @param value The original value to be clamped.
 * @param bound1 One bound to clamp the value.
 * @param bound2 The other bound to clamp the value.
 * @return The clamped value.
 */
template <typename T>
T Clamp(const T value, T bound1, T bound2) {
    if (bound1 > bound2) {
        std::swap(bound1, bound2);
    }

    if (value < bound1) {
        return bound1;
    } else if (value > bound2) {
        return bound2;
    }
    return value;
}


template <typename T, int M, int N, typename D>
void DenseToCSCMatrix(const Eigen::Matrix<T, M, N> &dense_matrix, std::vector<T> *data, std::vector<D> *indices, std::vector<D> *indptr)
{
    static constexpr double epsilon = 1e-9;
    int data_count = 0;
    for (int c = 0; c < dense_matrix.cols(); ++c) {
        indptr->emplace_back(data_count);
        for (int r = 0; r < dense_matrix.rows(); ++r) {
            if (std::fabs(dense_matrix(r, c)) < epsilon) {
                continue;
            }
            data->emplace_back(dense_matrix(r, c));
            ++data_count;
            indices->emplace_back(r);
        }
    }
    indptr->emplace_back(data_count);
}



}
}
}
