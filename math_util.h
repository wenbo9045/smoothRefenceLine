#pragma once
#include <array>
#include <vector>
#include "Eigen/Dense"

static const double MATH_FLAG_EPS = 1e-8;

constexpr double kMathEpsilon = 1e-10;

double correctAngle(double last, double cur);

double normalizeAngle(double angle);

double normalizeRadian(double radian);

void handleReferencePoints(std::vector<ReferencePoint_T> &refpoints);

std::vector<ReferencePoint_T> uniformReferencePoints(std::vector<ReferencePoint_T>& refpoints, double delta_s);

std::vector<ReferencePoint_T> ComputePathProfileWithXYCoordinate(const std::vector<Eigen::Vector2d> &xy_points, const std::vector<ReferencePoint_T> &refpoints);

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
