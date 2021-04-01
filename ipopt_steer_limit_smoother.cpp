//
// Created by luoyy on 2020/7/20.
//
#include "ipopt_steer_limit_smoother.h"

#include "../../common/math_util.h"
#include "Eigen/Dense"

#define tag_f 1
#define tag_g 2
#define tag_L 3

#define weight_center 1.0e2
#define weight_smooth 1.0e5
#define weight_distance 1.0e3
#define kappa_bound 0.2

using namespace xpilot::reference_line::math;

namespace xpilot {
namespace reference_line {
namespace smoother {

/* Constructor. */
SteerLimitIpoptInterface::SteerLimitIpoptInterface(const std::vector<ReferencePoint_T>& ref_points) {

  ref_points_ = ref_points;
  num_of_points_ = ref_points_.size();

  x_loc_bounds_.clear();
  y_loc_bounds_.clear();
  center_points_.clear();

  const double bound_offset_value = 1.1;  // half_ego_width
  for (const auto& ref_point : ref_points_) {
    auto left_bound = std::max(0.0, ref_point.left_boundary - bound_offset_value);
    auto right_bound = std::max(0.0, ref_point.right_boundary - bound_offset_value);
    auto x_left = ref_point.x + left_bound * std::cos(ref_point.theta + M_PI_2);
    auto y_left = ref_point.y + left_bound * std::sin(ref_point.theta + M_PI_2);
    auto x_right = ref_point.x + right_bound * std::cos(ref_point.theta - M_PI_2);
    auto y_right = ref_point.y + right_bound * std::sin(ref_point.theta - M_PI_2);

    x_loc_bounds_.emplace_back(std::make_pair(std::min(x_left, x_right), std::max(x_left, x_right)));
    y_loc_bounds_.emplace_back(std::make_pair(std::min(y_left, y_right), std::max(y_left, y_right)));
    center_points_.emplace_back(std::make_pair((x_left + x_right) / 2.0, (y_left + y_right) / 2.0));
  }

}

// nlp 问题求解
bool SteerLimitIpoptInterface::get_nlp_info(int& n, int& m, int& nnz_jac_g, int& nnz_h_lag,
                                            IndexStyleEnum& index_style) {
  n = 3 * num_of_points_;        // 优化变量个数
  m = 2 * (num_of_points_ - 1);  // 约束变量个数

  generate_tapes(n, m, &nnz_jac_g, &nnz_h_lag);
  // use the C style indexing (0-based)
  index_style = C_STYLE;
  return true;
}

bool SteerLimitIpoptInterface::get_bounds_info(int n, double* x_l, double* x_u, int m, double* g_l, double* g_u) {
  // x y 上下界
  for (int i = 0; i < num_of_points_; i++) {
    // x
    x_l[i] = -1.0e20;//x_loc_bounds_[i].first;
    x_u[i] = 1.0e20;//x_loc_bounds_[i].second;
    // y
    x_l[i + num_of_points_] = -1.0e20;//y_loc_bounds_[i].first;
    x_u[i + num_of_points_] = 1.0e20;//y_loc_bounds_[i].second;
    // theta
    x_l[i + 2 * num_of_points_] = -1.0e20;//
    x_u[i + 2 * num_of_points_] = 1.0e20;//
  }

  // dx*cos(theta) - dy*sin(theta)
  for (int i = 0; i < num_of_points_ - 1; i++) {
    g_l[i] = -1e-8;
    g_u[i] = 1e-8;
  }

  // kappa
  for (int i = num_of_points_ - 1; i < 2 * (num_of_points_ - 1); i++) {
    g_l[i] = -kappa_bound;
    g_u[i] = kappa_bound;
  }

  return true;
}

bool SteerLimitIpoptInterface::get_starting_point(int n, bool init_x, double* x, bool init_z, double* z_L, double* z_U,
                                                  int m, bool init_lambda, double* lambda) {
  // Here, we assume we only have starting values for x.
  for (int i = 0; i < num_of_points_; i++) {
    const auto ref_point = ref_points_[i];
    x[i] = ref_point.x;
    x[num_of_points_ + i] = ref_point.y;
    x[2 * num_of_points_ + i] = ref_point.theta;
  }

  return true;
}

bool SteerLimitIpoptInterface::eval_f(int n, const double* x, bool new_x, double& obj_value) {
  eval_obj(n, x, &obj_value);
  return true;
}

bool SteerLimitIpoptInterface::eval_grad_f(int n, const double* x, bool new_x, double* grad_f) {
  gradient(tag_f, n, x, grad_f);
  return true;
}

bool SteerLimitIpoptInterface::eval_g(int n, const double* x, bool new_x, int m, double* g) {
  eval_constraints(n, x, m, g);
  return true;
}

bool SteerLimitIpoptInterface::eval_jac_g(int n, const double* x, bool new_x, int m, int nele_jac, int* iRow, int* jCol,
                                          double* values) {
  if (values == nullptr) {
    // return the structure of the jacobian
    for (int idx = 0; idx < nnz_jac_; idx++) {
      iRow[idx] = rind_g_[idx];
      jCol[idx] = cind_g_[idx];
    }
  } else {
    // return the values of the jacobian of the constraints
    sparse_jac(tag_g, m, n, 1, x, &nnz_jac_, &rind_g_, &cind_g_, &jacval_, options_g_);
    for (int idx = 0; idx < nnz_jac_; idx++) {
      values[idx] = jacval_[idx];
    }
  }
  return true;
}

bool SteerLimitIpoptInterface::eval_h(int n, const double* x, bool new_x, double obj_factor, int m,
                                      const double* lambda, bool new_lambda, int nele_hess, int* iRow, int* jCol,
                                      double* values) {
  if (values == nullptr) {
    // return the structure. This is a symmetric matrix, fill the lower left triangle only.
    for (int idx = 0; idx < nnz_L_; idx++) {
      iRow[idx] = rind_L_[idx];
      jCol[idx] = cind_L_[idx];
    }
  } else {
    // return the values. This is a symmetric matrix, fill the lower left triangle only
    obj_lam_[0] = obj_factor;
    for (int idx = 0; idx < m; idx++) {
      obj_lam_[1 + idx] = lambda[idx];
    }
    set_param_vec(tag_L, m + 1, &obj_lam_[0]);
    sparse_hess(tag_L, n, 1, const_cast<double*>(x), &nnz_L_, &rind_L_, &cind_L_, &hessval_, options_L_);

    for (int idx = 0; idx < nnz_L_; idx++) {
      values[idx] = hessval_[idx];
    }
  }
  return true;
}

void SteerLimitIpoptInterface::finalize_solution(Ipopt::SolverReturn status, int n, const double* x, const double* z_L,
                                                 const double* z_U, int m, const double* g, const double* lambda,
                                                 double obj_value, const Ipopt::IpoptData* ip_data,
                                                 Ipopt::IpoptCalculatedQuantities* ip_cq) {
  std::vector<Eigen::Vector2d> xy_points;
  for (int i = 0; i < num_of_points_; ++i) {
    auto point = Eigen::Vector2d(x[i], x[num_of_points_ + i]);
    xy_points.emplace_back(point);
  }
  smooth_points_ = ComputePathProfileWithXYCoordinate(xy_points, ref_points_);

  free(rind_g_);
  free(cind_g_);
  free(rind_L_);
  free(cind_L_);
  free(jacval_);
  free(hessval_);
}

template <class T>
bool SteerLimitIpoptInterface::eval_obj(int n, const T* x, T* obj_value) {
  T sum_center = 0;
  for (int i = 0; i < num_of_points_; i++) {
    const auto ref_point = ref_points_[i];
    T dx = ref_point.x - x[i];                   // dx
    T dy = ref_point.y - x[num_of_points_ + i];  // dy

    const auto center_point = center_points_[i];
    dx = center_point.first - x[i];                    // dx
    dy = center_point.second - x[num_of_points_ + i];  // dy

    sum_center += pow(dx, 2) + pow(dy, 2);
  }

  T sum_dis = 0;
  for (int i = 0; i < num_of_points_ - 1; i++) {
    sum_dis += (pow(x[i + 1] - x[i], 2) + pow(x[i + 1 + num_of_points_] - x[i + num_of_points_], 2));
  }

  T sum_smooth = 0;
  for (int i = 1; i < num_of_points_ - 1; i++) {
    sum_smooth += (pow(x[i - 1] + x[i + 1] - 2 * x[i], 2) +
                   pow(x[i - 1 + num_of_points_] + x[i + 1 + num_of_points_] - 2 * x[i + num_of_points_], 2));
  }

  // 优化目标
  *obj_value = weight_center * sum_center + weight_distance * sum_dis + weight_smooth * sum_smooth;
  return true;
}

template <class T>
bool SteerLimitIpoptInterface::eval_constraints(int n, const T* x, int m, T* g) {
  int counter = 0;
  // sin(theta)*dx = cos(theta)*dy
  for (int i = 0; i < num_of_points_ - 1; ++i) {
    g[counter] = (x[i + 1] - x[i]) * sin(x[2 * num_of_points_ + i]) -
                 (x[num_of_points_ + i + 1] - x[num_of_points_ + i]) * cos(x[2 * num_of_points_ + i]);
    ++counter;
  }  // n-1

  // kappa = dtheta / ds
  for (int i = 0; i < num_of_points_ - 1; ++i) {
    T ds = sqrt(pow(x[i + 1] - x[i], 2) + pow(x[num_of_points_ + i + 1] - x[num_of_points_ + i], 2));
    g[counter] = (x[2 * num_of_points_ + i + 1] - x[2 * num_of_points_ + i]) / (ds + kMathEpsilon);
    ++counter;
  }  // n-1

  return true;
}

//***************    ADOL-C part ***********************************

void SteerLimitIpoptInterface::generate_tapes(int n, int m, int* nnz_jac_g, int* nnz_h_lag) {
  std::vector<double> xp(n, 0.0);
  std::vector<double> lamp(m, 0.0);
  std::vector<double> zl(m, 0.0);
  std::vector<double> zu(m, 0.0);
  std::vector<adouble> xa(n, 0.0);
  std::vector<adouble> g(m, 0.0);
  std::vector<double> lam(m, 0.0);

  double sig;
  adouble obj_value;

  double dummy = 0.0;
  obj_lam_.clear();
  obj_lam_.resize(m + 1, 0.0);
  get_starting_point(n, 1, &xp[0], 0, &zl[0], &zu[0], m, 0, &lamp[0]);

  // Trace on Objectives
  trace_on(tag_f);
  for (int idx = 0; idx < n; idx++) {
    xa[idx] <<= xp[idx];
  }
  eval_obj(n, &xa[0], &obj_value);
  obj_value >>= dummy;
  trace_off();

  // Trace on Jacobian
  trace_on(tag_g);
  for (int idx = 0; idx < n; idx++) {
    xa[idx] <<= xp[idx];
  }
  eval_constraints(n, &xa[0], m, &g[0]);
  for (int idx = 0; idx < m; idx++) {
    g[idx] >>= dummy;
  }
  trace_off();

  // Trace on Hessian
  trace_on(tag_L);
  for (int idx = 0; idx < n; idx++) {
    xa[idx] <<= xp[idx];
  }
  for (int idx = 0; idx < m; idx++) {
    lam[idx] = 1.0;
  }
  sig = 1.0;
  eval_obj(n, &xa[0], &obj_value);
  obj_value *= mkparam(sig);
  eval_constraints(n, &xa[0], m, &g[0]);
  for (int idx = 0; idx < m; idx++) {
    obj_value += g[idx] * mkparam(lam[idx]);
  }
  obj_value >>= dummy;
  trace_off();

  rind_g_ = nullptr;
  cind_g_ = nullptr;
  rind_L_ = nullptr;
  cind_L_ = nullptr;

  options_g_[0] = 0; /* sparsity pattern by index domains (default) */
  options_g_[1] = 0; /*                         safe mode (default) */
  options_g_[2] = 0;
  options_g_[3] = 0; /*                column compression (default) */

  jacval_ = nullptr;
  hessval_ = nullptr;

  sparse_jac(tag_g, m, n, 0, &xp[0], &nnz_jac_, &rind_g_, &cind_g_, &jacval_, options_g_);
  *nnz_jac_g = nnz_jac_;
  options_L_[0] = 0;
  options_L_[1] = 1;
  sparse_hess(tag_L, n, 0, &xp[0], &nnz_L_, &rind_L_, &cind_L_, &hessval_, options_L_);
  *nnz_h_lag = nnz_L_;
}

std::vector<ReferencePoint_T> SteerLimitIpoptInterface::send_smooth_points() { return this->smooth_points_; }

}  // namespace smoother
}  // namespace reference_line
}  // namespace xpilot