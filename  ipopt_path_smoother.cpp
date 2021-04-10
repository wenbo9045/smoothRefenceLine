//
// Created by luoyy on 2020/7/20.
//
#include "ipopt_path_smoother.h"

#include "math_util.h"

#define tag_f 1
#define tag_g 2
#define tag_L 3

#define weight_center 1.0e2
#define weight_smooth 1.0e5
#define weight_distance 1.0e3
#define weight_dkappa 1.0e4

#define weight_dl 1.0e1
#define weight_ddl 1.0e1
#define weight_dddl 1.0e1

#define kappa_bound 0.2

#define PATH_OPTIMIZATION_TYPE 4

using namespace xpilot::behavior_planning;

/* Constructor. */
PathSmoothIpoptInterface::PathSmoothIpoptInterface(const behavior_planning::MpInputMsgs& mp_input_msgs) {
  auto bp_output = mp_input_msgs.bp_output_msg();
  auto motion_policy = bp_output.motion_policies()[0];
  reference_line_ = std::vector<PathPoint>(motion_policy.reference_line_msg().path_point().begin(),
                                           motion_policy.reference_line_msg().path_point().end());

  origin_sl_points_ = std::vector<SLPoint>(motion_policy.sl_points().begin(), motion_policy.sl_points().end());
  std::vector<planning::PathPoint> path_points;
  matched_reference_points_.clear();
  origin_path_points_.clear();
  x_loc_bounds_.clear();
  y_loc_bounds_.clear();

  const double bound_offset_value = 1.1;  // half_ego_width
  for (const auto sl_point : origin_sl_points_) {
    auto matched_reference_point = math::PathMatcher::MatchToPath(reference_line_, sl_point.s());
    auto left_bound = std::max(0.0, matched_reference_point.left_width() - bound_offset_value);
    auto right_bound = std::max(0.0, matched_reference_point.right_width() - bound_offset_value);
    auto x_left = matched_reference_point.x() + left_bound * std::cos(matched_reference_point.theta() + M_PI_2);
    auto y_left = matched_reference_point.y() + left_bound * std::sin(matched_reference_point.theta() + M_PI_2);
    auto x_right = matched_reference_point.x() + right_bound * std::cos(matched_reference_point.theta() - M_PI_2);
    auto y_right = matched_reference_point.y() + right_bound * std::sin(matched_reference_point.theta() - M_PI_2);

    x_loc_bounds_.emplace_back(std::make_pair(std::min(x_left, x_right), std::max(x_left, x_right)));
    y_loc_bounds_.emplace_back(std::make_pair(std::min(y_left, y_right), std::max(y_left, y_right)));

    matched_reference_points_.emplace_back(matched_reference_point);
    std::array<double, 3> d_condition = {sl_point.l(), sl_point.dl(), sl_point.ddl()};
    PathPoint ego_point;
    planning::math::CartesianFrenetConverter::CartesianFrenetConverter::frenet_to_cartesian(
        matched_reference_point.x(), matched_reference_point.y(), matched_reference_point.theta(),
        matched_reference_point.kappa(), matched_reference_point.dkappa(), d_condition, &ego_point);
    origin_path_points_.emplace_back(ego_point);
  }
  num_of_points_ = origin_path_points_.size();
  XDLOG(INFO) << "num_of_points: " << num_of_points_;
  ds_ = origin_sl_points_[1].s() - origin_sl_points_.front().s();

  //handle max_passage
}

// nlp 问题求解
bool PathSmoothIpoptInterface::get_nlp_info(int& n, int& m, int& nnz_jac_g, int& nnz_h_lag,
                                            IndexStyleEnum& index_style) {
#if PATH_OPTIMIZATION_TYPE == 1
  n = 4 * num_of_points_;        //x, y, theta, kappa
  m = 2 * (num_of_points_ - 1);  // 约束变量个数
#elif PATH_OPTIMIZATION_TYPE == 2
  n = 3 * num_of_points_;        //x, y, theta
  m = 2 * (num_of_points_ - 1);  // 约束变量个数
#elif PATH_OPTIMIZATION_TYPE == 3
  n = 3 * num_of_points_;                         //l, dl, ddl
  m = num_of_points_ + 2 * (num_of_points_ - 1);  // 约束变量个数
#elif PATH_OPTIMIZATION_TYPE == 4
  n = num_of_points_;  //l
  m = num_of_points_;  // 约束变量个数
#endif

  generate_tapes(n, m, &nnz_jac_g, &nnz_h_lag);
  // use the C style indexing (0-based)
  index_style = C_STYLE;
  return true;
}

bool PathSmoothIpoptInterface::get_bounds_info(int n, double* x_l, double* x_u, int m, double* g_l, double* g_u) {
  for (int i = 0; i < num_of_points_; i++) {
#if PATH_OPTIMIZATION_TYPE == 1
    //x, y, theta, kappa
    if (i == 0) {
      x_l[i] = origin_path_points_.front().x() - 1e-3;
      x_u[i] = origin_path_points_.front().x() + 1e-3;
      x_l[i + num_of_points_] = origin_path_points_.front().y() - 1e-3;
      x_u[i + num_of_points_] = origin_path_points_.front().y() + 1e-3;
      x_l[i + 2 * num_of_points_] = origin_path_points_.front().theta() - 1e-3;
      x_u[i + 2 * num_of_points_] = origin_path_points_.front().theta() + 1e-3;
      x_l[i + 3 * num_of_points_] = origin_path_points_.front().kappa() - 1e-3;
      x_u[i + 3 * num_of_points_] = origin_path_points_.front().kappa() + 1e-3;
    } else {
      x_l[i] = -1.0e20;                   //x_loc_bounds_[i].first;//
      x_u[i] = 1.0e20;                    //x_loc_bounds_[i].second;//
      x_l[i + num_of_points_] = -1.0e20;  //y_loc_bounds_[i].first;//
      x_u[i + num_of_points_] = 1.0e20;   //y_loc_bounds_[i].second;//
      x_l[i + 2 * num_of_points_] = -1.0e20;
      x_u[i + 2 * num_of_points_] = 1.0e20;
      x_l[i + 3 * num_of_points_] = -kappa_bound;
      x_u[i + 3 * num_of_points_] = kappa_bound;
    }

#elif PATH_OPTIMIZATION_TYPE == 2
    //x, y, theta
    if (i == 0) {
      x_l[i] = origin_path_points_.front().x() - 1e-3;
      x_u[i] = origin_path_points_.front().x() + 1e-3;
      x_l[i + num_of_points_] = origin_path_points_.front().y() - 1e-3;
      x_u[i + num_of_points_] = origin_path_points_.front().y() + 1e-3;
      x_l[i + 2 * num_of_points_] = origin_path_points_.front().theta() - 1e-3;
      x_u[i + 2 * num_of_points_] = origin_path_points_.front().theta() + 1e-3;
    } else {
      x_l[i] = -1.0e20;
      x_u[i] = 1.0e20;
      x_l[i + num_of_points_] = -1.0e20;
      x_u[i + num_of_points_] = 1.0e20;
      x_l[i + 2 * num_of_points_] = -1.0e20;
      x_u[i + 2 * num_of_points_] = 1.0e20;
    }

#elif PATH_OPTIMIZATION_TYPE == 3
    //l, dl, ddl
    if (i == 0) {
      x_l[i] = origin_sl_points_.front().l() - 1e-3;
      x_u[i] = origin_sl_points_.front().l() + 1e-3;
      x_l[i + num_of_points_] = origin_sl_points_.front().dl() - 1e-3;
      x_u[i + num_of_points_] = origin_sl_points_.front().dl() + 1e-3;
      x_l[i + 2 * num_of_points_] = origin_sl_points_.front().ddl() - 1e-3;
      x_u[i + 2 * num_of_points_] = origin_sl_points_.front().ddl() + 1e-3;
    } else {
      x_l[i] = -1.0e20;
      x_u[i] = 1.0e20;
      x_l[i + num_of_points_] = -1.0e20;
      x_u[i + num_of_points_] = 1.0e20;
      x_l[i + 2 * num_of_points_] = -1.0e20;
      x_u[i + 2 * num_of_points_] = 1.0e20;
    }

#elif PATH_OPTIMIZATION_TYPE == 4
    //l
    if (i == 0) {
      x_l[i] = origin_sl_points_.front().l();
      x_u[i] = origin_sl_points_.front().l();
    } else {
      x_l[i] = -1.0e20;
      x_u[i] = 1.0e20;
    }
#endif
  }

#if PATH_OPTIMIZATION_TYPE == 1
  // dx*cos(theta) - dy*sin(theta)
  for (int i = 0; i < (num_of_points_ - 1); i++) {
    g_l[i] = -1e-8;
    g_u[i] = 1e-8;
  }
  // dtheta - kappa * ds
  for (int i = num_of_points_ - 1; i < 2 * (num_of_points_ - 1); i++) {
    g_l[i] = -1e-8;
    g_u[i] = 1e-8;
  }
#elif PATH_OPTIMIZATION_TYPE == 2
  // dx*cos(theta) - dy*sin(theta)
  for (int i = 0; i < (num_of_points_ - 1); i++) {
    g_l[i] = -1e-8;
    g_u[i] = 1e-8;
  }
  // kappa
  for (int i = num_of_points_ - 1; i < 2 * (num_of_points_ - 1); i++) {
    g_l[i] = -kappa_bound;
    g_u[i] = kappa_bound;
  }
#elif PATH_OPTIMIZATION_TYPE == 3
  // kappa
  for (int i = 0; i < num_of_points_; i++) {
    g_l[i] = -kappa_bound;
    g_u[i] = kappa_bound;
  }
  for (int i = num_of_points_; i < num_of_points_ + 2 * (num_of_points_ - 1); i++) {
    g_l[i] = -1e-8;
    g_u[i] = 1e-8;
  }
#elif PATH_OPTIMIZATION_TYPE == 4
  // kappa
  for (int i = 0; i < num_of_points_; i++) {
    g_l[i] = -kappa_bound;
    g_u[i] = kappa_bound;
  }
#endif

  return true;
}

bool PathSmoothIpoptInterface::get_starting_point(int n, bool init_x, double* x, bool init_z, double* z_L, double* z_U,
                                                  int m, bool init_lambda, double* lambda) {
  // Here, we assume we only have starting values for x.
  for (int i = 0; i < num_of_points_; i++) {
#if PATH_OPTIMIZATION_TYPE == 1
    const auto path_point = origin_path_points_[i];
    x[i] = path_point.x();
    x[num_of_points_ + i] = path_point.y();
    x[2 * num_of_points_ + i] = path_point.theta();
    x[3 * num_of_points_ + i] = path_point.kappa();
#elif PATH_OPTIMIZATION_TYPE == 2
    const auto path_point = origin_path_points_[i];
    x[i] = path_point.x();
    x[num_of_points_ + i] = path_point.y();
    x[2 * num_of_points_ + i] = path_point.theta();
#elif PATH_OPTIMIZATION_TYPE == 3
    const auto sl_point = origin_sl_points_[i];
    x[i] = sl_point.l();
    x[num_of_points_ + i] = sl_point.dl();
    x[2 * num_of_points_ + i] = sl_point.ddl();
#elif PATH_OPTIMIZATION_TYPE == 4
    const auto sl_point = origin_sl_points_[i];
    x[i] = sl_point.l();
#endif
  }

  return true;
}

bool PathSmoothIpoptInterface::eval_f(int n, const double* x, bool new_x, double& obj_value) {
  eval_obj(n, x, &obj_value);
  return true;
}

bool PathSmoothIpoptInterface::eval_grad_f(int n, const double* x, bool new_x, double* grad_f) {
  gradient(tag_f, n, x, grad_f);
  return true;
}

bool PathSmoothIpoptInterface::eval_g(int n, const double* x, bool new_x, int m, double* g) {
  eval_constraints(n, x, m, g);
  return true;
}

bool PathSmoothIpoptInterface::eval_jac_g(int n, const double* x, bool new_x, int m, int nele_jac, int* iRow, int* jCol,
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

bool PathSmoothIpoptInterface::eval_h(int n, const double* x, bool new_x, double obj_factor, int m,
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

void PathSmoothIpoptInterface::finalize_solution(Ipopt::SolverReturn status, int n, const double* x, const double* z_L,
                                                 const double* z_U, int m, const double* g, const double* lambda,
                                                 double obj_value, const Ipopt::IpoptData* ip_data,
                                                 Ipopt::IpoptCalculatedQuantities* ip_cq) {
  smooth_path_points_.clear();
#if PATH_OPTIMIZATION_TYPE == 1
  for (int i = 0; i < num_of_points_; ++i) {
    planning::PathPoint path_point;
    path_point.set_x(x[i]);
    path_point.set_y(x[num_of_points_ + i]);
    path_point.set_z(origin_path_points_[i].z());
    path_point.set_theta(x[2 * num_of_points_ + i]);
    path_point.set_kappa(x[3 * num_of_points_ + i]);

    smooth_path_points_.emplace_back(path_point);
  }
#elif PATH_OPTIMIZATION_TYPE == 2

  for (int i = 0; i < num_of_points_; ++i) {
    planning::PathPoint path_point;
    path_point.set_x(x[i]);
    path_point.set_y(x[num_of_points_ + i]);
    path_point.set_z(origin_path_points_[i].z());
    path_point.set_theta(x[2 * num_of_points_ + i]);
    smooth_path_points_.emplace_back(path_point);
  }
  double s = 0.0;
  smooth_path_points_.front().set_s(s);
  for (int i = 1; i < num_of_points_; ++i) {
    auto last_point = smooth_path_points_[i - 1];
    auto& cur_point = smooth_path_points_[i];
    s += std::hypot(cur_point.x() - last_point.x(), cur_point.y() - last_point.y());
    cur_point.set_s(s);
  }
  //kappa处理
  for (int i = 1; i < num_of_points_; i++) {
    auto& cur_point = smooth_path_points_[i];
    auto last_point = smooth_path_points_[i - 1];
    double dtheta = cur_point.theta() - last_point.theta();
    double ds = cur_point.s() - last_point.s();
    //kappa = dtheta / ds;
    double kappa = dtheta / ds;
    //    kappa = fabs(kappa) > kappa_bound ? (kappa / fabs(kappa)) * kappa_bound : kappa;
    cur_point.set_kappa(kappa);
  }
  smooth_path_points_.front().set_kappa(smooth_path_points_[1].kappa());  //front
#elif PATH_OPTIMIZATION_TYPE == 3
  for (int i = 0; i < num_of_points_; ++i) {
    planning::PathPoint path_point;
    auto matched_reference_point = matched_reference_points_[i];
    std::array<double, 3> d_condition = {x[i], x[num_of_points_ + i], x[2 * num_of_points_ + i]};
    planning::math::CartesianFrenetConverter::CartesianFrenetConverter::frenet_to_cartesian(
        matched_reference_point.x(), matched_reference_point.y(), matched_reference_point.theta(),
        matched_reference_point.kappa(), matched_reference_point.dkappa(), d_condition, &path_point);

    smooth_path_points_.emplace_back(path_point);
  }
#elif PATH_OPTIMIZATION_TYPE == 4
  std::vector<double> l_list;
  for (std::size_t i = 0; i < num_of_points_; ++i) {
    l_list.emplace_back(x[i]);
  }
  std::vector<double> dl_list;
  for (std::size_t i = 0; i < num_of_points_; ++i) {
    double dl = 0.0;
    if (i == num_of_points_ - 1) {
      dl = (x[i] - x[i - 1]) / ds_;
    } else {
      dl = (x[i + 1] - x[i]) / ds_;
    }
    dl_list.emplace_back(dl);
  }
  std::vector<double> ddl_list;
  for (std::size_t i = 0; i < num_of_points_; ++i) {
    double ddl = 0.0;
    if (i == num_of_points_ - 1) {
      ddl = (x[i] - 2 * x[i - 1] + x[i - 2]) / std::pow(ds_, 2.0);
    } else if (i == num_of_points_ - 2) {
      ddl = (x[i + 1] - 2 * x[i] + x[i - 1]) / std::pow(ds_, 2.0);
    } else {
      ddl = (x[i + 2] - 2 * x[i + 1] + x[i]) / std::pow(ds_, 2.0);
    }
    ddl_list.emplace_back(ddl);
  }

  for (int i = 0; i < num_of_points_; ++i) {
    planning::PathPoint path_point;
    auto matched_reference_point = matched_reference_points_[i];
    std::array<double, 3> d_condition = {l_list[i], dl_list[i], ddl_list[i]};
    planning::math::CartesianFrenetConverter::CartesianFrenetConverter::frenet_to_cartesian(
        matched_reference_point.x(), matched_reference_point.y(), matched_reference_point.theta(),
        matched_reference_point.kappa(), matched_reference_point.dkappa(), d_condition, &path_point);

    smooth_path_points_.emplace_back(path_point);
  }

#endif
  free(rind_g_);
  free(cind_g_);
  free(rind_L_);
  free(cind_L_);
  free(jacval_);
  free(hessval_);
}

template <class T>
bool PathSmoothIpoptInterface::eval_obj(int n, const T* x, T* obj_value) {
#if PATH_OPTIMIZATION_TYPE == 1
  T sum_center = 0;
  for (int i = 0; i < num_of_points_; i++) {
    const auto path_point = origin_path_points_[i];
    T dx = path_point.x() - x[i];                   // dx
    T dy = path_point.y() - x[num_of_points_ + i];  // dy
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

  T sum_dkappa = 0;
  for (int i = 0; i < num_of_points_ - 1; i++) {
    T dkappa = x[3 * num_of_points_ + i + 1] - x[3 * num_of_points_ + i];
    sum_dkappa += pow(dkappa, 2);
  }
  // 优化目标
  *obj_value =
      weight_center * sum_center + weight_distance * sum_dis + weight_smooth * sum_smooth + weight_dkappa * sum_dkappa;

#elif PATH_OPTIMIZATION_TYPE == 2
  T sum_center = 0;
  for (int i = 0; i < num_of_points_; i++) {
    const auto path_point = origin_path_points_[i];
    T dx = path_point.x() - x[i];                   // dx
    T dy = path_point.y() - x[num_of_points_ + i];  // dy
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

#elif PATH_OPTIMIZATION_TYPE == 3
  T sum_center = 0;
  for (int i = 0; i < num_of_points_; i++) {
    const auto sl_point = origin_sl_points_[i];
    // l_i - l_ref
    T dl = sl_point.l() - x[i];
    sum_center += pow(dl, 2);
  }

  T sum_dl = 0;
  for (int i = 0; i < num_of_points_; i++) {
    //dl
    sum_dl += pow(x[i + num_of_points_], 2);
  }

  T sum_ddl = 0;
  for (int i = 0; i < num_of_points_; i++) {
    //ddl
    sum_ddl += pow(x[i + 2 * num_of_points_], 2);
  }

  T sum_dddl = 0;
  for (int i = 0; i < num_of_points_ - 1; i++) {
    //ddl_i+1 - ddl_i
    sum_dddl += pow(x[i + 1 + 2 * num_of_points_] - x[i + 2 * num_of_points_], 2);
  }

  // 优化目标
  *obj_value = weight_center * sum_center + weight_dl * sum_dl + weight_ddl * sum_ddl + weight_dddl * sum_dddl;
#elif PATH_OPTIMIZATION_TYPE == 4
  T sum_center = 0;
  for (int i = 0; i < num_of_points_; i++) {
    const auto sl_point = origin_sl_points_[i];
    // l_i - l_ref
    sum_center += pow(sl_point.l() - x[i], 2);
  }

  T sum_dl = 0;
  for (int i = 0; i < num_of_points_ - 1; i++) {
    //dl = l_i+1 - l_i
    sum_dl += pow(x[i + 1] - x[i], 2);
  }

  T sum_ddl = 0;
  for (int i = 0; i < num_of_points_ - 2; i++) {
    //ddl = l_i+2 - 2*l_i+1 + l_i
    sum_ddl += pow(x[i + 2] - 2 * x[i + 1] + x[i], 2);
  }

  T sum_dddl = 0;
  for (int i = 0; i < num_of_points_ - 3; i++) {
    //dddl = l_i+3 - 3*l_i+2 + 3*l_i+1 - l_i
    sum_dddl += pow(x[i + 3] - 3 * x[i + 2] + 3 * x[i + 1] - x[i], 2);
  }

  // 优化目标
  *obj_value = weight_center * sum_center + weight_dl * sum_dl + weight_ddl * sum_ddl + weight_dddl * sum_dddl;

#endif

  return true;
}

template <class T>
T CalculateKappa(const double rkappa, const double rdkappa, const T l, const T dl, const T ddl) {
  T denominator = (dl * dl + (1 - l * rkappa) * (1 - l * rkappa));
  denominator = pow(denominator, 1.5);
  T numerator = rkappa + ddl - 2 * l * rkappa * rkappa - l * ddl * rkappa + l * l * rkappa * rkappa * rkappa +
                l * dl * rdkappa + 2 * dl * dl * rkappa;
  return numerator / (denominator + planning::math_util::kMathEpsilon);
}

template <class T>
bool PathSmoothIpoptInterface::eval_constraints(int n, const T* x, int m, T* g) {
  int counter = 0;

#if PATH_OPTIMIZATION_TYPE == 1
  //sin(theta)*dx = cos(theta)*dy
  for (int i = 0; i < num_of_points_ - 1; ++i) {
    g[counter] = (x[i + 1] - x[i]) * sin(x[2 * num_of_points_ + i]) -
                 (x[num_of_points_ + i + 1] - x[num_of_points_ + i]) * cos(x[2 * num_of_points_ + i]);
    ++counter;
  }  // n-1
  // dtheta = kappa * ds
  for (int i = 0; i < num_of_points_ - 1; ++i) {
    T ds = sqrt(pow(x[i + 1] - x[i], 2) + pow(x[num_of_points_ + i + 1] - x[num_of_points_ + i], 2));
    g[counter] = (x[2 * num_of_points_ + i + 1] - x[2 * num_of_points_ + i]) - x[3 * num_of_points_ + i] * ds;
    ++counter;
  }  // n-1
#elif PATH_OPTIMIZATION_TYPE == 2
  //sin(theta)*dx = cos(theta)*dy
  for (int i = 0; i < num_of_points_ - 1; ++i) {
    g[counter] = (x[i + 1] - x[i]) * sin(x[2 * num_of_points_ + i]) -
                 (x[num_of_points_ + i + 1] - x[num_of_points_ + i]) * cos(x[2 * num_of_points_ + i]);
    ++counter;
  }  // n-1
  // kappa = dtheta / ds
  for (int i = 0; i < num_of_points_ - 1; ++i) {
    T ds = sqrt(pow(x[i + 1] - x[i], 2) + pow(x[num_of_points_ + i + 1] - x[num_of_points_ + i], 2));
    g[counter] = (x[2 * num_of_points_ + i + 1] - x[2 * num_of_points_ + i]) / (ds + planning::math_util::kMathEpsilon);
    ++counter;
  }  // n-1
#elif PATH_OPTIMIZATION_TYPE == 3
  // kappa
  for (int i = 0; i < num_of_points_; ++i) {
    auto matched_point = matched_reference_points_[i];
    g[counter] = CalculateKappa(matched_point.kappa(), matched_point.dkappa(), x[i], x[num_of_points_ + i],
                                x[2 * num_of_points_ + i]);
    ++counter;
  }  // n

  //l_i+1 = l_i + dl * ds + 1/2 * ddl * ds^2 + 1/6 * (ddl_i+1 - ddl_i) * ds^3
  for (int i = 0; i < num_of_points_ - 1; i++) {
    g[counter] = x[i + 1] - x[i] - ds_ * x[i + num_of_points_] - 1.0 / 2.0 * pow(ds_, 2.0) * x[i + 2 * num_of_points_] -
                 1.0 / 6.0 * pow(ds_, 3.0) * (x[i + 1 + 2 * num_of_points_] - x[i + 2 * num_of_points_]);
    ++counter;
  }  // n-1

  //dl_i+1 - dl_i - ddl_i * ds - 1/2 (ddl_i+1 - ddl_i) * ds ^ 2
  for (int i = 0; i < num_of_points_ - 1; i++) {
    g[counter] = x[i + 1 + num_of_points_] - x[i + num_of_points_] - ds_ * x[i + 2 * num_of_points_] -
                 1.0 / 2.0 * pow(ds_, 2) * (x[i + 1 + 2 * num_of_points_] - x[i + 2 * num_of_points_]);
    ++counter;
  }  // n-1
#elif PATH_OPTIMIZATION_TYPE == 4
  // kappa
  for (int i = 0; i < num_of_points_; ++i) {
    auto matched_point = matched_reference_points_[i];
    //dl = (l_i+1 - l_i) / ds
    //ddl = (l_i+2 - 2*l_i+1 + l_i) / ds^2
    T l = x[i];
    T dl = 0.0;
    if (i == num_of_points_ - 1) {
      dl = (x[i] - x[i - 1]) / ds_;
    } else {
      dl = (x[i + 1] - x[i]) / ds_;
    }

    T ddl = 0.0;
    if (i == num_of_points_ - 1) {
      ddl = (x[i] - 2 * x[i - 1] + x[i - 2]) / std::pow(ds_, 2.0);
    } else if (i == num_of_points_ - 2) {
      ddl = (x[i + 1] - 2 * x[i] + x[i - 1]) / std::pow(ds_, 2.0);
    } else {
      ddl = (x[i + 2] - 2 * x[i + 1] + x[i]) / std::pow(ds_, 2.0);
    }

    g[counter] = CalculateKappa(matched_point.kappa(), matched_point.dkappa(), l, dl, ddl);
    ++counter;
  }  // n

#endif

  return true;
}

//***************    ADOL-C part ***********************************

void PathSmoothIpoptInterface::generate_tapes(int n, int m, int* nnz_jac_g, int* nnz_h_lag) {
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

void PathSmoothIpoptInterface::get_path_points(std::vector<planning::PathPoint>* origin_path_points,
                                               std::vector<planning::PathPoint>* smooth_path_points) {
  *origin_path_points = this->origin_path_points_;
  *smooth_path_points = this->smooth_path_points_;
}

