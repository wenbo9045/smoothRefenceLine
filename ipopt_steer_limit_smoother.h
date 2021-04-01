//
// Created by luoyy on 2020/7/20.
//

#pragma once

#include <tuple>

#include "adolc/adolc.h"
#include "adolc/adolc_sparse.h"
#include "adolc/adouble.h"
#include "coin/IpIpoptApplication.hpp"
#include "coin/IpIpoptCalculatedQuantities.hpp"
#include "coin/IpIpoptData.hpp"
#include "coin/IpOrigIpoptNLP.hpp"
#include "coin/IpSolveStatistics.hpp"
#include "coin/IpTNLP.hpp"
#include "coin/IpTNLPAdapter.hpp"
#include "coin/IpTypes.hpp"

class SteerLimitIpoptInterface : public Ipopt::TNLP {
 public:
  /** default constructor */
  SteerLimitIpoptInterface(const std::vector<ReferencePoint_T> &ref_points);;

  /** default destructor */
  virtual ~SteerLimitIpoptInterface() = default;

  /** Method to return some info about the nlp */
  bool get_nlp_info(int &n, int &m, int &nnz_jac_g, int &nnz_h_lag, IndexStyleEnum &index_style) override;

  /** Method to return the bounds for my problem */
  bool get_bounds_info(int n, double *x_l, double *x_u, int m, double *g_l, double *g_u) override;

  /** Method to return the starting point for the algorithm */
  bool get_starting_point(int n, bool init_x, double *x, bool init_z, double *z_L, double *z_U, int m, bool init_lambda,
                          double *lambda) override;

  /** Method to return the objective value */
  bool eval_f(int n, const double *x, bool new_x, double &obj_value) override;

  /** Method to return the gradient of the objective */
  bool eval_grad_f(int n, const double *x, bool new_x, double *grad_f) override;

  /** Method to return the constraint residuals */
  bool eval_g(int n, const double *x, bool new_x, int m, double *g) override;

  /** Method to return:
   *   1) The structure of the jacobian (if "values" is nullptr)
   *   2) The values of the jacobian (if "values" is not nullptr)
   */
  bool eval_jac_g(int n, const double *x, bool new_x, int m, int nele_jac, int *iRow, int *jCol,
                  double *values) override;

  /** Method to return:
   *   1) The structure of the hessian of the lagrangian (if "values" is
   * nullptr) 2) The values of the hessian of the lagrangian (if "values" is not
   * nullptr)
   */
  bool eval_h(int n, const double *x, bool new_x, double obj_factor, int m, const double *lambda, bool new_lambda,
              int nele_hess, int *iRow, int *jCol, double *values) override;

  /** @name Solution Methods */
  /** This method is called when the algorithm is complete so the TNLP can
   * store/write the solution */
  void finalize_solution(Ipopt::SolverReturn status, int n, const double *x, const double *z_L, const double *z_U,
                         int m, const double *g, const double *lambda, double obj_value,
                         const Ipopt::IpoptData *ip_data, Ipopt::IpoptCalculatedQuantities *ip_cq) override;

  //***************    start ADOL-C part ***********************************
  /** Template to return the objective value */
  template<class T>
  bool eval_obj(int n, const T *x, T *obj_value);

  /** Template to compute contraints */
  template <class T>
  bool eval_constraints(int n, const T *x, int m, T *g);

  /** Method to generate the required tapes */
  virtual void generate_tapes(int n, int m, int *nnz_jac_g, int *nnz_h_lag);
  //***************    end   ADOL-C part ***********************************

 private:
  /**@name Methods to block default compiler methods.
   * The compiler automatically generates the following three methods.
   *  Since the default compiler implementation is generally not what
   *  you want (for all but the most simple classes), we usually
   *  put the declarations of these methods in the private section
   *  and never implement them. This prevents the compiler from
   *  implementing an incorrect "default" behavior without us
   *  knowing. (See Scott Meyers book, "Effective C++")
   *
   */
  SteerLimitIpoptInterface(const SteerLimitIpoptInterface &);
  SteerLimitIpoptInterface &operator=(const SteerLimitIpoptInterface &);

  std::vector<double> obj_lam_;

  //** variables for sparsity exploitation
  unsigned int *rind_g_; /* row indices    */
  unsigned int *cind_g_; /* column indices */
  double *jacval_;       /* values         */
  unsigned int *rind_L_; /* row indices    */
  unsigned int *cind_L_; /* column indices */
  double *hessval_;      /* values */

  int nnz_jac_ = 0;
  int nnz_L_ = 0;
  int options_g_[4];
  int options_L_[4];

 private:
  std::vector<std::pair<double, double>> x_loc_bounds_;
  std::vector<std::pair<double, double>> y_loc_bounds_;
  std::vector<std::pair<double, double>> center_points_;

 private:
  int num_of_points_;
  std::vector<ReferencePoint_T> ref_points_;
  std::vector<ReferencePoint_T> smooth_points_;
};
