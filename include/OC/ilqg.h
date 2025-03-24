#ifndef __iLQG_core_h__
#define __iLQG_core_h__

#include "OC/model.hpp"
#include "OC/boxqp.h"
#include <memory>
#include <numeric>
#include <stdexcept>
#include <omp.h>

namespace OC {

class iLQG {
public:
  double dt = 0.01;      // time step
  int maxIter = 30;      // maximum number of iterations
  Eigen::MatrixXd xs;    // nx(T+1)
  Eigen::MatrixXd us;    // mxT

  iLQG(Model* p_dyn, Eigen::MatrixXd Q_k, Eigen::MatrixXd R_k, Eigen::MatrixXd Q_T);

  ~iLQG() = default;

  void generate_trajectory(Eigen::VectorXd& x_0, Eigen::VectorXd& x_d, Eigen::MatrixXd& u_0);

private:
  std::shared_ptr<Model> model;
  Eigen::MatrixXd Qk, Rk, QT;
  int T;        // number of state transitions
  static std::vector<double> Alpha;
  static constexpr double tolFun = 1e-7;
  static constexpr double tolGrad = 1e-4;
  double lambda = 1;
  double dlambda = 1;
  static constexpr double lambdaFactor = 1.6;
  static constexpr double lambdaMax = 1e10;
  static constexpr double lambdaMin = 1e-6;
  static constexpr double zMin = 0;

  Eigen::MatrixXd trace; // a trace of various convergence-related values [cost time_derivs time_forward time_backward]
  Eigen::VectorXd x0;
  Eigen::VectorXd xd;
  Eigen::MatrixXd xnew;  // nx(T+1)
  Eigen::MatrixXd unew;  // mxT
  double cost_s, cost_new;         // current cost

  // n = dims(state), m = dims(control)
  std::vector<Eigen::MatrixXd> fx; //nxnxT
  std::vector<Eigen::MatrixXd> fu; //nxmxT

  Eigen::Vector2d dV;
  Eigen::MatrixXd Vx; //nx(T+1)
  std::vector<Eigen::MatrixXd> Vxx; //nxnx(T+1)
  Eigen::MatrixXd k; //mxT
  std::vector<Eigen::MatrixXd> K; //mxnxT

  Eigen::VectorXd Qx, Qu, k_i;
  Eigen::MatrixXd Qxx, Qux, Quu, K_i, Qux_reg, QuuF;

  void init_traj();

  void forward_pass(const Eigen::MatrixXd& u);

  int backward_pass();
};

} // namespace OC

#endif
