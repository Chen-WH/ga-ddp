#include "OC/ilqg.h"
#include "OC/model.hpp"
#include "Utils/cubic_spline.hpp"

int main() {
  GA::robot GA_rbt("jaka");
  OC::iLQG* ilqg;
  OC::Model* rbt = new OC::Robot("jaka");

  // 初始化权重矩阵
  Eigen::MatrixXd Q_k = 10*Eigen::MatrixXd::Identity(rbt->x_dims, rbt->x_dims);
  Q_k.block(6, 6, 6, 6) = Eigen::MatrixXd::Zero(6, 6);
  Eigen::MatrixXd R_k = Eigen::MatrixXd::Zero(rbt->u_dims, rbt->u_dims);
  Eigen::MatrixXd Q_T = 100 * Eigen::MatrixXd::Identity(rbt->x_dims, rbt->x_dims);
  ilqg = new OC::iLQG(rbt, Q_k, R_k, Q_T);
  ilqg->dt = 0.016;

  // Define initial state
  // Set random seed for reproducibility
  std::srand(static_cast<unsigned int>(std::time(nullptr)));
  Eigen::VectorXd q0 = Eigen::VectorXd::Zero(rbt->u_dims);
  q0 << 0.5944, 2.3504, -0.1554, 0.2633, -0.8514, 2.0885;
  Eigen::VectorXd dq0 = Eigen::VectorXd::Zero(rbt->u_dims);
  Eigen::VectorXd x0(rbt->x_dims);
  x0 << q0, dq0;
  Eigen::VectorXd qd = Eigen::VectorXd::Random(rbt->u_dims);
  qd << 1.0345, 0.591165, 1.3817, 0.501313, -1.0589, 0.789955;
  Eigen::VectorXd dqd = Eigen::VectorXd::Zero(rbt->u_dims);
  Eigen::VectorXd xd(rbt->x_dims);
  xd << qd, dqd;

  // Make initialization for control sequence
  int T = 100;
  Eigen::VectorXd x(4);
  x << 0, 0.15*T*ilqg->dt, 0.85*T*ilqg->dt, T*ilqg->dt;
  Eigen::VectorXd y(4);
  Eigen::VectorXd u = Eigen::VectorXd::LinSpaced(T + 1, 0, T*ilqg->dt);
  Eigen::MatrixXd q = Eigen::MatrixXd::Zero(rbt->u_dims, T + 1);
  Eigen::MatrixXd dq = Eigen::MatrixXd::Zero(rbt->u_dims, T + 1);
  Eigen::MatrixXd ddq = Eigen::MatrixXd::Zero(rbt->u_dims, T + 1);
  Eigen::MatrixXd tau = Eigen::MatrixXd::Zero(rbt->u_dims, T);
  Eigen::VectorXd q_row_i(T + 1);
  Eigen::VectorXd dq_row_i(T + 1);
  Eigen::VectorXd ddq_row_i(T + 1);

  for (int i = 0; i < rbt->u_dims; ++i) {
    y << q0(i), 0.95*q0(i) + 0.05*qd(i), 0.05*q0(i) + 0.95*qd(i), qd(i);
    CubicSpline(x, y, u, dq0(i), dqd(i), q_row_i, dq_row_i, ddq_row_i);
    q.row(i) = q_row_i;
    dq.row(i) = dq_row_i;
    ddq.row(i) = ddq_row_i;
  }
  for (int i = 0; i < T; ++i) {
    tau.col(i) = GA_rbt.idyn(q.col(i), dq.col(i), ddq.col(i), Eigen::VectorXd::Zero(rbt->u_dims));
  }

  // Solve!
  std::cout << "Run iLQG!" << std::endl;
  auto start = std::chrono::system_clock::now();
  ilqg->generate_trajectory(x0, xd, tau);
  auto now = std::chrono::system_clock::now();
  long int elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(now - start).count();
  std::cout << "iLQG took: " << elapsed / 1000.0 << " seconds." << std::endl;

  return 0;
}