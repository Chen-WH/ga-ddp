#include "OC/ilqg.h"
#include "OC/model.hpp"

int main() {
  OC::iLQG* ilqg;
  Eigen::MatrixXd u0;
  Eigen::VectorXd x0;
  Eigen::VectorXd xd;

  // Make linear
  OC::Model* linear = new OC::Linear();

  // 初始化权重矩阵
  Eigen::MatrixXd Q_k = 0.01*Eigen::MatrixXd::Identity(linear->x_dims, linear->x_dims);
  Eigen::MatrixXd R_k = 0.001*Eigen::MatrixXd::Identity(linear->u_dims, linear->u_dims);
  Eigen::MatrixXd Q_T = 0.01*Eigen::MatrixXd::Identity(linear->x_dims, linear->x_dims);

  ilqg = new OC::iLQG(linear, Q_k, R_k, Q_T);
  ilqg->dt = 0.01;

  // Define initial state
  x0 = Eigen::VectorXd::Random(linear->x_dims);
  x0 << -0.3566, -1.5250, -0.6009, -0.5015, -0.3384, -0.2181, -1.8744, 1.8743, -0.1885, -1.3635;
  xd.resize(linear->x_dims);
  xd.setZero();

  // Make initialization for control sequence
  int T = 1000;
  u0 = 0.1*Eigen::MatrixXd::Ones(linear->u_dims, T);

  // Solve!
  std::cout << "Run iLQG!" << std::endl;
  auto start = std::chrono::system_clock::now();
  ilqg->generate_trajectory(x0, xd, u0);
  auto now = std::chrono::system_clock::now();
  long int elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(now - start).count();
  std::cout << "iLQG took: " << elapsed / 1000.0 << " seconds." << std::endl;

  return 0;
}