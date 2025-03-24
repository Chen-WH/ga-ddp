#include "OC/ilqg.h"

int main() {
  GA::robot GA_rbt("jaka");
  OC::iLQG* ilqg;
  OC::Model* rbt = new OC::Robot_kine("jaka");

  // 初始化权重矩阵
  Eigen::MatrixXd Q_k = Eigen::MatrixXd::Zero(rbt->x_dims, rbt->x_dims);
  Q_k.block(0, 0, 6, 6) = 10*Eigen::MatrixXd::Identity(6, 6);
  Eigen::MatrixXd R_k = Eigen::MatrixXd::Zero(rbt->u_dims, rbt->u_dims);
  Eigen::MatrixXd Q_T = Q_k;
  // Eigen::MatrixXd Q_T = 100 * Eigen::MatrixXd::Identity(rbt->x_dims, rbt->x_dims);
  ilqg = new OC::iLQG(rbt, Q_k, R_k, Q_T);
  ilqg->dt = 0.016;

  // Define initial state
  // Set random seed for reproducibility
  std::srand(static_cast<unsigned int>(std::time(nullptr)));
  Eigen::VectorXd q0 = Eigen::VectorXd::Zero(rbt->u_dims);
  q0 << 0.0, M_PI_2, 0.0, M_PI_2, 0.0, 0.0;
  Eigen::VectorXd M0 = GA_rbt.fkine(q0);

  Eigen::VectorXd qd = Eigen::VectorXd::Random(rbt->u_dims);
  qd << 1.0345, 0.591165, 1.3817, 0.501313, -1.0589, 1.4373;
  Eigen::VectorXd Md = GA_rbt.fkine(qd);
  
  Eigen::VectorXd x0(rbt->x_dims);
  Eigen::VectorXd xd(rbt->x_dims);
  x0 << GA::ga_log(GA::ga_prodM(M0, GA::ga_rev(Md))), q0;
  xd << Eigen::VectorXd::Zero(6), qd;

  // Make initialization for control sequence
  int T = 100;
  Eigen::MatrixXd u = ilqg->dt*Eigen::MatrixXd::Ones(rbt->u_dims, T);

  // Solve!
  std::cout << "Run iLQG!" << std::endl;
  auto start = std::chrono::system_clock::now();
  ilqg->generate_trajectory(x0, xd, u);
  auto now = std::chrono::system_clock::now();
  long int elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(now - start).count();
  std::cout << "iLQG took: " << elapsed / 1000.0 << " seconds." << std::endl;
  // Output the maximum absolute value of us
  std::cout << "Max absolute values of us:" << std::endl;
  for (int i = 0; i < 6; ++i) {
      std::cout << ilqg->us.row(i).cwiseAbs().maxCoeff() << std::endl;
  }
  std::cout << (ilqg->xs.col(T) - xd).head(6) << std::endl;

  return 0;
}