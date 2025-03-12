#ifndef OC_MODEL_HPP
#define OC_MODEL_HPP

#include <Eigen/Dense>
// #include <vector>

/*
  Borrowed from: https://github.com/kazuotani14/iLQR
*/

namespace OC {

class Model {
public:
  int x_dims, u_dims;

  Eigen::VectorXd u_min, u_max;

  virtual Eigen::VectorXd dynamics(const Eigen::VectorXd& x, const Eigen::VectorXd& u, const double dt) = 0;
  // vector<MatrixXd> for OpenMP
  virtual void dynamics_fo(const Eigen::VectorXd& x, const Eigen::VectorXd& u, const double dt, Eigen::MatrixXd& fx, Eigen::MatrixXd& fu) = 0;
  // full DDP: to be implemented
  // virtual void dynamics_so(const Eigen::VectorXd& x, const Eigen::VectorXd& u, const double dt, Eigen::MatrixXd& fx, Eigen::MatrixXd& fu) = 0;
};

} // namespace OC

#endif