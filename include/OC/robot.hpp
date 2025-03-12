#ifndef OC_ROBOT_HPP
#define OC_ROBOT_HPP

#include "OC/model.hpp"
#include "GA/GA_robot.hpp"
#include <cmath>

namespace OC {

class Robot : public Model {

public:
  GA::GA_robot GA_rbt;

  Robot(const std::string& robot_model) : GA_rbt(robot_model) {
    x_dims = 12;
    u_dims = 6;
    u_min = Eigen::VectorXd(6); u_max = Eigen::VectorXd(6);
    if (robot_model == "chin") {
      u_min << -150, -75, -75, -37.5, -37.5, -37.5;
      u_max << 150, 75, 75, 37.5, 37.5, 37.5;
    }
    else if (robot_model == "jaka") {
      u_min << -30, -120, -120, -30, -45, -15;
      u_max <<  30,  120,  120,  30,  45,  15;
    }
  }

  virtual Eigen::VectorXd dynamics(const Eigen::VectorXd& x, const Eigen::VectorXd& u, const double dt) override {
    Eigen::VectorXd f(x_dims);
    Eigen::VectorXd q = x.head(6), dq = x.tail(6);
    f << q + dt*dq, dq + dt*GA_rbt.fdyn(q, dq, u, Eigen::VectorXd::Zero(6));
    return f;
  }

  virtual void dynamics_fo(const Eigen::VectorXd& x, const Eigen::VectorXd& u, const double dt, Eigen::MatrixXd& fx, Eigen::MatrixXd& fu) override {
    fx.resize(x_dims, x_dims);
    fu.resize(x_dims, u_dims);
    Eigen::VectorXd q = x.head(6), dq = x.tail(6);
    Eigen::MatrixXd pddq_pq, pddq_pdq, pddq_ptau;
    GA_rbt.fdyn_fo(q, dq, u, Eigen::VectorXd::Zero(6), pddq_pq, pddq_pdq, pddq_ptau);
    fx << Eigen::MatrixXd::Identity(6, 6), dt*Eigen::MatrixXd::Identity(6, 6), dt*pddq_pq, Eigen::MatrixXd::Identity(6, 6) + dt*pddq_pdq;
    fu << Eigen::MatrixXd::Zero(6, 6), dt*pddq_ptau;
    return;
  }
};

} // namespace OC

#endif
