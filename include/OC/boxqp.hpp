#ifndef ILQR_BOXQP_HPP
#define ILQR_BOXQP_HPP

#include <Eigen/Dense>
#include <vector>
#include <iostream>
#include <chrono>

namespace OC {

struct BoxQPSolution {
    Eigen::VectorXd x;
    int result; 
    Eigen::LLT<Eigen::MatrixXd> llt;
    Eigen::VectorXd free;
};

BoxQPSolution boxQP(
    const Eigen::MatrixXd& H,
    const Eigen::VectorXd& g,
    const Eigen::VectorXd& lower,
    const Eigen::VectorXd& upper,
    const Eigen::VectorXd& x0
);

} // namespace OC

#endif
