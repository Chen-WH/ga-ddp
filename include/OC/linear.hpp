#ifndef OC_LINEAR_HPP
#define OC_LINEAR_HPP

#include "OC/model.hpp"
#include <cmath>

namespace OC {

class Linear : public Model {

public:

  Linear() {
    x_dims = 10;
    u_dims = 2;
    A = Eigen::MatrixXd::Random(x_dims, x_dims);
    A << 0.9978,    0.0051,    0.0057,   -0.0060,   -0.0158,    0.0313,    0.0504,   -0.0182,   -0.0004,  -0.0107,
   -0.0054,    0.9991,   -0.0139,   -0.0132,    0.0000,   -0.0150,    0.0047,   -0.0073,    0.0126,   -0.0308,
   -0.0044,    0.0133,    0.9991,    0.0099,    0.0095,   -0.0062,   -0.0262,   -0.0081,   -0.0063,   -0.0253,
    0.0080,    0.0123,   -0.0109,    0.9990,    0.0125,   -0.0254,   -0.0197,    0.0068,    0.0158,   -0.0105,
    0.0160,   -0.0003,   -0.0098,   -0.0128,    0.9995,    0.0142,   -0.0129,    0.0043,    0.0016,   -0.0118,
   -0.0311,    0.0156,    0.0064,    0.0253,   -0.0130,    0.9979,   -0.0041,   -0.0059,    0.0379,    0.0270,
   -0.0496,   -0.0038,    0.0265,    0.0200,    0.0143,    0.0007,    0.9974,    0.0154,    0.0238,    0.0232,
    0.0185,    0.0068,    0.0078,   -0.0077,   -0.0051,    0.0054,   -0.0150,    0.9989,    0.0352,   -0.0107,
    0.0020,   -0.0134,    0.0055,   -0.0167,   -0.0014,   -0.0376,   -0.0230,   -0.0354,    0.9982,    0.0010,
    0.0129,    0.0311,    0.0240,    0.0089,    0.0119,   -0.0273,   -0.0235,    0.0101,   -0.0018,    0.9983;
    B = Eigen::MatrixXd::Random(x_dims, u_dims);
    B << 0.0205,    1.2014,
   -0.8421,    0.4876,
    0.6305,   -0.2748,
    0.6958,    0.6668,
    0.3272,    0.5803,
   -0.8955,    1.3600,
   -0.7931,    0.4106,
    0.6148,   -0.9934,
    0.4061,    1.6993,
    0.0588,   -1.6746;

    u_min = Eigen::VectorXd(2); u_max = Eigen::VectorXd(2); 
    u_min << -0.6, -0.6;
    u_max << 0.6, 0.6;
  }

  virtual Eigen::VectorXd dynamics(const Eigen::VectorXd& x, const Eigen::VectorXd& u, const double dt) override {
    return A*x + dt*B*u;
  }

  virtual void dynamics_fo(const Eigen::VectorXd& x, const Eigen::VectorXd& u, const double dt, Eigen::MatrixXd& fx, Eigen::MatrixXd& fu) override {
    fx = A;
    fu = dt*B;
    return;
  }

private:
  Eigen::MatrixXd A, B;
};

} // namespace OC

#endif
