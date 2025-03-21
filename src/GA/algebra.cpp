#include "GA/algebra.h"

namespace GA {

Eigen::VectorXd ga_exp(const Eigen::VectorXd& B) {
  Eigen::VectorXd M(8);
  Eigen::VectorXd B_trans(3);
  Eigen::VectorXd B_rot(3);
  B_trans = -0.5*B.head(3);
  B_rot = -0.5*B.tail(3);
  double l = B_rot.squaredNorm();
  if (l < 1e-9) {
    // 纯平移
    M << 1, B_trans, 0, 0, 0, 0;
  } else {
    // 一般情况
    double m = B_trans.dot(B_rot);
    double a = sqrt(l);
    double c = cos(a);
    double s = sin(a) / a;
    double t = m / l * (c - s);
    M << c, s*B_trans + t*B_rot, s*B_rot, m*s;
  }
  return M;
}

Eigen::VectorXd ga_log(const Eigen::VectorXd& M) {
  Eigen::VectorXd B(6);
  double M0 = M(0);
  if (std::abs(M0 - 1) < 1e-9) {
    B << M.segment(1, 3), Eigen::Vector3d::Zero();
  } else {
    double a = 1 / (1 - M0*M0);    // inv squared length
    double b = acos(M0)*sqrt(a);   // rotation scale
    double c = a*M(7)*(1 - M0*b);  // translation scale
    B << c*M.segment(4, 3) + b*M.segment(1, 3), b*M.segment(4, 3);
  }
  B *= -2;
  return B;
}

Eigen::MatrixXd ga_AdM(const Eigen::VectorXd& M) {
  Eigen::MatrixXd Ad_M(6, 6);
  double M0 = M(0), M1 = M(1), M2 = M(2), M3 = M(3), M4 = M(4), M5 = M(5), M6 = M(6), M7 = M(7);
  double M0_sq = M0 * M0, M4_sq = M4 * M4, M5_sq = M5 * M5, M6_sq = M6 * M6;

  Ad_M(0, 0) = M0_sq + M4_sq - M5_sq - M6_sq;
  Ad_M(0, 1) = 2 * M0 * M6 + 2 * M4 * M5;
  Ad_M(0, 2) = 2 * M4 * M6 - 2 * M0 * M5;
  Ad_M(0, 3) = 2 * M1 * M4 - 2 * M0 * M7 - 2 * M2 * M5 - 2 * M3 * M6;
  Ad_M(0, 4) = 2 * M0 * M3 + 2 * M1 * M5 + 2 * M2 * M4 - 2 * M6 * M7;
  Ad_M(0, 5) = 2 * M1 * M6 - 2 * M0 * M2 + 2 * M3 * M4 + 2 * M5 * M7;

  Ad_M(1, 0) = 2 * M4 * M5 - 2 * M0 * M6;
  Ad_M(1, 1) = M0_sq - M4_sq + M5_sq - M6_sq;
  Ad_M(1, 2) = 2 * M0 * M4 + 2 * M5 * M6;
  Ad_M(1, 3) = 2 * M1 * M5 - 2 * M0 * M3 + 2 * M2 * M4 + 2 * M6 * M7;
  Ad_M(1, 4) = 2 * M2 * M5 - 2 * M0 * M7 - 2 * M1 * M4 - 2 * M3 * M6;
  Ad_M(1, 5) = 2 * M0 * M1 + 2 * M2 * M6 + 2 * M3 * M5 - 2 * M4 * M7;

  Ad_M(2, 0) = 2 * M0 * M5 + 2 * M4 * M6;
  Ad_M(2, 1) = 2 * M5 * M6 - 2 * M0 * M4;
  Ad_M(2, 2) = M0_sq - M4_sq - M5_sq + M6_sq;
  Ad_M(2, 3) = 2 * M0 * M2 + 2 * M1 * M6 + 2 * M3 * M4 - 2 * M5 * M7;
  Ad_M(2, 4) = 2 * M2 * M6 - 2 * M0 * M1 + 2 * M3 * M5 + 2 * M4 * M7;
  Ad_M(2, 5) = 2 * M3 * M6 - 2 * M0 * M7 - 2 * M2 * M5 - 2 * M1 * M4;

  Ad_M(3, 0) = 0;
  Ad_M(3, 1) = 0;
  Ad_M(3, 2) = 0;
  Ad_M(3, 3) = M0_sq + M4_sq - M5_sq - M6_sq;
  Ad_M(3, 4) = 2 * M0 * M6 + 2 * M4 * M5;
  Ad_M(3, 5) = 2 * M4 * M6 - 2 * M0 * M5;

  Ad_M(4, 0) = 0;
  Ad_M(4, 1) = 0;
  Ad_M(4, 2) = 0;
  Ad_M(4, 3) = 2 * M4 * M5 - 2 * M0 * M6;
  Ad_M(4, 4) = M0_sq - M4_sq + M5_sq - M6_sq;
  Ad_M(4, 5) = 2 * M0 * M4 + 2 * M5 * M6;

  Ad_M(5, 0) = 0;
  Ad_M(5, 1) = 0;
  Ad_M(5, 2) = 0;
  Ad_M(5, 3) = 2 * M0 * M5 + 2 * M4 * M6;
  Ad_M(5, 4) = 2 * M5 * M6 - 2 * M0 * M4;
  Ad_M(5, 5) = M0_sq - M4_sq - M5_sq + M6_sq;

  return Ad_M;
}

Eigen::MatrixXd ga_adB(const Eigen::VectorXd& B) {
  Eigen::MatrixXd ad_B = Eigen::MatrixXd::Zero(6, 6);
  ad_B(0, 1) = -B(5); ad_B(0, 2) =  B(4); ad_B(0, 4) = -B(2); ad_B(0, 5) =  B(1);
  ad_B(1, 0) =  B(5); ad_B(1, 2) = -B(3); ad_B(1, 3) =  B(2); ad_B(1, 5) = -B(0);
  ad_B(2, 0) = -B(4); ad_B(2, 1) =  B(3); ad_B(2, 3) = -B(1); ad_B(2, 4) =  B(0);
  ad_B(3, 4) = -B(5); ad_B(3, 5) =  B(4);
  ad_B(4, 3) =  B(5); ad_B(4, 5) = -B(3);
  ad_B(5, 3) = -B(4); ad_B(5, 4) =  B(3);
  return ad_B;
}

Eigen::MatrixXd ga_DHmotor(const Eigen::Vector4d& param, int type) {
  Eigen::MatrixXd M(8, 1);
  double d, theta, a, alpha;
  switch (type) {
    case 0: // DH 存储类型, param = [d, theta, a, alpha], M = exp(d*e03 + theta*e12)*exp(a*e01 + alpha*e23)
      d = 0.5*param(0);
      theta = 0.5*param(1);
      a = 0.5*param(2);
      alpha = 0.5*param(3);
      M(0) = cos(alpha) * cos(theta);
      M(1) =  d * sin(alpha) * sin(theta) - a * cos(alpha) * cos(theta);
      M(2) = -a * cos(alpha) * sin(theta) - d * sin(alpha) * cos(theta);
      M(3) =  a * sin(alpha) * sin(theta) - d * cos(alpha) * cos(theta);
      M(4) = -sin(alpha) * cos(theta);
      M(5) = -sin(alpha) * sin(theta);
      M(6) = -cos(alpha) * sin(theta);
      M(7) = a * sin(alpha) * cos(theta) + d * cos(alpha) * sin(theta);
      break;

    case 1: // MDH 存储类型, param = [a, alpha, d, theta], M = exp(a*e01 + alpha*e23)*exp(d*e03 + theta*e12)
      a = 0.5*param(0);
      alpha = 0.5*param(1);
      d = 0.5*param(2);
      theta = 0.5*param(3);
      M(0) = cos(alpha) * cos(theta);
      M(1) = d * sin(alpha) * sin(theta) - a * cos(alpha) * cos(theta);
      M(2) = a * cos(alpha) * sin(theta) + d * sin(alpha) * cos(theta);
      M(3) = a * sin(alpha) * sin(theta) - d * cos(alpha) * cos(theta);
      M(4) = -sin(alpha) * cos(theta);
      M(5) = sin(alpha) * sin(theta);
      M(6) = -cos(alpha) * sin(theta);
      M(7) = a * sin(alpha) * cos(theta) + d * cos(alpha) * sin(theta);
      break;

    default:
      std::cerr << "Unknown DH type!" << std::endl;
      break;
  }
  return M;
}

Eigen::VectorXd ga_dual(const Eigen::VectorXd& x) {
  int n = x.size();
  Eigen::VectorXd res;
  switch (n) {
    case 6:
      res.resize(6);
      res << -x(3), -x(4), -x(5), -x(0), -x(1), -x(2);
      break;
    
    case 8:
      res.resize(8);
      res << x(7), -x(4), -x(5), -x(6), -x(1), -x(2), -x(3), x(0);
      break;

    default:
      std::cerr << "Error: Invalid input size, expected 6 or 8!" << std::endl;
      break;
  }
  return res;
}

Eigen::RowVectorXd ga_metric(const Eigen::VectorXd& x) {
  Eigen::RowVectorXd res(6);
  res << x(3), x(4), x(5), x(0), x(1), x(2);
  return res;
}

Eigen::VectorXd ga_motor(const Eigen::VectorXd& L, double q, int type) {
  Eigen::VectorXd M(8);
  q = -0.5*q;
  switch (type) {
    case 1: {  // 旋转关节
      M << cos(q), sin(q) * L, 0;
      break;
    }
    case 2: {  // 平移关节
      M << 1, q * L.head(3), 0, 0, 0, 0;
      break;
    }
    default:
      std::cerr << "Error: Invalid type!" << std::endl;
      exit(1);
  }
  return M;
}

Eigen::VectorXd ga_prodM(const Eigen::VectorXd& M1, const Eigen::VectorXd& M2) {
  Eigen::VectorXd M(8);
  M(0) = M1(0)*M2(0) - M1(4)*M2(4) - M1(5)*M2(5) - M1(6)*M2(6);
  M(1) = M1(0)*M2(1) + M1(1)*M2(0) - M1(2)*M2(6) + M1(3)*M2(5) - M1(5)*M2(3) + M1(6)*M2(2) - M1(4)*M2(7) - M1(7)*M2(4);
  M(2) = M1(0)*M2(2) + M1(2)*M2(0) + M1(1)*M2(6) - M1(3)*M2(4) + M1(4)*M2(3) - M1(6)*M2(1) - M1(5)*M2(7) - M1(7)*M2(5);
  M(3) = M1(0)*M2(3) + M1(3)*M2(0) - M1(1)*M2(5) + M1(2)*M2(4) - M1(4)*M2(2) + M1(5)*M2(1) - M1(6)*M2(7) - M1(7)*M2(6);
  M(4) = M1(0)*M2(4) + M1(4)*M2(0) - M1(5)*M2(6) + M1(6)*M2(5);
  M(5) = M1(0)*M2(5) + M1(5)*M2(0) + M1(4)*M2(6) - M1(6)*M2(4);
  M(6) = M1(0)*M2(6) + M1(6)*M2(0) - M1(4)*M2(5) + M1(5)*M2(4);
  M(7) = M1(1)*M2(4) + M1(4)*M2(1) + M1(0)*M2(7) + M1(2)*M2(5) + M1(5)*M2(2) + M1(7)*M2(0) + M1(3)*M2(6) + M1(6)*M2(3);
  return M;
}

Eigen::VectorXd ga_rev(const Eigen::VectorXd& M) {
  Eigen::VectorXd iM(8);
  iM(0) = M(0);
  iM(7) = M(7);
  iM.segment(1, 6) = -M.segment(1, 6);
  return iM;
}

Eigen::VectorXd pga_Lie2(const Eigen::VectorXd& L, const Eigen::VectorXd& B) {
  Eigen::VectorXd dL(6);
  dL(0) = B(2)*L(4) - B(1)*L(5) - B(4)*L(2) + B(5)*L(1);
  dL(1) = B(0)*L(5) - B(2)*L(3) + B(3)*L(2) - B(5)*L(0);
  dL(2) = B(1)*L(3) - B(0)*L(4) - B(3)*L(1) + B(4)*L(0);
  dL(3) = B(5)*L(4) - B(4)*L(5);
  dL(4) = B(3)*L(5) - B(5)*L(3);
  dL(5) = B(4)*L(3) - B(3)*L(4);
  return dL;
}

Eigen::VectorXd pga_frame2(const Eigen::VectorXd& L, const Eigen::VectorXd& M) {
  Eigen::VectorXd Lp(6);

  Lp(0) = L(0)*(M(0)*M(0) + M(4)*M(4) - M(5)*M(5) - M(6)*M(6)) 
  + 2*( L(1)*M(0)*M(6) - L(2)*M(0)*M(5) + L(4)*M(0)*M(3) 
  - L(5)*M(0)*M(2) + L(3)*M(1)*M(4) + L(1)*M(4)*M(5) 
  - L(3)*M(0)*M(7) - L(3)*M(2)*M(5) + L(4)*M(1)*M(5) 
  + L(4)*M(2)*M(4) + L(2)*M(4)*M(6) - L(3)*M(3)*M(6) 
  + L(5)*M(1)*M(6) + L(5)*M(3)*M(4) - L(4)*M(6)*M(7) 
  + L(5)*M(5)*M(7));

  Lp(1) = L(1)*(M(0)*M(0) - M(4)*M(4) + M(5)*M(5) - M(6)*M(6)) 
  + 2*(-L(0)*M(0)*M(6) + L(2)*M(0)*M(4) - L(3)*M(0)*M(3) 
  + L(5)*M(0)*M(1) + L(0)*M(4)*M(5) + L(3)*M(1)*M(5) 
  + L(3)*M(2)*M(4) - L(4)*M(1)*M(4) - L(4)*M(0)*M(7) 
  + L(4)*M(2)*M(5) + L(2)*M(5)*M(6) - L(4)*M(3)*M(6) 
  + L(5)*M(2)*M(6) + L(5)*M(3)*M(5) + L(3)*M(6)*M(7) 
  - L(5)*M(4)*M(7));

  Lp(2) = L(2)*(M(0)*M(0) - M(4)*M(4) - M(5)*M(5) + M(6)*M(6)) 
  + 2*( L(0)*M(0)*M(5) - L(1)*M(0)*M(4) + L(3)*M(0)*M(2) 
  - L(4)*M(0)*M(1) + L(0)*M(4)*M(6) + L(3)*M(1)*M(6) 
  + L(3)*M(3)*M(4) - L(5)*M(1)*M(4) + L(1)*M(5)*M(6) 
  + L(4)*M(2)*M(6) + L(4)*M(3)*M(5) - L(5)*M(0)*M(7) 
  - L(5)*M(2)*M(5) + L(5)*M(3)*M(6) - L(3)*M(5)*M(7) 
  + L(4)*M(4)*M(7));

  Lp(3) = L(3)*(M(0)*M(0) + M(4)*M(4) - M(5)*M(5) - M(6)*M(6)) 
  + 2*(-L(5)*M(0)*M(5) + L(4)*M(0)*M(6) + L(4)*M(4)*M(5) 
  + L(5)*M(4)*M(6));

  Lp(4) = L(4)*(M(0)*M(0) - M(4)*M(4) + M(5)*M(5) - M(6)*M(6)) 
  + 2*( L(5)*M(0)*M(4) - L(3)*M(0)*M(6) + L(3)*M(4)*M(5) 
  + L(5)*M(5)*M(6));

  Lp(5) = L(5)*(M(0)*M(0) - M(4)*M(4) - M(5)*M(5) + M(6)*M(6)) 
  + 2*(-L(4)*M(0)*M(4) + L(3)*M(0)*M(5) + L(3)*M(4)*M(6) 
  + L(4)*M(5)*M(6));

  return Lp;
}

} // namespace GA