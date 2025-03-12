#ifndef GA_ALGEBRA_HPP
#define GA_ALGEBRA_HPP

#include <Eigen/Dense>
#include <cmath>
#include <iostream>

namespace GA {

/* Motor群公共元素 */

// GA_exp: 计算李代数 bivector B 对应的李群元素 M
Eigen::VectorXd GA_exp(const Eigen::VectorXd& B);

// GA_log: 计算李群元素 M 对应的李代数 bivector B
Eigen::VectorXd GA_log(const Eigen::VectorXd& M);

// GA_AdM: 计算伴随表示 Ad
Eigen::MatrixXd GA_AdM(const Eigen::VectorXd& M);

// GA_adB: 计算伴随表示 ad
Eigen::MatrixXd GA_adB(const Eigen::VectorXd& B);

// GA_DHmotor: 计算DH参数对应的Motor
Eigen::MatrixXd GA_DHmotor(const Eigen::Vector4d& param, int type);

// GA_dual: 计算对偶
Eigen::VectorXd GA_dual(const Eigen::VectorXd& x);

// GA_forque: 计算受力对应的forque
Eigen::VectorXd GA_forque(const Eigen::VectorXd& q, const Eigen::VectorXd& f);

// GA_metric: 计算李代数上的度规
Eigen::RowVectorXd GA_metric(const Eigen::VectorXd& x);

// GA_motor: 计算关节驱动量对应的Motor
Eigen::VectorXd GA_motor(const Eigen::VectorXd& L, double q, int type);

// GA_prodM: 计算两个Motor的乘积
Eigen::VectorXd GA_prodM(const Eigen::VectorXd& M1, const Eigen::VectorXd& M2);

// GA_Rev: 计算 motor 群的逆元
Eigen::VectorXd GA_Rev(const Eigen::VectorXd& M);

/* PGA元素 */

// PGA_Lie2: 计算李括号 dL = [L, B] = -(LB - BL)/2
Eigen::VectorXd PGA_Lie2(const Eigen::VectorXd& L, const Eigen::VectorXd& B);

// PGA_Frame2：通过刚体运动变换 Lp = M * L * iM
Eigen::VectorXd PGA_Frame2(const Eigen::VectorXd& L, const Eigen::VectorXd& M);

} // namespace GA

#endif // GA_ALGEBRA_HPP