#ifndef __GA_algebra_h__
#define __GA_algebra_h__

#include <Eigen/Dense>
#include <cmath>
#include <iostream>

namespace GA {

/* Motor群公共元素 */

// ga_exp: 计算李代数 bivector B 对应的李群元素 M
Eigen::VectorXd ga_exp(const Eigen::VectorXd& B);

// ga_log: 计算李群元素 M 对应的李代数 bivector B
Eigen::VectorXd ga_log(const Eigen::VectorXd& M);

// ga_AdM: 计算伴随表示 Ad
Eigen::MatrixXd ga_AdM(const Eigen::VectorXd& M);

// ga_adB: 计算伴随表示 ad
Eigen::MatrixXd ga_adB(const Eigen::VectorXd& B);

// ga_DHmotor: 计算DH参数对应的Motor
Eigen::MatrixXd ga_DHmotor(const Eigen::Vector4d& param, int type);

// ga_dual: 计算对偶
Eigen::VectorXd ga_dual(const Eigen::VectorXd& x);

// ga_metric: 计算李代数上的度规
Eigen::RowVectorXd ga_metric(const Eigen::VectorXd& x);

// ga_motor: 计算关节驱动量对应的Motor
Eigen::VectorXd ga_motor(const Eigen::VectorXd& L, double q, int type);

// ga_prodM: 计算两个Motor的乘积
Eigen::VectorXd ga_prodM(const Eigen::VectorXd& M1, const Eigen::VectorXd& M2);

// ga_rev: 计算 motor 群的逆元
Eigen::VectorXd ga_rev(const Eigen::VectorXd& M);

/* PGA元素 */

// pga_Lie2: 计算李括号 dL = [L, B] = -(LB - BL)/2
Eigen::VectorXd pga_Lie2(const Eigen::VectorXd& L, const Eigen::VectorXd& B);

// pga_frame2：通过刚体运动变换 Lp = M * L * iM
Eigen::VectorXd pga_frame2(const Eigen::VectorXd& L, const Eigen::VectorXd& M);

} // namespace GA

#endif