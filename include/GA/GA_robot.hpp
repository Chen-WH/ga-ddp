#ifndef GA_ROBOT_HPP
#define GA_ROBOT_HPP

#include <vector>
#include <string>
#include "GA/GA_algebra.hpp"

namespace GA {

Eigen::Matrix3d skew(const Eigen::Vector3d& v);

class GA_robot {
public:
    double g;                        // 重力加速度
    std::string model;               // 机器人型号

    // 构造函数
    GA_robot(const std::string& robot_model, double gravity = 9.81);

    // 析构函数
    ~GA_robot() = default;

    // 正运动学：计算末端执行器相对于基座的位姿
    Eigen::VectorXd fkine(const Eigen::VectorXd& q);

    // 逆动力学：计算机器人各关节的力矩
    Eigen::VectorXd idyn(const Eigen::VectorXd& q, const Eigen::VectorXd& dq, const Eigen::VectorXd& ddq, const Eigen::VectorXd& fext);
    
    // 正动力学：计算机器人各关节的加速度
    Eigen::VectorXd fdyn(const Eigen::VectorXd& q, const Eigen::VectorXd& dq, const Eigen::VectorXd& tau, const Eigen::VectorXd& fext);
    
    // 逆动力学一阶偏导：计算机器人各关节力矩对关节位置、速度、加速度的偏导
    void idyn_fo(const Eigen::VectorXd& q, const Eigen::VectorXd& dq, const Eigen::VectorXd& ddq, const Eigen::VectorXd& fext, Eigen::MatrixXd& ptau_pq, Eigen::MatrixXd& ptau_pdq, Eigen::MatrixXd& ptau_pddq);

    // 正动力学一阶偏导：计算机器人各关节加速度对关节位置、速度、力矩的偏导
    void fdyn_fo(const Eigen::VectorXd& q, const Eigen::VectorXd& dq, const Eigen::VectorXd& tau, const Eigen::VectorXd& fext, Eigen::MatrixXd& pddq_pq, Eigen::MatrixXd& pddq_pdq, Eigen::MatrixXd& pddq_ptau);

private:
    int n;                           // 机器人自由度
    Eigen::MatrixXd Mj;              // 8x(n+1) 初始连杆 motor
    Eigen::MatrixXd Lj;              // 6xn 关节轴线
    Eigen::VectorXi type;            // 1xn 关节类型
    std::vector<Eigen::MatrixXd> Ij; // 6x6xn 广义惯量矩阵
};

} // namespace GA

#endif // GA_ROBOT_HPP