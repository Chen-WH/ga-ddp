#include "GA/robotics.h"
#include <chrono>

int main() {
  GA::robot jaka("jaka");
  Eigen::VectorXd q(6), dq(6), ddq(6), fext(6);
  q << -1.2006, 1.4207, 1.7773, 1.2176, -3.0800, 2.1565;   // 设置关节角度
  dq << 4.2233, 2.7095, -4.5734, -1.2181, 2.0434, 2.2951;  // 设置关节速度
  ddq << -2.7572, -2.3095, 1.7303, -0.2251, 1.2372, -2.6356;  // 设置关节加速度
  fext << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;  // 设置外力

  Eigen::VectorXd M(8);
  M = jaka.fkine(q);
  std::cout << "Forward Kinematics:\n" << M << std::endl;

  Eigen::VectorXd tau = jaka.idyn(q, dq, ddq, fext);
  std::cout << "Inverse Dynamics:\n" << tau << std::endl;

  Eigen::VectorXd ddq_cal = jaka.fdyn(q, dq, tau, fext);
  std::cout << "Forward Dynamics:\n" << ddq_cal << std::endl;

  Eigen::MatrixXd ptau_pq, ptau_pdq, ptau_pddq;
  auto start_idyn_fo = std::chrono::high_resolution_clock::now();
  jaka.idyn_fo(q, dq, ddq, fext, ptau_pq, ptau_pdq, ptau_pddq);
  auto end_idyn_fo = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed_idyn_fo = end_idyn_fo - start_idyn_fo;
  std::cout << "Inverse Dynamics First Order Derivative ptau_pq:\n" << ptau_pq << std::endl;
  std::cout << "Inverse Dynamics First Order Derivative ptau_pdq:\n" << ptau_pdq << std::endl;
  std::cout << "Inverse Dynamics First Order Derivative ptau_pddq:\n" << ptau_pddq << std::endl;
  std::cout << "Time taken for idyn_fo: " << elapsed_idyn_fo.count() << " seconds" << std::endl;

  Eigen::MatrixXd pddq_pq, pddq_pdq, pddq_ptau;
  auto start_fdyn_fo = std::chrono::high_resolution_clock::now();
  jaka.fdyn_fo(q, dq, tau, fext, pddq_pq, pddq_pdq, pddq_ptau);
  auto end_fdyn_fo = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed_fdyn_fo = end_fdyn_fo - start_fdyn_fo;
  std::cout << "Forward Dynamics First Order Derivative pddq_pq:\n" << pddq_pq << std::endl;
  std::cout << "Forward Dynamics First Order Derivative pddq_pdq:\n" << pddq_pdq << std::endl;
  std::cout << "Forward Dynamics First Order Derivative pddq_ptau:\n" << pddq_ptau << std::endl;
  std::cout << "Time taken for fdyn_fo: " << elapsed_fdyn_fo.count() << " seconds" << std::endl;

  return 0;
}