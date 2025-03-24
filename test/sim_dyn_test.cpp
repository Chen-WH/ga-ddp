#include "OC/ilqg.h"
#include "Utils/cubic_spline.hpp"
#include "rclcpp/rclcpp.hpp"
#include "sensor_msgs/msg/joint_state.hpp"
#include "trajectory_msgs/msg/joint_trajectory.hpp"
#include "trajectory_msgs/msg/joint_trajectory_point.hpp"

Eigen::VectorXd q0 = Eigen::VectorXd::Zero(6);
Eigen::VectorXd dq0 = q0;
auto traj_msg = trajectory_msgs::msg::JointTrajectory();

void current_states_callback(const sensor_msgs::msg::JointState::SharedPtr msg) {
  // 更新关节状态
  traj_msg.header.stamp = msg->header.stamp;
  q0 = Eigen::Map<const Eigen::VectorXd>(msg->position.data(), msg->position.size());
  dq0 = Eigen::Map<const Eigen::VectorXd>(msg->velocity.data(), msg->velocity.size());
}

int main(int argc, char *argv[]) {
  rclcpp::init(argc, argv);
  auto node = rclcpp::Node::make_shared("sim_dyn_test");
  auto trajectory_publisher = node->create_publisher<trajectory_msgs::msg::JointTrajectory>("/joint_trajectory", 10);
  auto target_publisher = node->create_publisher<sensor_msgs::msg::JointState>("/desired_states", 10);
  auto states_subscriber = node->create_subscription<sensor_msgs::msg::JointState>("/current_states", 10, current_states_callback);

  // declare robot model and iLQG
  GA::robot GA_rbt("jaka");
  OC::iLQG* ilqg;
  OC::Model* rbt = new OC::Robot_dyn("jaka");

  // Initialize weight matrices
  Eigen::MatrixXd Q_k = Eigen::MatrixXd::Zero(rbt->x_dims, rbt->x_dims);
  Q_k.block(0, 0, 6, 6) = 10 * Eigen::MatrixXd::Identity(6, 6);
  Eigen::MatrixXd R_k = Eigen::MatrixXd::Zero(rbt->u_dims, rbt->u_dims);
  Eigen::MatrixXd Q_T = Q_k;
  ilqg = new OC::iLQG(rbt, Q_k, R_k, Q_T);

  // Define initial state
  Eigen::VectorXd q1(6), q2(6), q3(6), q4(6), qd(6);
  q1 << 0.0, M_PI_2, 0.0, M_PI_2, 0.0, 0.0;
  q2 << 1.0345, 0.591165, 1.3817, 0.501313, -1.0589, 1.4373;
  q3 << -0.6388, 1.0015, 0.2633, 1.5955, -1.3021, 0.789955;
  q4 << 0.5944, 2.3504, -0.1554, 0.2633, -0.8514, -1.0756;
  q0 = q1;
  qd = q2;

  std::srand(static_cast<unsigned int>(std::time(nullptr)));
  Eigen::VectorXd x0(rbt->x_dims);
  Eigen::VectorXd xd(rbt->x_dims);
  Eigen::VectorXd dqd = Eigen::VectorXd::Zero(rbt->u_dims);
  Eigen::VectorXd Md = GA_rbt.fkine(qd);
  xd << Eigen::VectorXd::Zero(6), qd, dqd;
  Eigen::VectorXd fext = Eigen::VectorXd::Zero(rbt->u_dims);

  // Make initialization for control sequence
  int T = 100;
  Eigen::VectorXd x(4);
  x << 0, 0.15*T*ilqg->dt, 0.85*T*ilqg->dt, T*ilqg->dt;
  Eigen::VectorXd y(4);
  Eigen::VectorXd u = Eigen::VectorXd::LinSpaced(T + 1, 0, T*ilqg->dt);
  Eigen::MatrixXd q = Eigen::MatrixXd::Zero(rbt->u_dims, T + 1);
  Eigen::MatrixXd dq = Eigen::MatrixXd::Zero(rbt->u_dims, T + 1);
  Eigen::MatrixXd ddq = Eigen::MatrixXd::Zero(rbt->u_dims, T + 1);
  Eigen::MatrixXd tau = Eigen::MatrixXd::Zero(rbt->u_dims, T);
  Eigen::VectorXd q_row_i(T + 1);
  Eigen::VectorXd dq_row_i(T + 1);
  Eigen::VectorXd ddq_row_i(T + 1);

  rclcpp::WallRate loop_rate(5);
  trajectory_msgs::msg::JointTrajectoryPoint point;
  sensor_msgs::msg::JointState desired_states;

  Eigen::VectorXd position = Eigen::VectorXd::Zero(rbt->u_dims);
  Eigen::VectorXd velocity = Eigen::VectorXd::Zero(rbt->u_dims);
  Eigen::VectorXd effort = Eigen::VectorXd::Zero(rbt->u_dims);
  Eigen::VectorXd error = Eigen::VectorXd::Zero(6);

  while (rclcpp::ok()) {
    Eigen::VectorXd M0 = GA_rbt.fkine(q0);
    x0 << GA::ga_log(GA::ga_prodM(M0, GA::ga_rev(Md))), q0, dq0;
    error = x0.head(6);
    if (error.norm() < 1e-3) {
        std::cout << "Target posture reached!" << std::endl;
    } else {
        std::cout << "Motor error: " << error.transpose() << std::endl;
    }
    
    // Generate trajectory
    for (int i = 0; i < rbt->u_dims; ++i) {
      y << q0(i), 0.95*q0(i) + 0.05*qd(i), 0.05*q0(i) + 0.95*qd(i), qd(i);
      CubicSpline(x, y, u, dq0(i), dqd(i), q_row_i, dq_row_i, ddq_row_i);
      q.row(i) = q_row_i;
      dq.row(i) = dq_row_i;
      ddq.row(i) = ddq_row_i;
    }
    for (int i = 0; i < T; ++i) {
      tau.col(i) = GA_rbt.idyn(q.col(i), dq.col(i), ddq.col(i), fext);
    }

    // Solve iLQG
    ilqg->generate_trajectory(x0, xd, tau);

    // Publish desired states
    sensor_msgs::msg::JointState desired_states;
    desired_states.position = std::vector<double>(qd.data(), qd.data() + rbt->u_dims);
    desired_states.velocity = std::vector<double>(dqd.data(), dqd.data() + rbt->u_dims);
    effort = GA_rbt.idyn(qd, dqd, fext, fext);
    desired_states.effort = std::vector<double>(effort.data(), effort.data() + rbt->u_dims);
    target_publisher->publish(desired_states);

    // Publish trajectory
    for (int i = 0; i < T; ++i) {
      position = ilqg->xs.col(i + 1).segment(6, rbt->u_dims);
      velocity = ilqg->xs.col(i + 1).tail(rbt->u_dims);
      effort = ilqg->us.col(i);
      point.time_from_start = rclcpp::Duration::from_seconds(i * ilqg->dt);
      point.positions = std::vector<double>(position.data(), position.data() + rbt->u_dims);
      point.velocities = std::vector<double>(velocity.data(), velocity.data() + rbt->u_dims);
      point.effort = std::vector<double>(effort.data(), effort.data() + rbt->u_dims);
      traj_msg.points.push_back(point);
    }
    trajectory_publisher->publish(traj_msg);
    traj_msg.points.clear();

    rclcpp::spin_some(node);
    loop_rate.sleep();
  }

  rclcpp::shutdown();
  return 0;
}