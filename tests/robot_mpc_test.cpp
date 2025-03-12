#include "OC/ilqr.hpp"
#include "OC/robot.hpp"
#include "rclcpp/rclcpp.hpp"
#include "sensor_msgs/msg/joint_state.hpp"
#include "trajectory_msgs/msg/joint_trajectory.hpp"
#include "trajectory_msgs/msg/joint_trajectory_point.hpp"

void CubicSpline(const Eigen::VectorXd& x, const Eigen::VectorXd& y, const Eigen::VectorXd& u, double dx0, double dxn, Eigen::VectorXd& v, Eigen::VectorXd& dv, Eigen::VectorXd& ddv) {
    size_t n = x.size();
    assert(y.size() == n);
    assert(n >= 2);

    Eigen::VectorXd dx(n - 1), dy(n - 1);
    for (size_t i = 0; i < n - 1; ++i) {
        dx[i] = x[i + 1] - x[i];
        dy[i] = y[i + 1] - x[i];
    }

    // Construct coefficient matrix A and right-hand side vector t
    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(n, n);
    Eigen::VectorXd t = Eigen::VectorXd::Zero(n);

    for (size_t i = 1; i < n - 1; ++i) {
        A(i, i - 1) = dx[i - 1];
        A(i, i) = 2 * (dx[i - 1] + dx[i]);
        A(i, i + 1) = dx[i];
        t(i) = 3 * (dy[i] / dx[i] - dy[i - 1] / dx[i - 1]);
    }

    A(0, 0) = 2 * dx[0];
    A(0, 1) = dx[0];
    t(0) = 3 * (dy[0] / dx[0] - dx0);

    A(n - 1, n - 2) = dx[n - 2];
    A(n - 1, n - 1) = 2 * dx[n - 2];
    t(n - 1) = 3 * (dxn - dy[n - 2] / dx[n - 2]);

    // Solve for c
    Eigen::VectorXd c = A.fullPivLu().solve(t);

    // Calculate b and d coefficients
    Eigen::VectorXd b(n - 1), d(n - 1);
    for (size_t i = 0; i < n - 1; ++i) {
        b[i] = dy[i] / dx[i] - dx[i] * (2 * c(i) + c(i + 1)) / 3.0;
        d[i] = (c(i + 1) - c(i)) / (3 * dx[i]);
    }
    size_t idx = 0;
    double e;
    for (size_t i = 0; i < u.size(); ++i) {
        if (u[i] > x[idx + 1]) {
            idx += 1;
        }
        e = u(i) - x(idx);
        v[i] = y[idx] + b[idx] * e + c(idx) * e * e + d[idx] * e * e * e;
        dv[i] = b[idx] + 2 * c(idx) * e + 3 * d[idx] * e * e;
        ddv[i] = 2 * c(idx) + 6 * d[idx] * e;
    }
    return;
}

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
    auto node = rclcpp::Node::make_shared("robot_ilqr_test");
    auto trajectory_publisher = node->create_publisher<trajectory_msgs::msg::JointTrajectory>("/joint_trajectory", 10);
    auto target_publisher = node->create_publisher<sensor_msgs::msg::JointState>("/desired_states", 10);
    auto states_subscriber = node->create_subscription<sensor_msgs::msg::JointState>("/current_states", 10, current_states_callback);
    
    // 初始化chin机器人
    GA::GA_robot GA_rbt("chin");
    OC::iLQR* ilqr;
    OC::Model* rbt = new OC::Robot("chin");

    // 初始化iLQR
    Eigen::MatrixXd Q_k = 0.1*Eigen::MatrixXd::Identity(rbt->x_dims, rbt->x_dims);
    Q_k.block(0, 0, 6, 6) = 10 * Eigen::MatrixXd::Identity(6, 6);
    Eigen::MatrixXd R_k = Eigen::MatrixXd::Zero(rbt->u_dims, rbt->u_dims);
    Eigen::MatrixXd Q_T = Q_k;
    ilqr = new OC::iLQR(rbt, Q_k, R_k, Q_T);

    // Define initial state
    Eigen::VectorXd x0(rbt->x_dims);
    std::srand(static_cast<unsigned int>(std::time(nullptr)));
    Eigen::VectorXd xd(rbt->x_dims);
    Eigen::VectorXd qd = Eigen::VectorXd::Zero(rbt->u_dims);
    //qd(1) += M_PI_2; qd(3) += M_PI_2;
    Eigen::VectorXd dqd = Eigen::VectorXd::Zero(rbt->u_dims);
    xd << qd, dqd;
    Eigen::VectorXd fext = Eigen::VectorXd::Zero(rbt->u_dims);

    // Make initialization for control sequence
    int T = 50;
    Eigen::VectorXd x(4);
    x << 0, 0.15*T*ilqr->dt, 0.85*T*ilqr->dt, T*ilqr->dt;
    Eigen::VectorXd y(4);
    Eigen::VectorXd u = Eigen::VectorXd::LinSpaced(T + 1, 0, T*ilqr->dt);
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

    double error = 0.0;

    while (rclcpp::ok()) {
        x0 << q0, dq0;
        error = (x0 - xd).norm();
        if (error < 1e-2) {
            std::cout << "Target posture reached, Generate next posture? (Press Enter to continue)";
            std::cin.ignore();
            qd = 0.7*Eigen::VectorXd::Random(rbt->u_dims);
            //qd(1) += M_PI_2; qd(3) += M_PI_2;
            xd << qd, dqd;
        } else {
            std::cout << "Angle error: " << error << std::endl;
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

        // Solve iLQR
        ilqr->generate_trajectory(x0, xd, tau);

        // Publish desired states
        sensor_msgs::msg::JointState desired_states;
        desired_states.position = std::vector<double>(qd.data(), qd.data() + 6);
        desired_states.velocity = std::vector<double>(dqd.data(), dqd.data() + 6);
        effort = GA_rbt.idyn(qd, dqd, fext, fext);
        desired_states.effort = std::vector<double>(effort.data(), effort.data() + 6);
        target_publisher->publish(desired_states);

        // Publish trajectory
        for (int i = 0; i < T; ++i) {
            position = ilqr->xs.col(i + 1).head(6);
            velocity = ilqr->xs.col(i + 1).tail(6);
            effort = ilqr->us.col(i);
            point.time_from_start = rclcpp::Duration::from_seconds(i * ilqr->dt);
            point.positions = std::vector<double>(position.data(), position.data() + 6);
            point.velocities = std::vector<double>(velocity.data(), velocity.data() + 6);
            point.effort = std::vector<double>(effort.data(), effort.data() + 6);
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