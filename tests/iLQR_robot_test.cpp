#include "OC/ilqr.hpp"
#include "OC/robot.hpp"

void CubicSpline(const Eigen::VectorXd& x, const Eigen::VectorXd& y, const Eigen::VectorXd& u, double dx0, double dxn, Eigen::VectorXd& v, Eigen::VectorXd& dv, Eigen::VectorXd& ddv) {
    size_t n = x.size();
    assert(y.size() == n);
    assert(n >= 2);

    Eigen::VectorXd dx(n - 1), dy(n - 1);
    for (size_t i = 0; i < n - 1; ++i) {
        dx[i] = x[i + 1] - x[i];
        dy[i] = y[i + 1] - y[i];
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

int main() {
    OC::iLQR* ilqr;
    // Make robot
    OC::Model* rbt = new OC::Robot("jaka");

    // 初始化权重矩阵
    Eigen::MatrixXd Q_k = 0.1*Eigen::MatrixXd::Identity(rbt->x_dims, rbt->x_dims);
    Q_k.block(6, 6, 6, 6) = Eigen::MatrixXd::Identity(6, 6);
    Eigen::MatrixXd R_k = Eigen::MatrixXd::Zero(rbt->u_dims, rbt->u_dims);
    Eigen::MatrixXd Q_T = 1000 * Eigen::MatrixXd::Identity(rbt->x_dims, rbt->x_dims);
    ilqr = new OC::iLQR(rbt, Q_k, R_k, Q_T);
    ilqr->dt = 0.016;

    // Define initial state
    // Set random seed for reproducibility
    std::srand(static_cast<unsigned int>(std::time(nullptr)));
    Eigen::VectorXd q0 = Eigen::VectorXd::Zero(rbt->u_dims);
    q0 << 0.5944, 2.3504, -0.1554, 0.2633, -0.8514, 2.0885;
    Eigen::VectorXd dq0 = Eigen::VectorXd::Zero(rbt->u_dims);
    Eigen::VectorXd x0(rbt->x_dims);
    x0 << q0, dq0;
    Eigen::VectorXd qd = Eigen::VectorXd::Random(rbt->u_dims);
    qd << 1.0345, 0.591165, 1.3817, 0.501313, -1.0589, 0.789955;
    Eigen::VectorXd dqd = Eigen::VectorXd::Zero(rbt->u_dims);
    Eigen::VectorXd xd(rbt->x_dims);
    xd << qd, dqd;

    // Make initialization for control sequence
    int T = 100;
    Eigen::VectorXd x(4);
    x << 0, 0.15*T*ilqr->dt, 0.85*T*ilqr->dt, T*ilqr->dt;
    Eigen::VectorXd y(4);
    Eigen::VectorXd u = Eigen::VectorXd::LinSpaced(T + 1, 0, T*ilqr->dt);
    Eigen::MatrixXd q = Eigen::MatrixXd::Zero(rbt->u_dims, T + 1);
    Eigen::MatrixXd dq = Eigen::MatrixXd::Zero(rbt->u_dims, T + 1);
    Eigen::MatrixXd ddq = Eigen::MatrixXd::Zero(rbt->u_dims, T + 1);
    Eigen::MatrixXd tau = Eigen::MatrixXd::Zero(rbt->u_dims, T);
    GA::GA_robot GA_rbt("jaka");
    Eigen::VectorXd q_row_i(T + 1);
    Eigen::VectorXd dq_row_i(T + 1);
    Eigen::VectorXd ddq_row_i(T + 1);

    for (int i = 0; i < rbt->u_dims; ++i) {
        y << q0(i), 0.95*q0(i) + 0.05*qd(i), 0.05*q0(i) + 0.95*qd(i), qd(i);
        CubicSpline(x, y, u, dq0(i), dqd(i), q_row_i, dq_row_i, ddq_row_i);
        q.row(i) = q_row_i;
        dq.row(i) = dq_row_i;
        ddq.row(i) = ddq_row_i;
    }
    for (int i = 0; i < T; ++i) {
        tau.col(i) = GA_rbt.idyn(q.col(i), dq.col(i), ddq.col(i), Eigen::VectorXd::Zero(rbt->u_dims));
    }

    // Solve!
    std::cout << "Run iLQR!" << std::endl;
    auto start = std::chrono::system_clock::now();
    ilqr->generate_trajectory(x0, xd, tau);
    auto now = std::chrono::system_clock::now();
    long int elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(now - start).count();
    std::cout << "iLQR took: " << elapsed / 1000.0 << " seconds." << std::endl;
    // Output the maximum absolute value of the last six rows of xs
    std::cout << "Max absolute values of the last six rows of xs:" << std::endl;
    for (int i = 6; i < 12; ++i) {
        std::cout << ilqr->xs.row(i).cwiseAbs().maxCoeff() << std::endl;
    }

    return 0;
}