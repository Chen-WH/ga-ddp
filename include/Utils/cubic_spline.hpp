#ifndef Utils_CubicSpline_HPP
#define Utils_CubicSpline_HPP

#include <Eigen/Dense>
#include <cmath>

void CubicSpline(const Eigen::VectorXd& x, const Eigen::VectorXd& y, const Eigen::VectorXd& u, double dx0, double dxn, Eigen::VectorXd& v, Eigen::VectorXd& dv, Eigen::VectorXd& ddv) {
    int n = x.size();
    assert(y.size() == n);
    assert(n >= 2);

    Eigen::VectorXd dx(n - 1), dy(n - 1);
    for (int i = 0; i < n - 1; ++i) {
        dx[i] = x[i + 1] - x[i];
        dy[i] = y[i + 1] - y[i];
    }

    // Construct coefficient matrix A and right-hand side vector t
    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(n, n);
    Eigen::VectorXd t = Eigen::VectorXd::Zero(n);

    for (int i = 1; i < n - 1; ++i) {
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
    for (int i = 0; i < n - 1; ++i) {
        b[i] = dy[i] / dx[i] - dx[i] * (2 * c(i) + c(i + 1)) / 3.0;
        d[i] = (c(i + 1) - c(i)) / (3 * dx[i]);
    }
    int idx = 0;
    double e;
    for (int i = 0; i < u.size(); ++i) {
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

#endif