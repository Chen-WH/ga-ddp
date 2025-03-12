#include "OC/boxqp.hpp"

int main() {
    // 构建示例数据
    /*
    int n = 500;
    Eigen::MatrixXd H = Eigen::MatrixXd::Random(n, n);
    H = H.transpose() * H; // 使H为正定矩阵
    Eigen::VectorXd g = Eigen::VectorXd::Random(n);
    Eigen::VectorXd lower = -Eigen::VectorXd::Ones(n);
    Eigen::VectorXd upper = Eigen::VectorXd::Ones(n);
    Eigen::VectorXd x0 = Eigen::VectorXd::Random(n);
    */
   
    Eigen::MatrixXd H(6, 6);
    H << 1.0000,   -0.0000,   -0.0000,   -0.0001,    0.0003,    0.0001,
   -0.0000,    1.0000,    0.0000,   -0.0000,    0.0000,   -0.0000,
   -0.0000,    0.0000,    1.0001,   -0.0001,    0.0001,   -0.0000,
   -0.0001,   -0.0000,   -0.0001,    1.0004,   -0.0008,   -0.0000,
    0.0003,    0.0000,    0.0001,   -0.0008,    1.0029,    0.0009,
    0.0001,   -0.0000,   -0.0000,   -0.0000,    0.0009,    1.0034;
    Eigen::VectorXd g(6);
    g << -0.1771,
    0.1681,
    0.2991,
   -0.3376,
   -0.8439,
   -1.1422;
    Eigen::VectorXd lower(6);
    lower << -100.6435,
  -46.5702,
  -51.1869,
  -24.8770,
  -25.0225,
  -25.0298;
    Eigen::VectorXd upper(6);
    upper << 99.3565,
   53.4298,
   48.8131,
   25.1230,
   24.9775,
   24.9702;
    Eigen::VectorXd x0(6);
    x0 << 0.1649,
   -0.1611,
   -0.2804,
    0.3844,
    0.7524,
    1.1800;
  
    // 开始计时
    auto start = std::chrono::high_resolution_clock::now();

    // 调用BoxQP
    OC::BoxQPSolution sol = OC::boxQP(H, g, lower, upper, x0);

    // 结束计时
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;

    // 输出结果
    std::cout << "Result: " << sol.result << std::endl;
    std::cout << "x: " << sol.x.transpose() << std::endl;
    std::cout << "BoxQP 用时: " << elapsed.count() << " 秒" << std::endl; // 输出运行时间

    return 0;
}
