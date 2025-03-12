#include "GA/GA_algebra.hpp"

int main() {
    Eigen::VectorXd B(6);
    B << 0.8147, 0.9058, 0.1270, 0.9134, 0.6324, 0.0975;
    //std::cout << "B given randomly:\n" << B << std::endl;

    // GA_exp
    Eigen::VectorXd M(8);
    M = GA::GA_exp(B);
    //std::cout << "M calculated by GA_exp:\n" << M << std::endl;

    // GA_log
    B = GA::GA_log(M);
    //std::cout << "B calculated by GA_log:\n" << B << std::endl;

    // GA_AdM
    Eigen::MatrixXd Ad_M(6, 6);
    Ad_M = GA::GA_AdM(M);
    //std::cout << "Ad_M calculated by GA_AdM:\n" << Ad_M << std::endl;

    // GA_adB
    Eigen::MatrixXd ad_B(6, 6);
    ad_B = GA::GA_adB(B);
    //std::cout << "ad_B calculated by GA_adB:\n" << ad_B << std::endl;

    // GA_DHmotor
    Eigen::Vector4d param(4);
    param << 0.8147, 0.9058, 0.1270, 0.9134;
    Eigen::MatrixXd Motor(8, 1);
    Motor = GA::GA_DHmotor(param, 1);
    std::cout << "Motor calculated by GA_DHmotor:\n" << Motor << std::endl;

    return 0;
}