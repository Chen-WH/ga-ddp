#include "GA/GA_algebra.hpp"

namespace GA {

// GA_exp: 计算李代数 bivector B 对应的李群元素 M
Eigen::VectorXd GA_exp(const Eigen::VectorXd& B) {
    // B 是6×1向量，M 是8×1向量
    Eigen::VectorXd M(8);
    Eigen::VectorXd B_trans(3);
    Eigen::VectorXd B_rot(3);

    B_trans = -0.5*B.head(3);
    B_rot = -0.5*B.tail(3);
    double l = B_rot.squaredNorm();
    
    if (l < 1e-9) {
        // 纯平移
        M << 1, B_trans, 0, 0, 0, 0;
    } else {
        // 一般情况
        double m = B_trans.dot(B_rot);
        double a = sqrt(l);
        double c = cos(a);
        double s = sin(a) / a;
        double t = m / l * (c - s);
        M << c, s*B_trans + t*B_rot, s*B_rot, m*s;
    }

    return M;
}

// GA_log: 计算李群元素 M 对应的李代数 bivector B
Eigen::VectorXd GA_log(const Eigen::VectorXd& M) {
    // M 是8×1向量，B 是6×1向量
    Eigen::VectorXd B(6);
    double M0 = M(0);

    if (std::abs(M0 - 1) < 1e-9) {
        B << M.segment(1, 3), Eigen::Vector3d::Zero();
    } else {
        double a = 1 / (1 - M0*M0);    // inv squared length
        double b = acos(M0)*sqrt(a);   // rotation scale
        double c = a*M(7)*(1 - M0*b);  // translation scale
        B << c*M.segment(4, 3) + b*M.segment(1, 3), b*M.segment(4, 3);
    }

    B *= -2;
    return B;
}

// GA_AdM: 计算伴随表示 Ad
Eigen::MatrixXd GA_AdM(const Eigen::VectorXd& M) {
    Eigen::MatrixXd Ad_M(6, 6);

    // 缓存常用的 M 向量的元素，避免重复访问
    double M0 = M(0), M1 = M(1), M2 = M(2), M3 = M(3), M4 = M(4), M5 = M(5), M6 = M(6), M7 = M(7);
    double M0_sq = M0 * M0, M4_sq = M4 * M4, M5_sq = M5 * M5, M6_sq = M6 * M6;

    // 填充 Ad_M 矩阵
    Ad_M(0, 0) = M0_sq + M4_sq - M5_sq - M6_sq;
    Ad_M(0, 1) = 2 * M0 * M6 + 2 * M4 * M5;
    Ad_M(0, 2) = 2 * M4 * M6 - 2 * M0 * M5;
    Ad_M(0, 3) = 2 * M1 * M4 - 2 * M0 * M7 - 2 * M2 * M5 - 2 * M3 * M6;
    Ad_M(0, 4) = 2 * M0 * M3 + 2 * M1 * M5 + 2 * M2 * M4 - 2 * M6 * M7;
    Ad_M(0, 5) = 2 * M1 * M6 - 2 * M0 * M2 + 2 * M3 * M4 + 2 * M5 * M7;

    Ad_M(1, 0) = 2 * M4 * M5 - 2 * M0 * M6;
    Ad_M(1, 1) = M0_sq - M4_sq + M5_sq - M6_sq;
    Ad_M(1, 2) = 2 * M0 * M4 + 2 * M5 * M6;
    Ad_M(1, 3) = 2 * M1 * M5 - 2 * M0 * M3 + 2 * M2 * M4 + 2 * M6 * M7;
    Ad_M(1, 4) = 2 * M2 * M5 - 2 * M0 * M7 - 2 * M1 * M4 - 2 * M3 * M6;
    Ad_M(1, 5) = 2 * M0 * M1 + 2 * M2 * M6 + 2 * M3 * M5 - 2 * M4 * M7;

    Ad_M(2, 0) = 2 * M0 * M5 + 2 * M4 * M6;
    Ad_M(2, 1) = 2 * M5 * M6 - 2 * M0 * M4;
    Ad_M(2, 2) = M0_sq - M4_sq - M5_sq + M6_sq;
    Ad_M(2, 3) = 2 * M0 * M2 + 2 * M1 * M6 + 2 * M3 * M4 - 2 * M5 * M7;
    Ad_M(2, 4) = 2 * M2 * M6 - 2 * M0 * M1 + 2 * M3 * M5 + 2 * M4 * M7;
    Ad_M(2, 5) = 2 * M3 * M6 - 2 * M0 * M7 - 2 * M2 * M5 - 2 * M1 * M4;

    Ad_M(3, 0) = 0;
    Ad_M(3, 1) = 0;
    Ad_M(3, 2) = 0;
    Ad_M(3, 3) = M0_sq + M4_sq - M5_sq - M6_sq;
    Ad_M(3, 4) = 2 * M0 * M6 + 2 * M4 * M5;
    Ad_M(3, 5) = 2 * M4 * M6 - 2 * M0 * M5;

    Ad_M(4, 0) = 0;
    Ad_M(4, 1) = 0;
    Ad_M(4, 2) = 0;
    Ad_M(4, 3) = 2 * M4 * M5 - 2 * M0 * M6;
    Ad_M(4, 4) = M0_sq - M4_sq + M5_sq - M6_sq;
    Ad_M(4, 5) = 2 * M0 * M4 + 2 * M5 * M6;

    Ad_M(5, 0) = 0;
    Ad_M(5, 1) = 0;
    Ad_M(5, 2) = 0;
    Ad_M(5, 3) = 2 * M0 * M5 + 2 * M4 * M6;
    Ad_M(5, 4) = 2 * M5 * M6 - 2 * M0 * M4;
    Ad_M(5, 5) = M0_sq - M4_sq - M5_sq + M6_sq;

    return Ad_M;
}

// GA_adB: 计算伴随表示 ad
Eigen::MatrixXd GA_adB(const Eigen::VectorXd& B) {
    Eigen::MatrixXd ad_B(6, 6);

    // 填充 ad_B 矩阵
    ad_B(0, 0) = 0;
    ad_B(0, 1) = -B(5);
    ad_B(0, 2) = B(4);
    ad_B(0, 3) = 0;
    ad_B(0, 4) = -B(2);
    ad_B(0, 5) = B(1);

    ad_B(1, 0) = B(5);
    ad_B(1, 1) = 0;
    ad_B(1, 2) = -B(3);
    ad_B(1, 3) = B(2);
    ad_B(1, 4) = 0;
    ad_B(1, 5) = -B(0);

    ad_B(2, 0) = -B(4);
    ad_B(2, 1) = B(3);
    ad_B(2, 2) = 0;
    ad_B(2, 3) = -B(1);
    ad_B(2, 4) = B(0);
    ad_B(2, 5) = 0;

    ad_B(3, 0) = 0;
    ad_B(3, 1) = 0;
    ad_B(3, 2) = 0;
    ad_B(3, 3) = 0;
    ad_B(3, 4) = -B(5);
    ad_B(3, 5) = B(4);

    ad_B(4, 0) = 0;
    ad_B(4, 1) = 0;
    ad_B(4, 2) = 0;
    ad_B(4, 3) = B(5);
    ad_B(4, 4) = 0;
    ad_B(4, 5) = -B(3);

    ad_B(5, 0) = 0;
    ad_B(5, 1) = 0;
    ad_B(5, 2) = 0;
    ad_B(5, 3) = -B(4);
    ad_B(5, 4) = B(3);
    ad_B(5, 5) = 0;

    return ad_B;
}

// GA_DHmotor: 计算DH参数对应的Motor
Eigen::MatrixXd GA_DHmotor(const Eigen::Vector4d& param, int type) {
    // DH 参数 (double类型的参数数组 param)
    Eigen::MatrixXd M(8, 1); // 8x1 矩阵

    double d, theta, a, alpha;
    switch (type) {
        case 0: // DH 存储类型, param = [d, theta, a, alpha], M = exp(d*e03 + theta*e12)*exp(a*e01 + alpha*e23)
            d = 0.5*param(0);
            theta = 0.5*param(1);
            a = 0.5*param(2);
            alpha = 0.5*param(3);

            // 计算 DH 变换矩阵
            M(0) = cos(alpha) * cos(theta);
            M(1) =  d * sin(alpha) * sin(theta) - a * cos(alpha) * cos(theta);
            M(2) = -a * cos(alpha) * sin(theta) - d * sin(alpha) * cos(theta);
            M(3) =  a * sin(alpha) * sin(theta) - d * cos(alpha) * cos(theta);
            M(4) = -sin(alpha) * cos(theta);
            M(5) = -sin(alpha) * sin(theta);
            M(6) = -cos(alpha) * sin(theta);
            M(7) = a * sin(alpha) * cos(theta) + d * cos(alpha) * sin(theta);
            break;

        case 1: // MDH 存储类型, param = [a, alpha, d, theta], M = exp(a*e01 + alpha*e23)*exp(d*e03 + theta*e12)
            a = 0.5*param(0);
            alpha = 0.5*param(1);
            d = 0.5*param(2);
            theta = 0.5*param(3);

            // 计算 MDH 变换矩阵
            M(0) = cos(alpha) * cos(theta);
            M(1) = d * sin(alpha) * sin(theta) - a * cos(alpha) * cos(theta);
            M(2) = a * cos(alpha) * sin(theta) + d * sin(alpha) * cos(theta);
            M(3) = a * sin(alpha) * sin(theta) - d * cos(alpha) * cos(theta);
            M(4) = -sin(alpha) * cos(theta);
            M(5) = sin(alpha) * sin(theta);
            M(6) = -cos(alpha) * sin(theta);
            M(7) = a * sin(alpha) * cos(theta) + d * cos(alpha) * sin(theta);
            break;

        default:
            std::cerr << "Unknown DH type!" << std::endl;
            break;
    }
    return M;
}

// GA_dual: 计算对偶
Eigen::VectorXd GA_dual(const Eigen::VectorXd& x) {
    int n = x.size();
    Eigen::VectorXd res;

    switch (n) {
        case 6:
            res.resize(6);
            res << -x(3), -x(4), -x(5), -x(0), -x(1), -x(2);
            break;
        
        case 8:
            res.resize(8);
            res << x(7), -x(4), -x(5), -x(6), -x(1), -x(2), -x(3), x(0);
            break;

        default:
            std::cerr << "Error: Invalid input size, expected 6 or 8!" << std::endl;
            break;
    }

    return res;
}

// GA_forque: 计算受力对应的forque
Eigen::VectorXd GA_forque(const Eigen::VectorXd& q, const Eigen::VectorXd& f) {
    // 确保输入是 3 维向量
    if (q.size() != 3 || f.size() != 3) {
        std::cerr << "Error: Both input vectors must be 3-dimensional!" << std::endl;
        exit(1);
    }

    Eigen::VectorXd F(6);  // 输出是 6 维向量

    F << q(1) * f(2) - q(2) * f(1),  // e01
        q(2) * f(0) - q(0) * f(2),   // e02
        q(0) * f(1) - q(1) * f(0),   // e03
        f(0),                        // e23
        f(1),                        // e31
        f(2);                        // e12

    return F;
}

// GA_metric: 计算李代数上的度规
Eigen::RowVectorXd GA_metric(const Eigen::VectorXd& x) {
    // 确保输入是 6 维向量
    if (x.size() != 6) {
        std::cerr << "Error: Input vector must be 6-dimensional!" << std::endl;
        exit(1);
    }

    Eigen::RowVectorXd res(6);  // 输出是 6 维向量
    // 根据给定的映射规则计算结果
    res << x(3), x(4), x(5), x(0), x(1), x(2);
    return res;
}

// GA_motor: 计算关节驱动量对应的Motor
Eigen::VectorXd GA_motor(const Eigen::VectorXd& L, double q, int type) {
    // 检查输入的L和q的维度
    if (L.size() != 6) {
        std::cerr << "Error: L must be a 6-dimensional vector!" << std::endl;
        exit(1);
    }

    Eigen::VectorXd M(8);  // 假设返回的矩阵大小为6
    q = -0.5*q;

    switch (type) {
        case 1: {  // 旋转关节
            M << cos(q), sin(q) * L, 0;
            break;
        }
        case 2: {  // 平移关节
            M << 1, q * L.head(3), 0, 0, 0, 0;
            break;
        }
        default:
            std::cerr << "Error: Invalid type!" << std::endl;
            exit(1);
    }

    return M;
}

// GA_prodM: 计算两个Motor的乘积
Eigen::VectorXd GA_prodM(const Eigen::VectorXd& M1, const Eigen::VectorXd& M2) {
    // 确保输入 M1 和 M2 都是长度为 8 的向量
    if (M1.size() != 8 || M2.size() != 8) {
        std::cerr << "Error: M1 and M2 must both be 8-dimensional vectors!" << std::endl;
        exit(1);
    }

    Eigen::VectorXd M(8);  // 返回的新变换矩阵大小为 8

    M(0) = M1(0)*M2(0) - M1(4)*M2(4) - M1(5)*M2(5) - M1(6)*M2(6);
    M(1) = M1(0)*M2(1) + M1(1)*M2(0) - M1(2)*M2(6) + M1(3)*M2(5) - M1(5)*M2(3) + M1(6)*M2(2) - M1(4)*M2(7) - M1(7)*M2(4);
    M(2) = M1(0)*M2(2) + M1(2)*M2(0) + M1(1)*M2(6) - M1(3)*M2(4) + M1(4)*M2(3) - M1(6)*M2(1) - M1(5)*M2(7) - M1(7)*M2(5);
    M(3) = M1(0)*M2(3) + M1(3)*M2(0) - M1(1)*M2(5) + M1(2)*M2(4) - M1(4)*M2(2) + M1(5)*M2(1) - M1(6)*M2(7) - M1(7)*M2(6);
    M(4) = M1(0)*M2(4) + M1(4)*M2(0) - M1(5)*M2(6) + M1(6)*M2(5);
    M(5) = M1(0)*M2(5) + M1(5)*M2(0) + M1(4)*M2(6) - M1(6)*M2(4);
    M(6) = M1(0)*M2(6) + M1(6)*M2(0) - M1(4)*M2(5) + M1(5)*M2(4);
    M(7) = M1(1)*M2(4) + M1(4)*M2(1) + M1(0)*M2(7) + M1(2)*M2(5) + M1(5)*M2(2) + M1(7)*M2(0) + M1(3)*M2(6) + M1(6)*M2(3);

    return M;
}

// GA_Rev: 计算 motor 群的逆元
Eigen::VectorXd GA_Rev(const Eigen::VectorXd& M) {
    // 确保输入 M 是长度为 8 的向量
    if (M.size() != 8) {
        std::cerr << "Error: M must be an 8-dimensional vector!" << std::endl;
        exit(1);
    }

    Eigen::VectorXd iM(8);  // 结果逆 motor
    iM(0) = M(0);
    iM(7) = M(7);

    // 反转其他元素的符号
    iM.segment(1, 6) = -M.segment(1, 6);

    return iM;
}

// PGA_Lie2: 计算李括号 dL = [L, B] = -(LB - BL)/2
Eigen::VectorXd PGA_Lie2(const Eigen::VectorXd& L, const Eigen::VectorXd& B) {
    // 确保输入向量 L 和 B 都是长度为 6 的向量
    if (L.size() != 6 || B.size() != 6) {
        std::cerr << "Error: Both L and B must be 6-dimensional vectors!" << std::endl;
        exit(1);
    }

    Eigen::VectorXd dL(6);  // 输出 bivector dL

    dL(0) = B(2)*L(4) - B(1)*L(5) - B(4)*L(2) + B(5)*L(1);
    dL(1) = B(0)*L(5) - B(2)*L(3) + B(3)*L(2) - B(5)*L(0);
    dL(2) = B(1)*L(3) - B(0)*L(4) - B(3)*L(1) + B(4)*L(0);
    dL(3) = B(5)*L(4) - B(4)*L(5);
    dL(4) = B(3)*L(5) - B(5)*L(3);
    dL(5) = B(4)*L(3) - B(3)*L(4);

    return dL;
}

// PGA_Frame2：通过刚体运动变换 Lp = M * L * iM
Eigen::VectorXd PGA_Frame2(const Eigen::VectorXd& L, const Eigen::VectorXd& M) {
    // 确保输入向量 L 和 M 的长度分别为 6 和 8
    if (L.size() != 6 || M.size() != 8) {
        std::cerr << "Error: L must be a 6-dimensional vector, M must be an 8-dimensional vector!" << std::endl;
        exit(1);
    }

    Eigen::VectorXd Lp(6);  // 输出 bivector Lp

    Lp(0) = L(0)*(M(0)*M(0) + M(4)*M(4) - M(5)*M(5) - M(6)*M(6)) 
            + 2*( L(1)*M(0)*M(6) - L(2)*M(0)*M(5) + L(4)*M(0)*M(3) 
            - L(5)*M(0)*M(2) + L(3)*M(1)*M(4) + L(1)*M(4)*M(5) 
            - L(3)*M(0)*M(7) - L(3)*M(2)*M(5) + L(4)*M(1)*M(5) 
            + L(4)*M(2)*M(4) + L(2)*M(4)*M(6) - L(3)*M(3)*M(6) 
            + L(5)*M(1)*M(6) + L(5)*M(3)*M(4) - L(4)*M(6)*M(7) 
            + L(5)*M(5)*M(7));

    Lp(1) = L(1)*(M(0)*M(0) - M(4)*M(4) + M(5)*M(5) - M(6)*M(6)) 
            + 2*(-L(0)*M(0)*M(6) + L(2)*M(0)*M(4) - L(3)*M(0)*M(3) 
            + L(5)*M(0)*M(1) + L(0)*M(4)*M(5) + L(3)*M(1)*M(5) 
            + L(3)*M(2)*M(4) - L(4)*M(1)*M(4) - L(4)*M(0)*M(7) 
            + L(4)*M(2)*M(5) + L(2)*M(5)*M(6) - L(4)*M(3)*M(6) 
            + L(5)*M(2)*M(6) + L(5)*M(3)*M(5) + L(3)*M(6)*M(7) 
            - L(5)*M(4)*M(7));

    Lp(2) = L(2)*(M(0)*M(0) - M(4)*M(4) - M(5)*M(5) + M(6)*M(6)) 
            + 2*( L(0)*M(0)*M(5) - L(1)*M(0)*M(4) + L(3)*M(0)*M(2) 
            - L(4)*M(0)*M(1) + L(0)*M(4)*M(6) + L(3)*M(1)*M(6) 
            + L(3)*M(3)*M(4) - L(5)*M(1)*M(4) + L(1)*M(5)*M(6) 
            + L(4)*M(2)*M(6) + L(4)*M(3)*M(5) - L(5)*M(0)*M(7) 
            - L(5)*M(2)*M(5) + L(5)*M(3)*M(6) - L(3)*M(5)*M(7) 
            + L(4)*M(4)*M(7));

    Lp(3) = L(3)*(M(0)*M(0) + M(4)*M(4) - M(5)*M(5) - M(6)*M(6)) 
            + 2*(-L(5)*M(0)*M(5) + L(4)*M(0)*M(6) + L(4)*M(4)*M(5) 
            + L(5)*M(4)*M(6));

    Lp(4) = L(4)*(M(0)*M(0) - M(4)*M(4) + M(5)*M(5) - M(6)*M(6)) 
            + 2*( L(5)*M(0)*M(4) - L(3)*M(0)*M(6) + L(3)*M(4)*M(5) 
            + L(5)*M(5)*M(6));

    Lp(5) = L(5)*(M(0)*M(0) - M(4)*M(4) - M(5)*M(5) + M(6)*M(6)) 
            + 2*(-L(4)*M(0)*M(4) + L(3)*M(0)*M(5) + L(3)*M(4)*M(6) 
            + L(4)*M(5)*M(6));

    return Lp;
}

} // namespace GA