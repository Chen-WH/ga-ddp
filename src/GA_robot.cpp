#include "GA/GA_robot.hpp"

namespace GA {

// 反对称矩阵
Eigen::Matrix3d skew(const Eigen::Vector3d& v) {
    Eigen::Matrix3d skew_matrix;
    skew_matrix <<  0,   -v(2),   v(1),
                    v(2),   0,    -v(0),
                    -v(1),   v(0),   0;
    return skew_matrix;
}

// 构造函数，初始化机器人参数
GA_robot::GA_robot(const std::string& robot_model, double gravity) : model(robot_model), g(gravity) {
    Eigen::MatrixXd Nc, rc;
    // 根据模型名称初始化机器人参数
    if (model == "chin") {
        n = 6;
        Mj.resize(8, n+1);
        Lj.resize(6, n);
        type.resize(n);
        Ij.resize(n);

        Mj.col(0) = GA_DHmotor(Eigen::Vector4d(    0,    M_PI,   -0.283,       0), 1);
        Mj.col(1) = GA_DHmotor(Eigen::Vector4d(    0,  M_PI_2,        0,  M_PI_2), 1);
        Mj.col(2) = GA_DHmotor(Eigen::Vector4d( -0.4,    M_PI,        0,       0), 1);
        Mj.col(3) = GA_DHmotor(Eigen::Vector4d(-0.36,    M_PI, -0.19875, -M_PI_2), 1);
        Mj.col(4) = GA_DHmotor(Eigen::Vector4d(    0, -M_PI_2, -0.16975,       0), 1);
        Mj.col(5) = GA_DHmotor(Eigen::Vector4d(    0,  M_PI_2, -0.14875,       0), 1);
        Mj.col(6) = GA_DHmotor(Eigen::Vector4d(    0,       0,        0,       0), 1);  // 末端执行器

        // 初始化 Lj 和 type
        for (int i = 0; i < n; ++i) {
            Lj.col(i) << 0.0, 0.0, 0.0, 0.0, 0.0, 1.0;  // 所有关节的轴线都在 z 轴上
            type(i) = 1;
        }

        // 计算惯量矩阵 Ij
        Nc.resize(n, 7);  // 惯量矩阵参数
        Nc << 5.675, 1.1632,     0, 0.001, 1.1632, 0, 0.0125,
              5.466,  0.014, 0.009,     0,  0.124, 0,  0.122,
              4.294,   0.01, -0.01,     0,  0.081, 0,  0.076,
              2.732,  0.005,     0,     0,  0.003, 0,  0.005,
              2.732,  0.005,     0,     0,  0.003, 0,  0.005,
               0.28,   0.01,     0,     0,   0.01, 0,   0.01;

        rc.resize(n, 3);  // 质量中心位置
        rc << 0,  0.043,  0.005,
         -0.298,      0, -0.126,
         -0.274,      0, -0.061,
              0, -0.055,  0.001,
              0,  0.055,  0.001,
              0,      0,  0.017;
    }
    else if (model == "jaka") {
        n = 6;
        Mj.resize(8, n+1);
        Lj.resize(6, n);
        type.resize(n);
        Ij.resize(n);

        Mj.col(0) = GA_DHmotor(Eigen::Vector4d(      0,       0, 0.10265, 0), 1);
        Mj.col(1) = GA_DHmotor(Eigen::Vector4d(      0,  M_PI_2,       0, 0), 1);
        Mj.col(2) = GA_DHmotor(Eigen::Vector4d(  0.595,       0,       0, 0), 1);
        Mj.col(3) = GA_DHmotor(Eigen::Vector4d( 0.5715,       0, -0.1315, 0), 1);
        Mj.col(4) = GA_DHmotor(Eigen::Vector4d(      0,  M_PI_2,   0.115, 0), 1);
        Mj.col(5) = GA_DHmotor(Eigen::Vector4d(      0, -M_PI_2,  0.1035, 0), 1);
        Mj.col(6) = GA_DHmotor(Eigen::Vector4d(      0,       0,       0, 0), 1);  // 末端执行器

        // 初始化 Lj 和 type
        for (int i = 0; i < n; ++i) {
            Lj.col(i) << 0.0, 0.0, 0.0, 0.0, 0.0, 1.0;  // 所有关节的轴线都在 z 轴���
            type(i) = 1;
        }

        // 计算惯量矩阵 Ij
        Nc.resize(n, 7);  // 惯量矩阵参数
        Nc <<  3.78,   0.0404051529004393, -5.35916937084468e-05,  8.56853276426298e-05,   0.0404457306741405,  -0.00039414179882763,  0.0343161527751278,
              12.44,    0.116456439495809,  -0.00156765863153774, -3.42571922960543e-06,     1.66880550591637, -1.42545474100985e-05,    1.65302926902866,
               5.06,   0.0212898342592158, -1.15196860083347e-05,  -0.00344346825692792,    0.497633875998665, -6.35439617896331e-07,   0.498513225896089,
              0.883,  0.00458464465961989,  2.63659533420518e-06, -1.52605286509264e-05,  0.00363879515047717, -3.19585101806877e-05, 0.00455871160423887,
                  1,  0.00702353102292446, -2.72522678679305e-06,  1.52464540194323e-05,  0.00413469714863856,  2.81695088352211e-05, 0.00699581168456812,
               0.22, 0.000597975009893646, -2.41634474044304e-05, -6.20059067536115e-06, 0.000647976470675928, -6.92746756534325e-07, 0.00106463342891045;

        rc.resize(n, 3);  // 质量中心位置
        rc << -0.000328431560018958,   0.00406089906856803,  -0.0250109305609668,
                  0.297499939118518, -1.25165127429216e-08,   -0.166069989366573,
                  0.294233782543355, -4.72319252706188e-06,  -0.0241861096678911,
               4.22017696188881e-06,   -0.0150343523848632,  0.00216593552372366,
               3.86594587120648e-06,   0.00426530078104992, -0.00184283112443297,
              -0.000786084663781494,  3.88801034911165e-05,  -0.0163246941858815;
    }
    else if (model == "panda") {
        n = 7;
        Mj.resize(8, n+1);
        Lj.resize(6, n);
        type.resize(n);
        Ij.resize(n);

        Mj.col(0) = GA_DHmotor(Eigen::Vector4d(      0,       0, 0.333, 0), 1);
        Mj.col(1) = GA_DHmotor(Eigen::Vector4d(      0, -M_PI_2,     0, 0), 1);
        Mj.col(2) = GA_DHmotor(Eigen::Vector4d(      0,  M_PI_2, 0.316, 0), 1);
        Mj.col(3) = GA_DHmotor(Eigen::Vector4d( 0.0825,  M_PI_2,     0, 0), 1);
        Mj.col(4) = GA_DHmotor(Eigen::Vector4d(-0.0825, -M_PI_2, 0.384, 0), 1);
        Mj.col(5) = GA_DHmotor(Eigen::Vector4d(      0,  M_PI_2,     0, 0), 1);
        Mj.col(6) = GA_DHmotor(Eigen::Vector4d(  0.088,  M_PI_2,     0, 0), 1);
        Mj.col(7) = GA_DHmotor(Eigen::Vector4d(      0,       0, 0.107, 0), 1);  // 末端执行器

        // 初始化 Lj 和 type
        for (int i = 0; i < n; ++i) {
            Lj.col(i) << 0.0, 0.0, 0.0, 0.0, 0.0, 1.0;  // 所有关节的轴线都在 z 轴���
            type(i) = 1;
        }

        // 计算惯量矩阵 Ij
        Nc.resize(n, 7);  // 惯量矩阵参数
        Nc << 4.970684,    0.70337, -0.00013900,   0.0067720,    0.70661,    0.019169,  0.0091170,
              0.646926,  0.0079620,  -3.9250e-3,  1.0254e-02, 2.8110e-02,  7.0400e-04, 2.5995e-02,
              3.228604, 3.7242e-02, -4.7610e-03, -1.1396e-02, 3.6155e-02, -1.2805e-02, 1.0830e-02,
              3.587895, 2.5853e-02,  7.7960e-03, -1.3320e-03, 1.9552e-02,  8.6410e-03, 2.8323e-02,
              1.225946, 3.5549e-02, -2.1170e-03, -4.0370e-03, 2.9474e-02,  2.2900e-04, 8.6270e-03,
              1.666555, 1.9640e-03,  1.0900e-04, -1.1580e-03, 4.3540e-03,  3.4100e-04, 5.4330e-03,
           7.35522e-01, 1.2516e-02, -4.2800e-04, -1.1960e-03, 1.0027e-02, -7.4100e-04, 4.8150e-03;

        rc.resize(n, 3);  // 质量中心位置
        rc << 0.003875,    0.002081,    -0.04762,
             -0.003141,    -0.02872,    0.003495,
                -0.274,           0,      -0.061,
            2.7518e-02,  3.9252e-02, -6.6502e-02,
            -5.317e-02, 1.04419e-01,  2.7454e-02,
           -1.1953e-02,  4.1065e-02, -3.8437e-02,
            6.0149e-02, -1.4117e-02, -1.0517e-02,
            1.0517e-02,  -4.252e-03,  6.1597e-02;
    }
    else {
        std::cerr << "Unsupported robot model: " << robot_model << std::endl;
    }
    // 计算广义惯量矩阵 Ij
    for (int i = 0; i < n; ++i) {
        Eigen::Matrix3d I;
        I << Nc(i, 1), Nc(i, 2), Nc(i, 3),
                Nc(i, 2), Nc(i, 4), Nc(i, 5),
                Nc(i, 3), Nc(i, 5), Nc(i, 6);
        Eigen::Matrix3d Rc = skew(rc.row(i));

        Eigen::MatrixXd Ij_full(6, 6);
        Ij_full.block(0, 0, 3, 3) = Nc(i, 0) * Rc;
        Ij_full.block(0, 3, 3, 3) = I - Nc(i, 0) * Rc * Rc;
        Ij_full.block(3, 0, 3, 3) = Nc(i, 0) * Eigen::Matrix3d::Identity();
        Ij_full.block(3, 3, 3, 3) = -Nc(i, 0) * Rc;

        Ij[i] = Ij_full;
    }
}

// 正运动学：计算末端执行器相对于基座的位姿
Eigen::VectorXd GA_robot::fkine(const Eigen::VectorXd& q) {
    Eigen::VectorXd M(8);
    M = Mj.col(0);
    for (int i = 0; i < n; ++i) {
        M = GA_prodM(GA_prodM(M, GA_motor(Lj.col(i), q(i), type(i))), Mj.col(i+1));
    }
    return M;
}

// 逆动力学：计算机器人各关节的力矩
Eigen::VectorXd GA_robot::idyn(const Eigen::VectorXd& q, const Eigen::VectorXd& dq, const Eigen::VectorXd& ddq, const Eigen::VectorXd& fext) {
    Eigen::MatrixXd M = Eigen::MatrixXd::Zero(8, n);
    Eigen::MatrixXd Mt = Eigen::MatrixXd::Zero(8, n);
    Eigen::MatrixXd V = Eigen::MatrixXd::Zero(6, n);
    Eigen::MatrixXd dV = Eigen::MatrixXd::Zero(6, n);
    Eigen::MatrixXd F = Eigen::MatrixXd::Zero(6, n);
    Eigen::VectorXd tau = Eigen::VectorXd::Zero(n);
    Eigen::VectorXd gravity(6);
    gravity << 0, 0, g, 0, 0, 0;

    // Forward iterations
    M.col(0) = GA_prodM(Mj.col(0), GA_motor(Lj.col(0), q(0), type(0)));
    Mt.col(0) = GA_Rev(M.col(0));
    V.col(0) = dq(0) * Lj.col(0);
    dV.col(0) = PGA_Frame2(gravity, Mt.col(0)) + ddq(0) * Lj.col(0);
    for (int i = 1; i < n; ++i) {
        M.col(i) = GA_prodM(Mj.col(i), GA_motor(Lj.col(i), q(i), type(i)));
        Mt.col(i) = GA_Rev(M.col(i));
        V.col(i) = PGA_Frame2(V.col(i - 1), Mt.col(i)) + dq(i) * Lj.col(i);
        dV.col(i) = PGA_Frame2(dV.col(i - 1), Mt.col(i)) + dq(i) * PGA_Lie2(V.col(i), Lj.col(i)) + ddq(i) * Lj.col(i);
    }

    // Backward iterations
    F.col(n - 1) = -PGA_Frame2(fext, Mj.col(n)) + Ij[n - 1] * dV.col(n - 1) - PGA_Lie2(Ij[n - 1] * V.col(n - 1), V.col(n - 1));
    tau(n - 1) = GA_metric(Lj.col(n - 1))*F.col(n - 1);
    for (int i = n - 2; i >= 0; --i) {
        F.col(i) = PGA_Frame2(F.col(i + 1), M.col(i + 1)) + Ij[i] * dV.col(i) - PGA_Lie2(Ij[i] * V.col(i), V.col(i));
        tau(i) = GA_metric(Lj.col(i))*F.col(i);
    }
    return tau;
}

// 正动力学：计算机器人各关节的加速度
Eigen::VectorXd GA_robot::fdyn(const Eigen::VectorXd& q, const Eigen::VectorXd& dq, const Eigen::VectorXd& tau, const Eigen::VectorXd& fext) {
    Eigen::MatrixXd M = Eigen::MatrixXd::Zero(8, n);
    Eigen::MatrixXd Mt = Eigen::MatrixXd::Zero(8, n);
    Eigen::MatrixXd V = Eigen::MatrixXd::Zero(6, n);
    Eigen::MatrixXd c = Eigen::MatrixXd::Zero(6, n);
    Eigen::MatrixXd b = Eigen::MatrixXd::Zero(6, n);
    Eigen::VectorXd d = Eigen::VectorXd::Zero(n);
    std::vector<Eigen::MatrixXd> hI(n, Eigen::MatrixXd::Zero(6, 6));
    std::vector<Eigen::MatrixXd> I = Ij;
    Eigen::MatrixXd hb = Eigen::MatrixXd::Zero(6, n);
    Eigen::MatrixXd dV = Eigen::MatrixXd::Zero(6, n);
    Eigen::VectorXd ddq = Eigen::VectorXd::Zero(n);
    Eigen::VectorXd gravity(6);
    gravity << 0, 0, g, 0, 0, 0;

    // Forward iterations 速度与偏置项
    M.col(0) = GA_prodM(Mj.col(0), GA_motor(Lj.col(0), q(0), type(0)));
    Mt.col(0) = GA_Rev(M.col(0));
    V.col(0) = dq(0) * Lj.col(0);
    b.col(0) = -PGA_Lie2(I[0] * V.col(0), V.col(0));
    for (int i = 1; i < n; ++i) {
        M.col(i) = GA_prodM(Mj.col(i), GA_motor(Lj.col(i), q(i), type(i)));
        Mt.col(i) = GA_Rev(M.col(i));
        V.col(i) = PGA_Frame2(V.col(i - 1), Mt.col(i)) + dq(i) * Lj.col(i);
        c.col(i) = dq(i) * PGA_Lie2(V.col(i), Lj.col(i));
        b.col(i) = -PGA_Lie2(I[i] * V.col(i), V.col(i));
    }

    // Backward iterations 铰接体惯性与偏置力
    for (int i = n - 1; i >= 1; --i) {
        d(i) = 1 / (GA_metric(Lj.col(i))*I[i] * Lj.col(i));
        hI[i] = I[i] - d(i) * I[i] * Lj.col(i) * GA_metric(Lj.col(i)) * I[i];
        hb.col(i) = b.col(i) + hI[i] * c.col(i) + d(i) * I[i] * Lj.col(i) * (tau(i) - GA_metric(Lj.col(i)) * b.col(i));
        I[i - 1] = I[i - 1] + GA_AdM(M.col(i)) * hI[i] * GA_AdM(Mt.col(i));
        b.col(i - 1) = b.col(i - 1) + PGA_Frame2(hb.col(i), M.col(i));
    }
    d(0) = 1 / (GA_metric(Lj.col(0)) * I[0] * Lj.col(0));

    // Forward iterations 加速度
    dV.col(0) = PGA_Frame2(gravity, Mt.col(0)) + c.col(0);
    ddq(0) = d(0) * (tau(0) - GA_metric(Lj.col(0)) * (b.col(0) + I[0] * dV.col(0)));
    dV.col(0) = dV.col(0) + ddq(0) * Lj.col(0);
    for (int i = 1; i < n; ++i) {
        dV.col(i) = PGA_Frame2(dV.col(i - 1), Mt.col(i)) + c.col(i);
        ddq(i) = d(i) * (tau(i) - GA_metric(Lj.col(i)) * (b.col(i) + I[i] * dV.col(i)));
        dV.col(i) = dV.col(i) + ddq(i) * Lj.col(i);
    }
    return ddq;
}

// 逆动力学一阶偏导：计算机器人各关节力矩对关节位置、速度、加速度的偏导
void GA_robot::idyn_fo(const Eigen::VectorXd& q, const Eigen::VectorXd& dq, const Eigen::VectorXd& ddq, const Eigen::VectorXd& fext, Eigen::MatrixXd& ptau_pq, Eigen::MatrixXd& ptau_pdq, Eigen::MatrixXd& ptau_pddq) {
    Eigen::MatrixXd M = Eigen::MatrixXd::Zero(8, n);
    Eigen::MatrixXd Mt = Eigen::MatrixXd::Zero(8, n);
    Eigen::MatrixXd V = Eigen::MatrixXd::Zero(6, n);
    std::vector<Eigen::MatrixXd> pV_pq(n, Eigen::MatrixXd::Zero(6, n));
    std::vector<Eigen::MatrixXd> pV_pdq(n, Eigen::MatrixXd::Zero(6, n));
    Eigen::MatrixXd dV = Eigen::MatrixXd::Zero(6, n);
    std::vector<Eigen::MatrixXd> pdV_pq(n, Eigen::MatrixXd::Zero(6, n));
    std::vector<Eigen::MatrixXd> pdV_pdq(n, Eigen::MatrixXd::Zero(6, n));
    std::vector<Eigen::MatrixXd> pdV_pddq(n, Eigen::MatrixXd::Zero(6, n));
    Eigen::MatrixXd F = Eigen::MatrixXd::Zero(6, n);
    std::vector<Eigen::MatrixXd> pF_pq(n, Eigen::MatrixXd::Zero(6, n));
    std::vector<Eigen::MatrixXd> pF_pdq(n, Eigen::MatrixXd::Zero(6, n));
    std::vector<Eigen::MatrixXd> pF_pddq(n, Eigen::MatrixXd::Zero(6, n));
    ptau_pq.resize(n, n);
    ptau_pdq.resize(n, n);
    ptau_pddq.resize(n, n);
    Eigen::VectorXd gravity(6);
    gravity << 0, 0, g, 0, 0, 0;

    // Forward iterations
    M.col(0) = GA_prodM(Mj.col(0), GA_motor(Lj.col(0), q(0), type(0)));
    Mt.col(0) = GA_Rev(M.col(0));
    V.col(0) = dq(0) * Lj.col(0);
    dV.col(0) = PGA_Frame2(gravity, Mt.col(0)) + ddq(0) * Lj.col(0);
    for (int i = 0; i < n; ++i) {
        pV_pq[0].col(i) = Eigen::VectorXd::Zero(6);
        pV_pdq[0].col(i) = Eigen::VectorXd::Zero(6);
        pdV_pq[0].col(i) = Eigen::VectorXd::Zero(6);
        pdV_pdq[0].col(i) = Eigen::VectorXd::Zero(6);
        pdV_pddq[0].col(i) = Eigen::VectorXd::Zero(6);
    }
    pV_pdq[0].col(0) = Lj.col(0);
    pdV_pddq[0].col(0) = Lj.col(0);
    for (int k = 1; k < n; ++k) {
        // Kinect pass
        M.col(k) = GA_prodM(Mj.col(k), GA_motor(Lj.col(k), q(k), type(k)));
        Mt.col(k) = GA_Rev(M.col(k));
        V.col(k) = PGA_Frame2(V.col(k - 1), Mt.col(k)) + dq(k) * Lj.col(k);
        dV.col(k) = PGA_Frame2(dV.col(k - 1), Mt.col(k)) + dq(k) * PGA_Lie2(V.col(k), Lj.col(k)) + ddq(k) * Lj.col(k);
        // first-order pass
        for (int i = 0; i < n; ++i) {
            if (i != k) {
                pV_pq[k].col(i) = PGA_Frame2(pV_pq[k - 1].col(i), Mt.col(k));
                pV_pdq[k].col(i) = PGA_Frame2(pV_pdq[k - 1].col(i), Mt.col(k));
                pdV_pq[k].col(i) = PGA_Frame2(pdV_pq[k - 1].col(i), Mt.col(k)) + PGA_Lie2(dq(k) * pV_pq[k].col(i), Lj.col(k));
                pdV_pdq[k].col(i) = PGA_Frame2(pdV_pdq[k - 1].col(i), Mt.col(k)) + PGA_Lie2(dq(k) * pV_pdq[k].col(i), Lj.col(k));
                pdV_pddq[k].col(i) = PGA_Frame2(pdV_pddq[k - 1].col(i), Mt.col(k));
            } else {
                pV_pq[k].col(i) = PGA_Frame2(pV_pq[k - 1].col(i), Mt.col(k)) + PGA_Lie2(V.col(k), Lj.col(k));
                pV_pdq[k].col(i) = PGA_Frame2(pV_pdq[k - 1].col(i), Mt.col(k)) + Lj.col(k);
                pdV_pq[k].col(i) = PGA_Frame2(pdV_pq[k - 1].col(i), Mt.col(k)) + PGA_Lie2(PGA_Frame2(dV.col(k - 1), Mt.col(k)) + dq(k) * pV_pq[k].col(i), Lj.col(k));
                pdV_pdq[k].col(i) = PGA_Frame2(pdV_pdq[k - 1].col(i), Mt.col(k)) + PGA_Lie2(dq(k) * pV_pdq[k].col(i) + V.col(k), Lj.col(k));
                pdV_pddq[k].col(i) = PGA_Frame2(pdV_pddq[k - 1].col(i), Mt.col(k)) + Lj.col(k);
            }
        }
    }

    // Backward iterations
    F.col(n - 1) = -PGA_Frame2(fext, Mj.col(n)) + Ij[n - 1] * dV.col(n - 1) - PGA_Lie2(Ij[n - 1] * V.col(n - 1), V.col(n - 1));
    for (int i = 0; i < n; ++i) {
        pF_pq[n - 1].col(i) = Ij[n - 1] * pdV_pq[n - 1].col(i) - PGA_Lie2(Ij[n - 1] * V.col(n - 1), pV_pq[n - 1].col(i)) - PGA_Lie2(Ij[n - 1] * pV_pq[n - 1].col(i), V.col(n - 1));
        pF_pdq[n - 1].col(i) = Ij[n - 1] * pdV_pdq[n - 1].col(i) - PGA_Lie2(Ij[n - 1] * V.col(n - 1), pV_pdq[n - 1].col(i)) - PGA_Lie2(Ij[n - 1] * pV_pdq[n - 1].col(i), V.col(n - 1));
        pF_pddq[n - 1].col(i) = Ij[n - 1] * pdV_pddq[n - 1].col(i);
        ptau_pq(n - 1, i) = GA_metric(Lj.col(n - 1)) * pF_pq[n - 1].col(i);
        ptau_pdq(n - 1, i) = GA_metric(Lj.col(n - 1)) * pF_pdq[n - 1].col(i);
        ptau_pddq(n - 1, i) = GA_metric(Lj.col(n - 1)) * pF_pddq[n - 1].col(i);
    }
    for (int k = n - 2; k >= 0; --k) {
        F.col(k) = PGA_Frame2(F.col(k + 1), M.col(k + 1)) + Ij[k] * dV.col(k) - PGA_Lie2(Ij[k] * V.col(k), V.col(k));
        for (int i = 0; i < n; ++i) {
            if (i != k + 1) {
                pF_pq[k].col(i) = PGA_Frame2(pF_pq[k + 1].col(i), M.col(k + 1)) + Ij[k] * pdV_pq[k].col(i) - PGA_Lie2(Ij[k] * V.col(k), pV_pq[k].col(i)) - PGA_Lie2(Ij[k] * pV_pq[k].col(i), V.col(k));
                pF_pdq[k].col(i) = PGA_Frame2(pF_pdq[k + 1].col(i), M.col(k + 1)) + Ij[k] * pdV_pdq[k].col(i) - PGA_Lie2(Ij[k] * V.col(k), pV_pdq[k].col(i)) - PGA_Lie2(Ij[k] * pV_pdq[k].col(i), V.col(k));
                pF_pddq[k].col(i) = PGA_Frame2(pF_pddq[k + 1].col(i), M.col(k + 1)) + Ij[k] * pdV_pddq[k].col(i);
            } else {
                pF_pq[k].col(i) = PGA_Frame2(pF_pq[k + 1].col(i) + PGA_Lie2(Lj.col(k + 1), F.col(k + 1)), M.col(k + 1)) + Ij[k] * pdV_pq[k].col(i) - PGA_Lie2(Ij[k] * V.col(k), pV_pq[k].col(i)) - PGA_Lie2(Ij[k] * pV_pq[k].col(i), V.col(k));
                pF_pdq[k].col(i) = PGA_Frame2(pF_pdq[k + 1].col(i), M.col(k + 1)) + Ij[k] * pdV_pdq[k].col(i) - PGA_Lie2(Ij[k] * V.col(k), pV_pdq[k].col(i)) - PGA_Lie2(Ij[k] * pV_pdq[k].col(i), V.col(k));
                pF_pddq[k].col(i) = PGA_Frame2(pF_pddq[k + 1].col(i), M.col(k + 1)) + Ij[k] * pdV_pddq[k].col(i);
            }
            ptau_pq(k, i) = GA_metric(Lj.col(k)) * pF_pq[k].col(i);
            ptau_pdq(k, i) = GA_metric(Lj.col(k)) * pF_pdq[k].col(i);
            ptau_pddq(k, i) = GA_metric(Lj.col(k)) * pF_pddq[k].col(i);
        }
    }
    return;
}

// 正动力学一阶偏导：计算机器人各关节加速度对关节位置、速度、力矩的偏导
void GA_robot::fdyn_fo(const Eigen::VectorXd& q, const Eigen::VectorXd& dq, const Eigen::VectorXd& tau, const Eigen::VectorXd& fext, Eigen::MatrixXd& pddq_pq, Eigen::MatrixXd& pddq_pdq, Eigen::MatrixXd& pddq_ptau) {
    Eigen::VectorXd ddq = GA_robot::fdyn(q, dq, tau, fext);
    Eigen::MatrixXd ptau_pq, ptau_pdq, ptau_pddq;
    GA_robot::idyn_fo(q, dq, ddq, fext, ptau_pq, ptau_pdq, ptau_pddq);
    pddq_ptau = ptau_pddq.inverse();
    pddq_pq = -pddq_ptau * ptau_pq;
    pddq_pdq = -pddq_ptau * ptau_pdq;
}

} // namespace GA