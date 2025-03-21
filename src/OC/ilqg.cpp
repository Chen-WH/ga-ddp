#include "OC/ilqg.h"
#include <algorithm>

namespace OC {

// 定义 sgn 函数
double sgn(double val) {
    return (val > 0) - (val < 0);
}

std::vector<double> OC::iLQG::Alpha = {1.0000, 0.5012, 0.2512, 0.1259, 0.0631, 0.0316, 0.0158, 0.0079, 0.0040, 0.0020, 0.0010};

// 定义静态常量变量
constexpr double iLQG::lambdaFactor;
constexpr double iLQG::lambdaMax;
constexpr double iLQG::lambdaMin;
constexpr double iLQG::zMin;

// Constructor
iLQG::iLQG(Model* p_dyn, Eigen::MatrixXd Q_k, Eigen::MatrixXd R_k, Eigen::MatrixXd Q_T) {
  model.reset(p_dyn);
  Qk = Q_k;
  Rk = R_k;
  QT = Q_T;
  omp_set_num_threads(8);
}

// Initialize trajectory with control sequence
void iLQG::init_traj() {
  cost_s = 0;
  xs.resize(model->x_dims, T + 1);
  xs.col(0) = x0;


  for (int t = 0; t < T; ++t) {
    us.col(t) = us.col(t).cwiseMin(model->u_max).cwiseMax(model->u_min); // clamp controls to limits
    cost_s += 0.5*((xs.col(t) - xd).transpose() * Qk * (xs.col(t) - xd) + us.col(t).transpose() * Rk * us.col(t)).value();
    xs.col(t + 1) = model->dynamics(xs.col(t), us.col(t), dt);
  }

  cost_s += 0.5*((xs.col(T) - xd).transpose() * QT * (xs.col(T) - xd)).value();

  trace.resize(4, maxIter);
  fx.resize(T);
  fu.resize(T);
  Vx.resize(model->x_dims, T + 1);
  Vxx.resize(T + 1);
  k.resize(model->u_dims, T);
  K.resize(T);
  lambda = 1;
  dlambda = 1;
  // std::cout << "Initial cost: " << cost_s << std::endl;
  return;
}

// Warm-start
void iLQG::generate_trajectory(Eigen::VectorXd& x_0, Eigen::VectorXd& x_d, Eigen::MatrixXd& u_0) {
  x0 = x_0;
  xd = x_d;
  T = u_0.cols();
  us = u_0;
  init_traj();

  // constants, timers, counters
  bool flgChange = true;
  double dcost = 0;
  double z = 0;
  double expected = 0;
  int diverge = 0; // flag for backward pass
  double g_norm, alpha;
  Eigen::MatrixXd u_ff;
  
  double t_differentiate = 0;
  double t_backward = 0;
  double t_forward = 0;
  auto all_start = std::chrono::high_resolution_clock::now();

  int iter;
  for (iter = 0; iter < maxIter; ++iter) {

    //--------------------------------------------------------------------------
    //STEP 1: Differentiate dynamics and cost along new trajectory

    auto start = std::chrono::high_resolution_clock::now();

    if (flgChange) {
      #pragma omp parallel for
      for (int t = 0; t < T; ++t) {
        model->dynamics_fo(xs.col(t), us.col(t), dt, fx[t], fu[t]);
      }
      flgChange = 0;
    }

    auto now = std::chrono::high_resolution_clock::now();
    trace(1, iter) = std::chrono::duration<double>(now - start).count();
    t_differentiate += trace(1, iter);

    //--------------------------------------------------------------------------
    // STEP 2: Backward pass, compute optimal control law and cost-to-go

    start = std::chrono::high_resolution_clock::now();

    bool backPassDone = false;
    while (!backPassDone) {
      diverge = backward_pass();

      if (diverge != 0) {
        dlambda = std::max(dlambda * lambdaFactor, lambdaFactor);
        lambda  = std::max(lambda * dlambda, lambdaMin);
        if (lambda > lambdaMax) break;
        continue;
      }
      backPassDone = true;
    }

    // check for termination due to small gradient
    g_norm = ((k.cwiseAbs().array() / (us.cwiseAbs().array() + 1.0)).colwise().maxCoeff()).mean();
    
    if (g_norm < tolGrad && lambda < 1e-5) {
      dlambda = std::min(dlambda / lambdaFactor, 1/lambdaFactor);
      lambda = lambda * dlambda * (lambda > lambdaMin);
      break;
    }

    now = std::chrono::high_resolution_clock::now();
    trace(2, iter) = std::chrono::duration<double>(now - start).count();
    t_backward += trace(2, iter);

    //--------------------------------------------------------------------------
    // STEP 3: Forward pass / line-search to find new control sequence, trajectory, cost

    start = std::chrono::high_resolution_clock::now();

    bool fwdPassDone = false;

    // serial backtracking line-search
    if (backPassDone) { 
      for (size_t i = 0; i < Alpha.size(); ++i) {
        alpha = Alpha[i];
        
        u_ff = us + alpha*k;
        forward_pass(u_ff);
        dcost    = cost_s - cost_new;
        expected = -alpha * (dV(0) + alpha*dV(1));

        if (expected > 0) {
          z = dcost/expected;
        }
        else {
          z = sgn(dcost); // 使用定义的 sgn 函数
          std::cout << "Warning: non-positive expected reduction. This should not occur" << std::endl; // 使用 std::cout
        }

        if(z > zMin) {
          fwdPassDone = true;
          break;
        }
      }

      if (!fwdPassDone) {
        // std::cout << "Forward pass failed" << std::endl;
        alpha = 0.0; // signals failure of forward pass
      }
    }

    now = std::chrono::high_resolution_clock::now();
    trace(3, iter) = std::chrono::duration<double>(now - start).count();
    t_forward += trace(3, iter);

    //--------------------------------------------------------------------------
    // STEP 4: accept step (or not), print status
    // 输出迭代信息
    // std::cout << "Iter: " << iter << " Cost: " << cost_s << " reduction: " << dcost << " expected: " << expected << " gradient: " << g_norm << " log10(lambda): " << log10(lambda) << std::endl;

    if (fwdPassDone) {

      // decrease lambda
      dlambda   = std::min(dlambda / lambdaFactor, 1/lambdaFactor);
      lambda    = lambda * dlambda * (lambda > lambdaMin);

      // accept changes
      us = unew;
      xs = xnew;
      cost_s = cost_new;
      trace(0, iter) = cost_s;
      flgChange = true;

      //terminate?
      if (dcost < tolFun) {
        std::cout << "\nSUCCESS: cost change < tolFun\n";
        break;
      }
    }
    else { // no cost improvement
    
      // increase lambda
      dlambda  = std::max(dlambda * lambdaFactor, lambdaFactor);
      lambda   = std::max(lambda * dlambda, lambdaMin);

      // terminate?
      if (lambda > lambdaMax) {
        std::cout << "\nEXIT: lambda > lambdaMax\n";
        break;
      }
    }

    if (iter == maxIter) std::cout << "\nEXIT: Maximum iterations reached.\n";
  } // end top-level for-loop

  auto all_end = std::chrono::high_resolution_clock::now();
  double t_total = std::chrono::duration<double>(all_end - all_start).count();
  
  std::cout << "Total    time: " << t_total << std::endl;
  //std::cout << "compute deriv: " << t_differentiate << std::endl;
  //std::cout << "backward pass: " << t_backward << std::endl;
  //std::cout << "forward  pass: " << t_forward << std::endl;
  //std::cout << "other   stuff: " << t_total - (t_differentiate + t_backward + t_forward) << std::endl;
  std::cout << "Final   error: " << (xs.col(T) - xd).norm() << std::endl;
  

} //generate_trajectory

// forward-pass (rollout)
void iLQG::forward_pass(const Eigen::MatrixXd& u) {
  cost_new = 0;
  xnew.resize(model->x_dims, T + 1);
  unew.resize(model->u_dims, T);
  xnew.col(0) = x0;
  
  for (int t = 0; t < T; ++t) {
    unew.col(t) = u.col(t) + K[t] * (xnew.col(t) - xs.col(t)); //apply LQR control gains after first iteration
    unew.col(t) = unew.col(t).cwiseMin(model->u_max).cwiseMax(model->u_min); // clamp controls to limits
    cost_new += 0.5*((xnew.col(t) - xd).transpose() * Qk * (xnew.col(t) - xd) + (unew.col(t).transpose() * Rk * unew.col(t))).value();
    xnew.col(t + 1) = model->dynamics(xnew.col(t), unew.col(t), dt);
  }

  cost_new += 0.5*((xnew.col(T) - xd).transpose() * QT * (xnew.col(T) - xd)).value();
  return;
}

// Perform the Ricatti-Mayne backward pass
int iLQG::backward_pass() {

  dV.setZero();

  //cost-to-go at end
  Vx.col(T) = QT*(xs.col(T) - xd);
  Vxx[T] = QT;

  for (int i = T - 1; i >= 0; --i) { // back up from end of trajectory
    Qx  = Qk * (xs.col(i) - xd) + fx[i].transpose() * Vx.col(i+1);
    Qu  = Rk * us.col(i) + fu[i].transpose() * Vx.col(i+1);
    Qxx = Qk + fx[i].transpose() * Vxx[i+1] * fx[i];
    Qux = fu[i].transpose() * Vxx[i+1] * fx[i];
    Quu = Rk + fu[i].transpose() * Vxx[i+1] * fu[i];

    // Similar to equations 10a and 10b in [Tassa 2012]. Note that regularization is different
    Qux_reg = fu[i].transpose() * Vxx[i+1] * fx[i];
    QuuF = Rk + lambda*Eigen::MatrixXd::Identity(model->u_dims,model->u_dims) + fu[i].transpose() * Vxx[i+1] * fu[i];

    BoxQPSolution res = boxQP(QuuF, Qu, model->u_min - us.col(i), model->u_max - us.col(i), k.col(std::min(i+1,T-1)));
    
    if(res.result < 1) return i;

    k_i = res.x;
    Eigen::VectorXd free = res.free;

    K_i = Eigen::MatrixXd::Zero(model->u_dims, model->x_dims);
    if (free.any()) {
      // 使用 Hfree 进行计算
      Eigen::MatrixXd Qux_reg_free(int(free.sum()), model->x_dims);
      int row_i = 0;
      for(int i = 0; i < model->u_dims; ++i){
          if (free(i)) {
              Qux_reg_free.row(row_i++) = Qux_reg.row(i);
          }
      }
      Eigen::MatrixXd Lfree = -res.llt.solve(Qux_reg_free); // 根据新的 Hfree 调整

      // 仅处理自由变量
      row_i = 0;
      for(int i = 0; i < model->u_dims; ++i) {
        if(free(i)) {
          K_i.row(i) = Lfree.row(row_i++);
        }
      }
    }

    // Update cost-to-go approximation. Equations 11 in [Tassa 2012]
    dV(0) += k_i.transpose()*Qu;
    dV(1) += 0.5*k_i.transpose()*Quu*k_i;

    Vx.col(i)  = Qx  + K_i.transpose()*Quu*k_i + K_i.transpose()*Qu + Qux.transpose()*k_i;
    Vxx[i] = Qxx + K_i.transpose()*Quu*K_i + K_i.transpose()*Qux + Qux.transpose()*K_i;
    Vxx[i] = 0.5 * (Vxx[i] + Vxx[i].transpose());

    // save controls/gains
    k.col(i) = k_i;
    K[i] = K_i;
  }

  return 0; // no divergence
}

} // namespace OC