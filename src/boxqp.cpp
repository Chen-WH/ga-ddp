#include "OC/boxqp.hpp"

namespace OC {

BoxQPSolution boxQP(
    const Eigen::MatrixXd& H,
    const Eigen::VectorXd& g,
    const Eigen::VectorXd& lower,
    const Eigen::VectorXd& upper,
    const Eigen::VectorXd& x0
) {
    BoxQPSolution sol;
    int n = H.rows();
    Eigen::VectorXd clamped = Eigen::VectorXd::Zero(n);
    Eigen::VectorXd free = Eigen::VectorXd::Ones(n);
    double oldvalue = 0.0;
    sol.result = 0;
    double gnorm = 0.0;
    int nfactor = 0;
    Eigen::LLT<Eigen::MatrixXd> llt;

    auto clamp = [&](const Eigen::VectorXd& vec) -> Eigen::VectorXd {
        return vec.cwiseMax(lower).cwiseMin(upper);
    };
    // initial state
    Eigen::VectorXd x = clamp(x0);
    for(int i = 0; i < n; ++i){ // 检查是否有非有限值
        if(!std::isfinite(x(i))){
            x(i) = (lower(i) + upper(i)) / 2;
        }
    }

    // options
    constexpr double maxIter = 100;
    constexpr double minGrad = 1e-8;
    constexpr double minRelImprove = 1e-8;
    constexpr double stepDec = 0.6;
    constexpr double minStep = 1e-22;
    constexpr double Armijo = 0.1;

    // initial objective value
    double value = x.dot(g) + 0.5 * x.dot(H * x);
    Eigen::VectorXd grad;
    Eigen::VectorXd old_clamped;
    bool factorize = true;
    Eigen::MatrixXd Hfree;
    Eigen::VectorXd search = Eigen::VectorXd::Zero(n);
    Eigen::VectorXd grad_clamped;
    Eigen::VectorXd grad_free;
    Eigen::VectorXd x_free;
    Eigen::VectorXd p_free;
    Eigen::VectorXd xc;
    // main loop
    for(int iter = 0; iter < maxIter; ++iter){
        if(sol.result != 0){
            break;
        }

        // check relative improvement
        if(iter > 0 && (oldvalue - value) < minRelImprove * std::abs(oldvalue)){
            sol.result = 4;
            break;
        }
        oldvalue = value;

        // get gradient
        grad = g + H * x;

        // find clamped dimensions
        old_clamped = clamped;
        clamped.setZero();
        for (int i = 0; i < n; ++i){
            if( (x(i) <= lower(i) + 1e-12) && grad(i) > 0 ){
                clamped(i) = 1;
            }
            if( (x(i) >= upper(i) - 1e-12) && grad(i) < 0 ){
                clamped(i) = 1;
            }
        }
        free = Eigen::VectorXd::Ones(n) - clamped;

        // check for all clamped
        if(free.sum() == 0){
            sol.result = 6;
            break;
        }

        // factorize if clamped has changed
        // 使用 .array() 进行逐元素比较
        if ( (clamped.array() != old_clamped.array()).any() ) {
            factorize = true;
        }

        // 取自由变量索引
        std::vector<int> free_indices;
        for (int i = 0; i < n; ++i) {
            if (free(i)) {
                free_indices.push_back(i);
            }
        }

        if (factorize){  // factorize if clamped has changed
            Hfree.resize(free_indices.size(), free_indices.size());
            for (int i = 0; i < free_indices.size(); ++i) {
                for(int j = 0; j < free_indices.size(); ++j){
                    Hfree(i, j) = H(free_indices[i], free_indices[j]);
                }
            }

            // Cholesky 分解
            llt.compute(Hfree);
            if (llt.info() != Eigen::Success) {
                sol.result = -1;
                break;
            }
            nfactor +=1;
        }

        // 检查梯度范数
        gnorm = grad.cwiseProduct(free).norm();
        if(gnorm < minGrad){
            sol.result = 5;
            break;
        }

        // 获取搜索方向
        grad_clamped = g + H * (x.cwiseProduct(clamped));
        search.setZero();

        // 构建 grad_free
        grad_free.resize(free_indices.size());
        x_free.resize(free_indices.size());
        for(int i = 0; i < free_indices.size(); ++i){
            grad_free(i) = grad_clamped(free_indices[i]);
            x_free(i) = x(free_indices[i]);
        }
        // 计算搜索方向自由部分:
        p_free = -llt.solve(grad_free) - x_free;
        // 填充搜索方向
        for(int i = 0;i < free_indices.size(); ++i){
            search(free_indices[i]) = p_free(i);
        }

        // check for descent direction
        double sdotg = search.dot(grad);
        if(sdotg >= 0){
            break;
        }

        // Armijo 线搜索
        double step = 1.0;
        int nstep = 0;
        xc = clamp(x + step * search);
        double vc = xc.dot(g) + 0.5 * xc.dot(H * xc);
        while( (vc - oldvalue)/(step * sdotg) < Armijo ){
            step *= stepDec;
            nstep += 1;
            xc = clamp(x + step * search);
            vc = xc.dot(g) + 0.5 * xc.dot(H * xc);
            if(step < minStep){
                sol.result =2;
                break;
            }
        }

        // 接受候选解
        x = xc;
        value = vc;
    }

    sol.x = x;
    sol.llt = llt;
    sol.free = free;

    return sol;
}

} // namespace OC
