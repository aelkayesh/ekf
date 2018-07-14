#include <iostream>
#include "tools.h"
#include <math.h>

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {

    int size = estimations.size();
    VectorXd sum(4);
    VectorXd rmse(4);
    rmse << 0 ,0 ,0 ,0;
	sum << 0, 0, 0, 0;

    if(size != ground_truth.size() || size == 0){
      cerr << "Error: Estimation and ground truth vectors sizes are not equal or zero.";
    } else{
        for (int i = 0; i < size; ++i) {
            VectorXd residual = estimations[i] - ground_truth[i];
            residual = residual.array().pow(2);
            sum += residual;
        }
        sum /= size;
        rmse << sum.array().sqrt();
    }
    return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {

  MatrixXd jacobian (3,4);
  jacobian <<   0,0,0,0,
                0,0,0,0,
                0,0,0,0;
  double px = x_state[0];
  double py = x_state[1];
  double vx = x_state[2];
  double vy = x_state[3];



  if (px == 0 && py == 0){
      cout << "Error: jacobian() - division by zero";
  }else{
      double denom1 = pow(px,2) + pow(py,2);
      double denom2 = sqrt(denom1);
      double denom3 = denom1 * denom2;

      jacobian <<   px/denom2,                  py/denom2,                  0 ,         0,
                    -py/denom1,                 px/denom1,                  0,          0,
                    py*((vx*py)-(vy*px))/denom3,  px*((vy*px)-(vx*py))/denom3,  px/denom2,  py/denom2;
  }
  return jacobian;
}

