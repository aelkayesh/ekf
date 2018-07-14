#include "kalman_filter.h"
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using namespace std;

// Please note that the Eigen library does not initialize 
// VectorXd or MatrixXd objects with zeros upon creation.

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}


void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in; //Object state 4X1
  P_ = P_in; //Object Covariance 4X1
  F_ = F_in; //State transition function 4X4
  H_ = H_in; //Measurement update function --> mapping predicted values to measurement
  R_ = R_in; //Measurement covariance matrix 
  Q_ = Q_in; //Process covariance matrix 4X4
}

void KalmanFilter::Predict() {
	x_ = F_ * x_ ;
	MatrixXd Ft = F_.transpose();
	P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  VectorXd y = z - H_ * x_;
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_* P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd K = P_ * Ht * Si;

  //new state
  x_ = x_ + (K *y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = ( I - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
	VectorXd zPred(3); //location Polar
	zPred(0) = sqrt((x_[0] * x_[0]) + (x_[1] * x_[1]));
	zPred(1) = atan2(x_[1], x_[0]);
    zPred(2) = ((x_[0] * x_[2]) + (x_[1] * x_[3])) / zPred(0);


	VectorXd y = z - zPred;
    //normalize phi
    if(y(1)< - M_PI ){
        y(1) += 2 * M_PI;
    }else if (y(1) >  M_PI ){
        y(1) -= 2 * M_PI;
    }
	//H is Hj
	MatrixXd Ht = H_.transpose();
	MatrixXd S = H_ * P_ * Ht + R_;
	MatrixXd Si = S.inverse();
	MatrixXd K = P_ * Ht * Si;

	//new state
	x_ = x_ + (K *y);
	long x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(x_size, x_size);
	P_ = (I - K * H_) * P_;
}


