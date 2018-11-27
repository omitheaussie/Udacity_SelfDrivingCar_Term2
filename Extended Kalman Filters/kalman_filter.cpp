#include "kalman_filter.h"
#include <iostream>
#define _USE_MATH_DEFINES // for C++  
#include <cmath>
#define LIM1 0.7
#define LIM2 0.6

using Eigen::MatrixXd;
using Eigen::VectorXd;
using namespace std;

// Please note that the Eigen library does not initialize 
// VectorXd or MatrixXd objects with zeros upon creation.

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  /**
  TODO:
    * predict the state
  */

	/*if (x_.rows() != F_.cols()) {
		std::cout << "x_ and F_ matrices mismatch for multiplication" << endl;
		std::cout << x_ << endl;
		std::cout << F_.cols() << endl;
	}
	if (P_.cols() != F_.rows()) {
		std::cout << "P_ and F_ matrices mismatch for multiplication" << endl;
		std::cout << F_.rows() << endl;
		std::cout << P_.cols() << endl;
	}*/
	x_ = F_ * x_;
	MatrixXd Ft = F_.transpose();
	P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Kalman Filter equations
  */

	VectorXd z_pred = H_ * x_;
	VectorXd y = z - z_pred;
	KF(y);
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Extended Kalman Filter equations
  */
	
	float rho = sqrt(x_[0] * x_[0] + x_[1] * x_[1]);
	if (fabs(rho) < 0.001) rho = 0.001;

	float px = x_[0];
	if (fabs(px) < 0.001) px = 0.001;
	float phi = atan2(x_[1] , px);

	float rhodot = (px * x_[2] + x_[1] * x_[3]) / rho;
	
	VectorXd z_pred = VectorXd(3); //h(x')
	z_pred << rho, phi, rhodot;
	
	VectorXd y = z - z_pred;
	
	KF(y);
}

void KalmanFilter::KF(const VectorXd &y) {
	MatrixXd Ht = H_.transpose();
	MatrixXd S = H_ * P_ * Ht + R_;
	MatrixXd Si = S.inverse();
	MatrixXd K = P_ * Ht * Si;
	// New state
	float curx = x_[0];
	float cury = x_[1];
	float curvx = x_[2];
	float curvy = x_[3];


	x_ = x_ + (K * y);

	// Noise removal
	if (fabs(x_[0] - curx) > LIM1) {
		if (x_[0] < 0) x_[0] = curx - LIM1;
		if (x_[0] > 0) x_[0] = curx + LIM1;
	}
	if (fabs(x_[1] - cury) > LIM1) {
		if (x_[1] < 0) x_[1] = cury - LIM1;
		if (x_[1] > 0) x_[1] = cury + LIM1;
	}
	if (fabs(x_[3] - curvy) > LIM2) {
		if (x_[3] < 0) x_[3] = curvy - LIM2;
		if (x_[3] > 0) x_[3] = curvy + LIM2;
	}
	if (fabs(x_[2] - curvx) > LIM1) {
		if (x_[2] < 0) x_[2] = curvx - LIM1;
		if (x_[2] > 0) x_[2] = curvx + LIM1;
	}
	int x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(x_size, x_size);
	P_ = (I - K * H_) * P_;
}
