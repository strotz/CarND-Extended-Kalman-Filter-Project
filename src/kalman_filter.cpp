#include "kalman_filter.h"

#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::cout;
using std::endl;

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
  x_ = F_ * x_;
  P_ = F_ * P_ * F_.transpose() + Q_;
}

void KalmanFilter::UpdateByError(const Eigen::VectorXd &y) {
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd K = P_ * Ht * Si;

  //new estimate
  x_ = x_ + K * y;
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
}

void KalmanFilter::Update(const VectorXd &z) {
  VectorXd z_prediction = H_ * x_;
  VectorXd y = z - z_prediction;
  UpdateByError(y);
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  // h(x)
  double px = x_(0);
  double py = x_(1);
  double vx = x_(2);
  double vy = x_(3);

  double hx1 = sqrt(px * px + py * py);

  if (hx1 < 0.0001)
  {
    cout << "UpdateEKF() - Devision By Zero" << endl;
  }

  double hx2 = atan2(py, px);
  double hx3 = (px * vx + py * vy) / hx1;

  VectorXd hx = VectorXd(3);
  hx << hx1, hx2, hx3;

  VectorXd y = z - hx;
  UpdateByError(y);
}
