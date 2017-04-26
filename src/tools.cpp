#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;
using std::cout;
using std::endl;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd>& estimations,
                              const vector<VectorXd>& ground_truth) {

  VectorXd rmse(4);
  rmse << 0, 0, 0, 0;

  if (estimations.size() == 0) {
    //cout << "CalculateRMSE () - Error - the estimation vector size should not be zero" << endl;
    return rmse;
  }

  if (estimations.size() != ground_truth.size()) {
    //cout << "CalculateRMSE () - Error - the estimation vector size should equal ground truth vector size" << endl;
    return rmse;
  }

  //accumulate squared residuals
  for (int i = 0; i < estimations.size(); ++i) {
    VectorXd temp = estimations[i] - ground_truth[i];
    rmse = rmse.array() + temp.array() * temp.array();
  }

  //calculate the mean
  rmse = rmse / estimations.size();

  //calculate the squared root
  rmse = rmse.array().sqrt();

  //return the result
  return rmse;

}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {

  MatrixXd Hj(3, 4);

  //recover state parameters
  float px = x_state(0);
  float py = x_state(1);
  float vx = x_state(2);
  float vy = x_state(3);

  //compute the Jacobian matrix
  float d = px * px + py * py;
  float d_1_2 = sqrt(d);
  float d_3_2= d * d_1_2;

  //check division by zero
  if (fabs(d) < 0.0001) {
    // cout << "CalculateJacobian () - Error - Division by Zero" << endl;
    return Hj;
  }

  Hj << px/d_1_2, py/d_1_2, 0, 0,
      -py/d, px / d, 0, 0,
      py*(vx*py - vy*px)/d_3_2, px*(vy*px - vx*py)/d_3_2, px/d_1_2, py/d_1_2;

  return Hj;
}
