#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
      0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
      0, 0.0009, 0,
      0, 0, 0.09;

  //TODO:  * Finish initializing the FusionEKF.

  ekf_.x_ = VectorXd(4); // state vector

  // TODO: is it correct initialization values?
  // TODO: * Create the covariance matrix.
  ekf_.P_ = MatrixXd(4,4); // state covariance matrix P
  ekf_.P_ << 1, 0, 0, 0,
      0, 1, 0, 0,
      0, 0, 1000, 0,
      0, 0, 0, 1000;

 // TODO:  * Set the process and measurement noises

  //set the acceleration noise components
  noise_ax = 9;
  noise_ay = 9;
}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

/*****************************************************************************
 *  Initialization
 ****************************************************************************/
void FusionEKF::init(const MeasurementPackage &measurement_pack) {

  // first measurement
  cout << "EKF: " << endl;
  ekf_.x_ << 1, 1, 1, 1;

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {

    float ro = measurement_pack.raw_measurements_(0);
    float phi = measurement_pack.raw_measurements_(1);
    ekf_.x_ << ro * cos(phi), ro * sin(phi), 0, 0;

  } else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {

    ekf_.x_ << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0, 0;
  }

  previous_timestamp_ = measurement_pack.timestamp_;
}

/*****************************************************************************
 *  Prediction
 ****************************************************************************/
void FusionEKF::predict(const MeasurementPackage &measurement_pack) {

  //compute the time elapsed between the current and previous measurements
  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;  //dt - expressed in seconds
  previous_timestamp_ = measurement_pack.timestamp_;

  ekf_.F_ << 1, 0, dt, 0,
      0, 1, 0, dt,
      0, 0, 1, 0,
      0, 0, 0, 1;

  //2. Set the process covariance matrix Q
  float dt2 = dt * dt;
  float dt3 = dt2 * dt / 2;
  float dt4 = dt3 * dt / 2;

  ekf_.Q_ = MatrixXd(4, 4);
  ekf_.Q_ << dt4 * noise_ax, 0, dt3 * noise_ax, 0,
      0, dt4 * noise_ay, 0, dt3 * noise_ay,
      dt3 * noise_ax, 0, dt2 * noise_ax, 0,
      0, dt3 * noise_ay, 0, dt2 * noise_ay;


  ekf_.Predict();
}


void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {

  if (!is_initialized_) {
    init(measurement_pack);

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  predict(measurement_pack);

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // TODO:  * Use the sensor type to perform the update step.
    // TODO:  * Update the state and covariance matrices.
  } else {
    // Laser updates
    ekf_.Update(measurement_pack.raw_measurements_);
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
