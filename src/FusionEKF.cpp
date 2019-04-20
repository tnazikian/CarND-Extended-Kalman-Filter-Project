#include "FusionEKF.h"
#include <iostream>
#include "Eigen/Dense"
#include "tools.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::cout;
using std::endl;
using std::vector;

/**
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

  /**
   * TODO: Finish initializing the FusionEKF.
   * TODO: Set the process and measurement noises
   */ 
  
  ekf_.P_ = MatrixXd(4, 4);
  ekf_.P_ << 1, 0, 0, 0,
            0, 1, 0, 0, 
            0, 0, 1000, 0,
            0, 0, 0, 1000;
  
  ekf_.F_ = MatrixXd(4, 4);
  ekf_.F_ << 1, 0, 1, 0,
             0, 1, 0, 1,
             0, 0, 1, 0,
             0, 0, 0, 1;
    
  H_laser_ << 1, 0, 0, 0,
              0, 1, 0, 0;
    
  ekf_.H_ = MatrixXd(4, 4);
  ekf_.H_ << 1, 0, 0, 0,
             0, 1, 0, 0,
             0, 0, 1, 0,
             0, 0, 0, 1;
}

/**
 * Destructor.
 */
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
  /**
   * Initialization
   */
  if (!is_initialized_) {
    /**
     * TODO: Initialize the state ekf_.x_ with the first measurement.
     * TODO: Create the covariance matrix.
     * You'll need to convert radar from polar to cartesian coordinates.
     */
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 1, 1;
    //initialize covm
    ekf_.Q_ = MatrixXd(4, 4);
    float px;
    float py;
    
    //check if radar or lidar
    //if radar, then get cartesian coords from polar measurements
    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      float rho = measurement_pack.raw_measurements_[0];
      float phi = measurement_pack.raw_measurements_[1];
      px = rho * cos(phi);
      py = rho * sin(phi);
    }
    //laser measurements
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      px = measurement_pack.raw_measurements_[0];
      py = measurement_pack.raw_measurements_[1];
    }
    //Check if values for px or py too small
    if (fabs(px) < 0.0001) {
      px = px * 1000;
    }
    
    if (fabs(py) < 0.0001) {
      py = py * 1000;
    }
    
    ekf_.x_ << px, py, 0, 0;
    previous_timestamp_ = measurement_pack.timestamp_;
    is_imitialized_ = true;
    return;
  }

  /**
   * Prediction
   */

  /**
   * TODO: Update the state transition matrix F according to the new elapsed time.
   * Time is measured in seconds.
   * TODO: Update the process noise covariance matrix.
   * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */

  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
  previous_timestamp_ = measurement_pack.timestamp_;
  ekf_.F_(0, 2) = dt;
  ekf_.F_(1, 3) = dt;
  
  float noise_ax = 9;
  float noise_ay = 9;
  
  float dt_2 = dt * dt;
  float dt_3 = dt_2 * dt;
  float dt_4 = dt_3 * dt;

  // set the process covariance matrix Q
  ekf_.Q_ = MatrixXd(4, 4);
  ekf_.Q_ <<  dt_4/4*noise_ax, 0, dt_3/2*noise_ax, 0,
              0, dt_4/4*noise_ay, 0, dt_3/2*noise_ay,
              dt_3/2*noise_ax, 0, dt_2*noise_ax, 0,
              0, dt_3/2*noise_ay, 0, dt_2*noise_ay;

  ekf_.Predict();

  /**
   * Update
   */

  /**
   * TODO:
   * - Use the sensor type to perform the update step.
   * - Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // TODO: Radar updates
    //convert polar measurements to cartesian
    float phi;
    float rho;
    float rhodot;
    
    float px = ekf_.x_[0];
    float py = ekf_.x_[1];
    float vx = ekf_.x_[2];
    float vy = ekf_.x_[3];
    
    if (fabs(px) < 0.0001 or fabs(py) < 0.0001) {
      if (fabs(px) < 0.0001) {
        px = px * 1000;
      }
      if (fabs(py) < 0.0001) {
        py = py * 1000;
      }
      rho = sqrt(px*px + py*py);
      //phi = atan2(py, px);
      phi = 0;
      //rhodot = (px*vx + py*vy) /rho;
      rhodot = 0;
    }
    else {
      rho = sqrt(px*px + py*py);
      phi = atan2(py, px); 
      rhodot = (px*vx + py*vy) / rho;
    }
    
    //calculate jacobian and store polar coordinates for update
    ekf_.hx_ << rho, phi, rhodot;
    Hj_ = tools.CalculateJacobian(ekf_.x_);
    
    //check jacobian
    if (Hj_.isZero(0)) {
      cout << "Jacobian zero";
      return;
    }
    ekf_.H_ = Hj;
    ekf_.R_ = R_radar_;
    //Update using extended kalman filter
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);   

  } else {
    // TODO: Laser updates
    ekf_.H_ = H_laser_;
    ekf_.R_ = R_laser_;
    ekf_.Update(measurement_pack.raw_measurements_);
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}