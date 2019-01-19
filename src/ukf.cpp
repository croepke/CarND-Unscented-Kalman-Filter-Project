#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 30;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 30;

  /**
   * DO NOT MODIFY measurement noise values below.
   * These are provided by the sensor manufacturer.
   */

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;

  /**
   * End DO NOT MODIFY section for measurement noise values
   */

  /**
   * TODO: Complete the initialization. See ukf.h for other member properties.
   * Hint: one or more values initialized above might be wildly off...
   */

   is_initialized_ = false;
   n_x_ = 5;
   n_aug_ = 7;
   Xsig_pred_ = MatrixXd(n_aug_, n_aug_);
   time_us_= 0;
   weights_ = VectorXd(2*n_aug_+1);
   lambda_ = 3-n_aug_;
}

UKF::~UKF() {}

void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Make sure you switch between lidar and radar
   * measurements.
   */
}

void UKF::Prediction(double delta_t) {
  /**
   * TODO: Complete this function! Estimate the object's location.
   * Modify the state vector, x_. Predict sigma points, the state,
   * and the state covariance matrix.
   */

   // Augmented Covariance Matrix
   MatrixXd P_aug_ = MatrixXd(n_aug_, n_aug_);
   P_aug_.topLeftCorner(n_x_, n_x_) = P_;
   P_aug_(5,5) = std_a_*std_a_;
   P_aug_(6,6) = std_yawdd_*std_yawdd_;

   // Square root matrix A
   MatrixXd A_ = P_aug_.llt().matrixL();

   // Augmented state x
   VectorXd x_aug_ = VectorXd(n_aug_);
   x_aug_.topRows(n_x_) = x_;

   // Generate Sigma Points
   Xsig_pred_.col(0) = x_;
   for (unsigned int i=1; i < n_aug_+1; ++i) {
     Xsig_pred_.col(i) = x_aug_+sqrt(lambda_+n_aug_)*A_;
     Xsig_pred_.col(i+n_aug_) = x_aug_-sqrt(lambda_+n_aug_)*A_;
   }

   // Predict Sigma Points into Measurement Space
   for (unsigned int i=0; i<2*n_aug_+1; ++i) {
     float v = Xsig_pred_(2, i);
     float yaw = Xsig_pred_(3, i);
     float yaw_rate = Xsig_pred_(4, i);
     float acc_noise = Xsig_pred_(5, i);
     float yaw_noise = Xsig_pred_(6, i);

     VectorXd predicted_state = VectorXd(5);
     float delta_t_sq = delta_t*delta_t;

     predicted_state(2) = 0;
     predicted_state(3) = yaw_rate*delta_t;
     predicted_state(4) = 0;

     if (abs(yaw_rate) < 0.0001) {
       predicted_state(0) = v*cos(yaw)*delta_t;
       predicted_state(1) = v*sin(yaw)*delta_t;
     }
     else {
       predicted_state(0) = (v/yaw_rate)*(sin(yaw+yaw_rate*delta_t)-sin(yaw));
       predicted_state(1) = (v/yaw_rate)*(-cos(yaw+yaw_rate*delta_t)+cos(yaw));
     }

     // Add noise terms separately
     predicted_state(0) += (0.5)*(delta_t_sq)*cos(yaw)*acc_noise;
     predicted_state(1) += (0.5)*(delta_t_sq)*sin(yaw)*acc_noise;
     predicted_state(2) += delta_t*acc_noise;
     predicted_state(3) += (0.5)*(delta_t_sq)*yaw_noise;
     predicted_state(4) += delta_t*yaw_noise;
   }

   // Calculate new mean and covariance
   VectorXd x_new = VectorXd(n_x_);
   MatrixXd P_new = MatrixXd(n_x_, n_x_);
   VectorXd weights = VectorXd(2*n_aug_+1);
   weights(0) = (lambda_/(lambda_+n_aug_));
   weights.middleRows(1, n_x_-1).fill(1/(2*(lambda_+n_aug_)));

   for (unsigned int i=0; i<2*n_aug_+1; ++i) {
     x_new += weights(i)*Xsig_pred_.col(i);
   }

   for (unsigned int i=0; i<2*n_aug_+1; ++i) {
     P_new += weights(i)*(Xsig_pred_.col(i)-x_new)*(Xsig_pred_.col(i)-x_new).transpose();
   }

   P_ = P_new;
   x_ = x_new;
}

void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use lidar data to update the belief
   * about the object's position. Modify the state vector, x_, and
   * covariance, P_.
   * You can also calculate the lidar NIS, if desired.
   */
}

void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use radar data to update the belief
   * about the object's position. Modify the state vector, x_, and
   * covariance, P_.
   * You can also calculate the radar NIS, if desired.
   */
}
