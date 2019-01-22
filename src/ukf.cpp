#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>
#include <math.h>

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
  std_a_ = 2;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.3;

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
   time_us_ = 0.0;
   n_x_ = 5;
   n_aug_ = 7;
   Xsig_pred_ = MatrixXd(n_x_, 2*n_aug_+1);
   weights_ = VectorXd(2*n_aug_+1);
   lambda_ = 3-n_aug_;

   // Set weights
   double weight_0 = lambda_/(lambda_+n_aug_);
   weights_(0) = weight_0;
   for (int i=1; i<2*n_aug_+1; ++i) {  // 2n+1 weights
     double weight = 0.5/(n_aug_+lambda_);
     weights_(i) = weight;
   }
}

UKF::~UKF() {}

void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Make sure you switch between lidar and radar
   * measurements.
   */

   // Should the current sensor measurement be ignored?
   if( (use_radar_ && meas_package.sensor_type_ == MeasurementPackage::RADAR) ||
       (use_laser_ && meas_package.sensor_type_ == MeasurementPackage::LASER )) {

     if(!is_initialized_) {
       P_ << MatrixXd::Identity(n_x_,n_x_);
       P_(0,0) = 0.15;
       P_(1,1) = 0.15;
       x_ << 1,1,1,1,0.1;
       time_us_= meas_package.timestamp_;

       if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
         double rho = meas_package.raw_measurements_[0];
         double phi = meas_package.raw_measurements_[1];
         x_(0) = rho * cos(phi);
         x_(1) = rho * sin(phi);
       }
       else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
         x_(0) = meas_package.raw_measurements_[0];
         x_(1) = meas_package.raw_measurements_[1];
       }
       is_initialized_ = true;
       return;
     }

     // Predict the state
     float dt = (meas_package.timestamp_-time_us_)/1000000.0;
     time_us_ = meas_package.timestamp_;
     Prediction(dt);

     // Update the state
     if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
       UpdateRadar(meas_package);
     }
     else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
       UpdateLidar(meas_package);
     }
  }
}

void UKF::Prediction(double delta_t) {
  /**
   * TODO: Complete this function! Estimate the object's location.
   * Modify the state vector, x_. Predict sigma points, the state,
   * and the state covariance matrix.
   */

   // Augmented Sigma Matrix
   MatrixXd Xsig_aug = MatrixXd(n_aug_, 2*n_aug_+1);

   // Augmented Covariance Matrix
   MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
   P_aug.fill(0.0);
   P_aug.topLeftCorner(n_x_, n_x_) = P_;
   P_aug(5,5) = std_a_*std_a_;
   P_aug(6,6) = std_yawdd_*std_yawdd_;

   // Square root matrix A
   MatrixXd A = P_aug.llt().matrixL();

   // Augmented state x
   VectorXd x_aug = VectorXd(n_aug_);
   x_aug.head(n_x_) = x_;
   x_aug(5) = 0;
   x_aug(6) = 0;

   // Generate Sigma Points
   Xsig_aug.col(0) = x_aug;
   for (unsigned int i=0; i < n_aug_; ++i) {
     Xsig_aug.col(i+1) = x_aug+sqrt(lambda_+n_aug_)*A.col(i);
     Xsig_aug.col(i+n_aug_+1) = x_aug-sqrt(lambda_+n_aug_)*A.col(i);
   }
   // Predict Sigma Points
   for (unsigned int i=0; i<2*n_aug_+1; ++i) {

     double px = Xsig_aug(0,i);
     double py = Xsig_aug(1,i);
     double v = Xsig_aug(2,i);
     double yaw = Xsig_aug(3,i);
     double yaw_rate = Xsig_aug(4,i);
     double acc_noise = Xsig_aug(5,i);
     double yaw_noise = Xsig_aug(6,i);

     VectorXd predicted_state = VectorXd(5);
     predicted_state << px, py, v, yaw, yaw_rate;
     double delta_t_sq = delta_t*delta_t;

     predicted_state(3) += yaw_rate*delta_t;

     if (fabs(yaw_rate) < 0.001) {
       predicted_state(0) += v*cos(yaw)*delta_t;
       predicted_state(1) += v*sin(yaw)*delta_t;
     }
     else {
       predicted_state(0) += (v/yaw_rate)*(sin(yaw+yaw_rate*delta_t)-sin(yaw));
       predicted_state(1) += (v/yaw_rate)*(-cos(yaw+yaw_rate*delta_t)+cos(yaw));
     }

     // Add noise terms separately
     predicted_state(0) += (0.5)*(delta_t_sq)*cos(yaw)*acc_noise;
     predicted_state(1) += (0.5)*(delta_t_sq)*sin(yaw)*acc_noise;
     predicted_state(2) += delta_t*acc_noise;
     predicted_state(3) += (0.5)*(delta_t_sq)*yaw_noise;
     predicted_state(4) += delta_t*yaw_noise;
     Xsig_pred_.col(i) = predicted_state;
   }

   // Calculate new mean and covariance
   VectorXd x_new = VectorXd(n_x_);
   x_new.fill(0.0);
   MatrixXd P_new = MatrixXd(n_x_, n_x_);
   P_new.fill(0.0);

   for (unsigned int i=0; i<2*n_aug_+1; ++i) {
     x_new += weights_(i)*Xsig_pred_.col(i);
   }

   for (unsigned int i=0; i<2*n_aug_+1; ++i) {
     VectorXd x_diff = Xsig_pred_.col(i)-x_new;
     while(x_diff(3)>M_PI) x_diff(3)-=2.*M_PI;
     while(x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;
     P_new += weights_(i)*x_diff*x_diff.transpose();
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
   // Projection Matrix
   MatrixXd H = MatrixXd(2, n_x_);
   H << 1,0,0,0,0,
        0,1,0,0,0;

   // Measurement covariance matrix
   MatrixXd R = MatrixXd(2, 2);
   R << std_laspx_*std_laspx_,0,
        0,std_laspy_*std_laspy_;

   VectorXd y = meas_package.raw_measurements_ - H*x_;
   MatrixXd S = H*P_*H.transpose()+R;
   MatrixXd K = P_*H.transpose()*S.inverse();
   MatrixXd I = MatrixXd::Identity(5,5);

   x_ = x_ + K*y;
   P_ = (I-K*H)*P_;
}

void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use radar data to update the belief
   * about the object's position. Modify the state vector, x_, and
   * covariance, P_.
   * You can also calculate the radar NIS, if desired.
   */

   // Define measurement covariance
   MatrixXd R = MatrixXd(3,3);
   R << std_radr_*std_radr_,0,0,
        0,std_radphi_*std_radphi_,0,
        0,0,std_radrd_*std_radrd_;

   // Project sigma points into measurement space
   MatrixXd Zsig = MatrixXd(3, 2*n_aug_+1);
   for (unsigned int i=0; i < 2*n_aug_+1; ++i) {
     // Strap into variables for easier calculation
     double px = Xsig_pred_(0,i);
     double py = Xsig_pred_(1,i);
     double v = Xsig_pred_(2,i);
     double yaw = Xsig_pred_(3,i);
     // calculate vx and vy
     double vx = cos(yaw)*v;
     double vy = sin(yaw)*v;

     Zsig(0,i) = sqrt(px*px+py*py);
     Zsig(1,i) = atan2(py,px);
     Zsig(2,i) = (px*vx+py*vy)/sqrt(px*px+py*py);
   }

   // Calculate new mean and covariancce
   VectorXd z = VectorXd(3);
   z.fill(0.0);
   MatrixXd S = MatrixXd(3,3);
   S.fill(0.0);

   for (unsigned int i=0; i<2*n_aug_+1; ++i) {
     z += weights_(i)*Zsig.col(i);
   }

   for (unsigned int i=0; i<2*n_aug_+1; ++i) {
     VectorXd z_diff = Zsig.col(i)-z;
     while(z_diff(1)>M_PI) z_diff(1)-=2.*M_PI;
     while(z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;
     S += weights_(i)*z_diff*z_diff.transpose();
   }
   S += R;

   // Calculate cross-correlation matrix T
   MatrixXd T = MatrixXd(5,3);
   T.fill(0.0);
   for (unsigned int i=0; i<2*n_aug_+1; ++i) {
     VectorXd z_diff = Zsig.col(i)-z;
     while(z_diff(1)>M_PI) z_diff(1)-=2.*M_PI;
     while(z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

     VectorXd x_diff = Xsig_pred_.col(i)-x_;
     while(x_diff(3)>M_PI) x_diff(3)-=2.*M_PI;
     while(x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

     T+= weights_(i)*x_diff*z_diff.transpose();
   }

   // Calculate Kalman gain
   MatrixXd K = T*S.inverse();
   // Make predictions
   VectorXd z_diff = meas_package.raw_measurements_-z;
   while(z_diff(1)>M_PI) z_diff(1)-=2.*M_PI;
   while(z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;
   x_ = x_+K*z_diff;
   P_ = P_-K*S*K.transpose();
}
