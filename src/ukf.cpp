#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
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

  //DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
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
  //DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.

  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */

  is_initialized_ = false;
  time_us_ = 0;

  ///* State dimension
  n_x_ = 5;

  ///* Augmented state dimension
  n_aug_ = 7;

  ///* Sigma point spreading parameter
  lambda_ = 3 - n_aug_;

  ///* Augmented sigma point prediction matrix
  Xsig_pred_ = MatrixXd(n_x_, 2*n_aug_+1);

  ///* Weights of sigma points
  for (int i = 1; i < n_aug_*2+1; i++) {
    double weight = 1 / (2 * (lambda_ + n_aug_));
    weights_(i) = weight;
  }

}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */
  if (!is_initialized_) {

    if (meas_package.sensor_type == MeasurementPackage::RADAR) {
      float rho = meas_package.raw_measurements_[0];
      float phi = meas_package.raw_measurements_[1];
      float rho_dot = meas_package.raw_measurements_[2];

      // convert to cartesian coordinates
      x_ << rho*cos(phi), rho*sin(phi), 0, 0, 0;

    } else if (meas_package.sensor_type == MeasurementPackage::LIDAR) {
      x_ << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], 0, 0, 0;
    }

    time_us_ = meas_package.timestamp_;
    is_initialized_ = true;

    P_ = P_.setIdentity();

    return;
  }

  /***************
    * Prediction
  /***************/
  float delta_t = (meas_package.timestamp_ - time_us_) / 1000000.0;
  Prediction(delta_t);


}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */

  // * Generate Sigma Points
  MatrixXd Q = MatrixXd(2, 2); // Process Noise (acceleration & turning acceleration)
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2*n_aug_+1);

  VectorXd x_aug = VectorXd(7)
  x_aug.head(5) = x_;

  Q << std_a_*std_a_, 0
       0, std_yawdd_*std_yawdd_;

  P_aug.topLeftCorner(n_x_, n_x_) = P_;
  P_aug.bottomRightCorner(2, 2) = Q;

  // Create Square root matrix
  MatrixXd A = P_aug.llt().matrixL();

  // Create Augmented sigma Points
  // 1st col = mean
  x_aug(5) = 0.0;
  x_aug(6) = 0.0;
  Xsig_aug.col(0) = x_aug;

  for (int i = 1; i < n_aug_; i++) {
    Xsig_aug.col(i) = x_aug + sqrt(lambda_ + n_aug_) * A.col(i-1);
    Xsig_aug.col(i+7) = x_aug - sqrt(lambda_ + n_aug_) * A.col(i-1);
  }

  // * Predict Sigma Points
  Xsig_pred_ = Xsig_aug.block(0,0,5,15);

  for (int i = 0; i < (2*n_aug_+1); i++) {
    double px = Xsig_aug(0, i);
    double py = Xsig_aug(1, i);
    double vel = Xsig_aug(2, i);
    double rho = Xsig_aug(3, i);
    double rho_dot = Xsig_aug(4, i);

    double v_noise = Xsig_aug(5, i);
    double r_noise = Xsig_aug(6, i);

    if (rho_dot > 0.001) {
      Xsig_pred_(0, i) = px + (vel/rho_dot)*(sin(rho + rho_dot*delta_t) - sin(rho)) + 0.5*(delta_t*delta_t)*cos(rho)*v_noise;
      Xsig_pred_(1, i) = py + (vel/rho_dot)*(-cos(rho + rho_dot*delta_t) + cos(rho)) + 0.5*(delta_t*delta_t)*sin(rho)*v_noise;
      Xsig_pred_(2, i) = vel + 0 + delta_t*v_noise;
      Xsig_pred_(3, i) = rho + rho_dot*delta_t + 0.5*(delta_t*delta_t)*r_noise;
      Xsig_pred_(4, i) = rho_dot + 0 + delta_t*r_noise;
    } else {
      Xsig_pred_(0, i) = px + vel*cos(rho)*delta_t + 0.5*(delta_t*delta_t)*cos(rho)*v_noise;
      Xsig_pred_(1, i) = py + vel*sin(rad)*delta_t + 0.5*(delta_t*delta_t)*sin(rho)*v_noise;
      Xsig_pred_(2, i) = vel + 0 + delta_t*v_noise;
      Xsig_pred_(3, i) = rho + rho*delta_t + 0.5*(delta_t*delta_t)*r_noise;
      Xsig_pred_(4, i) = rho_dot + 0 + delta_t*r_noise;
    }
  }

  // * Predict Mean & Covariance

  //predict state mean
  for (int i = 0; i < n_aug_*2+1; i++) {
    x_ = x_ + weights_(i)*Xsig_pred_.col(i)
  }

  //predict state covariance
  for (int i = 0; i < n_aug_*2+1; i++) {
    VectorXd x_diff = Xsig_pred_.col(i) - x;

    // Normalize the angle -M_PI < angle < M_PI
    while (x_diff(3) > M_PI) x_diff(3) -= 2*M_PI;
    while (x_diff(3) < -M_PI) x_diff(3) += 2*M_PI;

    P_ = P_ + weights_(i) * x_diff * x_diff.transpose();
  }
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */

  //measurement dimension (r, phi, r_dot)
  int n_z = 3;

  //measurement covariance matrix
  MatrixXd S = MatrixXd(n_z, n_z);

  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2*n_aug_+1);

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);

  Zsig.fill(0.0);
  for (int i = 0; i < n_aug_*2+1; i++) {
    VectorXd s_pt = Xsig_pred_.col(i);
    double px = s_pt(0);
    double py = s_pt(1);
    double v = s_pt(2);
    double yaw = s_pt(3);
    double yaw_dot = s_pt(4);

    double v1 = cos(yaw)*v;
    double v2 = sin(yaw)*v;

    Zsig(0, i) = sqrt((px*px) + (py*py));
    Zsig(1, i) = atan2(py, px);
    Zsig(2, i) = (px*v1 + py*v2) / sqrt((px*px) + (py*py));
  }

  //calculate mean predicted measurement
  z_pred.fill(0.0);
  for (int i = 0; i < n_aug_*2+1; i++) {
    z_pred = z_pred + (weights_(i) * Zsig.col(i));
  }

  //calculate innovation convariance matrix S
  MatrixXd R = MatrixXd(n_z, n_z);
  R.fill(0.0);
  R(0, 0) = std_radr_ * std_radr_;
  R(1, 1) = std_radphi_ * std_radphi_;
  R(2, 2) = std_radrd_ * std_radrd_;

  S.fill(0.0);
  for (int i = 0; i < n_aug_*2+1; i++) {
    VectorXd Z_diff = Zsig.col(i) - z_pred;

    while (Z_diff(1) > M_PI) Z_diff(1) -= 2*M_PI;
    while (Z_diff(1) < -M_PI) Z_diff(1) += 2*M_PI;

    S = S + (weights_(i) * Z_diff * Z_diff.transpose());
  }

  S = S + R;

  //create matrix for cross correlation (Tc)
  MatrixXd Tc = MatrixXd(n_x_, n_z);

  Tc.fill(0.0);
  for (int i = 0; i < n_aug_*2+1; i++) {
    VectorXd X_diff = Xsig_pred_.col(i) - x_;

    // Angle Normalization
    while (X_diff(3) > M_PI) X_diff(3) -= 2*M_PI;
    while (X_diff(3) < -M_PI) X_diff(3) += 2*M_PI;

    VectorXd Z_diff = Zsig.col(i) - z_pred;

    // Angle Normalization
    while (Z_diff(1) > M_PI) Z_diff(1) -= 2*M_PI;
    while (Z_diff(1) < -M_PI) Z_diff(1) += 2*M_PI;

    Tc = Tc + weights_(i) * X_diff * Z_diff.transpose();
  }

  // Calculate the Kalman gain
  MatrixXd K = MatrixXd(5, 3);
  K = Tc * S.inverse();

  //update state mean and covariance matrix
  VectorXd Z_diff = meas_package.raw_measurements_ - z_pred;

  //Angle Normalization
  while (Z_diff(1) > M_PI) Z_diff(1) -= 2*M_PI;
  while (Z_diff(1) < -M_PI) Z_diff(1) += 2*M_PI;

  x_ = x_ + K*Z_diff;
  P_ = P - K * S * K.transpose();
}
