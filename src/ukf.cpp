#include "ukf.h"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

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
  std_a_ = 0.5; // 30;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.15; // 30;

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
  TODO:
  Complete the initialization. See ukf.h for other member properties.
  Hint: one or more values initialized above might be wildly off...
  */
  // State dimension
  n_x_ = 5;

  // Augmented state dimension
  n_aug_ = 7;
  n_sigma_ = 2*n_aug_ + 1;

  // Sigma point spreading parameter
  lambda_ = 3 - n_aug_;

  // Limit for division by zero
  epsilon_ = 1e-12;

  // Weights of the sigma points
  weights_ = VectorXd(2*n_aug_ + 1);

  // Predicted sigma points matrix
  Xsig_pred_ = MatrixXd(n_x_, n_sigma_);

  // Not initialized yet
  is_initialized_ = false;

  NIS_radar_ = INFINITY;
  NIS_laser_ = INFINITY;
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

  // Initialize
  if (!is_initialized_) {
    Initialize(meas_package);

    is_initialized_ = true;
    return;
  }

  // Predict
  double dt_seconds = (meas_package.timestamp_ - time_us_)/1000000.0;
  time_us_ = meas_package.timestamp_;
  Prediction(dt_seconds);

  // Update
  if (use_radar_ && meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    UpdateRadar(meas_package);
  } else if (use_laser_ && meas_package.sensor_type_ == MeasurementPackage::LASER) {
    UpdateLidar(meas_package);
  }
}

void UKF::Initialize(MeasurementPackage meas_package) {

  cout << "Initializing UKF..." << endl;

  // Initialize weights
  weights_(0) = lambda_ / (lambda_ + n_aug_);
  for (int i = 1; i < n_sigma_; ++i) {
      weights_(i) = 1.0 / (2 * (lambda_ + n_aug_));
  }

  // Take the first measurement
  float rho;
  float phi;
  float rhod;
  float px;
  float py;
  float v;
  float yaw;
  float yawd;

  if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      cout << "First measurement is RADAR" << endl;

      rho = meas_package.raw_measurements_(0);
      phi = meas_package.raw_measurements_(1);
      rhod = meas_package.raw_measurements_(2); // unused
      px = rho * cos(phi);
      py = rho * sin(phi);
      v = 0;
      yaw = 0;
      yawd = 0;
    } else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      cout << "First measurement is LASER" << endl;

      px = meas_package.raw_measurements_(0);
      py = meas_package.raw_measurements_(1);
      v = 0;
      yaw = 0;
      yawd = 0;
    }

  // Assign a small value to px and py when both are zero
  if (abs(px) <= epsilon_ && abs(py) <= epsilon_) {
      px = epsilon_;
      py = epsilon_;
    }

  // Initialize state and state covariance
  x_ << px, py, v, yaw, yawd;
  P_.fill(0);
  P_.diagonal().setOnes();

  // Initial values
  cout << "Initial state x:" << endl << x_ << endl;
  cout << "Initial covariance state P:" << endl << P_ << endl;

  // Initialize timestamp
  time_us_ = meas_package.timestamp_;
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

  // Generate sigma points (augmented)
  MatrixXd Xsig_aug_ = GenerateSigmaPoints();

  // Predict sigma points
  PredictSigmaPoints(delta_t, Xsig_aug_);

  // Predict state mean and state covariance
  PredictStateMean();
  PredictStateCovariance();
}

void UKF::PredictStateCovariance() {
  P_.fill(0.0);
  for (int i = 0; i < n_sigma_; i++) {

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;

    // angle normalization
    while (x_diff(3) >  M_PI) x_diff(3) -= 2.*M_PI;
    while (x_diff(3) < -M_PI) x_diff(3) += 2.*M_PI;

    P_ = P_ + weights_(i) * x_diff * x_diff.transpose();
  }
}

void UKF::PredictStateMean() {
  x_.fill(0.0);
  for (int i = 0; i < n_sigma_; i++) {
    x_ += weights_(i) * Xsig_pred_.col(i);
  }
}

void UKF::PredictSigmaPoints(double delta_t, const MatrixXd &Xsig_aug_) {
  for (int i = 0; i < n_sigma_; i++) {
    double px = Xsig_aug_(0, i);
    double py = Xsig_aug_(1, i);
    double v = Xsig_aug_(2, i);
    double yaw = Xsig_aug_(3, i);
    double yawd = Xsig_aug_(4, i);
    double nu_a = Xsig_aug_(5, i);
    double nu_yawdd = Xsig_aug_(6, i);

    // predicted state values
    double px_p, py_p;
    double v_p = v;
    double yaw_p = yaw + yawd*delta_t;
    double yawd_p = yawd;

    // avoid division by zero
    if (fabs(yawd) > epsilon_) {
      px_p = px + (v/yawd) * (sin(yaw_p) - sin(yaw));
      py_p = py + (v/yawd) * (cos(yaw) - cos(yaw_p));
    } else {
      px_p = px + v*delta_t*cos(yaw);
      py_p = py + v*delta_t*sin(yaw);
    }

    // add noise
    px_p = px_p + 0.5*nu_a*delta_t*delta_t*cos(yaw);
    py_p = py_p + 0.5*nu_a*delta_t*delta_t*sin(yaw);
    v_p = v_p + nu_a*delta_t;
    yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
    yawd_p = yawd_p + nu_yawdd*delta_t;

    // write predicted sigma point into right column
    Xsig_pred_(0, i) = px_p;
    Xsig_pred_(1, i) = py_p;
    Xsig_pred_(2, i) = v_p;
    Xsig_pred_(3, i) = yaw_p;
    Xsig_pred_(4, i) = yawd_p;
  }
}

MatrixXd UKF::GenerateSigmaPoints() const {
  VectorXd x_aug_ = VectorXd(n_aug_);
  MatrixXd P_aug_ = MatrixXd(n_aug_, n_aug_);
  MatrixXd Xsig_aug_ = MatrixXd(n_aug_, n_sigma_);

  x_aug_ << x_, 0, 0;

  P_aug_.fill(0);
  P_aug_.block(0, 0, n_x_, n_x_) << P_;
  P_aug_.block(n_x_, n_x_, 2, 2) << pow(std_a_, 2),                  0,
                                                 0, pow(std_yawdd_, 2);

  MatrixXd A = P_aug_.llt().matrixL();

  Xsig_aug_.col(0) = x_aug_;
  for (int i = 0; i < n_aug_; i++) {
    Xsig_aug_.col(i+1)          = x_aug_ + sqrt(lambda_ + n_aug_) * A.col(i);
    Xsig_aug_.col(i + 1 + n_aug_) = x_aug_ - sqrt(lambda_ + n_aug_) * A.col(i);
  }
  return Xsig_aug_;
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
  // ------
  // Update
  // ------

  //set measurement dimension, lidar can measure px and py
  int n_z = 2;

  // Predict measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0);

  MatrixXd Zsig = MatrixXd(n_z, n_sigma_);
  Zsig.fill(0);

  for (int i = 0; i < n_sigma_; i++) {
    // extract values for better readability
    double px = Xsig_pred_(0, i);
    double py = Xsig_pred_(1, i);

    Zsig(0, i) = px;
    Zsig(1, i) = py;
 }

  for (int i = 0; i < n_sigma_; i++) {
    z_pred = z_pred + weights_(i)*Zsig.col(i);
  }

  // measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z, n_z);

  S.fill(0.0);
  for (int i = 0; i < n_sigma_; i++) {
    // residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    // angle normalization
    while (z_diff(1) >  M_PI) z_diff(1) -= 2.*M_PI;
    while (z_diff(1) < -M_PI) z_diff(1) += 2.*M_PI;

    S = S + weights_(i)*z_diff*z_diff.transpose();
  }

  // add measurement noise covariance matrix
  MatrixXd R = MatrixXd(n_z, n_z);

  R << std_laspx_*std_laspx_,                     0,
                           0, std_laspy_*std_laspy_;

  S = S + R;

  // Update state

  // create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);

  //calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < n_sigma_; i++) {

    // residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    // angle normalization
    while (z_diff(1) >  M_PI) z_diff(1) -= 2.*M_PI;
    while (z_diff(1) < -M_PI) z_diff(1) += 2.*M_PI;

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    while (x_diff(3) >  M_PI) x_diff(3) -= 2.*M_PI;
    while (x_diff(3) < -M_PI) x_diff(3) += 2.*M_PI;

    Tc = Tc + weights_(i)*x_diff*z_diff.transpose();
  }

  // Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  // residual
  VectorXd z = VectorXd(n_z);
  z = meas_package.raw_measurements_;
  VectorXd z_diff = z - z_pred;

  // angle normalization
  while (z_diff(1) >  M_PI) z_diff(1) -= 2.*M_PI;
  while (z_diff(1) < -M_PI) z_diff(1) += 2.*M_PI;

  // update state mean and covariance matrix
  x_ = x_ + K*z_diff;
  P_ = P_ - K*S*K.transpose();

  // Laser NIS
  NIS_laser_ = z_diff.transpose()*S.inverse()*z_diff;
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
  // ------
  // Update
  // ------

  // Predict measurement

  // create matrix for sigma points in measurement space
  int n_z = 3;
  MatrixXd Zsig = MatrixXd(n_z, n_sigma_);
  Zsig.fill(0);

  // transform sigma points into measurement space
  for (int i = 0; i < n_sigma_; i++) {

    // extract values for better readability
    double px = Xsig_pred_(0, i);
    double py = Xsig_pred_(1, i);
    double v = Xsig_pred_(2, i);
    double yaw = Xsig_pred_(3, i);

    // handle division by zero
    double r = sqrt(px*px + py*py);
    if (r < epsilon_) { r = epsilon_; }

    // measurement model
    Zsig(0, i) = r;                                 // r
    Zsig(1, i) = atan2(py, px);                     // phi
    Zsig(2, i) = (px*cos(yaw) + py*sin(yaw))*v / r; // r_dot
  }

  // mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);

  z_pred.fill(0.0);
  for (int i = 0; i < n_sigma_; i++) {
    z_pred = z_pred + weights_(i) * Zsig.col(i);
  }

  // measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z, n_z);

  S.fill(0.0);
  for (int i = 0; i < n_sigma_; i++) {
    // residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    // angle normalization
    while (z_diff(1) >  M_PI) z_diff(1) -= 2.*M_PI;
    while (z_diff(1) < -M_PI) z_diff(1) += 2.*M_PI;

    S = S + weights_(i)*z_diff*z_diff.transpose();
  }

  // add measurement noise covariance matrix
  MatrixXd R = MatrixXd(n_z,n_z);

  R << std_radr_*std_radr_,                       0,                     0,
                         0, std_radphi_*std_radphi_,                     0,
                         0,                       0, std_radrd_*std_radrd_;

  S = S + R;

  // Update state

  // create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);

  // calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < n_sigma_; i++) {

    // residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    // angle normalization
    while (z_diff(1) >  M_PI) z_diff(1) -= 2.*M_PI;
    while (z_diff(1) < -M_PI) z_diff(1) += 2.*M_PI;

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    // angle normalization
    while (x_diff(3) >  M_PI) x_diff(3) -= 2.*M_PI;
    while (x_diff(3) < -M_PI) x_diff(3) += 2.*M_PI;

    Tc = Tc + weights_(i)*x_diff*z_diff.transpose();
  }

  // Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  // residual
  VectorXd z = VectorXd(n_z);
  z = meas_package.raw_measurements_;
  VectorXd z_diff = z - z_pred;

  // angle normalization
  while (z_diff(1) >  M_PI) z_diff(1) -= 2.*M_PI;
  while (z_diff(1) < -M_PI) z_diff(1) += 2.*M_PI;

  // update state mean and covariance matrix
  x_ = x_ + K*z_diff;
  P_ = P_ - K*S*K.transpose();

  // Radar NIS
  NIS_radar_ = z_diff.transpose()*S.inverse()*z_diff;
}
