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
	/*
	// Open NIS data files
	NISvals_radar_.open("C:/Users/Omkar/CarND-Term2-Unscented-Kalman-Filter-Project-master/postprocess_NIS/NISvals_radar.txt", ios::out);
	NISvals_laser_.open("C:/Users/Omkar/CarND-Term2-Unscented-Kalman-Filter-Project-master/postprocess_NIS/NISvals_laser.txt", ios::out);

	// Check for errors opening the files
	if (!NISvals_radar_.is_open())
	{
		cout << "Error opening NISvals_radar.txt" << endl;
		exit(1);
	}

	if (!NISvals_laser_.is_open())
	{
		cout << "Error opening NISvals_laser.txt" << endl;
		exit(1);
	}
	*/
	//set state dimension
	n_x_ = 5;

	//set augmented dimension
	n_aug_ = 7;
	//set measurement dimension, radar can measure r, phi, and r_dot
	n_radar_ = 3;

	n_laser_ = 2;

	//define spreading parameter
	lambda_ = 3 - n_aug_;

	// if this is false, laser measurements will be ignored (except during init)
	use_laser_ = true;

	// if this is false, radar measurements will be ignored (except during init)
	use_radar_ = true;

	// initial state vector
	x_ = VectorXd(5);

	// initial covariance matrix
	P_ = MatrixXd(5,5);

	// Process noise standard deviation longitudinal acceleration in m/s^2
	std_a_ = 3.0;

	// Process noise standard deviation yaw acceleration in rad/s^2
	std_yawdd_ = 1.0;

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
	//create vector for weights
	
	x_aug_ = VectorXd(n_aug_);

	deltax_ = VectorXd(n_aug_);

	Xsig_aug_ = MatrixXd(n_aug_, 2 * n_aug_ + 1);

	Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);

	P_aug_ = MatrixXd(n_aug_, n_aug_);

	L_ = MatrixXd(n_aug_, n_aug_);

	weights_ = VectorXd(2 * n_aug_ + 1);
	weights_(0) = lambda_ / (lambda_ + n_aug_);
	for (int i = 1; i<2 * n_aug_ + 1; i++)
		weights_(i) = 0.5 / (n_aug_ + lambda_);

	// Variables used for radar update
	z_pred_radar_ = VectorXd(n_radar_);
	deltaz_radar_ = VectorXd(n_radar_);
	Zsig_radar_ = MatrixXd(n_radar_, 2 * n_aug_ + 1);
	// Radar measurement noise covariance matrix is constant/persistent
	R_radar_ = MatrixXd(n_radar_, n_radar_);
	R_radar_.fill(0.);
	R_radar_(0, 0) = std_radr_*std_radr_;
	R_radar_(1, 1) = std_radphi_*std_radphi_;
	R_radar_(2, 2) = std_radrd_*std_radrd_;
	S_radar_ = MatrixXd(n_radar_, n_radar_);
	Tc_radar_ = MatrixXd(n_x_, n_radar_);
	K_radar_ = MatrixXd(n_x_, n_radar_);

	// Variables used for laser update
	z_pred_laser_ = VectorXd(n_laser_);
	deltaz_laser_ = VectorXd(n_laser_);
	Zsig_laser_ = MatrixXd(n_laser_, 2 * n_aug_ + 1);
	// Laser measurement noise covariance matrix is constant/persistent
	R_laser_ = MatrixXd(n_laser_, n_laser_);
	R_laser_.fill(0.);
	R_laser_(0, 0) = std_laspx_*std_laspx_;
	R_laser_(1, 1) = std_laspy_*std_laspy_;
	S_laser_ = MatrixXd(n_laser_, n_laser_);
	Tc_laser_ = MatrixXd(n_x_, n_laser_);
	K_laser_ = MatrixXd(n_x_, n_laser_);
	/**
	TODO:

	Complete the initialization. See ukf.h for other member properties.

	Hint: one or more values initialized above might be wildly off...
	*/
}

UKF::~UKF() {
	/*
	NISvals_radar_.close();
	NISvals_laser_.close();
	*/
}

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
		// first measurement
		previous_timestamp_ = meas_package.timestamp_;
		if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
			// Initialise state
			x_ << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], 0., 0., 0.;

		}
		else if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
			/**
			Initialize state.
			*/
			float px;
			float py;
			px = meas_package.raw_measurements_[0] * cos(meas_package.raw_measurements_[1]);
			py = meas_package.raw_measurements_[0] * sin(meas_package.raw_measurements_[1]);
	
			x_ << px, py, 0., 0., 0.;

		}
		P_.fill(0.);
		P_(0, 0) = 1.;
		P_(1, 1) = 1.;
		P_(2, 2) = 1.;
		P_(3, 3) = 1.;
		P_(4, 4) = 1.;
		// done initializing, no need to predict or update
		is_initialized_ = true;
		return;
	}
	float dt = (meas_package.timestamp_ - previous_timestamp_) / 1000000.0;	//dt - expressed in seconds
	previous_timestamp_ = meas_package.timestamp_;

	if(dt > 0.0001) Prediction(dt);

	if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
		UpdateLidar(meas_package);
	}
	else if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
		UpdateRadar(meas_package);
	}
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
	// Augmentated state

	//create augmented mean state
	x_aug_.fill(0.0);
	x_aug_.head(n_x_) = x_;
	x_aug_(5) = 0;
	x_aug_(6) = 0;
	//create augmented covariance matrix
	P_aug_.fill(0.0);
	P_aug_.block<5, 5>(0, 0) = P_;
	MatrixXd Q = MatrixXd(2, 2);
	Q.fill(0.0);
	Q << std_a_*std_a_, 0, 0, std_yawdd_*std_yawdd_;
	P_aug_.block<2, 2>(5, 5) = Q;
	//create square root matrix
	double sqrt_l_p_nx = sqrt(lambda_ + n_aug_);
	MatrixXd A = P_aug_.llt().matrixL();
	MatrixXd sqrt_term = sqrt_l_p_nx*A;
	//create augmented sigma points
	Xsig_aug_.fill(0.0);
	Xsig_aug_.block<5, 1>(0, 0) = x_;
	for (int i = 1; i<n_aug_ + 1; i++) {
		Xsig_aug_.block<7, 1>(0, i) = x_aug_ + sqrt_term.block<7, 1>(0, i - 1);
	}
	for (int i = 1; i<n_aug_ + 1; i++) {
		Xsig_aug_.block<7, 1>(0, i + 7) = x_aug_ - sqrt_term.block<7, 1>(0, i - 1);
	}
	// Predict Sigma points
	for (int i = 0; i< 2 * n_aug_ + 1; i++)
	{
		//extract values for better readability
		double p_x = Xsig_aug_(0, i);
		double p_y = Xsig_aug_(1, i);
		double v = Xsig_aug_(2, i);
		double yaw = Xsig_aug_(3, i);
		double yawd = Xsig_aug_(4, i);
		double nu_a = Xsig_aug_(5, i);
		double nu_yawdd = Xsig_aug_(6, i);
		//predicted state values
		double px_p, py_p;

		//avoid division by zero
		if (fabs(yawd) > 0.001) {
			px_p = p_x + v / yawd * (sin(yaw + yawd*delta_t) - sin(yaw));
			py_p = p_y + v / yawd * (cos(yaw) - cos(yaw + yawd*delta_t));
		}
		else {
			px_p = p_x + v*delta_t*cos(yaw) + 0.5f*delta_t*cos(yaw)*nu_a;
			py_p = p_y + v*delta_t*sin(yaw) + 0.5f*delta_t*sin(yaw)*nu_a;
		}

		double v_p = v;
		double yaw_p = yaw + yawd*delta_t + 0.5f*delta_t*nu_yawdd;
		double yawd_p = yawd;

		//add noise
		px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
		py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
		v_p = v_p + nu_a*delta_t;

		yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
		yawd_p = yawd_p + nu_yawdd*delta_t;

		//write predicted sigma point into right column
		Xsig_pred_(0, i) = px_p;
		Xsig_pred_(1, i) = py_p;
		Xsig_pred_(2, i) = v_p;
		Xsig_pred_(3, i) = yaw_p;
		Xsig_pred_(4, i) = yawd_p;
	}
	// Predict mean and covariance
	//predicted state mean
	x_.fill(0.0);
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points
		x_ = x_ + weights_(i) * Xsig_pred_.col(i);
	}

	//predicted state covariance matrix
	P_.fill(0.0);
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points

											   // state difference
		VectorXd x_diff;
		x_diff.fill(0.0);
		x_diff = Xsig_pred_.col(i) - x_;
		//angle normalization
		while (x_diff(3) > M_PI) {
			x_diff(3) -= 2.*M_PI;
		}
		while (x_diff(3) < -M_PI) {
			x_diff(3) += 2.*M_PI;
		}

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
	// Transform sigma points into measurement space
	for (int i = 0; i < 2 * n_aug_ + 1; i++)
	{
		Zsig_laser_(0, i) = Xsig_pred_(0, i);
		Zsig_laser_(1, i) = Xsig_pred_(1, i);
	}

	//calculate mean predicted measurement
	z_pred_laser_.fill(0.);
	for (int i = 0; i < 2 * n_aug_ + 1; i++)
		z_pred_laser_ = z_pred_laser_ + weights_(i)*Zsig_laser_.col(i);

	//calculate measurement covariance matrix S_laser_
	S_laser_.fill(0.);
	for (int i = 0; i < 2 * n_aug_ + 1; i++)
	{
		deltaz_laser_ = Zsig_laser_.col(i) - z_pred_laser_;
		S_laser_ = S_laser_ + weights_(i)*deltaz_laser_*deltaz_laser_.transpose();
	}

	S_laser_ = S_laser_ + R_laser_;

	//calculate cross correlation matrix
	Tc_laser_.fill(0.);
	for (int i = 0; i < 2 * n_aug_ + 1; i++)
	{
		deltax_ = Xsig_pred_.col(i) - x_;
		deltaz_laser_ = Zsig_laser_.col(i) - z_pred_laser_;
		Tc_laser_ = Tc_laser_ + weights_(i)*deltax_*deltaz_laser_.transpose();
	}

	//calculate K_laser_alman gain K_laser_;
	K_laser_ = Tc_laser_*S_laser_.inverse();

	//update state mean and covariance matrix
	deltaz_laser_ = meas_package.raw_measurements_ - z_pred_laser_;

	x_ = x_ + K_laser_*deltaz_laser_;
	P_ = P_ - K_laser_*S_laser_*K_laser_.transpose();

	// NIS Laser
	NIS_laser_ = deltaz_laser_.transpose()*S_laser_.inverse()*deltaz_laser_;
	NISvals_laser_ << NIS_laser_ << endl;
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
	//transform sigma points into measurement space
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

											   // extract values for better readibility
		double p_x = Xsig_pred_(0, i);
		double p_y = Xsig_pred_(1, i);
		double v = Xsig_pred_(2, i);
		double yaw = Xsig_pred_(3, i);

		double v1 = cos(yaw)*v;
		double v2 = sin(yaw)*v;

		// measurement model
		Zsig_radar_(0, i) = sqrt(p_x*p_x + p_y*p_y);                        //r
		Zsig_radar_(1, i) = atan2(p_y, p_x);                                 //phi
		Zsig_radar_(2, i) = (p_x*v1 + p_y*v2) / sqrt(p_x*p_x + p_y*p_y);   //r_dot
	}

	//mean predicted measurement
	z_pred_radar_.fill(0.0);
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {
		z_pred_radar_ = z_pred_radar_ + weights_(i) * Zsig_radar_.col(i);
	}

	//innovation covariance matrix S
	S_radar_.fill(0.0);
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
											   //residual
		deltaz_radar_ = Zsig_radar_.col(i) - z_pred_radar_;

		//angle normalization
		while (deltaz_radar_(1)> M_PI) deltaz_radar_(1) -= 2.*M_PI;
		while (deltaz_radar_(1)<-M_PI) deltaz_radar_(1) += 2.*M_PI;

		S_radar_ = S_radar_ + weights_(i) * deltaz_radar_ * deltaz_radar_.transpose();
	}

	//add measurement noise covariance matrix
	S_radar_ = S_radar_ + R_radar_;


	/*
	UKF Update
	*/

	//calculate cross correlation matrix
	Tc_radar_.fill(0.0);
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

											   //residual
		deltax_ = Xsig_pred_.col(i) - x_;
		deltaz_radar_ = Zsig_radar_.col(i) - z_pred_radar_;
		//angle normalization
		while (deltaz_radar_(1)> M_PI) deltaz_radar_(1) -= 2.*M_PI;
		while (deltaz_radar_(1)<-M_PI) deltaz_radar_(1) += 2.*M_PI;

		// state difference
		//angle normalization
		while (deltax_(3)> M_PI) deltax_(3) -= 2.*M_PI;
		while (deltax_(3)<-M_PI) deltax_(3) += 2.*M_PI;

		Tc_radar_ = Tc_radar_ + weights_(i) * deltax_ * deltaz_radar_.transpose();
	}

	//Kalman gain K;
	K_radar_ = Tc_radar_*S_radar_.inverse();

	//residual
	deltaz_radar_ = meas_package.raw_measurements_ - z_pred_radar_;

	//angle normalization
	while (deltaz_radar_(1)> M_PI) deltaz_radar_(1) -= 2.*M_PI;
	while (deltaz_radar_(1)<-M_PI) deltaz_radar_(1) += 2.*M_PI;

	//update state mean and covariance matrix
	x_ = x_ + K_radar_*deltaz_radar_;
	P_ = P_ - K_radar_*S_radar_*K_radar_.transpose();

	//NIS Radar
	NIS_radar_ = deltaz_radar_.transpose()*S_radar_.inverse()*deltaz_radar_;
	NISvals_radar_ << NIS_radar_ << endl;;
}
