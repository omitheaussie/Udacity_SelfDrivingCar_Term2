#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>
#define SMLNMBR 0.0001

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

  /**
  TODO:
    * Finish initializing the FusionEKF.
    * Set the process (Q) and measurement noises
  */
  H_laser_ <<	1, 0, 0, 0,
				0, 1, 0, 0;
  Hj_ << 1, 1, 0, 0,
	  1, 1, 0, 0,
	  1, 1, 1, 1;

  ekf_.F_ = MatrixXd(4, 4);
  ekf_.F_ << 1, 0, 1, 0,
	  0, 1, 0, 1,
	  0, 0, 1, 0,
	  0, 0, 0, 1;

  //ekf_.x_ = VectorXd(4);
  ekf_.Q_ = MatrixXd(4, 4);
  ekf_.P_ = MatrixXd(4, 4);
  ekf_.P_ << 1, 0, 0, 0,
	  0, 1, 0, 0,
	  0, 0, 1000, 0,
	  0, 0, 0, 1000;
  noise_ax = 9.0;
  noise_ay = 9.0;

}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {


  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    /**
    TODO:
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */
    // first measurement
    //cout << "EKF: " << endl;
	previous_timestamp_ = measurement_pack.timestamp_;
	
    if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
		ekf_.x_ = VectorXd(4);
		ekf_.x_ << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0, 0;
		
		ekf_.curdiffx = 0;
		ekf_.curdiffy = 0;
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Initialize state.
      */
		ekf_.x_= VectorXd(4);
		float px;
		float py;
		float vx;
		float vy;

		px = measurement_pack.raw_measurements_[0] * cos(measurement_pack.raw_measurements_[1]);
		py = measurement_pack.raw_measurements_[0] * sin(measurement_pack.raw_measurements_[1]);
		vx = measurement_pack.raw_measurements_[2] * cos(measurement_pack.raw_measurements_[1]);
		vy = measurement_pack.raw_measurements_[2] * sin(measurement_pack.raw_measurements_[1]);

		ekf_.x_ << px, py , vx, vy;

		ekf_.curdiffx = 0;
		ekf_.curdiffy = 0;
    }
    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }
  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;	//dt - expressed in seconds
  previous_timestamp_ = measurement_pack.timestamp_;
  float dt_2 = dt * dt;
  float dt_3 = dt_2 * dt;
  float dt_4 = dt_3 * dt;

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /**
   TODO:
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
     * Update the process noise covariance matrix.
     * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */
  	ekf_.F_(0, 2) = dt;
	ekf_.F_(1, 3) = dt;

	ekf_.Q_ << dt_4 / 4 * noise_ax, 0, dt_3 / 2 * noise_ax, 0,
		0, dt_4 / 4 * noise_ay, 0, dt_3 / 2 * noise_ay,
		dt_3 / 2 * noise_ax, 0, dt_2*noise_ax, 0,
		0, dt_3 / 2 * noise_ay, 0, dt_2*noise_ay;

	ekf_.Predict();
 

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
   TODO:
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */
  // Check for small ekf_x_
  if(fabs(ekf_.x_[0]) < SMLNMBR) ekf_.x_[0] = SMLNMBR;
  if (fabs(ekf_.x_[1]) < SMLNMBR) ekf_.x_[1] = SMLNMBR;

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {

    // Radar updates
	  Tools tools;
	  ekf_.H_ = tools.CalculateJacobian(ekf_.x_);
	  //cout << "R" << " " << measurement_pack.raw_measurements_[0] << " " << measurement_pack.raw_measurements_[1] << " " << measurement_pack.raw_measurements_[2] << " " << atan2(ekf_.x_[1], ekf_.x_[0]) << endl;
	  //std::cout << ekf_.H_ << endl;
	  ekf_.R_ = R_radar_;
	  
	  ekf_.UpdateEKF(measurement_pack.raw_measurements_);
	  //cout << "R" << " " << ekf_.x_[0] << " " << ekf_.x_[1] << " " << ekf_.x_[2] << " " << ekf_.x_[3] << endl;
  } else {
    // Laser updates
	  ekf_.H_ = H_laser_;
	  ekf_.R_ = R_laser_;
	  //cout << "L" << " " << measurement_pack.raw_measurements_[0] << " " << measurement_pack.raw_measurements_[1] << endl;
	  ekf_.Update(measurement_pack.raw_measurements_);
	  //cout << "L" << " " << ekf_.x_[0] << " " << ekf_.x_[1] << " " << ekf_.x_[2] << " " << ekf_.x_[3] << endl;
  }

  //print the output
  
  //cout << "P_ = " << ekf_.P_ << endl;
}
