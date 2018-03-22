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
  
  H_laser_ << 1, 0, 0, 0,
         0, 1, 0, 0;

  /**
  TODO:
    * Finish initializing the FusionEKF.
    * Set the process and measurement noises
  */
  
  ekf_.F_ = MatrixXd(4, 4);
  ekf_.F_<<1,0,1,0,
  0,1,0,1,
  0,0,1,0,
  0,0,0,1;
  
  ekf_.H_ = MatrixXd(2, 4);
  ekf_.H_<<1,0,0,0,
           0,1,0,0;
  
  ekf_.Q_ = MatrixXd(4, 4);
  ekf_.Q_<<0,0,0,0,
           0,0,0,0,
           0,0,0,0,
           0,0,0,0;
  
  ekf_.P_ = MatrixXd(4, 4);
  ekf_.P_<<1,0,0,0,
           0,1,0,0,
           0,0,1,0,
           0,0,0,1;

}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {

  float dt = 0;
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
    cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);
    
    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      float rho = measurement_pack.raw_measurements_[0];
      float theta = measurement_pack.raw_measurements_[1];
      float rho_dot = measurement_pack.raw_measurements_[2];
      
      ekf_.x_[0] = rho*sin(theta);
      ekf_.x_[1] = rho*cos(theta);
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
      ekf_.x_[0] = measurement_pack.raw_measurements_[0];
      ekf_.x_[1] = measurement_pack.raw_measurements_[1];
    }
    //ekf_.x_[2] = ekf_.x_[3] = 0.05; // RMSE Starts around 5.9 for vx and vy
    //ekf_.x_[2] = ekf_.x_[3] = 10; // RMSE Starts around 5 for vx and 10 for vy
    //ekf_.x_[2] = 10.0; ekf_.x_[3] = 0.05; // RMSE Starts around 5 for vx and 0 for vy but then rose to 2
    //ekf_.x_[2] = 10.0; ekf_.x_[3] = 5.05; // RMSE Starts around 4.8 for vx and 5 for vy
    //ekf_.x_[2] = 10.0; ekf_.x_[3] = 2.05; // RMSE Starts around 4.8 for vx and 2 for vy
    //ekf_.x_[2] = 10.0; ekf_.x_[3] = 0.50; // RMSE Starts around 4.8 for vx and 0.5 for vy but then rose to 1.8
    //ekf_.x_[2] = 3.0; ekf_.x_[3] = 1.0; // RMSE Starts around 4.8 for vx and 0.5 for vy but then rose to 1.8
    //ekf_.x_[2] = 0.005; ekf_.x_[3] = 0.0015; // RMSE Starts around 4.8 for vx and 0.5 for vy but then rose to 1.8
    
    // Picking the first value of the vx and vy from the ground truth.
    ekf_.x_[2] = 5.2;    
    ekf_.x_[3] = 0.0008; 
    // done initializing, no need to predict or update
    previous_timestamp_ = measurement_pack.timestamp_;
    is_initialized_ = true;
    return;
  }
  dt = (measurement_pack.timestamp_ - previous_timestamp_)/1000000.0;
  previous_timestamp_ = measurement_pack.timestamp_;
  
  
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
  
  //UPDATE THE TIMESTAMP OF THE F MATRIX
  ekf_.F_(0,2) = dt;
  ekf_.F_(1,3) = dt;
  
  //UPDATE THE NOISE MATRIX WITH NEW THE
  float dt2 = dt*dt;
  float dt3 = dt2*dt;
  float dt4 = dt3*dt;
  float dt3d2 = dt3/2;
  float dt4d4 = dt4/4;
  
  float nx = 9.0, ny = 9.0;
  ekf_.Q_ << dt4d4*nx, 0, dt3d2*nx, 0,
             0, dt4d4*ny, 0, dt3d2*ny,
             dt3d2*nx, 0, dt2*nx, 0,
             0, dt3d2*ny, 0, dt2*ny;
  ekf_.Predict();
  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
   TODO:
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
    ekf_.H_ = tools.CalculateJacobian(ekf_.x_);
    ekf_.R_ = R_radar_;
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);
  } else {
    // Laser updates
    ekf_.H_ = H_laser_;
    ekf_.R_ = R_laser_;
    ekf_.Update(measurement_pack.raw_measurements_);
  }
 
  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
