#pragma once


#include <vector>
#include <Eigen\Dense>
#include "sensors.h"

class ExtendedKalmanFilter {
public:
	ExtendedKalmanFilter() :init(false), state{ Eigen::Vector3d::Zero() } {};
	Eigen::VectorXd getState()const ;
	Eigen::MatrixXd getCovariance() const ;
	void setState(const Eigen::VectorXd& new_state) ;
	void setCovariance(const Eigen::MatrixXd& new_covariance) ;

	// Important functions fot the Kalman Filter including pysical Model and Prediction and Update
	void predictionStep(double dt);
	void predictionStep(sensorMeas::AccelMeas accMeas, sensorMeas::GyroMeas gyroMeas, sensorMeas::MagMeas magMeas, double dt);
	void updateStep(sensorMeas::GPSMeas gpsMeas, double dt);

private:
	bool init;
	Eigen::VectorXd state;
	Eigen::MatrixXd covariance;
};
