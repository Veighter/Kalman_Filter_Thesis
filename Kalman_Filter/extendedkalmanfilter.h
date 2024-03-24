#pragma once


#include <vector>
#include <Eigen\Dense>
#include "sensors.h"


namespace calibration {
	struct Calibration_Params {
		Eigen::Vector3d bias;
		Eigen::Matrix3d theta;
	};

	/// <summary>
	/// The struct contains all calibrations parameters for the sensors of the used IMU according to their sensor error model
	/// </summary>
	struct IMU_Calibration {
		Calibration_Params accelCali;
		Calibration_Params gyroCali;
		Calibration_Params magCali;
	};
}

class ExtendedKalmanFilter {
public:
	ExtendedKalmanFilter() :init(false),calibrated(false){};
	bool isInitialised() const { return init; }
	Eigen::VectorXd getState() const {
		return state;
	}
	Eigen::MatrixXd getCovariance() const { return covariance; }
	void setState(const Eigen::VectorXd& new_state) { state = new_state; }
	void setCovariance(const Eigen::MatrixXd& new_covariance) { covariance = new_covariance; }

	void setAccelBias(Eigen::Vector3d b) { calibration_params.accelCali.bias = b; }
	void setAccelTransformMatrix(Eigen::Matrix3d t) { calibration_params.accelCali.theta = t; }
	void setGyroBias(Eigen::Vector3d b) { calibration_params.gyroCali.bias = b; }
	void setGyroTransformMatrix(Eigen::Matrix3d t) { calibration_params.gyroCali.theta = t; }
	void setMagBias(Eigen::Vector3d b) { calibration_params.magCali.bias = b; }
	void setMagTransformMatrix(Eigen::Matrix3d t) { calibration_params.magCali.theta = t; }
	calibration::IMU_Calibration getCalibrationParams() { return calibration_params; }

	void setCoords(Eigen::Vector3d c) { coords = c; }
	Eigen::Vector3d getCoords() { return coords; }
	void setOrientation(Eigen::Quaternion<double> q) { orientation = q; }
	Eigen::Quaternion<double> getOrientation() { return orientation; }



	// Important functions fot the Kalman Filter including pysical Model and Prediction and Update
//	void predictionStep(double dt);
	void predictionStep(Eigen::Vector3d accMeas, Eigen::Vector3d gyroMeas, Eigen::Vector3d magMeas, double dt);
	void updateStep(sensorMeas::GPSMeas gpsMeas, double dt);

private:
	bool init, calibrated;
	Eigen::VectorXd state;
	Eigen::MatrixXd covariance;

	calibration::IMU_Calibration calibration_params;

	Eigen::Vector3d coords; // coordinates in mm
	Eigen::Quaternion<double> orientation;



};
