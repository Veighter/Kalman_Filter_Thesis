#pragma once


#include <vector>
#include <Eigen\Dense>
#include "coordTransformation.h"



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
	ExtendedKalmanFilter() :init(false), calibrated(false) {};
	bool isInitialised() const { return init; }
	Eigen::VectorXd getState() const {
		return state;
	}
	Eigen::MatrixXd getCovariance() const { return covariance; }
	Eigen::Vector3d getReferenceECEFPosition() const { return referenceECEFPosition; }
	void setState(const Eigen::VectorXd& new_state) { state = new_state; }
	void setCovariance(const Eigen::MatrixXd& new_covariance) { covariance = new_covariance; }
	void setReferenceGeodeticPosition(const Eigen::Vector3d& referenceGeodeticPosition) {
		setReferenceECEFPosition(referenceGeodeticPosition);
	}


	void setAccelBias(Eigen::Vector3d b) { calibration_params.accelCali.bias = b; }
	void setAccelTransformMatrix(Eigen::Matrix3d t) { calibration_params.accelCali.theta = t; }
	void setGyroBias(Eigen::Vector3d b) { calibration_params.gyroCali.bias = b; }
	void setGyroTransformMatrix(Eigen::Matrix3d t) { calibration_params.gyroCali.theta = t; }
	void setMagBias(Eigen::Vector3d b) { calibration_params.magCali.bias = b; }
	void setMagTransformMatrix(Eigen::Matrix3d t) { calibration_params.magCali.theta = t; }
	calibration::IMU_Calibration getCalibrationParams() { return calibration_params; }

	void setCoords(Eigen::Vector3d c) { coords = c; }
	Eigen::Vector3d getCoords() { return coords; }
	void setOrientation(Eigen::Quaternion<double> q) { orientation = q; }// Orientation from IMU to VIMU is static!!
	Eigen::Quaternion<double> getOrientation() { return orientation; }



	// Important functions fot the Kalman Filter including pysical Model and Prediction and Update
	void predictionStep(Eigen::Vector3d gyroMeas, double dt);
	void updateAcc(Eigen::Vector3d accMeas, double dt);
	void updateMag(Eigen::Vector3d magMeas, double dt);
	void updateGPS(Eigen::Vector3d gpsMeas, double dt, Eigen::Vector3d gpsVelocityInitial, Eigen::Quaternion<double> orientationInitial);

private:
	bool init, calibrated;
	Eigen::VectorXd state;
	Eigen::MatrixXd covariance;
	CoordTransformer coordTransformer;
	Eigen::Vector3d referenceECEFPosition;

	calibration::IMU_Calibration calibration_params;

	Eigen::Vector3d coords; // coordinates in m (maybe geodetic)
	Eigen::Quaternion<double> orientation;

	void setReferenceECEFPosition(const Eigen::Vector3d& referenceGeodeticPosition) {
		referenceECEFPosition = coordTransformer.geo_to_ecef(referenceGeodeticPosition);
	}



};
