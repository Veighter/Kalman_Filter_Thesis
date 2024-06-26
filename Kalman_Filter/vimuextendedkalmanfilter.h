#pragma once

#include <vector>
#include <Eigen\Dense>
#include "coordTransformation.h"

class VIMUExtendedKalmanFilter {
public:
	VIMUExtendedKalmanFilter() :init(false), calibrated(false) {};
	bool isInitialised() const { return init; }
	Eigen::VectorXd getState() const {
		return state;
	}
	Eigen::MatrixXd getCovariance() const { return covariance; }
	Eigen::Vector3d getReferenceECEFPosition() const { return referenceECEFPosition; }
	void setState(const Eigen::VectorXd& new_state) { state = new_state; }
	void initFinished() { init = true; }
	void setCovariance(const Eigen::MatrixXd& new_covariance) { covariance = new_covariance; }
	void setReferenceGeodeticPosition(Eigen::Vector3d& referenceGeodeticPosition) {
		setReferenceECEFPosition(referenceGeodeticPosition);
	}
	Eigen::Vector3i getInitMeasurementCount() { return initMeasurementCounter; }
	void updateMeasurementCount(uint8_t coloumn) { initMeasurementCounter[coloumn] += 1; }
	uint8_t getMinInitMeasurementCount() const { return minInitMeasurements; }


	void setAccelBias(Eigen::Vector3d b) { calibrationParams.accelCali.bias = b; }
	void setAccelTransformMatrix(Eigen::Matrix3d t) { calibrationParams.accelCali.theta = t; }
	void setGyroBias(Eigen::Vector3d b) { calibrationParams.gyroCali.bias = b; }
	void setGyroTransformMatrix(Eigen::Matrix3d t) { calibrationParams.gyroCali.theta = t; }
	void setMagBias(Eigen::Vector3d b) { calibrationParams.magCali.bias = b; }
	void setMagTransformMatrix(Eigen::Matrix3d t) { calibrationParams.magCali.theta = t; }
	calibration::IMU_Calibration getCalibrationParams() { return calibrationParams; }

	void setCoords(Eigen::Vector3d c) { coords = c / 1e3; }
	Eigen::Vector3d getCoords() { return coords; }
	void setOrientation(Eigen::Quaternion<double> q) { orientation = q; }// Orientation from IMU to VIMU is static!!
	Eigen::Quaternion<double> getOrientation() { return orientation; }
	Eigen::Quaternion<double> computeOrientation(double roll, double pitch, double yaw);



	// Important functions fot the Kalman Filter including pysical Model and Prediction and Update
	void predictionStep(Eigen::Vector3d gyroMeas, double dt);
	void updateAcc(Eigen::Vector3d accMeas, double dt);
	void updateMag(Eigen::Vector3d magMeas, double dt);
	void updateGPS(Eigen::Vector3d gpsMeas, double dt, Eigen::Vector3d gpsVelocityInitial = Eigen::Vector3d(), Eigen::Quaternion<double> orientationInitial = Eigen::Quaternion<double>{ 0,0,0,0 });

	bool isValid(Sensortype sensor);

private:
	bool init, calibrated;
	uint8_t minInitMeasurements = 10;
	Eigen::Vector3i initMeasurementCounter;	// Counter for initial Measurements for INIT Orientation, Position
	// [Mag, Acc, GPS]

	Eigen::VectorXd state;
	Eigen::MatrixXd covariance;
	CoordTransformer coordTransformer;
	Eigen::Vector3d referenceECEFPosition;
	Eigen::Matrix3d rotationMatrix;

	calibration::IMU_Calibration calibrationParams;

	Eigen::Vector3d coords; // coordinates in mm relativ to the VIMU in the construct
	Eigen::Quaternion<double> orientation;

	void setReferenceECEFPosition(Eigen::Vector3d& referenceGeodeticPosition) {
		setRotationMatrixECEF2NED(referenceGeodeticPosition);
		referenceECEFPosition = coordTransformer.geo_to_ecef(referenceGeodeticPosition);
	}

	Eigen::Vector3d computeECEF2NED(Eigen::Vector3d& geodeticPosition) {

		return rotationMatrix * (coordTransformer.geo_to_ecef(geodeticPosition) - referenceECEFPosition);
	}

	void setRotationMatrixECEF2NED(Eigen::Vector3d& geodeticPosition) {
		rotationMatrix = coordTransformer.ecef_to_ned_RotationMatrix(geodeticPosition);
	}
};
