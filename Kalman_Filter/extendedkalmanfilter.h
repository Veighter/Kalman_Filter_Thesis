#pragma once


#include <vector>
#include <Eigen\Dense>
#include "coordTransformation.h"
#include "Fusion.h"



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

enum Sensortype {
	gyroscope,
	accelerometer,
	magnetometer,
	gps
};

class ExtendedKalmanFilterBase {
public:
	ExtendedKalmanFilterBase() {};
	ExtendedKalmanFilterBase(FusionInit fg) :init(false), fusionInit(fg) {
		if (fusionInit == FusionInit::MAG) {
			initMeasurementCounter = Eigen::Vector2i{ 0,0 }; // init with MAG & ACC
		}
		if (fusionInit == FusionInit::MAGGPS) {
			initMeasurementCounter = Eigen::Vector3i{ 0,0,0 }; // init with MAG, ACC & GPS
		}
	};
	bool isInitialised() const { return init; }
	Eigen::VectorXd getState() const {
		return state;
	}
	Eigen::MatrixXd getCovariance() const { return covariance; }
	Eigen::Vector3d getReferenceECEFPosition() const { return referenceECEFPosition; }
	void setFusionInit(FusionInit f_init) {
		fusionInit = f_init;
		if (fusionInit == FusionInit::MAG) {
			initMeasurementCounter = Eigen::Vector2i{ 0,0 }; // init with MAG & ACC
		}
		if (fusionInit == FusionInit::MAGGPS) {
			initMeasurementCounter = Eigen::Vector3i{ 0,0,0 }; // init with MAG, ACC & GPS
		}
	}
	void setState(const Eigen::VectorXd& new_state) { state = new_state; }
	void initFinished() { init = true; }
	void setCovariance(const Eigen::MatrixXd& new_covariance) { covariance = new_covariance; }
	void setReferenceGeodeticPosition(Eigen::Vector3d& referenceGeodeticPosition) {
		setReferenceECEFPosition(referenceGeodeticPosition);
	}
	FusionInit getFusionConfig() { return fusionInit; }
	Eigen::VectorXi getInitMeasurementCounter() { return initMeasurementCounter; }
	void updateMeasurementCount(uint8_t coloumn) {
	initMeasurementCounter[coloumn] += 1; 
	}
	void setMeasurementCount(Eigen::VectorXi counts) { initMeasurementCounter = counts; }
	uint8_t getMinInitMeasurementCount() const { return minInitMeasurements; }

	

	

	// Orientation of the vehicle Frame to which the sensor refers
	Eigen::Quaternion<double> computeInitOrientation(Eigen::VectorXd state);
	Eigen::Quaternion<double> computeOrientation(double roll, double pitch, double yaw);
	

	void setReferenceECEFPosition(Eigen::Vector3d& referenceGeodeticPosition) {
		setRotationMatrixECEF2NED(referenceGeodeticPosition);
		referenceECEFPosition = coordTransformer.geo_to_ecef(referenceGeodeticPosition);
	}

	Eigen::Vector3d computeECEF2NED(Eigen::Vector3d& geodeticPosition) {
		return rotationMatrix * (coordTransformer.geo_to_ecef(geodeticPosition) - referenceECEFPosition);
	}

	Eigen::Vector3d computeNED2ECEFwithRef(Eigen::Vector3d nedPosition) {
		return rotationMatrix.inverse() * nedPosition;
	}

	void setRotationMatrixECEF2NED(Eigen::Vector3d& geodeticPosition) {
		rotationMatrix = coordTransformer.ecef_to_ned_RotationMatrix(geodeticPosition);
	}

	bool isValidMeasurement(Sensortype sensor, Eigen::Vector3d meas);
	
	bool initProcedureStart();
	bool initSamplingFinished();
	

private:
	bool init;
	FusionInit fusionInit;

	uint8_t minInitMeasurements = 10;
	Eigen::VectorXi initMeasurementCounter;	// Counter for initial Measurements for INIT Orientation, Position
	// [Mag, Acc] or [Mag, Acc, GPS]

	Eigen::VectorXd state;
	Eigen::MatrixXd covariance;
	CoordTransformer coordTransformer;
	Eigen::Vector3d referenceECEFPosition;
	Eigen::Matrix3d rotationMatrix;





};

class ExtendedKalmanFilter : public ExtendedKalmanFilterBase {
public:
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

	// Important functions fot the Kalman Filter including pysical Model and Prediction and Update
	void predictionStep(Eigen::Vector3d gyroMeas, double dt);
	void updateAcc(Eigen::Vector3d accMeas, double dt);
	void updateMag(Eigen::Vector3d magMeas, double dt);
	void updateGPS(Eigen::Vector3d gpsMeas, double dt, Eigen::Vector3d gpsVelocityInitial = Eigen::Vector3d(), Eigen::Quaternion<double> orientationInitial = Eigen::Quaternion<double>{ 0,0,0,0 });

	Eigen::Vector3d transform_Gyro(const Eigen::Vector3d& gyroMeas) {
		return getOrientation()._transformVector(gyroMeas);
	}

	Eigen::Vector3d transform_Mag(const Eigen::Vector3d& magMeas) {
		Eigen::Vector3d transformed_magMeas = magMeas;
		transformed_magMeas(1) *= -1;
		transformed_magMeas(2) *= -1;
		return getOrientation()._transformVector(transformed_magMeas);
	}

	Eigen::Vector3d transform_Acc(const Eigen::Vector3d& accMeas) {
		return getOrientation()._transformVector(accMeas);
	}

	void setInitialStateAndCovariance();


private:
	calibration::IMU_Calibration calibrationParams;

	Eigen::Vector3d coords; // coordinates in mm relativ to the VIMU in the construct
	Eigen::Quaternion<double> orientation; // orientation relativ to the vehicle Frame X (vorwaerts)  Y (rechts)  D (own)
};

class VIMUExtendedKalmanFilter : public ExtendedKalmanFilterBase {
public:
	VIMUExtendedKalmanFilter() : ExtendedKalmanFilterBase() {};

	// Konstruktor ohne Init, ueberpruefung ausstehend
	VIMUExtendedKalmanFilter(FusionInit fg) :numIMUs(0), VIMU_Orientations(std::vector<Eigen::Quaternion<double>>()), VIMU_Coords(std::vector<Eigen::Vector3d>()), ExtendedKalmanFilterBase(fg) {
		setState(Eigen::VectorXd::Zero(19));
		setMeasurementCount(Eigen::Vector3i{ getMinInitMeasurementCount(), getMinInitMeasurementCount(), getMinInitMeasurementCount() });
	} VIMUExtendedKalmanFilter(FusionInit fg,int numIMU, std::vector<Eigen::Quaternion<double>> IMU_Orientation, std::vector<Eigen::Vector3d> IMU_Coords) :numIMUs(numIMU), VIMU_Orientations(IMU_Orientation), VIMU_Coords(IMU_Coords), ExtendedKalmanFilterBase(fg) {
			};

	// VIMU has to calculate the average of the measurements
	void predictionStep(std::vector<Eigen::Vector3d> gyroMeas, double dt);
	void updateAcc(std::vector<Eigen::Vector3d> accMeas, double dt);
	void updateMag(std::vector<Eigen::Vector3d> magMeas, double dt);
	void updateGPS(std::vector<Eigen::Vector3d> gpsMeas, double dt, Eigen::Vector3d gpsVelocityInitial = Eigen::Vector3d(), Eigen::Quaternion<double> orientationInitial = Eigen::Quaternion<double>{ 0,0,0,0 });


	Eigen::Vector3d transform_Gyro(const Eigen::Vector3d& gyroMeas, const Eigen::Quaternion<double>& orientation) {
		return orientation._transformVector(gyroMeas);
	}

	Eigen::Vector3d transform_Mag(const Eigen::Vector3d& magMeas, const Eigen::Quaternion<double>& orientation) {
		Eigen::Vector3d transformed_magMeas = magMeas;
		transformed_magMeas(1) *= -1;
		transformed_magMeas(2) *= -1;
		return orientation._transformVector(transformed_magMeas);
	}


	Eigen::Vector3d transform_Acc(const Eigen::Vector3d& accMeas, const Eigen::Vector3d& psi_dot_vimu, const Eigen::Quaternion<double>& orientation, const Eigen::Vector3d& coords, Eigen::Vector3d psi_dot_dot_vimu = Eigen::Vector3d::Zero()) {
			return orientation._transformVector(accMeas);
			// Equation (2) of Data Fusion Algorithms for Multiple Inertial Measurement Units
	//	return orientation.conjugate()._transformVector(accMeas - orientation._transformVector(psi_dot_dot_vimu.cross(coords)) - orientation._transformVector(psi_dot_vimu.cross(psi_dot_vimu.cross(coords))));
	}

	void setInitialStateAndCovariance();

	void federtatedStateAVG(std::vector<Eigen::VectorXd> states);



private:
	int numIMUs;
	std::vector<Eigen::Quaternion<double>> VIMU_Orientations; // not neccessary
	std::vector<Eigen::Vector3d> VIMU_Coords; // not necessary
};
