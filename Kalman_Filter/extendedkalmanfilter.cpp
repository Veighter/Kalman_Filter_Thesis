#include "extendedkalmanfilter.h"
#include <iostream>

// For computing the Fehlerfortpflanzungsmatrix
constexpr double ACCEL_STD = 1.0;
constexpr double GYRO_STD = 0.1 / 180.0 * M_PI;
constexpr double MAG_STD = 1.0;
constexpr double GPS_POS_STD = 3.0;

// Initial standard derivation has to be high, bc we dont know exactly what they are
constexpr double INIT_VEL_STD = 10.0;
constexpr double INIT_ORIENTATION_VAR = 0.125;
constexpr double INIT_ORIENTATION_COVAR = 0.0003;

constexpr double g = 9.81;

// reference values to detect outliers simple way
constexpr double ref_lat = 51.0;
constexpr double ref_long = 7.0;


// magnetic-declination.com: Location is Haltern with 3 degree and 11 sec = 3.1833 degree
constexpr double declinationAngle = 3.1833 * M_PI / 180.0;



void ExtendedKalmanFilter::predictionStep(Eigen::Vector3d gyroMeas, double dt) {
	if (isInitialised()) {
		Eigen::VectorXd state = getState();
		Eigen::MatrixXd covariance = getCovariance(); // correlates to the gyro and accelerometer data due to orientation and position integration form the sensor data

		double x = state(0);
		double y = state(1);
		double z = state(2);
		double x_dot = state(3);
		double y_dot = state(4);
		double z_dot = state(5);
		double x_dot_dot = state(6);
		double y_dot_dot = state(7);
		double z_dot_dot = state(8);
		double q_0 = state(9);
		double q_1 = state(10);
		double q_2 = state(11);
		double q_3 = state(12);

		gyroMeas = transform_Gyro(gyroMeas);

		Eigen::Quaternion<double> orientation = Eigen::Quaternion<double>{ q_0,q_1,q_2,q_3 };

		// Conversion from ("body") XYD -> ("local") NED
		gyroMeas = orientation._transformVector(gyroMeas);

		double psi_x_dot = gyroMeas(0);
		double psi_y_dot = gyroMeas(1);
		double psi_z_dot = gyroMeas(2);


		// state = x,y,z,x_dot,y_dot,z_dot,x_dot_dot,y_dot_dot,z_dot_dot,psi_x_dot,psi_y_dot,q0,q1,q2,q3,B_x,B_y,B_z
		Eigen::MatrixXd F = Eigen::MatrixXd::Zero(13, 13);
		F.row(0) << 1, 0, 0, dt, 0, 0, 1. / 2 * dt * dt, 0, 0, 0, 0, 0, 0;
		F.row(1) << 0, 1, 0, 0, dt, 0, 0, 1. / 2 * dt * dt, 0, 0, 0, 0, 0;
		F.row(2) << 0, 0, 1, 0, 0, dt, 0, 0, 1. / 2 * dt * dt, 0, 0, 0, 0;
		F.row(3) << 0, 0, 0, 1, 0, 0, dt, 0, 0, 0, 0, 0, 0;
		F.row(4) << 0, 0, 0, 0, 1, 0, 0, dt, 0, 0, 0, 0, 0;
		F.row(5) << 0, 0, 0, 0, 0, 1, 0, 0, dt, 0, 0, 0, 0;
		F.row(6) << 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0;
		F.row(7) << 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0;
		F.row(8) << 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;
		F.row(9) << 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1. / 2 * dt * psi_x_dot, -1. / 2 * dt * psi_y_dot, -1. / 2 * dt * psi_z_dot;
		F.row(10) << 0, 0, 0, 0, 0, 0, 0, 0, 0, 1. / 2 * dt * psi_x_dot, 1, 1. / 2 * dt * psi_z_dot, -1. / 2 * dt * psi_y_dot;
		F.row(11) << 0, 0, 0, 0, 0, 0, 0, 0, 0, 1. / 2 * dt * psi_y_dot, -1. / 2 * dt * psi_z_dot, 1, 1. / 2 * dt * psi_x_dot;
		F.row(12) << 0, 0, 0, 0, 0, 0, 0, 0, 0, 1. / 2 * dt * psi_z_dot, 1. / 2 * dt * psi_y_dot, -1. / 2 * dt * psi_x_dot, 1;

		Eigen::MatrixXd Q = Eigen::MatrixXd::Zero(13, 13);

		// Only Orientation error is looked at now, other dependencies not known (bad)
		double gyro_var = GYRO_STD * GYRO_STD;
		Q(9, 9) = gyro_var;
		Q(10, 10) = gyro_var;
		Q(11, 11) = gyro_var;
		Q(12, 12) = gyro_var;



		state = F * state;
		covariance = F * covariance * F.transpose() + Q;

		Eigen::Quaternion<double> q = Eigen::Quaternion<double>{ state(9), state(10), state(11), state(12) };
		q.normalize();
		state(9) = q.w();
		state(10) = q.x();
		state(11) = q.y();
		state(12) = q.z();

		setState(state);
		setCovariance(covariance);
		//	std::cout << "GYRO IMU\n";
	}

}

void ExtendedKalmanFilter::updateAcc(Eigen::Vector3d accMeas, double dt) {

	if (!isInitialised()) {
		if (initSamplingFinished()) {
			setInitialStateAndCovariance();
		}
		else {
			Eigen::VectorXd state;
			if (initProcedureStart()) {
				state = Eigen::VectorXd::Zero(13);
			}
			else {
				state = getState();
			}
			// Transformation from Sensor Acc to vehicle Body Acc
			accMeas = transform_Acc(accMeas);

			state(6) += accMeas(0);
			state(7) += accMeas(1);
			state(8) += accMeas(2);

			updateMeasurementCount(0);

			setState(state);
			return;

		}
	}

	if (isInitialised()) {
		Eigen::VectorXd state = getState();
		Eigen::MatrixXd covariance = getCovariance();

		// Transformation from Sensor Acc to vehicle Body Acc
		accMeas = transform_Acc(accMeas);
		//Conversion from ("body") XYD -> ("local") NED
		Eigen::Quaternion<double> orientation = Eigen::Quaternion<double>{ state(9), state(10), state(11), state(12) };
		accMeas = orientation._transformVector(accMeas);

		Eigen::Vector3d gMeas = accMeas;
		gMeas.normalize();

		// Correction for the gravity	
		accMeas(2) += g;

		// Acc Correction, Measurement 
		Eigen::MatrixXd H = Eigen::MatrixXd::Zero(3, 13);
		Eigen::MatrixXd R = Eigen::MatrixXd::Zero(3, 3);

		double acc_var{};
		acc_var = ACCEL_STD * ACCEL_STD;
		R(0, 0) = acc_var;
		R(1, 1) = acc_var;
		R(2, 2) = acc_var;

		H.row(0) << 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0;
		H.row(1) << 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0;
		H.row(2) << 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;

		Eigen::Vector3d y = accMeas - H * state;
		Eigen::MatrixXd S = H * covariance * H.transpose() + R;
		Eigen::MatrixXd K = covariance * H.transpose() * S.inverse();

		state = state + K * y;
		covariance = (Eigen::MatrixXd::Identity(13, 13) - K * H) * covariance;



		// Orientation Correciton
		Eigen::Vector3d p_dot_dot_hat = Eigen::Vector3d::Zero();

		double q_0 = state(9);
		double q_1 = state(10);
		double q_2 = state(11);
		double q_3 = state(12);
		p_dot_dot_hat(0) = 2 * q_1 * q_3 - 2 * q_0 * q_2;
		p_dot_dot_hat(1) = 2 * q_2 * q_3 + 2 * q_0 * q_1;
		p_dot_dot_hat(2) = q_0 * q_0 - q_1 * q_1 - q_2 * q_2 + q_3 * q_3;

		p_dot_dot_hat = g * p_dot_dot_hat;
		y = gMeas - p_dot_dot_hat;

		// non linear modell, bc of direction cosine
		H.row(0) << 0, 0, 0, 0, 0, 0, 0, 0, 0, -2 * q_2, 2 * q_3, -2 * q_0, 2 * q_1;
		H.row(1) << 0, 0, 0, 0, 0, 0, 0, 0, 0, 2 * q_1, 2 * q_0, 2 * q_3, 2 * q_1;
		H.row(2) << 0, 0, 0, 0, 0, 0, 0, 0, 0, 2 * q_0, -2 * q_1, -2 * q_2, 2 * q_3;

		S = H * covariance * H.transpose() + R;
		K = covariance * H.transpose() * S.inverse();


		// Dont update the YAW Component with Accelerometer-Data
		Eigen::VectorXd state_error = Eigen::VectorXd::Zero(13);
		state_error = K * y;
		state_error(12) = 0;

		state = state + state_error;
		covariance = (Eigen::MatrixXd::Identity(13, 13) - K * H) * covariance;

		setState(state);
		setCovariance(covariance);
		//	std::cout << "ACC IMU\n";

	}

}

void ExtendedKalmanFilter::updateMag(Eigen::Vector3d magMeas, double dt) {

	if (!isInitialised()) {
		if (initSamplingFinished()) {
			setInitialStateAndCovariance();
		}
		else {
			Eigen::VectorXd state;
			if (initProcedureStart()) { // Init State
				state = Eigen::VectorXd::Zero(13);
			}
			else {
				state = getState();
			}
			// Transformation from Sensor Acc to vehicle Body Acc
			magMeas = transform_Mag(magMeas);

			// Use uninitialisied Quaternion state for placeholder
			state(9) += magMeas(0);
			state(10) += magMeas(1);
			state(11) += magMeas(2);

			updateMeasurementCount(1);

			setState(state);
			return;

		}
	}

	if (isInitialised()) {
		Eigen::VectorXd state = getState();
		Eigen::MatrixXd covariance = getCovariance();

		Eigen::Vector3d mag_hat = Eigen::Vector3d::Zero();
		// Richtung des Magnetfelds
		double q_0 = state(9);
		double q_1 = state(10);
		double q_2 = state(11);
		double q_3 = state(12);


		mag_hat(0) = q_0 * q_0 + q_1 * q_1 - q_2 * q_2 - q_3 * q_3;
		mag_hat(1) = 2 * q_1 * q_2 - 2 * q_0 * q_3;
		mag_hat(2) = 2 * q_1 * q_3 + 2 * q_0 * q_2;

		// Transformation from Sensor Acc to vehicle Body Acc
		magMeas = transform_Mag(magMeas);


		Eigen::Quaternion<double> orientation = Eigen::Quaternion<double>{ q_0,q_1,q_2,q_3 };
		magMeas = orientation._transformVector(magMeas);
		magMeas.normalize();
		Eigen::Vector3d y = magMeas - mag_hat;

		Eigen::MatrixXd H = Eigen::MatrixXd::Zero(3, 13);
		Eigen::MatrixXd R = Eigen::MatrixXd::Zero(3, 3);

		R(0, 0) = 1;
		R(1, 1) = 1;
		R(2, 2) = 1;

		H.row(0) << 0, 0, 0, 0, 0, 0, 0, 0, 0, 2 * q_0, 2 * q_1, -2 * q_2, -2 * q_3;
		H.row(1) << 0, 0, 0, 0, 0, 0, 0, 0, 0, -2 * q_3, 2 * q_2, 2 * q_1, -2 * q_0;
		H.row(2) << 0, 0, 0, 0, 0, 0, 0, 0, 0, 2 * q_2, 2 * q_3, 2 * q_0, 2 * q_1;

		Eigen::MatrixXd S = H * covariance * H.transpose() + R;
		Eigen::MatrixXd K = covariance * H.transpose() * S.inverse();

		// Dont update the Roll and Pitch Component with Magnetometer-Data
		Eigen::VectorXd state_error = Eigen::VectorXd::Zero(13);
		state_error = K * y;
		state_error(10) = 0;
		state_error(11) = 0;

		state = state + state_error;
		covariance = (Eigen::MatrixXd::Identity(13, 13) - K * H) * covariance;

		setState(state);
		setCovariance(covariance);
		//	std::cout << "MAG\n";

	}

}

/// <summary>
/// Method for the Measurement Update of the navigational Object state with the GPS 
/// </summary>
/// <param name="gpsMeas"></param>GPS given in LLA
/// <param name="dt"></param>
void ExtendedKalmanFilter::updateGPS(Eigen::Vector3d gpsMeas, double dt, Eigen::Vector3d gpsVelocityInitial, Eigen::Quaternion<double> orientationInitial) {

	if (!isInitialised()) {

		if (initSamplingFinished()) {
			Eigen::VectorXd state;
			if (initProcedureStart()) {
				state = Eigen::VectorXd::Zero(13);
			}
			else {
				state = getState();
			}
			state(0) += gpsMeas(0);
			state(1) += gpsMeas(1);
			state(2) += 0;

			setState(state);
			updateMeasurementCount(2); // Update property in Measurement Count for GPS measurement
		}
		else {
			Eigen::VectorXd state = getState();



			// Mean (better Median) of the measurements
			gpsMeas(0) += state(0);
			gpsMeas(1) += state(1);
			gpsMeas(2) += 0;

			gpsMeas = gpsMeas / (getInitMeasurementCounter()(2) + 1); // mean of the init GPS Positions (11th measurement for GPS)

			setReferenceGeodeticPosition(gpsMeas);

			setInitialStateAndCovariance();
		}
	}
	else {
		Eigen::VectorXd state = getState();
		Eigen::MatrixXd covariance = getCovariance();

		// GPS in NED locally Coordinates at Position of the (V)IMU
		gpsMeas = computeECEF2NED(gpsMeas) + getCoords();



		Eigen::MatrixXd H = Eigen::MatrixXd::Zero(3, 13);
		Eigen::MatrixXd R = Eigen::MatrixXd::Zero(3, 3);

		double gps_var = GPS_POS_STD * GPS_POS_STD;

		R(0, 0) = gps_var;
		R(1, 1) = gps_var;
		R(2, 2) = gps_var;

		H.row(0) << 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
		H.row(1) << 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
		H.row(2) << 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;

		Eigen::Vector3d y = gpsMeas - H * state;
		Eigen::MatrixXd S = H * covariance * H.transpose() + R;
		Eigen::MatrixXd K = covariance * H.transpose() * S.inverse();

		state = state + K * y;
		covariance = (Eigen::MatrixXd::Identity(13, 13) - K * H) * covariance;

		setState(state);
		setCovariance(covariance);
		//	std::cout << "GPS IMU\n";
	}

}

// VIMU has to calculate the average of the measurements
void VIMUExtendedKalmanFilter::predictionStep(std::vector<Eigen::Vector3d> gyroMeas, double dt) {

	if (isInitialised()) {
		Eigen::VectorXd state = getState();
		Eigen::MatrixXd covariance = getCovariance(); // correlates to the gyro and accelerometer data due to orientation and position integration form the sensor data

		double x = state(0);
		double y = state(1);
		double z = state(2);
		double x_dot = state(3);
		double y_dot = state(4);
		double z_dot = state(5);
		double x_dot_dot = state(6);
		double y_dot_dot = state(7);
		double z_dot_dot = state(8);
		double q_0 = state(9);
		double q_1 = state(10);
		double q_2 = state(11);
		double q_3 = state(12);
		double psi_x_dot_old = state(13);
		double psi_y_dot_old = state(14);
		double psi_z_dot_old = state(15);




		// Gyro AVG over measurements
		// Counts the number of imu influenced in the current measurement, ignores outliers
		int num_IMUs_avg = 0;
		Eigen::Vector3d gyroMeasAvg = Eigen::Vector3d::Zero();

		// 1. transform
		for (int i = 0; i < gyroMeas.size(); i++) {
			Eigen::Vector3d transformed_gyroMeas = transform_Gyro(gyroMeas[i], VIMU_Orientations[i]);
			if (isValidMeasurement(gyroscope, transformed_gyroMeas)) {
				gyroMeasAvg += transformed_gyroMeas;
				num_IMUs_avg++;
			}
		}
		// 2. AVG
		gyroMeasAvg /= num_IMUs_avg;

		Eigen::Vector3d gyroMeas = gyroMeasAvg;

		// Conversion from ("body") XYD -> ("local") NED
		Eigen::Quaternion<double> orientation = Eigen::Quaternion<double>{ q_0,q_1,q_2,q_3 };
		gyroMeas = orientation._transformVector(gyroMeas);

		double psi_x_dot = gyroMeas(0);
		double psi_y_dot = gyroMeas(1);
		double psi_z_dot = gyroMeas(2);

		// state = x,y,z,x_dot,y_dot,z_dot,x_dot_dot,y_dot_dot,z_dot_dot,psi_x_dot,psi_y_dot,q0,q1,q2,q3,B_x,B_y,B_z
		Eigen::MatrixXd F = Eigen::MatrixXd::Zero(19, 19);
		F.row(0) << 1, 0, 0, dt, 0, 0, 1. / 2 * dt * dt, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
		F.row(1) << 0, 1, 0, 0, dt, 0, 0, 1. / 2 * dt * dt, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
		F.row(2) << 0, 0, 1, 0, 0, dt, 0, 0, 1. / 2 * dt * dt, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
		F.row(3) << 0, 0, 0, 1, 0, 0, dt, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
		F.row(4) << 0, 0, 0, 0, 1, 0, 0, dt, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
		F.row(5) << 0, 0, 0, 0, 0, 1, 0, 0, dt, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
		F.row(6) << 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
		F.row(7) << 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
		F.row(8) << 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
		F.row(9) << 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1. / 2 * dt * psi_x_dot, -1. / 2 * dt * psi_y_dot, -1. / 2 * dt * psi_z_dot, 0, 0, 0, 0, 0, 0;
		F.row(10) << 0, 0, 0, 0, 0, 0, 0, 0, 0, 1. / 2 * dt * psi_x_dot, 1, 1. / 2 * dt * psi_z_dot, -1. / 2 * dt * psi_y_dot, 0, 0, 0, 0, 0, 0;
		F.row(11) << 0, 0, 0, 0, 0, 0, 0, 0, 0, 1. / 2 * dt * psi_y_dot, -1. / 2 * dt * psi_z_dot, 1, 1. / 2 * dt * psi_x_dot, 0, 0, 0, 0, 0, 0;
		F.row(12) << 0, 0, 0, 0, 0, 0, 0, 0, 0, 1. / 2 * dt * psi_z_dot, 1. / 2 * dt * psi_y_dot, -1. / 2 * dt * psi_x_dot, 1, 0, 0, 0, 0, 0, 0;
		F.row(13) << 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
		F.row(14) << 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
		F.row(15) << 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
		F.row(16) << 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
		F.row(17) << 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
		F.row(18) << 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;

		Eigen::MatrixXd Q = Eigen::MatrixXd::Zero(19, 19);

		// Only Orientation error is looked at now, other dependencies not known (bad)
		double gyro_var = GYRO_STD * GYRO_STD;
		Q(9, 9) = gyro_var;
		Q(10, 10) = gyro_var;
		Q(11, 11) = gyro_var;
		Q(12, 12) = gyro_var;
		Q(13, 13) = gyro_var;
		Q(14, 14) = gyro_var;
		Q(15, 15) = gyro_var;
		Q(16, 16) = gyro_var;
		Q(17, 17) = gyro_var;
		Q(18, 18) = gyro_var;

		state = F * state;
		covariance = F * covariance * F.transpose() + Q;

		Eigen::Quaternion<double> q = Eigen::Quaternion<double>{ state(9), state(10), state(11), state(12) };
		q.normalize();
		state(9) = q.w();
		state(10) = q.x();
		state(11) = q.y();
		state(12) = q.z();

		// update avg angular rate and angular acceleration
		state(13) = psi_x_dot;
		state(14) = psi_y_dot;
		state(15) = psi_z_dot;
		state(16) = (psi_x_dot - psi_x_dot_old) / dt;
		state(17) = (psi_y_dot - psi_y_dot_old) / dt;
		state(18) = (psi_z_dot - psi_z_dot_old) / dt;

		setState(state);
		setCovariance(covariance);
		//	std::cout << "GYRO\n";
	}
};


void VIMUExtendedKalmanFilter::updateAcc(std::vector<Eigen::Vector3d> accMeas, double dt) {

	if (!isInitialised()) {
		if (initSamplingFinished()) {
			setInitialStateAndCovariance();
		}
		else {
			Eigen::VectorXd state;
			if (initProcedureStart()) {
				state = Eigen::VectorXd::Zero(19);
			}
			else {
				state = getState();
			}

			// Acc AVG over measurements
			// Counts the number of imu influenced in the current measurement, ignores outliers
			int num_IMUs_avg = 0;
			Eigen::Vector3d accMeasAvg = Eigen::Vector3d::Zero();

			// 1. transform
			for (int i = 0; i < accMeas.size(); i++) {
				Eigen::Vector3d transformed_accMeas = transform_Acc(accMeas[i], state.segment<3>(13), VIMU_Orientations[i], VIMU_Coords[i], state.segment<3>(16));
				if (isValidMeasurement(accelerometer, transformed_accMeas)) {
					accMeasAvg += transformed_accMeas;
					num_IMUs_avg++;
				}
			}
			// 2. AVG
			accMeasAvg /= num_IMUs_avg;

			state(6) += accMeasAvg(0);
			state(7) += accMeasAvg(1);
			state(8) += accMeasAvg(2);

			updateMeasurementCount(0);

			setState(state);
			return;
		}
	}

	if (isInitialised()) {
		Eigen::VectorXd state = getState();
		Eigen::MatrixXd covariance = getCovariance();

		// Acc AVG over measurements
		// Counts the number of imu influenced in the current measurement, ignores outliers
		int num_IMUs_avg = 0;
		Eigen::Vector3d accMeasAvg = Eigen::Vector3d::Zero();

		// 1. transform
		for (int i = 0; i < accMeas.size(); i++) {
			Eigen::Vector3d transformed_accMeas = transform_Acc(accMeas[i], state.segment<3>(13), VIMU_Orientations[i], VIMU_Coords[i], state.segment<3>(16));
			if (isValidMeasurement(accelerometer, transformed_accMeas)) {
				accMeasAvg += transformed_accMeas;
				num_IMUs_avg++;
			}
		}
		// 2. AVG
		accMeasAvg /= num_IMUs_avg;

		Eigen::Vector3d accMeas = accMeasAvg;

		//Conversion from ("body") XYD -> ("local") NED
		Eigen::Quaternion<double> orientation = Eigen::Quaternion<double>{ state(9), state(10), state(11), state(12) };
		accMeas = orientation._transformVector(accMeas);

		Eigen::Vector3d gMeas = accMeas;
		gMeas.normalize();

		// Correction for the gravity	
		accMeas(2) += g;



		// Acc Correction, Measurement noise noch hinzufuegen [R]
		Eigen::MatrixXd H = Eigen::MatrixXd::Zero(3, 19);
		Eigen::MatrixXd R = Eigen::MatrixXd::Zero(3, 3);

		double acc_var{};
		acc_var = ACCEL_STD * ACCEL_STD;
		R(0, 0) = acc_var;
		R(1, 1) = acc_var;
		R(2, 2) = acc_var;

		H.row(0) << 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
		H.row(1) << 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
		H.row(2) << 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;

		Eigen::Vector3d y = accMeas - H * state;
		Eigen::MatrixXd S = H * covariance * H.transpose() + R;
		Eigen::MatrixXd K = covariance * H.transpose() * S.inverse();

		state = state + K * y;
		covariance = (Eigen::MatrixXd::Identity(19, 19) - K * H) * covariance;


		// Orientation Correciton
		Eigen::Vector3d p_dot_dot_hat = Eigen::Vector3d::Zero();

		double q_0 = state(9);
		double q_1 = state(10);
		double q_2 = state(11);
		double q_3 = state(12);
		p_dot_dot_hat(0) = 2 * q_1 * q_3 - 2 * q_0 * q_2;
		p_dot_dot_hat(1) = 2 * q_2 * q_3 + 2 * q_0 * q_1;
		p_dot_dot_hat(2) = q_0 * q_0 - q_1 * q_1 - q_2 * q_2 + q_3 * q_3;

		p_dot_dot_hat = g * p_dot_dot_hat;
		y = gMeas - p_dot_dot_hat;

		// non linear modell, bc of direction cosine
		H.row(0) << 0, 0, 0, 0, 0, 0, 0, 0, 0, -2 * q_2, 2 * q_3, -2 * q_0, 2 * q_1, 0, 0, 0, 0, 0, 0;
		H.row(1) << 0, 0, 0, 0, 0, 0, 0, 0, 0, 2 * q_1, 2 * q_0, 2 * q_3, 2 * q_1, 0, 0, 0, 0, 0, 0;
		H.row(2) << 0, 0, 0, 0, 0, 0, 0, 0, 0, 2 * q_0, -2 * q_1, -2 * q_2, 2 * q_3, 0, 0, 0, 0, 0, 0;

		S = H * covariance * H.transpose() + R;
		K = covariance * H.transpose() * S.inverse();


		// Dont update the YAW Component with Accelerometer-Data
		Eigen::VectorXd state_error = Eigen::VectorXd::Zero(19);
		state_error = K * y;
		state_error(12) = 0;

		state = state + state_error;
		covariance = (Eigen::MatrixXd::Identity(19, 19) - K * H) * covariance;

		setState(state);
		setCovariance(covariance);
		//	std::cout << "ACC\n";

	}


};
void VIMUExtendedKalmanFilter::updateMag(std::vector<Eigen::Vector3d> magMeas, double dt) {
	if (!isInitialised()) {
		if (initSamplingFinished()) {
			setInitialStateAndCovariance();
		}
		else {
			Eigen::VectorXd state;
			if (initProcedureStart()) { // Init State
				state = Eigen::VectorXd::Zero(19);
			}
			else {
				state = getState();
			}


			// MAG AVG over measurements
			// Counts the number of imu influenced in the current measurement, ignores outliers
			int num_IMUs_avg = 0;
			Eigen::Vector3d magMeasAvg = Eigen::Vector3d::Zero();

			// 1. transform
			for (int i = 0; i < magMeas.size(); i++) {
				Eigen::Vector3d transformed_magMeas = transform_Mag(magMeas[i], VIMU_Orientations[i]);
				if (isValidMeasurement(magnetometer, transformed_magMeas)) {
					magMeasAvg += transformed_magMeas;
					num_IMUs_avg++;
				}
			}
			// 2. AVG
			magMeasAvg /= num_IMUs_avg;

			// Use uninitialisied Quaternion state for placeholder
			state(9) += magMeasAvg(0);
			state(10) += magMeasAvg(1);
			state(11) += magMeasAvg(2);

			updateMeasurementCount(1);

			setState(state);
			return;
		}
	}

	if (isInitialised()) {
		Eigen::VectorXd state = getState();
		Eigen::MatrixXd covariance = getCovariance();

		Eigen::Vector3d mag_hat = Eigen::Vector3d::Zero();
		// Richtung des Magnetfelds
		double q_0 = state(9);
		double q_1 = state(10);
		double q_2 = state(11);
		double q_3 = state(12);


		// Magnet is in north direction, only x directed 

		mag_hat(0) = q_0 * q_0 + q_1 * q_1 - q_2 * q_2 - q_3 * q_3;
		mag_hat(1) = 2 * q_1 * q_2 - 2 * q_0 * q_3;
		mag_hat(2) = 2 * q_1 * q_3 + 2 * q_0 * q_2;


		// MAG AVG over measurements
		// Counts the number of imu influenced in the current measurement, ignores outliers
		int num_IMUs_avg = 0;
		Eigen::Vector3d magMeasAvg = Eigen::Vector3d::Zero();

		// 1. transform
		for (int i = 0; i < magMeas.size(); i++) {
			Eigen::Vector3d transformed_magMeas = transform_Mag(magMeas[i], VIMU_Orientations[i]);
			if (isValidMeasurement(gyroscope, transformed_magMeas)) {
				magMeasAvg += transformed_magMeas;
				num_IMUs_avg++;
			}
		}
		// 2. AVG
		magMeasAvg /= num_IMUs_avg;

		Eigen::Quaternion<double> orientation = Eigen::Quaternion<double>{ q_0,q_1,q_2,q_3 };
		Eigen::Vector3d	magMeas = orientation._transformVector(magMeasAvg);
		magMeas.normalize();
		Eigen::Vector3d y = magMeas - mag_hat;

		Eigen::MatrixXd H = Eigen::MatrixXd::Zero(3, 19);
		Eigen::MatrixXd R = Eigen::MatrixXd::Zero(3, 3);

		R(0, 0) = 1;
		R(1, 1) = 1;
		R(2, 2) = 1;

		H.row(0) << 0, 0, 0, 0, 0, 0, 0, 0, 0, 2 * q_3, 2 * q_2, 2 * q_1, 2 * q_2, 0, 0, 0, 0, 0, 0;
		H.row(1) << 0, 0, 0, 0, 0, 0, 0, 0, 0, 2 * q_0, -2 * q_1, -2 * q_2, -2 * q_3, 0, 0, 0, 0, 0, 0;
		H.row(2) << 0, 0, 0, 0, 0, 0, 0, 0, 0, -2 * q_1, -2 * q_0, 2 * q_3, 2 * q_2, 0, 0, 0, 0, 0, 0;

		Eigen::MatrixXd S = H * covariance * H.transpose() + R;
		Eigen::MatrixXd K = covariance * H.transpose() * S.inverse();

		// Dont update the Roll and Pitch Component with Magnetometer-Data
		Eigen::VectorXd state_error = Eigen::VectorXd::Zero(19);
		state_error = K * y;
		state_error(10) = 0;
		state_error(11) = 0;

		state = state + state_error;
		covariance = (Eigen::MatrixXd::Identity(19, 19) - K * H) * covariance;

		setState(state);
		setCovariance(covariance);
		//	std::cout << "MAG\n";

	}

};
void VIMUExtendedKalmanFilter::updateGPS(std::vector<Eigen::Vector3d> gpsMeas, double dt, Eigen::Vector3d gpsVelocityInitial, Eigen::Quaternion<double> orientationInitial) {


	if (!isInitialised()) {
		if (!initSamplingFinished()) {
			Eigen::VectorXd state;
			if (initProcedureStart()) {
				state = Eigen::VectorXd::Zero(19);
			}
			else {
				state = getState();
			}
			// GPS AVG over measurements
			// Counts the number of imu influenced in the current measurement, ignores outliers
			int num_IMUs_avg = 0;
			Eigen::Vector3d gpsMeasAvg = Eigen::Vector3d::Zero();

			// 1. transform
			for (int i = 0; i < gpsMeas.size(); i++) {
				if (isValidMeasurement(gps, gpsMeas[i])) {
					gpsMeasAvg += gpsMeas[i];
					num_IMUs_avg++;
				}
			}
			// 2. AVG
			gpsMeasAvg /= num_IMUs_avg;

			state(0) += gpsMeasAvg(0);
			state(1) += gpsMeasAvg(1);
			state(2) += 0;

			setState(state);
			updateMeasurementCount(2); // Update property in Measurement Count for GPS measurement
		}
		else {
			Eigen::VectorXd state = getState();

			// GPS AVG over measurements
			// Counts the number of imu influenced in the current measurement, ignores outliers
			int num_IMUs_avg = 0;
			Eigen::Vector3d gpsMeasAvg = Eigen::Vector3d::Zero();

			// 1. transform
			for (int i = 0; i < gpsMeas.size(); i++) {
				if (isValidMeasurement(gps, gpsMeas[i])) {
					gpsMeasAvg += gpsMeas[i];
					num_IMUs_avg++;
				}
			}
			// 2. AVG
			gpsMeasAvg /= num_IMUs_avg;

			Eigen::Vector3d gpsMeas = gpsMeasAvg;


			// Mean (better Median) of the measurements
			gpsMeas(0) += state(0);
			gpsMeas(1) += state(1);
			gpsMeas(2) += 0;

			gpsMeas = gpsMeas / (getInitMeasurementCounter()(2) + 1); // mean of the init GPS Positions (11th measurement for GPS)

			setReferenceGeodeticPosition(gpsMeas);

			// Init of Orientation, State and Covariance
			setInitialStateAndCovariance();

		}
	}
	else {
		Eigen::VectorXd state = getState();
		Eigen::MatrixXd covariance = getCovariance();

		// GPS AVG over measurements
		// Counts the number of imu influenced in the current measurement, ignores outliers
		int num_IMUs_avg = 0;
		Eigen::Vector3d gpsMeasAvg = Eigen::Vector3d::Zero();

		// 1. transform
		for (int i = 0; i < gpsMeas.size(); i++) {
			if (isValidMeasurement(gps, gpsMeas[i])) {
				gpsMeasAvg += gpsMeas[i];
				num_IMUs_avg++;
			}
		}
		// 2. AVG
		gpsMeasAvg /= num_IMUs_avg;

		Eigen::Vector3d gpsMeas = gpsMeasAvg;


		// GPS in NED locally Coordinates at Position of the (V)IMU
		gpsMeas = computeECEF2NED(gpsMeas);



		Eigen::MatrixXd H = Eigen::MatrixXd::Zero(3, 19);
		Eigen::MatrixXd R = Eigen::MatrixXd::Zero(3, 3);

		double gps_var = GPS_POS_STD * GPS_POS_STD;

		R(0, 0) = gps_var;
		R(1, 1) = gps_var;
		R(2, 2) = gps_var;

		H.row(0) << 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
		H.row(1) << 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
		H.row(2) << 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;

		Eigen::Vector3d y = gpsMeas - H * state;
		Eigen::MatrixXd S = H * covariance * H.transpose() + R;
		Eigen::MatrixXd K = covariance * H.transpose() * S.inverse();

		state = state + K * y;
		covariance = (Eigen::MatrixXd::Identity(19, 19) - K * H) * covariance;

		setState(state);
		setCovariance(covariance);
		//	std::cout << "GPS\n";
	}


};


void VIMUExtendedKalmanFilter::federtatedStateAVG(std::vector<Eigen::VectorXd> states) {

	Eigen::VectorXd stateVIMU = Eigen::VectorXd::Zero(states[0].rows());

	for (Eigen::VectorXd state : states) {
		stateVIMU += state;
	}

	stateVIMU /= states.size();
	setState(stateVIMU);

};

void ExtendedKalmanFilter::setInitialStateAndCovariance() {
	Eigen::VectorXd state = getState();

	// No VIMU, INS from local IMU not centered of GPS Signal
	Eigen::Vector3d coords = getCoords();
	// set Reference to zero + the Coords of the (V)IMU in the construct
	state(0) = 0 + coords(0);
	state(1) = 0 + coords(1);
	state(2) = 0 + coords(2);

	// Init Velocity zero
	state(3) = 0;
	state(4) = 0;
	state(5) = 0;

	Eigen::Quaternion<double> initOrientation = computeInitOrientation(state);

	state(9) = initOrientation.w();
	state(10) = initOrientation.x();
	state(11) = initOrientation.y();
	state(12) = initOrientation.z();

	Eigen::Vector3d initAcc = Eigen::Vector3d{ state(6), state(7), state(8) };
	initAcc = initAcc / getInitMeasurementCounter()(0);
	// init the acceleration with the first Orientation from ("Body") XYD -> ("local") NED
	initAcc = initOrientation._transformVector(initAcc);

	state(6) = initAcc(0);
	state(7) = initAcc(1);
	state(8) = initAcc(2) + g;



	//	state = x, y, z, x_dot, y_dot, z_dot, x_dot_dot, y_dot_dot, z_dot_dot, q0, q1, q2, q3
		// Fist Init of Covariane
	Eigen::MatrixXd covariance = Eigen::MatrixXd::Zero(13, 13);

	double gps_var{}, vel_var{}, acc_var{};
	gps_var = GPS_POS_STD * GPS_POS_STD;
	vel_var = INIT_VEL_STD * INIT_VEL_STD;
	acc_var = ACCEL_STD * ACCEL_STD;


	covariance(0, 0) = gps_var;
	covariance(1, 1) = gps_var;
	covariance(2, 2) = gps_var;
	covariance(3, 3) = vel_var;
	covariance(4, 4) = vel_var;
	covariance(5, 5) = vel_var;
	covariance(6, 6) = acc_var;
	covariance(7, 7) = acc_var;
	covariance(8, 8) = acc_var;

	// Use P0 from "A double stage kalman filter for orientation and Tracking with an IMU..."
	covariance.row(9) << 0, 0, 0, 0, 0, 0, 0, 0, 0, INIT_ORIENTATION_VAR, INIT_ORIENTATION_COVAR, INIT_ORIENTATION_COVAR, INIT_ORIENTATION_COVAR;
	covariance.row(10) << 0, 0, 0, 0, 0, 0, 0, 0, 0, INIT_ORIENTATION_COVAR, INIT_ORIENTATION_VAR, INIT_ORIENTATION_COVAR, INIT_ORIENTATION_COVAR;
	covariance.row(11) << 0, 0, 0, 0, 0, 0, 0, 0, 0, INIT_ORIENTATION_COVAR, INIT_ORIENTATION_COVAR, INIT_ORIENTATION_VAR, INIT_ORIENTATION_COVAR;
	covariance.row(12) << 0, 0, 0, 0, 0, 0, 0, 0, 0, INIT_ORIENTATION_COVAR, INIT_ORIENTATION_COVAR, INIT_ORIENTATION_COVAR, INIT_ORIENTATION_VAR;

	setState(state);
	setCovariance(covariance);
	initFinished();

}

void VIMUExtendedKalmanFilter::setInitialStateAndCovariance() {
	Eigen::VectorXd state = getState();
	// set Reference Position to zero
	state(0) = 0;
	state(1) = 0;
	state(2) = 0;

	// Init Velocity zero
	state(3) = 0;
	state(4) = 0;
	state(5) = 0;

	std::cout << state << std::endl;

	Eigen::Quaternion<double> initOrientation = computeInitOrientation(state);

	state(9) = initOrientation.w();
	state(10) = initOrientation.x();
	state(11) = initOrientation.y();
	state(12) = initOrientation.z();


	Eigen::Vector3d initAcc = Eigen::Vector3d{ state(6), state(7), state(8) };
	initAcc = initAcc / getInitMeasurementCounter()(0);
	// init the acceleration with the first Orientation from ("Body") XYD -> ("local") NED
	initAcc = initOrientation._transformVector(initAcc);

	state(6) = initAcc(0);
	state(7) = initAcc(1);
	state(8) = initAcc(2) + g;

	//	state = x, y, z, x_dot, y_dot, z_dot, x_dot_dot, y_dot_dot, z_dot_dot, q0, q1, q2, q3, x_psi_dot, y_psi_dot, z_psi_dot, x_psi_dot_dot, y_psi_dot_dot, z_psi_dot_dot
// Fist Init of Covariane
	Eigen::MatrixXd covariance = Eigen::MatrixXd::Zero(19, 19);

	double gps_var{}, vel_var{}, acc_var{}, gyro_var{};
	gps_var = GPS_POS_STD * GPS_POS_STD;
	vel_var = INIT_VEL_STD * INIT_VEL_STD;
	acc_var = ACCEL_STD * ACCEL_STD;
	gyro_var = GYRO_STD * GYRO_STD;

	covariance(0, 0) = gps_var;
	covariance(1, 1) = gps_var;
	covariance(2, 2) = gps_var;
	covariance(3, 3) = vel_var;
	covariance(4, 4) = vel_var;
	covariance(5, 5) = vel_var;
	covariance(6, 6) = acc_var;
	covariance(7, 7) = acc_var;
	covariance(8, 8) = acc_var;
	covariance(13, 13) = gyro_var;
	covariance(14, 14) = gyro_var;
	covariance(15, 15) = gyro_var;
	covariance(16, 16) = gyro_var;
	covariance(17, 17) = gyro_var;
	covariance(18, 18) = gyro_var;

	// Use P0 from "A double stage kalman filter for orientation and Tracking with an IMU..."
	covariance.row(9) << 0, 0, 0, 0, 0, 0, 0, 0, 0, INIT_ORIENTATION_VAR, INIT_ORIENTATION_COVAR, INIT_ORIENTATION_COVAR, INIT_ORIENTATION_COVAR, 0, 0, 0, 0, 0, 0;
	covariance.row(10) << 0, 0, 0, 0, 0, 0, 0, 0, 0, INIT_ORIENTATION_COVAR, INIT_ORIENTATION_VAR, INIT_ORIENTATION_COVAR, INIT_ORIENTATION_COVAR, 0, 0, 0, 0, 0, 0;
	covariance.row(11) << 0, 0, 0, 0, 0, 0, 0, 0, 0, INIT_ORIENTATION_COVAR, INIT_ORIENTATION_COVAR, INIT_ORIENTATION_VAR, INIT_ORIENTATION_COVAR, 0, 0, 0, 0, 0, 0;
	covariance.row(12) << 0, 0, 0, 0, 0, 0, 0, 0, 0, INIT_ORIENTATION_COVAR, INIT_ORIENTATION_COVAR, INIT_ORIENTATION_COVAR, INIT_ORIENTATION_VAR, 0, 0, 0, 0, 0, 0;

	std::cout << state << std::endl;

	setState(state);
	setCovariance(covariance);
	initFinished();
}


Eigen::Quaternion<double> ExtendedKalmanFilterBase::computeInitOrientation(Eigen::VectorXd state) {
	// Only is valid, if the sensor is stationary!!!
				// Mean over acceleration init measurements
	Eigen::Vector3d initAcc = Eigen::Vector3d{ state(6), state(7), state(8) };
	initAcc = initAcc / getInitMeasurementCounter()(0);

	// Determine Pitch and Roll from Accelerometer -> Elevation angle, has to be nearly zero
	// compute pitch
	double pitch = std::asin(initAcc(0) / g);

	// compute roll -> turn/ bank angle, has also to be nearly zero
	double param = initAcc(1) / initAcc(2);

	double roll = std::atan(param);//m_y,m_z -> atan does the job [-90,90] is adequat

	// Mean over magnetometer init measurements [with "Angle of Magnetic Declination" - Angle between magnetic and true north]
	Eigen::Vector3d initMag = Eigen::Vector3d{ state(9), state(10), state(11) };
	initMag = initMag / getInitMeasurementCounter()(1);

	// Determine Yaw from magnetometer
	double yaw = std::atan2(initMag(1), initMag(0)); //m_y/m_x ich

	yaw += declinationAngle;

	return computeOrientation(yaw, pitch, roll);
}

// sequence is ZYX -> Yaw, Pitch, Roll : Is it adaequat?? 1 of 12 possibilities
// Quaternion & Rotation Sequences p.167
Eigen::Quaternion<double> ExtendedKalmanFilterBase::computeOrientation(double yaw, double pitch, double roll) {


	double y_c = std::cos(yaw / 2);
	double y_s = std::sin(yaw / 2);
	double p_c = std::cos(pitch / 2);
	double p_s = std::sin(pitch / 2);
	double r_c = std::cos(roll / 2);
	double r_s = std::sin(roll / 2);

	Eigen::Quaternion<double> q = Eigen::Quaternion<double>::Identity();
	q.w() = y_c * p_c * r_c + y_s * p_s * r_s;
	q.x() = y_c * p_c * r_s - y_s * p_s * r_c;
	q.y() = y_c * p_s * r_c + y_s * p_c * r_s;
	q.z() = y_s * p_c * r_c - y_c * p_s * r_s;

	q.normalize();

	return q;
}

// Mahalanobis Distance check for incoming Measurement, if its to wide outlying -> not what we expect, primaryly GPS
// now outliers just, when it is too unrealistic
bool ExtendedKalmanFilterBase::isValidMeasurement(Sensortype sensor, Eigen::Vector3d meas) {
	switch (sensor)
	{
	case gyroscope:
		if (meas.norm() >= 0.52) {//  30 Graddrehung in einer Achse
			return false;
		}
		else {
			return true;
		}
	case accelerometer:
		if (meas.norm() >= 20) {
			return false;
		}
		else {
			return true;
		}
	case magnetometer:
		return true;
	case gps:
		if (abs(meas(0) - ref_lat) > 1) {
			return false;
		}
		else if (abs(meas(1) - ref_long) > 1) {
			return false;
		}
		else {
			return true;
		}
	};
}

bool ExtendedKalmanFilterBase::initSamplingFinished() {

	for (int i = 0; i < initMeasurementCounter.size(); i++) {
		if (initMeasurementCounter(i) < minInitMeasurements) {
			return false;
		}
	}
	return true;

}

bool ExtendedKalmanFilterBase::initProcedureStart() {
	for (int i = 0; i < initMeasurementCounter.size(); i++) {
		if (initMeasurementCounter(i) != 0) {
			return false;
		}
	}
	return true;
}

