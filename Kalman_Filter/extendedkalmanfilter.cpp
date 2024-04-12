#include "extendedkalmanfilter.h"
#include <iostream>

// For computing the Fehlerfortpflanzungsmatrix
constexpr double ACCEL_STD = 1.0;
constexpr double GYRO_STD = 0.1 / 180.0 * M_PI;
constexpr double MAG_STD = 1.0;
constexpr double GPS_POS_STD =3.0;

// Initial standard derivation has to be high, bc we dont know exactly what they are
constexpr double INIT_VEL_STD = 10.0;
constexpr double INIT_ORIENTATION_VAR = 0.125;
constexpr double INIT_ORIENTATION_COVAR = 0.0003;

constexpr double g = 9.81;

// magnetic-declination.com: Location is Halter with 3 degree and 11 sec = 3.1833 degree
constexpr double declinationAngle = 3.1833 * M_PI / 180.0;



// Prediciton Step with implicit Gyro inputs, decreases number of states https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=6316172
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

		Eigen::Quaternion<double> orientation = Eigen::Quaternion<double>{ q_0,q_1,q_2,q_3 };
	
		// Conversion from ("body") XYD -> ("local") NED
		gyroMeas = orientation._transformVector(gyroMeas);

		double psi_x_dot = gyroMeas(0);
		double psi_y_dot = gyroMeas(1);
		double psi_z_dot = gyroMeas(2);
	


		//// Predict position
		//state(0) = x + dt * x_dot + 1. / 2 * dt * dt * x_dot_dot;
		//state(1) = y + dt * y_dot + 1. / 2 * dt * dt * y_dot_dot;
		//state(2) = z + dt * z_dot + 1. / 2 * dt * dt * z_dot_dot;

		//// Predict velocity
		//state(3) = x_dot + dt * x_dot_dot;
		//state(4) = y_dot + dt * y_dot_dot;
		//state(5) = z_dot + dt * z_dot_dot;

		//// Update Acceleration with measurements
		//state(6) = x_dot_dot;
		//state(7) = y_dot_dot;
		//state(8) = z_dot_dot;

		//// Update Angular velocity
		//state(9) = psi_x_dot;
		//state(10) = psi_y_dot;
		//state(11) = psi_z_dot;

		//// Predict Orientation
		//Eigen::Quaternion<double> q = Eigen::Quaternion<double>{ q_0 + 1. / 2 * (-psi_x_dot * q_1 - psi_y_dot * q_2 - psi_z_dot * q_3) * dt, q_1 + 1. / 2 * (psi_x_dot * q_0 + psi_z_dot * q_2 - psi_y_dot * q_3) * dt, q_2 + 1. / 2 * (psi_y_dot * q_0 - psi_z_dot * q_1 + psi_x_dot * q_3) * dt, q_3 + 1. / 2 * (psi_z_dot * q_0 + psi_y_dot * q_1 - psi_x_dot * q_2) * dt
		//};
		//q.normalize();
		//state(12) = q.w();
		//state(13) = q.x();
		//state(14) = q.y();
		//state(15) = q.z();


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
		std::cout << "GYRO\n";
	}
	
}

// ist das noch ein Prediction Step oder ist das schon ein Update?? Evtl Funktion anpassen wie geht das hier?????


// Das was geschaetzt werden soll ist eigentlich nur Position, Geschwindigkeit und Orientierung des zu navigierenden Objekts. 
// Anpassen des Zustandsvektors ist also von Noeten!!
void ExtendedKalmanFilter::updateAcc(Eigen::Vector3d accMeas, double dt) {

	if (isInitialised()) {
		Eigen::VectorXd state = getState();
		Eigen::MatrixXd covariance = getCovariance();
		
		//Conversion from ("body") XYD -> ("local") NED
		Eigen::Quaternion<double> orientation = Eigen::Quaternion<double>{ state(9), state(10), state(11), state(12) };
		accMeas = orientation._transformVector(accMeas);

		Eigen::Vector3d gMeas = accMeas;
		gMeas.normalize();

		// Correction for the gravity	
		accMeas(2) += g;
		
		

		// Acc Correction, Measurement noise noch hinzufuegen [R]
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


		/*
		* Rechne ich nun so weiter oder kann ich das auch zu einem nichtlinearen Modell machen?
		*/

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
		std::cout << "ACC\n";

	}
	else {
		Eigen::Vector3i measCount = getInitMeasurementCount();
		Eigen::VectorXd state;
		if (measCount(0) == 0 && measCount(1) == 0 && measCount(2) == 0) {
			state = Eigen::VectorXd::Zero(13);
		}
		else {
			state = getState();
		}
		state(6) += accMeas(0);
		state(7) += accMeas(1);
		state(8) += accMeas(2);



		updateMeasurementCount(0);

		setState(state);
	}
	
}

void ExtendedKalmanFilter::updateMag(Eigen::Vector3d magMeas, double dt) {

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

		// Correction seems a bit of, mag_hat has to be in x-axis of the ned! But what does that mean for now?
		mag_hat(0) = q_0 * q_0 + q_1 * q_1 - q_2 * q_2 - q_3 * q_3;
		mag_hat(1) = 2 * q_1 * q_2 - 2 * q_0 * q_3;
		mag_hat(2) = 2 * q_1 * q_3 + 2 * q_0 * q_2; 

		// where do i put the declination angle? ("Angle between magnetic and true north")

		Eigen::Quaternion<double> orientation = Eigen::Quaternion<double>{ q_0,q_1,q_2,q_3 };
		magMeas = orientation._transformVector(magMeas);
		magMeas.normalize();
		Eigen::Vector3d y = magMeas - mag_hat;

		Eigen::MatrixXd H = Eigen::MatrixXd::Zero(3, 13);
		Eigen::MatrixXd R = Eigen::MatrixXd::Zero(3, 3);

		R(0, 0) = 1;
		R(1, 1) = 1;
		R(2, 2) = 1;

		H.row(0) << 0, 0, 0, 0, 0, 0, 0, 0, 0, 2 * q_3, 2 * q_2, 2 * q_1, 2 * q_2;
		H.row(1) << 0, 0, 0, 0, 0, 0, 0, 0, 0, 2 * q_0, -2 * q_1, -2 * q_2, -2 * q_3;
		H.row(2) << 0, 0, 0, 0, 0, 0, 0, 0, 0, -2 * q_1, -2 * q_0, 2 * q_3, 2 * q_2;

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
		std::cout << "MAG\n";

	}
	else {

		Eigen::Vector3i measCount = getInitMeasurementCount();
		Eigen::VectorXd state;
		if (measCount(0) == 0 && measCount(1) == 0 && measCount(2) == 0) { // Init State
			state = Eigen::VectorXd::Zero(13);
		}
		else {
			state = getState();
		}// Use uninitialisied Quaternion state for placeholder
		state(9) += magMeas(0);
		state(10) += magMeas(1);
		state(11) += magMeas(2);

		updateMeasurementCount(1);

		setState(state);
	}
	
}

// die Initialisierung muss hier drin geschehen, da ich die Position brauche!!
// Anfagnsorientierung, wie bestimme ich diese, muss ja wissen wo mein NED hinzeigt (Aus Accelerometer und Magnetometerwert grob bestimmen wie 
// oben in den Update schritten, einmal 10 Measurements den Durchschnitt nehmen fuer eine Halbwegs gute Initialisierung

/// <summary>
/// Method for the Measurement Update of the navigational Object state with the GPS 
/// </summary>
/// <param name="gpsMeas"></param>GPS given in LLA
/// <param name="dt"></param>
void ExtendedKalmanFilter::updateGPS(Eigen::Vector3d gpsMeas, double dt, Eigen::Vector3d gpsVelocityInitial, Eigen::Quaternion<double> orientationInitial) {

	if (!isInitialised()) {
		Eigen::Vector3i initMeasurementCount = getInitMeasurementCount(); // ist das hier auch mit 0 initialisiert??
		uint8_t initMeasurements = getMinInitMeasurementCount();

		if (initMeasurementCount(0) < initMeasurements || initMeasurementCount(1) < initMeasurements || initMeasurementCount(2) < initMeasurements) {
			Eigen::VectorXd state;
			if (initMeasurementCount(0) == 0 && initMeasurementCount(1) == 0 && initMeasurementCount(2) == 0) {
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

			gpsMeas = gpsMeas / (getInitMeasurementCount()(2)+1); // mean of the init GPS Positions (11th measurement for GPS)

			setReferenceGeodeticPosition(gpsMeas);

			// set Reference to zero
			state(0) = 0;
			state(1) = 0;
			state(2) = 0;

			if (false) { // for velocity when added
				state(3) = gpsVelocityInitial(0);
				state(4) = gpsVelocityInitial(1);
				state(5) = gpsVelocityInitial(2);
			}
			else {
				state(3) = 0;
				state(4) = 0;
				state(5) = 0;
			}
			// Only is valid, if the sensor is stationary!!!
			// Mean over acceleration init measurements
			Eigen::Vector3d initAcc = Eigen::Vector3d{ state(6), state(7), state(8) };
			initAcc = initAcc / initMeasurementCount(0);

			// Determine Pitch and Roll from Accelerometer -> Elevation angle, has to be nearly zero
			// compute pitch
			double pitch = std::asin(initAcc(0) / g);

			// compute roll -> turn/ bank angle, has also to be nearly zero
			double param = initAcc(1) / initAcc(2);
		
			double roll = std::atan(param);//m_y,m_z -> atan does the job [-90,90] is adequat

			// Mean over magnetometer init measurements [with "Angle of Magnetic Declination" - Angle between magnetic and true north]
			Eigen::Vector3d initMag = Eigen::Vector3d{ state(9), state(10), state(11) };
			initMag = initMag / initMeasurementCount(1);

			// Determine Yaw from magnetometer
			double yaw = std::atan2(initMag(1), initMag(0)); //m_y/m_x ich
	
			
			yaw += declinationAngle;

			Eigen::Quaternion<double> initOrientation = computeOrientation(yaw, pitch, roll);

			state(9) = initOrientation.w();
			state(10) = initOrientation.x();
			state(11) = initOrientation.y();
			state(12) = initOrientation.z();
			
			// init the acceleration with the first Orientation from ("Body") XYD -> ("local") NED
			initAcc = initOrientation._transformVector(initAcc);

			state(6) = initAcc(0);	
			state(7) = initAcc(1);
			state(8) = initAcc(2)+g;

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
	}
	else {
		Eigen::VectorXd state = getState();
		Eigen::MatrixXd covariance = getCovariance();

		// GPS in NED locally Coordinates
		gpsMeas = computeECEF2NED(gpsMeas);



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
		std::cout << "GPS\n";
	}

}
// sequence is ZYX -> Yaw, Pitch, Roll : Is it adaequat?? 1 of 12 possibilities
// Quaternion & Rotation Sequences p.167
Eigen::Quaternion<double> ExtendedKalmanFilter::computeOrientation(double yaw, double pitch, double roll) {
	/*Eigen::Matrix3d roll = Eigen::Matrix3d::Zero();
	roll << 1, 0, 0, 0, std::cos(roll), -std::sin(roll), 0, std::sin(roll), std::cos(roll);

	Eigen::Matrix3d pitch = Eigen::Matrix3d::Zero();
	pitch << std::cos(pitch), 0, std::sin(pitch), 0, 1, 0, -std::sin(pitch), 0, std::cos(pitch);

	Eigen::Matrix3d yaw = Eigen::Matrix3d::Zero();
	yaw << std::cos(yaw), -std::sin(yaw), 0, std::sin(yaw), std::cos(yaw), 0, 0, 0, 1;*/


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

// Maha.... Distance check for incoming Measurement, if its to wide outlying -> not what we expect, primaryly GPS
bool ExtendedKalmanFilter::isValid(Sensortype sensor) {
	return false;
}


















//F(0, 0) = 1;
//F(0, 3) = dt;
//F(0, 6)1. / 2 * dt * dt * x_dot;
//F(1, 1) = 1;
//F(1, 4) = dt;
//F(1, 7) = 1. / 2 * dt * dt * y_dot_dot;
//F(2, 2) = 1;
//F(2, 5) = dt;
//F(2, 8) = 1. / 2 * dt * dt * z_dot_dot;
//
//
//
//F(3, 3) = 1;
//F(3, 6) = dt * x_dot_dot;
//F(4, 4) = 1;
//F(4, 7) = dt * y_dot_dot;
//F(5, 5) = 1;
//F(5, 8) = dt * z_dot_dot;
//
//
//
//F(6, 6) = x_dot_dot;
//F(7, 7) = y_dot_dot;
//F(8, 8) = z_dot_dot;
//
//
//F(9, 9) = psi_x_dot;
//F(10, 10) = psi_y_dot;
//F(11, 11) = psi_z_dot;
//
//
//
//F(12, 12);
//F(12, 13);
//F(12, 14);
//F(12, 15);
//F(13, 12);
//F(13, 13);
//F(13, 14);
//F(13, 15);
//F(14, 12);
//F(14, 13);
//F(14, 14);
//F(14, 15);
//F(15, 12);
//F(15, 13);
//F(15, 14);
//F(15, 15);
//
//F(16, 16);
//F(17, 17);
//F(18, 18);
//
//state = F * state;