#include "extendedkalmanfilter.h"

// For computing the Fehlerfortpflanzungsmatrix
constexpr double ACCEL_STD = 1.0;
constexpr double GYRO_STD = 0.01 / 180.0 * M_PI;
constexpr double INIT_VEL_STD = 10.0;
constexpr double INIT_PSI_STD = 45.0 / 180.0 * M_PI;
constexpr double GPS_POS_STD = 3.0;

// state = p, p_dot, p_dot_dot, q

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
		double psi_x_dot = gyroMeas(0);
		double psi_y_dot = gyroMeas(1);
		double psi_z_dot = gyroMeas(2);
		double q_0 = state(9);
		double q_1 = state(10);
		double q_2 = state(11);
		double q_3 = state(12);


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
		F.row(8) << 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1. / 2 * dt * psi_x_dot, -1. / 2 * dt * psi_y_dot, -1. / 2 * dt * psi_z_dot;
		F.row(9) << 0, 0, 0, 0, 0, 0, 0, 0, 0, 1. / 2 * dt * psi_x_dot, 1, 1. / 2 * dt * psi_z_dot, -1. / 2 * dt * psi_y_dot;
		F.row(10) << 0, 0, 0, 0, 0, 0, 0, 0, 0, 1. / 2 * dt * psi_y_dot, -1. / 2 * dt * psi_z_dot, 1, 1. / 2 * dt * psi_x_dot;
		F.row(11) << 0, 0, 0, 0, 0, 0, 0, 0, 0, 1. / 2 * dt * psi_z_dot, 1. / 2 * dt * psi_y_dot, -1. / 2 * dt * psi_x_dot, 1;

		state = F * state;
		covariance = F * covariance * F.trace();

		Eigen::Quaternion<double> q = Eigen::Quaternion<double>{ state(9), state(10), state(11), state(12) };
		q.normalize();
		state(9) = q.w();
		state(10) = q.x();
		state(11) = q.y();
		state(12) = q.z();

		setState(state);

	}
}

// ist das noch ein Prediction Step oder ist das schon ein Update?? Evtl Funktion anpassen wie geht das hier?????


// Das was geschaetzt werden soll ist eigentlich nur Position, Geschwindigkeit und Orientierung des zu navigierenden Objekts. 
// Anpassen des Zustandsvektors ist also von Noeten!!
void ExtendedKalmanFilter::updateAcc(Eigen::Vector3d accMeas, double dt) {

	if (isInitialised()) {
		Eigen::VectorXd state = getState();
		Eigen::MatrixXd covariance = getCovariance();

		Eigen::Vector3d gMeas = accMeas;
		gMeas.normalize();

		// Acc Correction, Measurement noise noch hinzufuegen [R]
		Eigen::MatrixXd H = Eigen::MatrixXd::Zero(3, 13);
		Eigen::MatrixXd R = Eigen::MatrixXd::Zero(3, 3);
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
		Eigen::Vector3d p_dot_dot_hat = Eigen::Vector3d::Zero(3);

		double q_0 = state(9);
		double q_1 = state(10);
		double q_2 = state(11);
		double q_3 = state(12);
		p_dot_dot_hat(0) = 2 * q_1 * q_3 - 2 * q_0 * q_2;
		p_dot_dot_hat(1) = 2 * q_2 * q_3 + 2 * q_0 * q_1;
		p_dot_dot_hat(2) = q_0 * q_0 - q_1 * q_1 - q_2 * q_2 + q_3 * q_3;

		p_dot_dot_hat = 9.81 * p_dot_dot_hat;
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

	}
}

void ExtendedKalmanFilter::updateMag(Eigen::Vector3d magMeas, double dt) {

	if (isInitialised()) {
		Eigen::VectorXd state = getState();
		Eigen::MatrixXd covariance = getCovariance();

	}
}

// die Initialisierung muss hier drin geschehen, da ich die Position brauche!!
// Anfagnsorientierung, wie bestimme ich diese, muss ja wissen wo mein NED hinzeigt (Aus Accelerometer und Magnetometerwert grob bestimmen wie 
// oben in den Update schritten

/// <summary>
/// Method for the Measurement Update of the navigational Object state with the GPS 
/// </summary>
/// <param name="gpsMeas"></param>GPS given in LLA
/// <param name="dt"></param>
void ExtendedKalmanFilter::updateGPS(Eigen::Vector3d gpsMeas, double dt, Eigen::Vector3d gpsVelocityInitial = Eigen::Vector3d(), Eigen::Quaternion<double> orientationInitial = Eigen::Quaternion<double>{ 0,0,0,0 }) {

	if (!isInitialised()) {
		Eigen::VectorXd state = Eigen::VectorXd::Zero(13);
		Eigen::MatrixXd covariance = Eigen::MatrixXd::Zero(13, 13);


		state(0) = gpsMeas(0); // -> Transform to NED Coordinates needed
		state(1) = gpsMeas(1);
		state(2) = 0;

		if (gpsVelocityInitial.size() == 3) {
			state(3) = gpsVelocityInitial(0);
			state(4) = gpsVelocityInitial(1);
			state(5) = gpsVelocityInitial(2);
		}
		else {
			state(3) = 0;
			state(4) = 0;
			state(5) = 0;
		}
		
		state(6) = 0;		// acceleration
		state(7) = 0;
		state(8) = 0;

		// orientation
		if (orientationInitial.norm() > 0) { // >0 reicht, wegen rundungsfehlern kann die Norm auch groesser 1 sein, obwohl das normalisiert uebernommen werden soltle
			state(9) = orientationInitial.w();
			state(10) = orientationInitial.x();
			state(11) = orientationInitial.y();
			state(12) = orientationInitial.z();

		}
		else {
			state(9) = 0;					
			state(10) = 0;
			state(11) = 0;
			state(12) = 0;
		}
		





		setState(state);
		setCovariance(covariance);
	}

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