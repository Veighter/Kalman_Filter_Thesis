#include "extendedkalmanfilter.h"

// For computing the Fehlerfortpflanzungsmatrix
constexpr double ACCEL_STD = 1.0;
constexpr double GYRO_STD = 0.01 / 180.0 * M_PI;
constexpr double INIT_VEL_STD = 10.0;
constexpr double INIT_PSI_STD = 45.0 / 180.0 * M_PI;
constexpr double GPS_POS_STD = 3.0;


void ExtendedKalmanFilter::predictionStep(double dt) {
	if (isInitialised()) {
		Eigen::VectorXd state = getState();
		Eigen::MatrixXd covariance = getCovariance();

		double x = state(0);
		double y = state(1);
		double z = state(2);
		double x_dot = state(3);
		double y_dot = state(4);
		double z_dot = state(5);
		double x_dot_dot = state(6);
		double y_dot_dot = state(7);
		double z_dot_dot = state(8);
		double psi_x_dot = state(9);
		double psi_y_dot = state(10);
		double psi_z_dot = state(11);
		double q_0 = state(12);
		double q_1 = state(13);
		double q_2 = state(14);
		double q_3 = state(15);


		// Predict position
		state(0) = x + dt * x_dot + 1. / 2 * dt * dt * x_dot_dot;
		state(1) = y + dt * y_dot + 1. / 2 * dt * dt * y_dot_dot;
		state(2) = z + dt * z_dot + 1. / 2 * dt * dt * z_dot_dot;

		// Predict velocity
		state(3) = x_dot + dt * x_dot_dot;
		state(4) = y_dot + dt * y_dot_dot;
		state(5) = z_dot + dt * z_dot_dot;

		// Update Acceleration with measurements
		state(6) = x_dot_dot;
		state(7) = y_dot_dot;
		state(8) = z_dot_dot;

		// Update Angular velocity
		state(9) = psi_x_dot;
		state(10) = psi_y_dot;
		state(11) = psi_z_dot;

		// Predict Orientation
		Eigen::Quaternion<double> q = Eigen::Quaternion<double>{ q_0 + 1. / 2 * (-psi_x_dot * q_1 - psi_y_dot * q_2 - psi_z_dot * q_3) * dt, q_1 + 1. / 2 * (psi_x_dot * q_0 + psi_z_dot * q_2 - psi_y_dot * q_3) * dt, q_2 + 1. / 2 * (psi_y_dot * q_0 - psi_z_dot * q_1 + psi_x_dot * q_3) * dt, q_3 + 1. / 2 * (psi_z_dot * q_0 + psi_y_dot * q_1 - psi_x_dot * q_2) * dt
		};
		q.normalize();
		state(12) = q.w();
		state(13) = q.x();
		state(14) = q.y();
		state(15) = q.z();


		// state = x,y,z,x_dot,y_dot,z_dot,x_dot_dot,y_dot_dot,z_dot_dot,psi_x_dot,psi_y_dot,q0,q1,q2,q3,B_x,B_y,B_z
		Eigen::MatrixXd F = Eigen::MatrixXd::Zero(16, 16);
		F.row(0) << 1, 0, 0, dt, 0, 0, 1. / 2 * dt * dt, 0, 0, 0, 0, 0, 0, 0, 0, 0;
		F.row(1) << 0, 1, 0, 0, dt, 0, 0, 1. / 2 * dt * dt, 0, 0, 0, 0, 0, 0, 0, 0;
		F.row(2) << 0, 0, 1, 0, 0, dt, 0, 0, 1. / 2 * dt * dt, 0, 0, 0, 0, 0, 0, 0;
		F.row(3) << 0, 0, 0, 1, 0, 0, dt, 0, 0, 0, 0, 0, 0, 0, 0, 0;
		F.row(4) << 0, 0, 0, 0, 1, 0, 0, dt, 0, 0, 0, 0, 0, 0, 0, 0;
		F.row(5) << 0, 0, 0, 0, 0, 1, 0, 0, dt, 0, 0, 0, 0, 0, 0, 0;
		F.row(6) << 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0;
		F.row(7) << 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0;
		F.row(8) << 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0;
		F.row(9) << 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0;
		F.row(10) << 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0;
		F.row(11) << 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;
		F.row(12) << 0, 0, 0, 0, 0, 0, 0, 0, 0, -1. / 2 * dt * q.x(), -1. / 2 * dt * q.y(), -1. / 2 * dt * q.z(), 1, -1. / 2 * dt * psi_x_dot, -1. / 2 * dt * psi_y_dot, -1. / 2 * dt * psi_z_dot;
		F.row(13) << 0, 0, 0, 0, 0, 0, 0, 0, 0, 1. / 2 * dt * q.w(), 1. / 2 * dt * q.y(), -1. / 2 * dt * q.z(), 1. / 2 * dt * psi_x_dot, 1, 1. / 2 * dt * psi_z_dot, -1. / 2 * dt * psi_y_dot;
		F.row(14) << 0, 0, 0, 0, 0, 0, 0, 0, 0, 1. / 2 * dt * q.z(), 1. / 2 * dt * q.w(), -1. / 2 * dt * q.x(), 1. / 2 * dt * psi_y_dot, -1. / 2 * dt * psi_z_dot, 1. / 2 * dt * psi_x_dot;
		F.row(15) << 0, 0, 0, 0, 0, 0, 0, 0, 0, -1. / 2 * dt * q.y(), 1. / 2 * dt * q.x(), 1. / 2 * dt * q.w(), 1. / 2 * dt * psi_z_dot, 1. / 2 * dt * psi_y_dot, -1. / 2 * dt * psi_x_dot, 1;

		setState(state);

	}
}

// ist das noch ein Prediction Step oder ist das schon ein Update?? Evtl Funktion anpassen wie geht das hier?????


// Das was geschaetzt werden soll ist eigentlich nur Position, Geschwindigkeit und Orientierung des zu navigierenden Objekts. 
// Anpassen des Zustandsvektors ist also von Noeten!!
void ExtendedKalmanFilter::predictionStep(Eigen::Vector3d accMeas, Eigen::Vector3d gyroMeas, Eigen::Vector3d magMeas, double dt) {
	
	if(isInitialised()) {
		Eigen::VectorXd state = getState();
		Eigen::MatrixXd covariance = getCovariance();

	}
}

void ExtendedKalmanFilter::predictionStep(Eigen::Vector3d gyroMeas, double dt) {

	if (isInitialised()) {
		Eigen::VectorXd state = getState();
		Eigen::MatrixXd covariance = getCovariance();

	}
}

void ExtendedKalmanFilter::predictionStep(Eigen::Vector3d magMeas, double dt) {

	if (isInitialised()) {
		Eigen::VectorXd state = getState();
		Eigen::MatrixXd covariance = getCovariance();

	}
}

// die Initialisierung muss hier drin geschehen, da ich die Position brauche!!
// Anfagnsorientierung, wie bestimme ich diese, muss ja wissen wo mein NED hinzeigt
void ExtendedKalmanFilter::updateStep(Eigen::Vector3d gpsMeas, double dt) {

	if (!isInitialised()) {
		Eigen::VectorXd state = Eigen::VectorXd::Zero(16);
		Eigen::MatrixXd covariance = Eigen::MatrixXd::Zero(16, 16);


		state(0) = gpsMeas(0); // -> Transform to NED Coordinates needed
		state(1) = gpsMeas(1);
		state(2) = gpsMeas(2);
		state(3) = 0;
		state(4) = 0;
		state(5) = 0;
		state(6) = 0;
		state(7) = 0;
		state(8) = 0;
		state(9) = 0; 
		state(10) = 0;
		state(11) = 0;
		state(12) = 0;
		state(13) = 0;
		state(14) = 0;
		state(15) = 0;

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