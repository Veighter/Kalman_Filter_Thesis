#include "extendedkalmanfilter.h"

// For computing the Fehlerfortpflanzungsmatrix
constexpr double ACCEL_STD = 1.0;
constexpr double GYRO_STD = 0.01 / 180.0 * M_PI;
constexpr double INIT_VEL_STD = 10.0;
constexpr double INIT_PSI_STD = 45.0 / 180.0 * M_PI;
constexpr double GPS_POS_STD = 3.0;


void ExtendedKalmanFilter::predictionStep(double dt) {
	if (isInitialised()) {
	

	}
}

// ist das noch ein Prediction Step oder ist das schon ein Update?? Evtl Funktion anpassen wie geht das hier?????
void ExtendedKalmanFilter::predictionStep(Eigen::Vector3d accMeas, Eigen::Vector3d gyroMeas, Eigen::Vector3d magMeas, double dt) {
	if (!isInitialised()) {
		Eigen::VectorXd state = Eigen::VectorXd::Zero(19);

		Eigen::Quaternion<double> orientation = getOrientation();

		state(0) = 0;
		state(1) = 0;
		state(2) = 0;
		state(3) = 0;
		state(4) = 0;
		state(5) = 0;
		state(6) = accMeas(0);
		state(7) = accMeas(1);
		state(8) = accMeas(2);
		state(9) = gyroMeas(0);
		state(10) = gyroMeas(1);
		state(11) = gyroMeas(2);
		state(12) = orientation.w;
		state(13) = orientation.x;
		state(14) = orientation.y;
		state(15) = orientation.z;
		state(16) = magMeas(0);
		state(17) = magMeas(1);
		state(18) = magMeas(2);
	}
	else {
		Eigen::VectorXd state = getState();
		Eigen::MatrixXd covariance = getCovariance();

		double x = state(0);
		double y = state(1);
		double z = state(2);
		double x_dot = state(3);
		double y_dot = state(4);
		double z_dot = state(5);
		double x_dot_dot = accMeas(0);
		double y_dot_dot = accMeas(1);
		double z_dot_dot = accMeas(2);
		double psi_x_dot = gyroMeas(0);
		double psi_y_dot = gyroMeas(1);
		double psi_z_dot = gyroMeas(2);
		double q_0 = state(12);
		double q_1 = state(13);
		double q_2 = state(14);
		double q_3 = state(15);
		double B_x = magMeas(0);
		double B_y = masMeas(1);
		double B_z = magMeas(2);

		// Predict position
		state(0) = x + dt * x_dot + 1. / 2 * dt * dt * x_dot_dot;
		state(1) = y + dt * y_dot + 1. / 2 * dt * dt * y_dot_dot;
		state(2) = z + dt * z_dot + 1. / 2 * dt * dt * z_dot_dot;

		// Predict velocity
		state(3) = x_dot + dt * x_dot_dot;
		state(4) = y_dot + dt * y_dot_dot;
		state(5) = z_dot + dt * z_dot_dot;

		state(6) = x_dot_dot;
		state(7) = y_dot_dot;
		state(8) = z_dot_dot;

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

		state(16) = B_x;
		state(17) = B_y;
		state(18) = B_z;
		
		setState(state);
	}
}

// die Initialisierung muss hier drin geschehen, da ich die Position brauche!!
void ExtendedKalmanFilter::updateStep(sensorMeas::GPSMeas gpsMeas, double dt) {
}