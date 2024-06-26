﻿#define _USE_MATH_DEFINES

#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <math.h>
#include "sensors.h"
#include "extendedkalmanfilter.h"
#include "Configuration.h"

// 8851 Daten sind da drin vorhanden
// In welcher Einheit sind die Sekunden angegeben? Laut Abgabe in Microsekunden
// alle 2 sekunden kommt ein GPS rein
using namespace sensorMeas;
int IMU_DATA_ROWS = 0;
int GPS_DATA_ROWS = 0;
constexpr int num_GPSs = 3;
constexpr size_t num_IMUs = 1;
constexpr double g = 9.81;

double time_constant = 1e6;


struct INS_state {
	Eigen::VectorXd state;
	Eigen::MatrixXd covariance;
}fusion_center;

struct IMU_Data {
	Eigen::Vector3d accelMeas{};
	Eigen::Vector3d magMeas{};
	Eigen::Vector3d gyroMeas{};
};




/// <summary>
/// The struct INS contains the following members (called INS, because the struct can instanciate an ekf in the second fusion domain (estimation fusion))
/// imu_port: connected imu_port of the imu connected to the ins
/// timeData: vector containing the time stamps of the data sampled
/// accMeas: vector containing the Accelerometer Measurements in IMU reference Frame XYZ at position of IMU (mg/s^2)
/// gyroMeas: vector containing the Gyroskop Measurements in IMU reference Frame XYZ at position of IMU (rad/s)
/// magMeas: vector containing the Magnetometer Measurements in IMU reference Frame XYZ at position of IMU (mgauss)
/// ekf: Extended Kalman Filter for the INS 
/// calib: calibration parameters (Vector and Matrix) for the measurement vectors according to their sensor error models (Ax-b resp. A(x-b))
/// </summary>
struct INS
{
	int imu_port{}; // connected port on the multiplexer

	// IMU specific Data
	std::vector<double> timeDataIMU{ std::vector<double>() };
	std::vector<IMU_Data> imuData{ std::vector<IMU_Data>() };

	// GPS specific Data
	std::vector<double> timeDataGPS{ std::vector<double>() };
	std::vector<Eigen::Vector3d> GPSData{ std::vector<Eigen::Vector3d>() };

	ExtendedKalmanFilter ekf{};

} ins_1, ins_2, ins_3, ins_4;


/// <summary>
/// The struct represents the Center Of Mass INS systems [VIMU] combining all 4 INS of the used IMUs in one INS
/// inss: INS-Array of the INS to the given IMUs in the 4 positions connected to the construct
/// ekf: Extended Kalman-Filter in the center of mass for interpreting the (calibrated) and averaged raw sensor data in the Center Of Mass
/// </summary>
struct CM_INS {
	std::vector<double> timeDataIMU{ std::vector<double>() };
	std::vector<IMU_Data> centralized_imuData{ std::vector<IMU_Data>() };

	std::vector<double> timeDataGPS{ std::vector<double>() };

	// 2 D Vector of the 3 GPS
	std::vector<std::vector<Eigen::Vector3d>> GPSData{ num_GPSs,std::vector<Eigen::Vector3d>() };

	INS* inss[num_IMUs]{};
	ExtendedKalmanFilter ekf{};
	VIMUExtendedKalmanFilter vekf{};

}cm_INS;


Eigen::Vector3d transform_Gyro(const Eigen::Vector3d& gyroMeas, const Eigen::Quaternion<double>& orientation) {
	return orientation._transformVector(gyroMeas);
}

Eigen::Vector3d transform_Mag(const Eigen::Vector3d& magMeas, const Eigen::Quaternion<double>& orientation) {
	Eigen::Vector3d transformed_magMeas = magMeas;
	transformed_magMeas(1) *= -1;
	transformed_magMeas(2) *= -1;
	return orientation._transformVector(transformed_magMeas);
}

Eigen::Vector3d transform_Accel(int row, const INS& ins, const Eigen::Quaternion<double>& orientation) {
	Eigen::Vector3d psi_dot_dot;
	psi_dot_dot.setZero();

	int row_before = row - 1;

	if (1 <= row) {
		psi_dot_dot = (ins.imuData[row].accelMeas - ins.imuData[row_before].accelMeas) / (ins.timeDataIMU[row] - ins.timeDataIMU[row_before]);
	}

	return orientation._transformVector(ins.imuData[row].accelMeas);
	// Equation (2) of Data Fusion Algorithms for Multiple Inertial Measurement Units
	// return orientation._transformVector(imuData.accelMeas) - orientation._transformVector(psi_dot_dot.cross(ins.ekf.getCoords())) - orientation._transformVector(imuData.gyroMeas.cross(imuData.gyroMeas.cross(ins.ekf.getCoords())));
}


/// <summary>
/// Transforms the measurements of the INS IMU Information at position distinct from Center Of Mass to VIMU Values in Center of Mass
/// </summary>
/// <param name="ins"></param>
IMU_Data imu_2_vimu(INS& ins, int row) {
	Eigen::Quaternion<double> orientation = ins.ekf.getOrientation();

	IMU_Data transformed_Data{ IMU_Data{} };

	transformed_Data.gyroMeas = transform_Gyro(ins.imuData[row].gyroMeas, orientation);
	transformed_Data.magMeas = transform_Mag(ins.imuData[row].magMeas, orientation);
	transformed_Data.accelMeas = transform_Accel(row, ins, orientation);

	return transformed_Data;
}

Eigen::Vector4d computeDerivationQuaternion(const Eigen::Vector4d& q, const Eigen::Vector3d gyroMeas) {
	double q_0_dot = 1.0 / 2 * (-gyroMeas.x() * q(1) - gyroMeas.y() * q(2) - gyroMeas.z() * q(3));
	double q_1_dot = 1.0 / 2 * (gyroMeas.x() * q(0) + gyroMeas.z() * q(2) - gyroMeas.y() * q(3));
	double q_2_dot = 1.0 / 2 * (gyroMeas.y() * q(0) - gyroMeas.z() * q(1) + gyroMeas.x() * q(3));
	double q_3_dot = 1.0 / 2 * (gyroMeas.z() * q(0) + gyroMeas.y() * q(1) - gyroMeas.x() * q(2));

	return Eigen::Vector4d{ q_0_dot, q_1_dot, q_2_dot, q_3_dot };
}

Eigen::Quaternion<double> rungeKutta4thOrder(const Eigen::Quaternion<double>& q_prev, const Eigen::Vector3d& gyroMeas, double dt) {
	const double c1 = 0, c2 = 1.0 / 2, c3 = 1.0 / 2, c4 = 1;
	const double a21 = 1.0 / 2, a31 = 0, a32 = 1.0 / 2, a41 = 0, a42 = 0, a43 = 0;

	// Transform q_k to Vector4d for operations
	Eigen::Vector4d q_k{ q_prev.w(), q_prev.x(), q_prev.y(), q_prev.z() };

	// k1
	Eigen::Vector4d k1 = computeDerivationQuaternion(q_k, gyroMeas);

	// k2
	Eigen::Vector4d q_2 = q_k + dt * a21 * k1;
	Eigen::Vector4d k2 = computeDerivationQuaternion(q_2, gyroMeas);

	// k3
	Eigen::Vector4d q_3 = q_k + dt * (a31 * k1 + a32 * k2);
	Eigen::Vector4d k3 = computeDerivationQuaternion(q_3, gyroMeas);

	// k4
	Eigen::Vector4d q_4 = q_k + dt * (a41 * k1 + a42 * k2 + a43 * k3);
	Eigen::Vector4d k4 = computeDerivationQuaternion(q_4, gyroMeas);

	q_k = q_k + dt * 1.0 / 6 * (k1 + 2 * k2 + 2 * k3 + k4);

	return Eigen::Quaternion<double> {q_k(0), q_k(1), q_k(2), q_k(3)};
}

INS_state naiveFusion(const std::vector<INS_state>& localTracks);

void validateIMUData(CM_INS& centralized_ins, FusionInit fusion_init, std::string configName) {

	std::ofstream state_writer{};
	Eigen::VectorXd state;
	std::string path_time, path_position, path_orientation;

	double dt{};

	centralized_ins.timeDataIMU = centralized_ins.inss[0]->timeDataIMU;
	centralized_ins.vekf = VIMUExtendedKalmanFilter();

	for (INS* ins : centralized_ins.inss) {
		ins->ekf.setFusionInit(fusion_init);
	}

	int row_imu = 1, row_before = 0;

	for (INS* ins : centralized_ins.inss) {
		int row_imu = 1, row_before = 0;
		bool init = false;

		std::string file_appendix = std::to_string(ins->imu_port);


		//path_time = "C:/dev/Thesis/Kalman_Filter_Thesis/Kalman_Filter/Datalogs/time_integration_IMU_" + file_appendix + ".txt";
		//path_position = "C:/dev/Thesis/Kalman_Filter_Thesis/Kalman_Filter/Datalogs/position_integration_IMU_" + file_appendix + ".txt";
		//path_orientation = "C:/dev/Thesis/Kalman_Filter_Thesis/Kalman_Filter/Datalogs/orientation_integration_IMU_" + file_appendix + ".txt";

		//path_time = "C:/dev/Thesis/Kalman_Filter_Thesis/Kalman_Filter/Datalogs/Viereck Wohnzimmer/time_integration_IMU_" + file_appendix + ".txt";
		//path_position = "C:/dev/Thesis/Kalman_Filter_Thesis/Kalman_Filter/Datalogs/Viereck Wohnzimmer/position_integration_IMU_" + file_appendix + ".txt";
		//path_orientation = "C:/dev/Thesis/Kalman_Filter_Thesis/Kalman_Filter/Datalogs/Viereck Wohnzimmer/orientation_integration_IMU_" + file_appendix + ".txt";


		path_time = "C:/dev/Thesis/Kalman_Filter_Thesis/Kalman_Filter/Datalogs/Rotationen/time_integration_IMU_" + file_appendix + "_" + configName + ".txt";
		path_position = "C:/dev/Thesis/Kalman_Filter_Thesis/Kalman_Filter/Datalogs/Rotationen/position_integration_IMU_" + file_appendix + "_" + configName + ".txt";
		path_orientation = "C:/dev/Thesis/Kalman_Filter_Thesis/Kalman_Filter/Datalogs/Rotationen/orientation_integration_IMU_" + file_appendix + "_" + configName + ".txt";



		std::remove(path_time.c_str());
		std::remove(path_position.c_str());
		std::remove(path_orientation.c_str());


		while (row_imu < IMU_DATA_ROWS) {

			dt = centralized_ins.timeDataIMU[row_imu] - centralized_ins.timeDataIMU[row_imu - (int)1];
			if (!ins->ekf.isInitialised()) {
				ins->ekf.updateAcc(ins->imuData[row_imu].accelMeas, dt);
				ins->ekf.updateMag(ins->imuData[row_imu].magMeas, dt);
			}
			else {
				Eigen::VectorXd state = ins->ekf.getState();

				double x = state(0);
				double y = state(1);
				double z = state(2);
				double q_0 = state(9);
				double q_1 = state(10);
				double q_2 = state(11);
				double q_3 = state(12);

				// Orientation to the local tangent plane NED
				Eigen::Quaternion<double> orientation_vehicle = Eigen::Quaternion<double>{ q_0,q_1,q_2,q_3 };

				// transformation in origin of vehicle frame
				Eigen::Vector3d accMeas = ins->ekf.transform_Acc(ins->imuData[row_imu].accelMeas);
				if (ins->ekf.isValidMeasurement(accelerometer, accMeas)) {
					accMeas = orientation_vehicle._transformVector(accMeas);

					x = x + accMeas(0) * dt + 1.0 / 2.0 * accMeas(0) * dt * dt;
					y = y + accMeas(1) * dt + 1.0 / 2.0 * accMeas(1) * dt * dt;
					accMeas(2) += g;
					z = z + accMeas(2) * dt + 1.0 / 2.0 * accMeas(2) * dt * dt;

					//// Orientation bestimmen
					//double pitch = std::asin(accMeas(0) / g);

					//// compute roll -> turn/ bank angle, has also to be nearly zero
					//double param = accMeas(1) / accMeas(2);

					//double roll = std::atan(param);//m_y,m_z -> atan does the job [-90,90] is adequat

					//Eigen::Vector3d magMeas = ins->ekf.transform_Mag(ins->imuData[row_imu].magMeas);
					//if (ins->ekf.isValidMeasurement(magnetometer, magMeas)) {
					//	magMeas = orientation_vehicle._transformVector(magMeas);
					//	double yaw = std::atan2(magMeas(1), magMeas(0));
					//	
					//	orientation_vehicle = ins->ekf.computeOrientation(yaw, pitch, roll);
					//	orientation_vehicle.normalize();

					//}
				}

				// update with the gyro data

				Eigen::Vector3d gyroMeas = ins->ekf.transform_Gyro(ins->imuData[row_imu].gyroMeas);

				if (ins->ekf.isValidMeasurement(gyroscope, gyroMeas)) {
					gyroMeas = orientation_vehicle._transformVector(gyroMeas);

					double psi_x_dot = gyroMeas(0);
					double psi_y_dot = gyroMeas(1);
					double psi_z_dot = gyroMeas(2);

					/*Eigen::Quaternion<double> gyro_Quat = Eigen::Quaternion<double>{ 0, psi_x_dot, psi_y_dot, psi_z_dot };
					Eigen::Quaternion<double> placeholder_Quat = Eigen::Quaternion<double>{ 0,0,0,0 };

					placeholder_Quat = gyro_Quat * orientation_vehicle;
					placeholder_Quat = Eigen::Quaternion<double>{dt* 1.0 / 2 * placeholder_Quat.w(),dt* 1.0 / 2 * placeholder_Quat.x(),dt*1.0 / 2 * placeholder_Quat.y(),dt* 1.0 / 2 * placeholder_Quat.z() };

					orientation_vehicle = Eigen::Quaternion<double>{ orientation_vehicle.w() + placeholder_Quat.w(),orientation_vehicle.x() + placeholder_Quat.x(),orientation_vehicle.y() + placeholder_Quat.y(),orientation_vehicle.z() + placeholder_Quat.z() };*/

					///
					/// Runge Kutta Integration Algorithm -> Rk4n (Robust Methode IMU Calibration Page ??)
					///

					orientation_vehicle = rungeKutta4thOrder(orientation_vehicle, gyroMeas, dt);

					/// 3.: Mit Euler Winkeln ZYX-Sequenz
				//	Eigen::Vector3d euler = orientation_vehicle.toRotationMatrix().eulerAngles(2,1,0);


					orientation_vehicle.normalize();


				}


				// update with mag for orientation spaeter

				// Save output to file
				state_writer.open(path_time, std::ios_base::app);
				state_writer << centralized_ins.timeDataIMU[row_imu] << "\n";
				state_writer.close();

				state_writer.open(path_position, std::ios_base::app);
				state_writer << x << ", " << y << ", " << z << "\n";
				state_writer.close();

				state_writer.open(path_orientation, std::ios_base::app);
				state_writer << orientation_vehicle.w() << "," << orientation_vehicle.x() << "," << orientation_vehicle.y() << "," << orientation_vehicle.z() << "\n";
				state_writer.close();


				// save the double integrated data in the state
				state(0) = x;
				state(1) = y;
				state(2) = z;
				state(9) = orientation_vehicle.w();
				state(10) = orientation_vehicle.x();
				state(11) = orientation_vehicle.y();
				state(12) = orientation_vehicle.z();
				ins->ekf.setState(state);
			}

			row_imu++;
		}
	}
}



/*
* Raw Fusion
*/

void multiple_imu_fusion_raw(CM_INS& centralized_ins, std::string file_appendix) {

	int row_imu = 1; 
	std::ofstream state_writer{};
	std::string path_time, path_position, path_orientation;
	Eigen::VectorXd state;

	std::vector<Eigen::Quaternion<double>> VIMU_Orientations = std::vector<Eigen::Quaternion<double>>();
	std::vector<Eigen::Vector3d> VIMU_Coords = std::vector<Eigen::Vector3d>();

	centralized_ins.timeDataIMU = centralized_ins.inss[0]->timeDataIMU;

	for (INS* ins : centralized_ins.inss) {
		VIMU_Orientations.push_back(ins->ekf.getOrientation());
		VIMU_Coords.push_back(ins->ekf.getCoords());
	}

	centralized_ins.vekf = VIMUExtendedKalmanFilter(FusionInit::MAG, num_IMUs, VIMU_Orientations, VIMU_Coords);

	double dt_imu{};
	std::vector<Eigen::Vector3d> meas = std::vector<Eigen::Vector3d>();
	while (row_imu < IMU_DATA_ROWS) {

		// Prediction Step
		dt_imu = centralized_ins.timeDataIMU[row_imu] - centralized_ins.timeDataIMU[row_imu - (int)1];
		// Append gyro meas of each imu to vector
		for (INS* ins : centralized_ins.inss) {
			meas.push_back(ins->imuData[row_imu].gyroMeas);
		}

		centralized_ins.vekf.predictionStep(meas, dt_imu);

		// Update Step
			// Append acc meas of each imu to vector
		meas.clear();
		for (INS* ins : centralized_ins.inss) {
			meas.push_back(ins->imuData[row_imu].accelMeas);
		}

		centralized_ins.vekf.updateAcc(meas, dt_imu);

		// MAG only used for initialisation process
		if (!centralized_ins.vekf.isInitialised()) {

			// Append mag meas of each imu to vector
		meas.clear();
		for (INS* ins : centralized_ins.inss) {
			meas.push_back(ins->imuData[row_imu].magMeas);
		}

		centralized_ins.vekf.updateMag(meas, dt_imu);
			}

		if (centralized_ins.vekf.isInitialised()) {



			path_time = "C:/dev/Thesis/Kalman_Filter_Thesis/Kalman_Filter/Datalogs/time_vimu_" + file_appendix + ".txt";

			state_writer.open(path_time, std::ios_base::app);
			state_writer << centralized_ins.timeDataIMU[row_imu] << "\n";
			state_writer.close();

			state = centralized_ins.vekf.getState();
			//	std::cout << state << std::endl;

			path_position = "C:/dev/Thesis/Kalman_Filter_Thesis/Kalman_Filter/Datalogs/position_xyz_vimu_" + file_appendix + ".txt";
			state_writer.open(path_position, std::ios_base::app);
			//Eigen::Vector3d ecefPoint = centralized_ins.vekf.computeNED2ECEFwithRef(Eigen::Vector3d{ state(0), state(1), state(2) });
			//state_writer << ecefPoint(0) << ", " << ecefPoint(1) << ", " << ecefPoint(2) << "\n";
			state_writer << state(0) << "," << state(1) << "," << state(2) << "\n";
			state_writer.close();

			path_orientation = "C:/dev/Thesis/Kalman_Filter_Thesis/Kalman_Filter/Datalogs/orientation_vimu_" + file_appendix + ".txt";
			state_writer.open(path_orientation, std::ios_base::app);

			state_writer << state(9) << ", " << state(10) << ", " << state(11) << "," << state(12) << "\n";
			state_writer.close();




			//if (state(0) > 1000 || state(1) > 1000||state(2)>1000) {

			//	if(row_imu<530){
			//			std::cout << "Zeit: " << centralized_ins.timeDataIMU[row_imu] << std::endl;
			//	std::cout << "State:\n " << centralized_ins.vekf.getState() << std::endl;// << ", Covariance: " << cm_INS.ekf.getCovariance() << std::endl;
			//}
		}
		meas.clear();
		row_imu++;
	}


};

void multiple_imu_fusion_raw_gps(CM_INS& centralized_ins, std::string file_appendix) {

	int row_imu = 1; // wie kann man das noch eleganter loesen, auch wegen unten den Zeiten etc.
	int row_gps = 0;

	std::ofstream state_writer{};
	std::string path_time, path_position;

	Eigen::VectorXd state;

	std::vector<Eigen::Quaternion<double>> VIMU_Orientations = std::vector<Eigen::Quaternion<double>>();
	std::vector<Eigen::Vector3d> VIMU_Coords = std::vector<Eigen::Vector3d>();

	centralized_ins.timeDataIMU = centralized_ins.inss[0]->timeDataIMU;

	for (INS* ins : centralized_ins.inss) {
		VIMU_Orientations.push_back(ins->ekf.getOrientation());
		VIMU_Coords.push_back(ins->ekf.getCoords());
	}

	centralized_ins.vekf = VIMUExtendedKalmanFilter(FusionInit::MAGGPS, num_IMUs, VIMU_Orientations, VIMU_Coords);

	double dt_imu{}, dt_gps{};
	std::vector<Eigen::Vector3d> meas = std::vector<Eigen::Vector3d>();
	while (row_gps < GPS_DATA_ROWS && row_imu < IMU_DATA_ROWS) {

		// Prediction Step
		dt_imu = centralized_ins.timeDataIMU[row_imu] - centralized_ins.timeDataIMU[row_imu - (int)1];
		// Append gyro meas of each imu to vector
		for (INS* ins : centralized_ins.inss) {
			meas.push_back(ins->imuData[row_imu].gyroMeas);
		}

		centralized_ins.vekf.predictionStep(meas, dt_imu);

		// Update Step
		// Append acc meas of each imu to vector
		meas.clear();
		for (INS* ins : centralized_ins.inss) {
			meas.push_back(ins->imuData[row_imu].accelMeas);
		}

		centralized_ins.vekf.updateAcc(meas, dt_imu);

		// MAG only used for initialisation process
		if (!centralized_ins.vekf.isInitialised()) {

			// Append mag meas of each imu to vector
			meas.clear();
			for (INS* ins : centralized_ins.inss) {
				meas.push_back(ins->imuData[row_imu].magMeas);
			}

			centralized_ins.vekf.updateMag(meas, dt_imu);
		}

		int row_before = row_imu - 1;
		if (centralized_ins.timeDataIMU[row_before] < centralized_ins.timeDataGPS[row_gps] && centralized_ins.timeDataGPS[row_gps] < centralized_ins.timeDataIMU[row_imu]) {// was ist wenn es mehrere GPS Messungen in der Zeit der Aquirierung gibt??

			// Important for velocioty (if added)
			if (row_gps > 0) {
				dt_gps = centralized_ins.timeDataGPS[row_gps] - centralized_ins.timeDataGPS[row_gps - (int)1];
			}
			else { dt_gps = 0; }

			// Append gps meas to vector
			meas.clear();
			for (std::vector<Eigen::Vector3d> gpsMeas : centralized_ins.GPSData) {
				meas.push_back(gpsMeas[row_gps]);
			}

			centralized_ins.vekf.updateGPS(meas, dt_gps);
			row_gps++;

		}

		if (centralized_ins.vekf.isInitialised()) {
			path_time = "C:/dev/Thesis/Kalman_Filter_Thesis/Kalman_Filter/Datalogs/time_vimu_" + file_appendix + ".txt";

			state_writer.open(path_time, std::ios_base::app);
			state_writer << centralized_ins.timeDataIMU[row_imu] << "\n";
			state_writer.close();

			state = centralized_ins.vekf.getState();

			path_position = "C:/dev/Thesis/Kalman_Filter_Thesis/Kalman_Filter/Datalogs/position_xyz_vimu_" + file_appendix + ".txt";
			state_writer.open(path_position, std::ios_base::app);
			state_writer << state(0) << "," << state(1) << "," << state(2) << "\n";
			state_writer.close();

			if (state(0) > 1000 || state(1) > 1000) {

				std::cout << "Zeit: " << centralized_ins.timeDataIMU[row_imu] << std::endl;
				std::cout << "State:\n " << centralized_ins.vekf.getState() << std::endl;// << ", Covariance: " << cm_INS.ekf.getCovariance() << std::endl;
			}
		}
		meas.clear();
		row_imu++;
		row_before++;
	}

}


// Naive Fusion like http://fusion.isif.org/proceedings/fusion01CD/fusion/searchengine/pdf/WeA22.pdf
INS_state naiveFusion(const std::vector<INS_state>& localTracks) {

	// init covariance of system 
	Eigen::MatrixXd inverseTotalCovariance = Eigen::MatrixXd::Zero(localTracks[0].covariance.rows(), localTracks[0].covariance.cols());

	// compute the covariance
	for (const INS_state& localTrack : localTracks) {
		inverseTotalCovariance += localTrack.covariance.inverse();
	}

	// compute inverse covariance 
	Eigen::MatrixXd totalCovariance = inverseTotalCovariance.inverse();

	// compute the combined state 
	Eigen::VectorXd totalState = Eigen::VectorXd::Zero(localTracks[0].state.rows(), localTracks[0].state.cols());

	for (const INS_state& localTrack : localTracks) {
		totalState += localTrack.covariance.inverse() * localTrack.state;
	}

	totalState = totalCovariance * totalState;

	INS_state fusion_center;
	fusion_center.covariance = totalCovariance;
	fusion_center.state = totalState;

	return fusion_center;
}


/*
* Estimation Fusion
*/

void multiple_imu_fusion_estimation(CM_INS& centralized_ins, std::string file_appendix, FusionInit fusion_init) {

	double sampleTime = 1000000; // 1 s
	std::ofstream state_writer{};
	Eigen::VectorXd state;
	std::string path_time, path_position;

	double dt_imu{};

	centralized_ins.timeDataIMU = centralized_ins.inss[0]->timeDataIMU;
	centralized_ins.vekf = VIMUExtendedKalmanFilter();

	// With what sensors does our fusion start??
	for (INS* ins : centralized_ins.inss) {
		ins->ekf.setFusionInit(fusion_init);
	}

	// jeder einzelne Sensor schaetzt seinen Zustand in sich selbst (nach transformieren in den Rahmen der VIMU)
	// dann beim eintreffen der GPS wird das Mittel der Zustaende genommen
	// zustand aller imus wird zu dem der VIMU\

	int row_imu = 1, row_before = 0;
	bool init = false;

	std::vector<Eigen::Vector3d> meas = std::vector<Eigen::Vector3d>();

	while (row_imu < IMU_DATA_ROWS) {

		dt_imu = centralized_ins.timeDataIMU[row_imu] - centralized_ins.timeDataIMU[row_imu - (int)1];

		// Update the IMU extended Kalman filter with the values rotated in the virtual IMU frame in the Center of Mass
		for (INS* ins : centralized_ins.inss) {

			//dt_imu = ins->timeDataIMU[row_imu] - ins->timeDataIMU[row_imu - 1];
			ins->ekf.predictionStep(ins->imuData[row_imu].gyroMeas, dt_imu);

			ins->ekf.updateAcc(ins->imuData[row_imu].accelMeas, dt_imu);

			// MAG only used for initialisation process
			if (!ins->ekf.isInitialised()) {
				ins->ekf.updateMag(ins->imuData[row_imu].magMeas, dt_imu);
			}
		}

		row_before = row_imu - 1;

		// durchlaufen der init abfrage der ins
		if (!init) {
			init = true;
			for (INS* ins : centralized_ins.inss) {
				init = init && ins->ekf.isInitialised();
			}
		}


		if (init) {

			std::vector<INS_state> ins_states = std::vector<INS_state>();
			INS_state ins_state;
			for (INS* ins : centralized_ins.inss) {
				ins_state.state = ins->ekf.getState();
				ins_state.covariance = ins->ekf.getCovariance();
				ins_states.push_back(ins_state);
			}
			fusion_center = naiveFusion(ins_states);

			/*	state_writer.open("C:/dev/Thesis/Kalman_Filter_Thesis/Kalman_Filter/time_federated.txt", std::ios_base::app);
				state_writer << centralized_ins.timeDataIMU[row_imu] << "\n";
				state_writer.close();
				*/


			path_time = "C:/dev/Thesis/Kalman_Filter_Thesis/Kalman_Filter/Datalogs/time_estimation_fusion_" + file_appendix + ".txt";

			state_writer.open(path_time, std::ios_base::app);
			state_writer << centralized_ins.timeDataIMU[row_imu] << "\n";
			state_writer.close();

			state = fusion_center.state;

			path_position = "C:/dev/Thesis/Kalman_Filter_Thesis/Kalman_Filter/Datalogs/position_xyz_estimation_fusion_" + file_appendix + ".txt";
			state_writer.open(path_position, std::ios_base::app);
			state_writer << state(0) << "," << state(1) << "," << state(2) << "\n";
			state_writer.close();
			//if (state(0) > 1000 || state(1) > 1000) {

			//	std::cout << "Zeit: " << centralized_ins.timeDataIMU[row_imu] << std::endl;
			//	std::cout << "State:\n " << centralized_ins.vekf.getState() << std::endl;// << ", Covariance: " << cm_INS.ekf.getCovariance() << std::endl;
			//}

		}
		row_imu++;
	}

};

void multiple_imu_fusion_estimation_gps(CM_INS& centralized_ins, std::string file_appendix, FusionInit fusion_init) {
	double sampleTime = 1000000; // 1 ms
	std::ofstream state_writer{};
	Eigen::VectorXd state;
	std::string path_time, path_position;



	double dt_imu{}, dt_gps{};

	centralized_ins.timeDataIMU = centralized_ins.inss[0]->timeDataIMU;
	centralized_ins.vekf = VIMUExtendedKalmanFilter();

	// With what sensors does our fusion start??
	for (INS* ins : centralized_ins.inss) {
		ins->ekf.setFusionInit(fusion_init);
	}

	// jeder einzelne Sensor schaetzt seinen Zustand in sich selbst (nach transformieren in den Rahmen der VIMU)
	// dann beim eintreffen der GPS wird das Mittel der Zustaende genommen
	// zustand aller imus wird zu dem der VIMU\

	int row_imu = 1, row_gps = 0, row_before = 0;
	bool init = false;

	std::vector<Eigen::Vector3d> meas = std::vector<Eigen::Vector3d>();

	while (row_gps < GPS_DATA_ROWS && row_imu < IMU_DATA_ROWS) {

		dt_imu = centralized_ins.timeDataIMU[row_imu] - centralized_ins.timeDataIMU[row_imu - (int)1];

		// Update the IMU extended Kalman filter with the values rotated in the virtual IMU frame in the Center of Mass
		for (INS* ins : centralized_ins.inss) {

			//dt_imu = ins->timeDataIMU[row_imu] - ins->timeDataIMU[row_imu - 1];
			ins->ekf.predictionStep(ins->imuData[row_imu].gyroMeas, dt_imu);

			ins->ekf.updateAcc(ins->imuData[row_imu].accelMeas, dt_imu);

			// MAG only used for initialisation process
			if (!ins->ekf.isInitialised()) {
				ins->ekf.updateMag(ins->imuData[row_imu].magMeas, dt_imu);
			}
		}

		row_before = row_imu - 1;

		// durchlaufen der init abfrage der ins
		if (!init) {
			init = true;
			for (INS* ins : centralized_ins.inss) {
				init = init && ins->ekf.isInitialised();
			}
		}


		if (centralized_ins.timeDataIMU[row_before] < centralized_ins.timeDataGPS[row_gps] && centralized_ins.timeDataGPS[row_gps] < centralized_ins.timeDataIMU[row_imu]) {

			if (row_gps > 0) {
				dt_gps = centralized_ins.timeDataGPS[row_gps] - centralized_ins.timeDataGPS[row_gps - (int)1];
			}
			else { dt_gps = 0; }

			// Append gps meas to vector
			meas.clear();
			for (std::vector<Eigen::Vector3d> gpsMeas : centralized_ins.GPSData) {
				meas.push_back(gpsMeas[row_gps]);
			}



			/*
			* Second try with the "naive" Fusion in a Fusion Center, that is not a EKF
			* Just a Struct that Holds a covariance and an imu_state, not dependend on the state prior
			* https://www.intechopen.com/chapters/57985
			*/
			Eigen::Vector3d gpsMeasAvg = Eigen::Vector3d::Zero();
			int num_IMUs_avg = 0;
			// 1. transform
			for (int i = 0; i < meas.size(); i++) {
				if (centralized_ins.vekf.isValidMeasurement(gps, meas[i])) {
					gpsMeasAvg += meas[i];
					num_IMUs_avg++;
				}
			}
			gpsMeasAvg /= num_IMUs_avg;

			for (INS* ins : centralized_ins.inss)
			{
				ins->ekf.updateGPS(gpsMeasAvg, dt_gps);

			}



			row_gps++;
		}


		if (init) {

			std::vector<INS_state> ins_states = std::vector<INS_state>();
			INS_state ins_state;
			for (INS* ins : centralized_ins.inss) {
				ins_state.state = ins->ekf.getState();
				ins_state.covariance = ins->ekf.getCovariance();
				ins_states.push_back(ins_state);
			}
			fusion_center = naiveFusion(ins_states);

			/*	state_writer.open("C:/dev/Thesis/Kalman_Filter_Thesis/Kalman_Filter/time_federated.txt", std::ios_base::app);
				state_writer << centralized_ins.timeDataIMU[row_imu] << "\n";
				state_writer.close();*/

				//	state = centralized_ins.vekf.getState();


			path_time = "C:/dev/Thesis/Kalman_Filter_Thesis/Kalman_Filter/Datalogs/time_estimation_fusion_" + file_appendix + ".txt";

			state_writer.open(path_time, std::ios_base::app);
			state_writer << centralized_ins.timeDataIMU[row_imu] << "\n";
			state_writer.close();

			state = fusion_center.state;

			path_position = "C:/dev/Thesis/Kalman_Filter_Thesis/Kalman_Filter/Datalogs/position_xyz_estimation_fusion_" + file_appendix + ".txt";
			state_writer.open(path_position, std::ios_base::app);
			state_writer << state(0) << "," << state(1) << "," << state(2) << "\n";
			state_writer.close();
			if (state(0) > 1000 || state(1) > 1000) {

				std::cout << "Zeit: " << centralized_ins.timeDataIMU[row_imu] << std::endl;
				std::cout << "State:\n " << centralized_ins.vekf.getState() << std::endl;// << ", Covariance: " << cm_INS.ekf.getCovariance() << std::endl;
			}

		}
		row_imu++;
	}
}



/// <summary>
/// Transforms the measurements of the INS delivered from the IMU given the sensor error Model of the IMU sensors
/// </summary>
/// <param name="ins"></param> ins carrying the sensor data which has to be calibrated
void get_calibrated_meas(INS& ins) {
	calibration::IMU_Calibration calibration_Params = ins.ekf.getCalibrationParams();
	calibration::Calibration_Params accel_Params = calibration_Params.accelCali;
	calibration::Calibration_Params gyro_params = calibration_Params.gyroCali;
	calibration::Calibration_Params mag_params = calibration_Params.magCali;

	for (IMU_Data& data : ins.imuData) {
		//data.accelMeas = (accel_Params.theta * data.accelMeas - accel_Params.bias) * 9.81 / 1e3; // m/s2
		//data.gyroMeas = (gyro_params.theta * data.gyroMeas - gyro_params.bias) * M_PI / 180.0; // rad/s
		data.magMeas = mag_params.theta * (data.magMeas - mag_params.bias);

		data.accelMeas = (data.accelMeas) * 9.81 / 1e3;
		data.gyroMeas = (gyro_params.theta * data.gyroMeas - gyro_params.bias) * M_PI / 180.0;
	}
}

void init_local_INS(Orientation& orientation) {
	// Define INS1
	ins_1.imu_port = 0;
	ins_1.imuData = std::vector<IMU_Data>();
	ins_1.timeDataIMU = std::vector<double>();
	ins_1.GPSData = std::vector<Eigen::Vector3d>();
	ins_1.timeDataGPS = std::vector<double>();
	ins_1.ekf.setCoords(Eigen::Vector3d{ 76.7441, -1.6672, -53.0125 });
	ins_1.ekf.setAccelBias(Eigen::Vector3d{ -4.31342161, -22.0591438, 29.0506018 });
	ins_1.ekf.setAccelTransformMatrix((Eigen::Matrix3d() << 0.997624911, 0.00501776681, 0.0211610225, -0.00811466326, 0.986648117, 0.136514105, -0.0214393877, -0.138505947, 0.985038735).finished());
	ins_1.ekf.setMagBias(Eigen::Vector3d
		{ -21.155646, 15.08731, 60.443188 });
	ins_1.ekf.setMagTransformMatrix((Eigen::Matrix3d() << 0.945738, -0.003007, -0.006546, -0.003007, 0.922916, 0.007009, -0.006546, 0.007009, 0.949019).finished());
	ins_1.ekf.setGyroBias(Eigen::Vector3d{ 0.01302428, 0.00608158, 0.00806136 });
	ins_1.ekf.setGyroTransformMatrix((Eigen::Matrix3d() << 0.92552028, -0.03281577, 0.33111104, 0.04988341, 0.88751086, 0.42732689, 0.13176533, -0.15967813, 0.80570043).finished());

	// Define INS2
	ins_2.imu_port = 1;
	ins_2.imuData = std::vector<IMU_Data>();
	ins_2.timeDataIMU = std::vector<double>();
	ins_2.GPSData = std::vector<Eigen::Vector3d>();
	ins_2.timeDataGPS = std::vector<double>();

	ins_2.ekf.setCoords(Eigen::Vector3d{ 1.6521 ,  75.3151 ,  54.5905 });
	ins_2.ekf.setAccelBias(Eigen::Vector3d{ -3.1927911, -19.8014002, 5.25052353 });
	ins_2.ekf.setAccelTransformMatrix((Eigen::Matrix3d() << 0.998175208, -0.0131904022, 0.00489315879, 0.0138542229, 0.998597102, 0.0123811444, -0.00724020321, -0.0177614771, 0.993377592).finished());
	ins_2.ekf.setMagBias(Eigen::Vector3d{ 6.847327, -22.258473, 54.19993 });
	ins_2.ekf.setMagTransformMatrix((Eigen::Matrix3d() << 1.264232, -0.044928, -0.002404, -0.044928, 1.225985, -0.006535, -0.002404, -0.006535, 1.244639).finished());
	ins_2.ekf.setGyroBias(Eigen::Vector3d{ 0.00381464, 0.01862675, 0.00931548 });
	ins_2.ekf.setGyroTransformMatrix((Eigen::Matrix3d() << 0.91231546, 0.1936244, -0.01587249, -0.04097905, 0.95112852, 0.0414368, -0.01217961, 0.2924943, 0.90286951).finished());


	// Define INS3
	ins_3.imu_port = 6;
	ins_3.imuData = std::vector<IMU_Data>();
	ins_3.timeDataIMU = std::vector<double>();
	ins_3.GPSData = std::vector<Eigen::Vector3d>();
	ins_3.timeDataGPS = std::vector<double>();

	ins_3.ekf.setCoords(Eigen::Vector3d{ -76.5838, -1.6672, -53.0125 });
	ins_3.ekf.setAccelBias(Eigen::Vector3d{ -3.00372421, -8.129569, 16.655453 });
	ins_3.ekf.setAccelTransformMatrix((Eigen::Matrix3d() << 0.997494151, -0.0224596029, 0.0299220177, 0.0216129065, 0.996283148, -0.0106842197, -0.032514178, 0.00961013661, 0.993912779).finished());
	ins_3.ekf.setMagBias(Eigen::Vector3d{ -13.293332, -8.99686, 0.35697 });
	ins_3.ekf.setMagTransformMatrix((Eigen::Matrix3d() << 1.225777, 0.032457, 0.029723, 0.032457, 1.188017, -0.026329, 0.029723, -0.026329, 1.139707).finished());
	ins_3.ekf.setGyroBias(Eigen::Vector3d{ -0.03347582, 0.05826998, 0.0033403 });
	ins_3.ekf.setGyroTransformMatrix((Eigen::Matrix3d() << 1.03047692, 0.11318723, 0.05332105, -0.04702399, 0.73067087, -0.11064558, 0.09328256, -0.4519223, 0.57949588).finished());


	// Define INS4
	ins_4.imu_port = 7;
	ins_4.imuData = std::vector<IMU_Data>();
	ins_4.timeDataIMU = std::vector<double>();
	ins_4.GPSData = std::vector<Eigen::Vector3d>();
	ins_4.timeDataGPS = std::vector<double>();

	ins_4.ekf.setCoords(Eigen::Vector3d{ 1.6521 , -75.3151, 54.5905 }); // given in mm, conversion to m in ekf
	ins_4.ekf.setAccelBias(Eigen::Vector3d{ 0.513195483, -6.38354307, 18.6818155 });
	ins_4.ekf.setAccelTransformMatrix((Eigen::Matrix3d() << 0.996643923, 0.00345847616, -0.0105455711, -0.0028747351, 0.997351685, 0.0171765314, 0.0104712047, -0.0120770725, 0.996643184).finished());
	ins_4.ekf.setMagBias(Eigen::Vector3d{ 6.961981, -44.798494, 36.018804 });
	ins_4.ekf.setMagTransformMatrix((Eigen::Matrix3d() << 0.919869, -0.063969, 0.036467, -0.063969, 1.003396, 0.002549, 0.036467, 0.002549, 1.048945).finished());
	ins_4.ekf.setGyroBias(Eigen::Vector3d{ -0.04707186, 0.00766569, 0.00236872 });
	ins_4.ekf.setGyroTransformMatrix((Eigen::Matrix3d() << 0.93908668, -0.0466796, -0.20330397, -0.13423415, 0.67620875, -0.33500868, 0.08244915, -0.28759273, 0.49361518).finished());

	// set Orientation
	orientation.o_imu_0.normalize();
	orientation.o_imu_1.normalize();
	orientation.o_imu_6.normalize();
	orientation.o_imu_7.normalize();

	ins_1.ekf.setOrientation(orientation.o_imu_0);
	ins_2.ekf.setOrientation(orientation.o_imu_1);
	ins_3.ekf.setOrientation(orientation.o_imu_6);
	ins_4.ekf.setOrientation(orientation.o_imu_7);
}

int read_imu_data(std::string data_IMU_path) {

	for (INS* ins : cm_INS.inss)
	{
		std::stringstream filepath{};
		filepath << data_IMU_path << "/IMU_" << ins->imu_port << "_data.txt";
		std::ifstream data_IMU{ filepath.str() };

		if (!data_IMU.is_open())
		{
			std::cerr << "Fehler beim Öffnen der Datei!" << std::endl;
			return 1;
		}

		double time{};
		IMU_Data imuData{};
		Eigen::Vector3d accelMeas{};
		Eigen::Vector3d magMeas{};
		Eigen::Vector3d gyroMeas{};
		IMU_DATA_ROWS = 0;

		for (std::string values; std::getline(data_IMU, values);)
		{

			std::istringstream iss(values);

			iss >> time;

			// Time in seconds
			time /= time_constant;
			ins->timeDataIMU.push_back(time);

			iss >> accelMeas[0] >> accelMeas[1] >> accelMeas[2];

			iss >> gyroMeas[0] >> gyroMeas[1] >> gyroMeas[2];

			iss >> magMeas[0] >> magMeas[1] >> magMeas[2];



			imuData.accelMeas = accelMeas;
			imuData.gyroMeas = gyroMeas;
			imuData.magMeas = magMeas;

			ins->imuData.push_back(imuData);

			IMU_DATA_ROWS++;

		}
	}
	return 0;

}

int read_gps_data(std::string data_GPS_path) {
	// Format: Time [us], Latitude [deg], Longitude Degree[deg]

	for (int gps_number = 0; gps_number < num_GPSs; gps_number++) {
		std::stringstream filepath{};
		filepath << data_GPS_path << "/GPS_" << gps_number + 1 << "_data.txt";
		std::ifstream data_GPS{ filepath.str() };

		if (!data_GPS.is_open())
		{
			std::cout << filepath.str() << std::endl;
			std::cerr << "Fehler beim Öffnen der Datei!" << std::endl;
			return 1;
		}

		double time;
		Eigen::Vector3d GPSData;
		GPS_DATA_ROWS = 0;

		for (std::string values; std::getline(data_GPS, values);)
		{
			std::istringstream iss(values);

			iss >> time;
			time /= time_constant;
			cm_INS.timeDataGPS.push_back(time);

			// Latitude, Longitude
			iss >> GPSData(0) >> GPSData(1);

			// Altitude
			GPSData(2) = 0;

			cm_INS.GPSData[gps_number].push_back(GPSData);
			GPS_DATA_ROWS++;
		}

	}
	return 1;
}


void fuse(Configuration& configuration) {
	FusionMethod fusion_method = configuration.getFusion_Method();

	init_local_INS(configuration.getOrientations());

	// Multi-IMU
//	cm_INS.inss[0] = &ins_1;
//	cm_INS.inss[1] = &ins_2;
	//cm_INS.inss[2] = &ins_3;
	//cm_INS.inss[3] = &ins_4;

	// Single IMU
	cm_INS.inss[0] = &ins_4;

	read_imu_data(configuration.getIMU_Data_Path());

	if (fusion_method == FusionMethod::Raw_GPS || fusion_method == FusionMethod::Federated_GPS) {
		read_gps_data(configuration.getGPS_Data_Path());
	}

	for (INS* ins : cm_INS.inss) {
		get_calibrated_meas(*ins);
	}

	//validateIMUData(cm_INS, configuration.getFusion_Init(), configuration.getName());
	//return;

	switch (fusion_method) {
	case FusionMethod::Raw:
		multiple_imu_fusion_raw(cm_INS, configuration.getName());
		break;
	case FusionMethod::Raw_GPS:
		multiple_imu_fusion_raw_gps(cm_INS, configuration.getName());
		break;
	case FusionMethod::Federated:
		multiple_imu_fusion_estimation(cm_INS, configuration.getName(), configuration.getFusion_Init());
		break;
	case FusionMethod::Federated_GPS:
		multiple_imu_fusion_estimation_gps(cm_INS, configuration.getName(), configuration.getFusion_Init());
		break;
	}

}


int main()
{

	// Orientation of the IMUs in the previous thesis
	Orientation orientation_bochum{};
	orientation_bochum.o_imu_0 = Eigen::Quaternion<double>{ -0.23 , -0.769, -0.444, -0.398 };
	orientation_bochum.o_imu_1 = Eigen::Quaternion<double>{ -0.23, 0.119, 0.444, 0.858 };
	orientation_bochum.o_imu_6 = Eigen::Quaternion<double>{ -0.39807, -0.44399, 0.76902, 0.22983 };
	orientation_bochum.o_imu_7 = Eigen::Quaternion<double>{ 0.85781 ,0.44404, -0.11898 , 0.22985 };

	// Orientation of the IMUs in Test Data Aquisition from 18.4.2024
	Orientation orientation_default{};
	orientation_default.o_imu_0 = Eigen::Quaternion<double>{ 0.23, 0.769, 0.444, 0.398 };
	orientation_default.o_imu_1 = Eigen::Quaternion<double>{ -0.858,0.444,-0.119,-0.23 };
	orientation_default.o_imu_6 = Eigen::Quaternion<double>{ -0.23, 0.769, 0.444, -0.398 };
	orientation_default.o_imu_7 = Eigen::Quaternion<double>{ -0.23, -0.119, -0.444, 0.858 };

	/*
	* Format of the Input Data
	* whitespace as delimiter (easier reading with inputstream)
	*
	* Time [ms],ACC_X [mg],ACC_Y [mg],ACC_Z [mg],GYRO_X [dps],GYRO_Y [dps],GYRO_Z [dps],MAG_X [uT],MAG_Y [uT],MAG_Z [uT]
	*/

	/*
	* Create Configs for the different test data
	*/
	// Testdata Bochum
	Configuration bochum_raw_gps = Configuration("bochum_raw_gps", orientation_bochum, "C:/Users/veigh/Desktop/Bachelor-Arbeit/Code/old_data_from_sd", FusionMethod::Raw_GPS, FusionInit::MAGGPS, "C:/Users/veigh/Desktop/Bachelor-Arbeit/Code/old_data_from_sd");

	Configuration bochum_federated_gps = Configuration("bochum_federated_gps", orientation_bochum, "C:/Users/veigh/Desktop/Bachelor-Arbeit/Code/old_data_from_sd", FusionMethod::Federated_GPS, FusionInit::MAGGPS, "C:/Users/veigh/Desktop/Bachelor-Arbeit/Code/old_data_from_sd");

	//	fuse(bochum_raw_gps);
		//fuse(bochum_federated_gps);


		/*
		* Create Configs for Froemern Drive 18.4
		*/
	Configuration froemern_raw = Configuration("froemern_raw", orientation_default, "C:/Users/veigh/Desktop/Bachelor-Arbeit/Code/Datalogs Veit/Froemern Fahrt/18_4_2024/Konstrukt Daten", FusionInit::MAG);

	Configuration froemern_federated = Configuration("fromern_federated", orientation_default, "C:/Users/veigh/Desktop/Bachelor-Arbeit/Code/Datalogs Veit/Froemern Fahrt/18_4_2024/Konstrukt Daten", FusionMethod::Federated, FusionInit::MAG);


	//	fuse(froemern_raw);
	//	fuse(froemern_federated);


			/*
		* Create Configs for Billmerich Rundgang 20.4	*/
	Configuration billmerich_raw = Configuration("billmerich_raw", orientation_default, "C:/Users/veigh/Desktop/Bachelor-Arbeit/Code/Datalogs Veit/Spaziergang Billmericher Dorfstrasse/Konstrukt Daten", FusionInit::MAG);

	Configuration billmerich_federated = Configuration("billmerich_federated", orientation_default, "C:/Users/veigh/Desktop/Bachelor-Arbeit/Code/Datalogs Veit/Spaziergang Billmericher Dorfstrasse/Konstrukt Daten", FusionMethod::Federated, FusionInit::MAG);

	//fuse(billmerich_raw);
	//	fuse(billmerich_federated);


	/*
	* Config Vierecktest
	*/

	Configuration wohnzimmer_raw = Configuration("wohnzimmer_raw", orientation_default, "C:/Users/veigh/Desktop/Bachelor-Arbeit/Code/Datalogs Veit/Viereck Wohnzimmer/Konstrukt Daten", FusionInit::MAG);

	//	fuse(wohnzimmer_raw);

		/*
		* Configs Rotation
		*/

	Configuration rotation_raw_x = Configuration("rotation_raw_x", orientation_default, "C:/Users/veigh/Desktop/Bachelor-Arbeit/Code/Datalogs Veit/Rotationen/Rotation um X/Konstrukt Daten", FusionInit::MAG);

	fuse(rotation_raw_x);


	Configuration rotation_raw_y = Configuration("rotation_raw_y", orientation_default, "C:/Users/veigh/Desktop/Bachelor-Arbeit/Code/Datalogs Veit/Rotationen/Rotation um Y/Konstrukt Daten", FusionInit::MAG);

	fuse(rotation_raw_y);


	Configuration rotation_raw_z = Configuration("rotation_raw_z", orientation_default, "C:/Users/veigh/Desktop/Bachelor-Arbeit/Code/Datalogs Veit/Rotationen/Rotation um Z/Konstrukt Daten", FusionInit::MAG);

	fuse(rotation_raw_z);




	return 0;



}

