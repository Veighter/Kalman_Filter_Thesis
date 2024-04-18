#pragma once

#include <iostream>
#include <string>
#include <vector>

struct Orientation {
	Eigen::Quaternion<double> o_imu_0{};
	Eigen::Quaternion<double> o_imu_1{};
	Eigen::Quaternion<double> o_imu_6{};
	Eigen::Quaternion<double> o_imu_7{};
};


enum class FusionMethod {
	None, // no fusion, but just integration of the values of the IMUs -> to be implemented
	Raw,
	Raw_GPS,
	Federated, 
	Federated_GPS
};

class Configuration {
public:
	Configuration(std::string name, Orientation orientations_imu, std::string imu_data_path ) :c_name(name), c_orientations_imu(orientations_imu), c_imu_data_path(imu_data_path), c_gps_data_path(""), c_fusion_method(FusionMethod::Raw) {};
	Configuration(std::string name, Orientation orientations_imu, std::string imu_data_path, FusionMethod fusion_method) :c_name(name), c_orientations_imu(orientations_imu), c_imu_data_path(imu_data_path), c_gps_data_path(""), c_fusion_method(fusion_method) {};
	Configuration(std::string name, Orientation orientations_imu, std::string imu_data_path, std::string gps_data_path) :c_name(name), c_orientations_imu(orientations_imu), c_imu_data_path(imu_data_path), c_gps_data_path(gps_data_path), c_fusion_method(FusionMethod::Raw) {};
	Configuration(std::string name, Orientation orientations_imu, std::string imu_data_path, std::string gps_data_path, FusionMethod fusion_method) :c_name(name), c_orientations_imu(orientations_imu), c_imu_data_path(imu_data_path), c_gps_data_path(gps_data_path), c_fusion_method(fusion_method) {};


	std::string getName() {	return c_name;}
	std::string getIMU_Data_Path() { return c_imu_data_path; }
	std::string getGPS_Data_Path() { return c_gps_data_path; }
	FusionMethod getFusion_Method() { return c_fusion_method; }
	Orientation& getOrientations() { return c_orientations_imu; }

private:
	std::string c_name;
	std::string c_imu_data_path;
	std::string c_gps_data_path;

	FusionMethod c_fusion_method;
	Orientation c_orientations_imu;
};
