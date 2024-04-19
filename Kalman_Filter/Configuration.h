#pragma once

#include <iostream>
#include <string>
#include <vector>
#include "Fusion.h"

struct Orientation {
	Eigen::Quaternion<double> o_imu_0{};
	Eigen::Quaternion<double> o_imu_1{};
	Eigen::Quaternion<double> o_imu_6{};
	Eigen::Quaternion<double> o_imu_7{};
};


class Configuration {
public:
	Configuration(std::string name, Orientation orientations_imu, std::string imu_data_path, FusionInit fusion_init) :c_name(name), c_orientations_imu(orientations_imu), c_imu_data_path(imu_data_path), c_gps_data_path(""), c_fusion_method(FusionMethod::Raw), c_fusion_init(fusion_init) {};
	Configuration(std::string name, Orientation orientations_imu, std::string imu_data_path, FusionMethod fusion_method, FusionInit fusion_init) :c_name(name), c_orientations_imu(orientations_imu), c_imu_data_path(imu_data_path), c_gps_data_path(""), c_fusion_method(fusion_method), c_fusion_init(fusion_init) {};
	Configuration(std::string name, Orientation orientations_imu, std::string imu_data_path, FusionInit fusion_init, std::string gps_data_path) :c_name(name), c_orientations_imu(orientations_imu), c_imu_data_path(imu_data_path), c_gps_data_path(gps_data_path), c_fusion_method(FusionMethod::Raw), c_fusion_init(fusion_init) {};
	Configuration(std::string name, Orientation orientations_imu, std::string imu_data_path, FusionInit fusion_init, std::string gps_data_path, FusionMethod fusion_method) :c_name(name), c_orientations_imu(orientations_imu), c_imu_data_path(imu_data_path), c_gps_data_path(gps_data_path), c_fusion_method(fusion_method), c_fusion_init(fusion_init) {};


	std::string getName() {	return c_name;}
	std::string getIMU_Data_Path() { return c_imu_data_path; }
	std::string getGPS_Data_Path() { return c_gps_data_path; }
	FusionMethod getFusion_Method() { return c_fusion_method; }
	FusionInit getFusion_Init() { return c_fusion_init; }
	Orientation& getOrientations() { return c_orientations_imu; }
	

private:
	std::string c_name;
	std::string c_imu_data_path;
	std::string c_gps_data_path;

	FusionMethod c_fusion_method;
	FusionInit c_fusion_init;
	Orientation c_orientations_imu;
};
