#pragma once

enum class FusionMethod {
	None, // no fusion, but just integration of the values of the IMUs -> to be implemented
	Raw,
	Raw_GPS,
	Federated,
	Federated_GPS
};

/// <summary>
/// Fusion Configurations, priori knowledge of all available Sensors for the initialisation
/// </summary>
enum class FusionInit {
	MAG, // Magnetometer, Accelerometer, Gyroskop
	MAGGPS // Magnetometer, Accelerometer, Gyroskop, GPS
};
