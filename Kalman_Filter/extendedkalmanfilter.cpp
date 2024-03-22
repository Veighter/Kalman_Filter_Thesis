#include "extendedkalmanfilter.h"

// For computing the Fehlerfortpflanzungsmatrix
constexpr double ACCEL_STD = 1.0;
constexpr double GYRO_STD = 0.01 / 180.0 * M_PI;
constexpr double INIT_VEL_STD = 10.0;
constexpr double INIT_PSI_STD = 45.0 / 180.0 * M_PI;
constexpr double GPS_POS_STD = 3.0;


void ExtendedKalmanFilter::predictionStep(double dt) {
}

void ExtendedKalmanFilter::predictionStep(sensorMeas::AccelMeas accMeas, sensorMeas::GyroMeas gyroMeas, sensorMeas::MagMeas magMeas, double dt) {
}

void ExtendedKalmanFilter::updateStep(sensorMeas::GPSMeas gpsMeas, double dt) {
}