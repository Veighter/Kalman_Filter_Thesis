#define _USE_MATH_DEFINES

#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <math.h>
#include "sensors.h"
#include "extendedkalmanfilter.h"
// 8851 Daten sind da drin vorhanden
// In welcher Einheit sind die Sekunden angegeben? Laut Abgabe in Microsekunden
// alle 2 sekunden kommt ein GPS rein
using namespace sensorMeas;
constexpr int IMU_DATA_ROWS = 8851;
constexpr int GPS_DATA_ROWS = 872;

// Was ist mit den GPS Daten, wie lege ich die an??
struct IMU_Data {
    Eigen::Vector3d accelMeas{};
    Eigen::Vector3d magMeas{};
    Eigen::Vector3d gyroMeas{};
};

/// <summary>
/// Calibration Params for one Sensor Model
/// </summary>
struct Calibration_Params {
    Eigen::Vector3d bias{};
    Eigen::Matrix3d theta{};
};

/// <summary>
/// The struct contains all calibrations parameters for the sensors of the used IMU according to their sensor error model
/// </summary>
struct IMU_Calibration {
    Calibration_Params accelCali{};
    Calibration_Params gyroCali{};
    Calibration_Params magCali{};
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
    int imu_port{}; // connect port on the multiplexer
    Eigen::Vector3d coords{}; // coordinats in mm
    Eigen::Quaternion<float> quat{};
    // IMU Data delivered
    std::vector<double> timeData{ std::vector<double>(IMU_DATA_ROWS) };
    std::vector<IMU_Data> imuData{ std::vector<IMU_Data>(IMU_DATA_ROWS) };
    ExtendedKalmanFilter ekf{};
    IMU_Calibration calib{};
    // quaternion or rotation matrix for transforming in the reference Frame
} ins_1, ins_2, ins_3, ins_4;


/// <summary>
/// The struct represents the Center Of Mass INS systems [VIMU] combining all 4 INS of the used IMUs in one INS
/// inss: INS-Array of the INS to the given IMUs in the 4 positions connected to the construct
/// ekf: Extended Kalman-Filter in the center of mass for interpreting the (calibrated) and averaged raw sensor data in the Center Of Mass
/// </summary>
struct CM_INS {
    std::vector<double> timeData{ std::vector<double>(IMU_DATA_ROWS) };
    INS inss[4]{};
    ExtendedKalmanFilter ekf{};
    
}global_INS;


Eigen::Vector3d transformToVector(double& data_1, double& data_2, double& data_3) {
      Eigen::Vector3d data = Eigen::Vector3d::Zero();
      data << data_1, data_2, data_3;
      return data;
}

/// <summary>
/// Berechnung der Multi IMU Fusion des Tetreader Konstrukts
/// </summary>
/// <param name="ins_1"></param> ins_1, carrying the information of the imu at port 0 and position 1, not the center of mass
/// <param name="ins_2"></param> ins_2, carrying the information of the imu at port 1 and position 2, not the center of mass
/// <param name="ins_3"></param> ins_3, carrying the information of the imu at port 6 and position 3, not the center of mass
/// <param name="ins_4"></param> ins_4, carrying the information of the imu at port 7 and position 4, not the center of mass
/// <param name="gpsmeas_1"></param> 
/// <param name="gpsmeas_2"></param>
/// <param name="gpsmeas_3"></param>
void multiple_imu_fusion(INS& ins_1, INS& ins_2, INS& ins_3, INS& ins_4, GPSMeas gpsmeas_1, GPSMeas gpsmeas_2, GPSMeas gpsmeas_3) {

}

/// <summary>
/// Transforms the measurements of the INS IMU Information at position distinct from Center Of Mass to VIMU Values in Center of Mass
/// </summary>
/// <param name="ins"></param>
void imu_2_vimu(INS& ins) {
    


   }


/// <summary>
/// Transforms the measurements of the INS delivered from the IMU given the sensor error Model of the IMU sensors
/// </summary>
/// <param name="ins"></param> ins carrying the sensor data which has to be calibrated
void get_calibrated_meas(INS& ins) {
    Calibration_Params accel_Params = ins.calib.accelCali;
    Calibration_Params gyro_params = ins.calib.gyroCali;
    Calibration_Params mag_params = ins.calib.magCali;

    for (IMU_Data& data : ins.imuData) {
        data.accelMeas = accel_Params.theta * data.accelMeas - accel_Params.bias;
        data.gyroMeas = gyro_params.theta * data.gyroMeas - gyro_params.bias;
        data.magMeas = mag_params.theta * (data.magMeas - mag_params.bias);
    }

}

int main()
{   
    // Flag fuer die Estimation Fusion
    constexpr bool estimation_fusion = false;

    // read in of datafiles with the coloum structure:
    // Time [us]	ACC_X [mg]	ACC_Y [mg]	ACC_Z [mg]	GYRO_X [dps]	GYRO_Y [dps]	GYRO_Z [dps]	MAG_X [uT]	MAG_Y [uT]	MAG_Z [uT]
    std::string data_IMU_path{ "C:/Users/veigh/Desktop/Bachelor-Arbeit/Code/old_data_from_sd/IMU_" };


    // Define INS1
    ins_1.imu_port = 0;
    ins_1.coords << 6.7441, 1.5375, 53.0125;
    ins_1.quat << -0.398 , 0.444 ,- 0.769 , 0.23;
    ins_1.calib.accelCali.bias << -4.31342161, -22.0591438, 29.0506018;
    ins_1.calib.accelCali.theta << 0.997624911, 0.00501776681, 0.0211610225, -0.00811466326, 0.986648117, 0.136514105, -0.0214393877, -0.138505947, 0.985038735;
    ins_1.calib.magCali.bias << -21.155646, 15.08731, 60.443188;
    ins_1.calib.magCali.theta << 0.945738, -0.003007, -0.006546, -0.003007, 0.922916, 0.007009, -0.006546, 0.007009, 0.949019;
    ins_1.calib.gyroCali.bias << 0.01302428, 0.00608158, 0.00806136;
    ins_1.calib.gyroCali.theta << 0.92552028, -0.03281577, 0.33111104, 0.04988341, 0.88751086, 0.42732689, 0.13176533, -0.15967813, 0.80570043;

    // Define INS2
    ins_2.imu_port = 1;
    ins_2.coords << 1.6521 ,- 75.4449 ,- 54.5905;
    ins_2.quat << -0.23,0.119,0.444, 0.858;
    ins_2.calib.accelCali.bias << -3.1927911, -19.8014002, 5.25052353;
    ins_2.calib.accelCali.theta << 0.998175208, -0.0131904022, 0.00489315879, 0.0138542229, 0.998597102, 0.0123811444, -0.00724020321, -0.0177614771, 0.993377592;
    ins_2.calib.magCali.bias << 6.847327, -22.258473, 54.19993;
    ins_2.calib.magCali.theta << 1.264232, -0.044928, -0.002404, -0.044928, 1.225985, -0.006535, -0.002404, -0.006535, 1.244639;
    ins_2.calib.gyroCali.bias << 0.00381464, 0.01862675, 0.00931548;
    ins_2.calib.gyroCali.theta << 0.91231546, 0.1936244, -0.01587249, -0.04097905, 0.95112852, 0.0414368, -0.01217961, 0.2924943, 0.90286951;


    // Define INS3
    ins_3.imu_port = 6;
    ins_3.coords << -76.5838, 1.5375, 53.0125;
    ins_3.quat << -0.23, 0.769, 0.444, - 0.398;
    ins_3.calib.accelCali.bias << -3.00372421, -8.129569, 16.655453;
    ins_3.calib.accelCali.theta << 0.997494151, -0.0224596029, 0.0299220177, 0.0216129065, 0.996283148, -0.0106842197, -0.032514178, 0.00961013661, 0.993912779;
    ins_3.calib.magCali.bias << -13.293332, -8.99686, 0.35697;
    ins_3.calib.magCali.theta << 1.225777, 0.032457, 0.029723, 0.032457, 1.188017, -0.026329, 0.029723, -0.026329, 1.139707;
    ins_3.calib.gyroCali.bias << -0.03347582, 0.05826998, 0.0033403;
    ins_3.calib.gyroCali.theta << 1.03047692, 0.11318723, 0.05332105, -0.04702399, 0.73067087, -0.11064558, 0.09328256, -0.4519223, 0.57949588;


    // Define INS4
    ins_4.imu_port = 7;
    ins_4.coords << 1.6521, 75.3151, -54.5905;
    ins_4.quat << -0.23, - 0.119, - 0.444, 0.858;
    ins_4.calib.accelCali.bias << 0.513195483, -6.38354307, 18.6818155;
    ins_4.calib.accelCali.theta << 0.996643923, 0.00345847616, -0.0105455711, -0.0028747351, 0.997351685, 0.0171765314, 0.0104712047, -0.0120770725, 0.996643184;
    ins_4.calib.magCali.bias << 6.961981, -44.798494, 36.018804;
    ins_4.calib.magCali.theta << 0.919869, -0.063969, 0.036467, -0.063969, 1.003396, 0.002549, 0.036467, 0.002549, 1.048945;
    ins_4.calib.gyroCali.bias << -0.04707186, 0.00766569, 0.00236872;
    ins_4.calib.gyroCali.theta << 0.93908668, -0.0466796, -0.20330397, -0.13423415, 0.67620875, -0.33500868, 0.08244915, -0.28759273, 0.49361518;


      
    INS* inss[4];
    inss[0] = &ins_1;
    inss[1] = &ins_2;
    inss[2] = &ins_3;
    inss[3] = &ins_4;

    int coloumn = 0; // Index for the array of values
    int row = 0;
    size_t pos = 0; // position of the delimiter in string
       

    // Einlesen der Messdaten in die IMUs 
    for (INS* ins : inss)
    {
        std::stringstream filepath{};
        filepath << data_IMU_path << ins->imu_port << "/IMU_"<< ins->imu_port<< "_data.txt";
        std::ifstream data_IMU{ filepath.str() };

        if (!data_IMU.is_open())
        {
            std::cerr << "Fehler beim Öffnen der Datei!" << std::endl;
            //return 1; // Rückgabe eines Fehlercodes
        }

        for (std::string values; std::getline(data_IMU, values) && row < IMU_DATA_ROWS; row++)
        {

            std::istringstream iss(values);

            iss >> ins->timeData[row];

            iss >> ins->imuData[row].accelMeas[0] >> ins->imuData[row].accelMeas[1] >> ins->imuData[row].accelMeas[2];

            iss >> ins->imuData[row].gyroMeas[0] >> ins->imuData[row].gyroMeas[1] >> ins->imuData[row].gyroMeas[2];

            iss >> ins->imuData[row].magMeas[0] >> ins->imuData[row].magMeas[1] >> ins->imuData[row].magMeas[2];
           
        }
        row = 0;

    }
  
    for (INS* ins : inss) {
        get_calibrated_meas(*ins);
    }
    // Domain Fusion -> Raw Data
    if (!estimation_fusion) {
        // Transformation in Punkt
        

    }
    // Estimation Fusion
    if (estimation_fusion) {

    }
    return 0;
}

/*coloumn = 0;
           while (coloumn < 10)
           {
               if (coloumn != 9)
               {
                   delimiter = "\t";
               }
               else
               {
                   delimiter = "\r";
               }
               pos = values.find(delimiter);

               data = values.substr(0, pos);
               float number = std::stof(data);
               if (coloumn == 0)
               {
                   ins.timeData[row] = std::stod(data);
               }
               if (coloumn == 1)
               {
                   ins.accMeas[row].x_dot_dot = std::stod(data);
               }
               if (coloumn == 2)
               {
                   ins.accMeas[row].y_dot_dot = std::stod(data);
               }
               if (coloumn == 3)
               {
                   ins.accMeas[row].z_dot_dot = std::stod(data);
               }
               if (coloumn == 4)
               {
                   ins.gyroMeas[row].psi_dot_x = std::stod(data) * M_PI / 180.0;
               }
               if (coloumn == 5)
               {
                   ins.gyroMeas[row].psi_dot_y = std::stod(data) * M_PI / 180.0;
               }
               if (coloumn == 6)
               {
                   ins.gyroMeas[row].psi_dot_z = std::stod(data) * M_PI / 180.0;
               }
               if (coloumn == 7)
               {
                   ins.magMeas[row].B_x = std::stod(data);
               }
               if (coloumn == 8)
               {
                   ins.magMeas[row].B_y = std::stod(data);
               }
               if (coloumn == 9)
               {
                   ins.magMeas[row].B_z = std::stod(data);
               }
               values.erase(0, pos + delimiter.length());
               coloumn++;
           }*/


