#pragma once

namespace sensorMeas
{

    struct AccelMeas
    {
        double x_dot_dot{};
        double y_dot_dot{};
        double z_dot_dot{};
    };
    struct GyroMeas
    {
        double psi_dot_x{};
        double psi_dot_y{};
        double psi_dot_z{};
    };
    struct MagMeas
    {
        double B_x{};
        double B_y{};
        double B_z{};
    };
    struct GPSMeas
    {
        double x{};
        double y{};
        double z{};
    };

}

