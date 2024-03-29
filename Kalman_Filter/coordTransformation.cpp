#include "coordTransformation.h"

Eigen::Vector3d coordTransformation::ecef_to_geo(const Eigen::Vector3d& ecef) {
    Eigen::Vector3d geo;   // Results go here (Lat, Lon, Altitude)
    x = ecef[0];
    y = ecef[1];
    z = ecef[2];
    zp = std::abs(z);
    w2 = x * x + y * y;
    w = std::sqrt(w2);
    r2 = w2 + z * z;
    r = std::sqrt(r2);
    geo[1] = std::atan2(y, x);       // Lon (final)
    s2 = z * z / r2;
    c2 = w2 / r2;
    u = a2 / r;
    v = a3 - a4 / r;
    if (c2 > 0.3) {
        s = (zp / r) * (1.0 + c2 * (a1 + u + s2 * v) / r);
        geo[0] = std::asin(s);      // Lat
        ss = s * s;
        c = std::sqrt(1.0 - ss);
    }
    else {
        c = (w / r) * (1.0 - s2 * (a5 - u - c2 * v) / r);
        geo[0] = std::acos(c);      // Lat
        ss = 1.0 - c * c;
        s = std::sqrt(ss);
    }
    g = 1.0 - e2 * ss;
    rg = a / std::sqrt(g);
    rf = a6 * rg;
    u = w - rg * c;
    v = zp - rf * s;
    f = c * u + s * v;
    m = c * v - s * u;
    p = m / (rf / g + f);
    geo[0] = geo[0] + p;      // Lat
    geo[2] = f + m * p / 2.0; // Altitude
    if (z < 0.0) {
        geo[0] *= -1.0;     // Lat
    }
    return geo;    // Return Lat, Lon, Altitude in that order
}

Eigen::Vector3d coordTransformation::geo_to_ecef(const Eigen::Vector3d& geo) {
    Eigen::Vector3d ecef;  // Results go here (x, y, z)
    lat = geo[0];
    lon = geo[1];
    alt = geo[2];
    n = a / std::sqrt(1 - e2 * std::sin(lat) * std::sin(lat));
    ecef[0] = (n + alt) * std::cos(lat) * std::cos(lon);    // ECEF x
    ecef[1] = (n + alt) * std::cos(lat) * std::sin(lon);    // ECEF y
    ecef[2] = (n * (1 - e2) + alt) * std::sin(lat);         // ECEF z
    return ecef;     // Return x, y, z in ECEF
}