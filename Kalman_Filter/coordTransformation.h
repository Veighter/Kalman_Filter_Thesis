#pragma once

#include <Eigen/Core>
#include <cmath>

/// <summary>
/// Class to transform between geodetic, ECEF and NED Coordinates
/// first two methods from https://danceswithcode.net/engineeringnotes/geodetic_to_ecef/geodetic_to_ecef.html
/// </summary>
class CoordTransformer {
private:
	double a = 6378137.0;                    // WGS-84 semi-major axis
	double e2 = 6.6943799901377997e-3;       // WGS-84 first eccentricity squared
	double a1 = 4.2697672707157535e+4;       // a1 = a*e2
	double a2 = 1.8230912546075455e+9;       // a2 = a1*a1
	double a3 = 1.4291722289812413e+2;       // a3 = a1*e2/2
	double a4 = 4.5577281365188637e+9;       // a4 = 2.5*a2
	double a5 = 4.2840589930055659e+4;       // a5 = a1+a3
	double a6 = 9.9330562000986220e-1;       // a6 = 1-e2

	double zp, w2, w, r2, r, s2, c2, s, c, ss;
	double g, rg, rf, u, v, m, f, p, x, y, z;
	double n, lat, lon, alt;

public:
	// Convert Earth-Centered-Earth-Fixed (ECEF) to lat, Lon, Altitude
	// Input is a three element array containing x, y, z in meters
	// Returned array contains lat and lon in radians, and altitude in meters
	Eigen::Vector3d ecef_to_geo(const Eigen::Vector3d& ecef) {
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



	// Convert Lat, Lon, Altitude to Earth-Centered-Earth-Fixed (ECEF)
	// Input is a three element array containing lat, lon (rads) and alt (m)
	// Returned array contains x, y, z in meters
	Eigen::Vector3d geo_to_ecef(const Eigen::Vector3d& geo) {
		Eigen::Vector3d ecef = Eigen::Vector3d::Zero();  // Results go here (x, y, z)
		lat = geo[0] * M_PI / 180.0;
		lon = geo[1] * M_PI / 180.0;
		alt = geo[2];
		n = a / std::sqrt(1 - e2 * std::sin(lat) * std::sin(lat));
		ecef[0] = (n + alt) * std::cos(lat) * std::cos(lon);    // ECEF x
		ecef[1] = (n + alt) * std::cos(lat) * std::sin(lon);    // ECEF y
		ecef[2] = (n * (1 - e2) + alt) * std::sin(lat);         // ECEF z
		return ecef;     // Return x, y, z in ECEF
	}

	// Transformation out of Jekeli
	Eigen::Matrix3d ecef_to_ned_RotationMatrix(Eigen::Vector3d& geo) {
		double phi = geo[0] * M_PI / 180.0; // Latitude in radians
		double lambda = geo[1] * M_PI / 180.0; // Longitude in radians

		Eigen::Matrix3d rotation_matrix;
		rotation_matrix << -std::sin(phi) * std::cos(lambda), -std::sin(phi) * std::sin(lambda), std::cos(phi),
			-std::sin(lambda), std::cos(lambda), 0,
			-std::cos(phi) * std::cos(lambda), -std::cos(phi) * std::sin(lambda), -std::sin(phi);

		return rotation_matrix;
	}

	
};

