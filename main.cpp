/*
 * Copyright (C) 2025 H. A. Guler
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#include <iostream>
#include <string>
#include <array>
#include <vector>
#include <algorithm>
#include <ctime>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <unordered_map>
#include <algorithm>

#ifdef _WIN32
#include <windows.h>
#else
#include <dirent.h>
#include <sys/stat.h>
#endif

extern "C"
{
#include "SpiceUsr.h"
}

// nlohmann/json is licensed under the terms of MIT License
// see https://github.com/nlohmann/json/blob/develop/LICENSE.MIT
#include "nlohmann/json.hpp"
using json = nlohmann::json;

// kilometers per astronomic unit
double AU = 149597870.7;

// I wrote these backwards once and spent hours 
// trying to debug what was going wrong...
double rad2deg(double x)
{
    return x * 180 / pi_c();
}

double deg2rad(double x)
{
    return x * pi_c() / 180;
}

double clamp(double v, double lo, double hi)
{
    if (lo < v && v < hi)
    {
        return v;
    }
    else if (v < lo)
    {
        return lo;
    }

    return hi;
}

// --- --- --- STRUCT AND CLASS DEFINITIONS --- --- --- [START]
struct Vec3 // extremely self-explanatory
{
    double x, y, z;

    Vec3()
    {
        x = 0;
        y = 0;
        z = 0;
    }

    Vec3(double xp, double yp, double zp)
    {
        x = xp;
        y = yp;
        z = zp;
    }

    Vec3(std::vector<double> vec)
    {
        x = vec[0];
        y = vec[1];
        z = vec[2];
    }

    Vec3(std::array<double, 3Ui64> vec)
    {
        x = vec[0];
        y = vec[1];
        z = vec[2];
    }

    Vec3(double RA, double DEC) // this is equivalent to a spherical2cartezian() function
    {
        x = cos(deg2rad(DEC)) * cos(deg2rad(RA));
        y = cos(deg2rad(DEC)) * sin(deg2rad(RA));
        z = sin(deg2rad(DEC));
    }

    Vec3 operator+(const Vec3& other) const
    {
        return { x + other.x, y + other.y, z + other.z };
    }

    Vec3 operator-(const Vec3& other) const
    {
        return { x - other.x, y - other.y, z - other.z };
    }

    Vec3 operator*(double scalar) const
    {
        return { x * scalar, y * scalar, z * scalar };
    }

    Vec3 operator/(double scalar) const
    {
        return { x / scalar, y / scalar, z / scalar };
    }

    Vec3& operator+=(const Vec3& other)
    {
        x += other.x; y += other.y; z += other.z;
        return *this;
    }

    Vec3 cross(const Vec3& other)
    {
        return Vec3(y * other.z - z * other.y,
            z * other.x - x * other.z,
            x * other.y - y * other.x);
    }

    double dot(const Vec3& other)
    {
        return x * other.x + y * other.y + z * other.z;
    }

    Vec3 normalized()
    {
        return Vec3(x, y, z) / Vec3(x, y, z).mag();
    }

    double mag()
    {
        return sqrt(x * x + y * y + z * z);
    }

    void printout()
    {
        std::cout << "Vec3(" << x << ", " << y << ", " << z << ")\n";
    }
};

// represents a physical gravitational field generating body (Sun + planets usually)
// used for orbit propagation
struct Body
{
    std::string name;
    double GM;
    Vec3 pos;

    Body(std::string namep, double GMp)
    {
        name = namep;
        GM = GMp;
        pos = Vec3();
    }
};

// represents an object which does not generate a notable gravitational field and is
// propagated via Yoshida8 instead of retrieving through SPICE
struct MinorPlanet
{
    Vec3 pos;
    Vec3 vel;

    MinorPlanet(Vec3 p, Vec3 v)
    {
        pos = p;
        vel = v;
    }
};

// ground observatories with MPC obscodes
class Observatory {
public:
    std::string name;
    std::string code;
    double l; // long
    double s; // sin
    double c; // cos

    Observatory()
    {
        name = "None";
        code = "---";
        l = 0;
        s = 0;
        c = 1;
    }

    Observatory(const std::string& name, const std::string& code,
        double l, double c, double s) : name(name), code(code), l(l), s(s), c(c) {}

    std::string toString() const
    {
        return code + " (" + name + ")";
    }
};

// Obs80 (MPC 80-column format) observations
class Observation
{
public:
    std::string obs_str;
    std::string perm;
    std::string prov;
    SpiceDouble et;
    double RA;
    double DEC;
    double mag;
    std::string obs_code;

    Observation(const std::string& line)
    {
        if (line.length() < 79)
        {
            throw std::runtime_error("Invalid observation line length.");
        }

        obs_str = line;

        perm = line.substr(0, 5);
        prov = line.substr(5, 7);
        std::string date_str = line.substr(15, 16);
        std::string RA_str = line.substr(32, 11);
        std::string DEC_str = line.substr(44, 11);
        std::string mag_str = line.substr(65, 4);
        obs_code = line.substr(77, 3);

        int year, month;
        double day_frac;
        std::istringstream ds(date_str);
        ds >> year >> month >> day_frac;

        int day_int = static_cast<int>(day_frac);
        double decimal_day = day_frac - day_int;
        double secs = decimal_day * 86400.0;

        char utc_str[64];
        std::snprintf(utc_str, sizeof(utc_str),
            "%04d-%02d-%02dT00:00:00", year, month, day_int);

        SpiceDouble base_et;
        str2et_c(utc_str, &base_et);
        et = base_et + secs;

        std::istringstream ra_stream(RA_str);
        double ra_h, ra_m, ra_s;
        ra_stream >> ra_h >> ra_m >> ra_s;
        RA = ra_h * 15.0 + ra_m * 15.0 / 60.0 + ra_s * 15.0 / 3600.0;

        // Parse DEC
        int dec_sign = (line[44] == '-') ? -1 : 1;
        std::istringstream dec_stream(DEC_str);
        double dec_d, dec_m, dec_s;
        dec_stream >> dec_d >> dec_m >> dec_s;
        DEC = dec_sign * (std::abs(dec_d) + dec_m / 60.0 + dec_s / 3600.0);

        // Parse magnitude
        try {
            mag = std::stod(mag_str);
        }
        catch (...) {
            mag = 0.0;
        }
    }

    std::string toString() const {
        return obs_str;
    }
};

// Return type: [body_index][position_or_velocity][x/y/z]
using StateMatrix = std::array<std::array<std::array<double, 3>, 2>, 9>; // I don't even know why I did it this way
// --- --- --- STRUCT AND CLASS DEFINITIONS --- --- --- [END]

// --- --- --- MISC UTILS --- --- --- [START]
void loadAllKernels(const std::string& directory)
{
#ifdef _WIN32
    std::string search_path = directory + "\\*.*";
    WIN32_FIND_DATAA fd;
    HANDLE hFind = ::FindFirstFileA(search_path.c_str(), &fd);

    if (hFind == INVALID_HANDLE_VALUE) {
        std::cerr << "Unable to open directory: " << directory << '\n';
        return;
    }

    do {
        if (!(fd.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY)) {
            std::string filepath = directory + "\\" + fd.cFileName;
            furnsh_c(filepath.c_str());
            // std::cout << "Loaded kernel: " << filepath << '\n';
        }
    } while (::FindNextFileA(hFind, &fd));

    ::FindClose(hFind);
#else
    DIR* dir = opendir(directory.c_str());
    if (!dir) {
        std::cerr << "Unable to open directory: " << directory << '\n';
        return;
    }

    struct dirent* entry;
    while ((entry = readdir(dir)) != NULL) {
        std::string filename = entry->d_name;

        // Skip "." and ".."
        if (filename == "." || filename == "..") {
            continue;
        }

        std::string filepath = directory + "/" + filename;

        struct stat path_stat;
        if (stat(filepath.c_str(), &path_stat) == 0 && S_ISREG(path_stat.st_mode)) {
            furnsh_c(filepath.c_str());
        }
    }

    closedir(dir);
#endif
}

std::string trim(const std::string& str) // trim trailing/leading whitespace, nothing fancy, stole it from somewhere
{
    size_t start = str.find_first_not_of(" \t\r\n");
    size_t end = str.find_last_not_of(" \t\r\n");
    return (start == std::string::npos) ? "" : str.substr(start, end - start + 1);
}

int digitsBeforeDecimal(double value) {
    if (value == 0.0) return 1;
    value = std::fabs(value); // Make it positive
    return static_cast<int>(std::floor(std::log10(value))) + 1;
}

double etToJD(SpiceDouble et)
{
    return 2451545.0 + et / 86400.0;
}

SpiceDouble JDToEt(double JD)
{
    return (JD - 2451545.0) * 86400.0;
}

std::string alignDecimal(const std::string& prefix, double value, int precision = 5) {
    std::ostringstream oss;

    // Width for just the number, not including the prefix
    const int number_field_width = 4 + 1 + precision; // 4 digits before decimal, 1 dot, rest after
    oss << prefix;

    // Format number with fixed precision and right-align in number_field_width
    std::ostringstream numstream;
    numstream << std::fixed << std::setprecision(precision) << value;
    std::string num_str = numstream.str();

    int pad_spaces = number_field_width - static_cast<int>(num_str.length());
    if (pad_spaces > 0) {
        oss << std::string(pad_spaces, ' ');
    }

    oss << num_str;
    return oss.str();
}

// this is the fun and simpler part of the equation
std::vector<double> stateVector2Kepler(Vec3 r_equ, Vec3 v_equ, SpiceDouble et, double mu = 1.3271244004193938E+11)
{
    // convert from equatorial to ecliptic frame
    double rot[3][3];
    pxform_c("J2000", "ECLIPJ2000", et, rot);

    double r_equ_arr[3] = { r_equ.x, r_equ.y, r_equ.z };
    double v_equ_arr[3] = { v_equ.x, v_equ.y, v_equ.z };

    double r_arr[3] = { 0, 0, 0 };
    double v_arr[3] = { 0, 0, 0 };

    mxv_c(rot, r_equ_arr, r_arr);
    mxv_c(rot, v_equ_arr, v_arr);

    Vec3 r = Vec3(r_arr[0], r_arr[1], r_arr[2]);
    Vec3 v = Vec3(v_arr[0], v_arr[1], v_arr[2]);

    // now calculate orbital params
    double r_mag = r.mag();
    double v_mag = v.mag();

    Vec3 h = r.cross(v);
    double h_mag = h.mag();

    double inclination = rad2deg(acos(h.z / h_mag));

    Vec3 k = Vec3(0, 0, 1);
    Vec3 n = k.cross(h);
    double n_mag = n.mag();

    double omega;
    double eccentricity;
    double arg_periapsis;
    double true_anomaly;
    double sma;
    double mean_anomaly;

    if (n_mag != 0)
    {
        omega = rad2deg(acos(n.x / n_mag));
        if (n.y < 0)
        {
            omega = 360 - omega;
        }
    }
    else
    {
        omega = 0;
    }

    Vec3 e_vec = (v.cross(h) - r * mu / r_mag) * (1 / mu);
    eccentricity = e_vec.mag();

    if (n_mag != 0)
    {
        if (eccentricity != 0)
        {
            arg_periapsis = rad2deg(acos(n.dot(e_vec) / (n_mag * eccentricity)));
            if (e_vec.z < 0)
            {
                arg_periapsis = 360 - arg_periapsis;
            }
        }
        else
        {
            arg_periapsis = 0;
        }
    }
    else
    {
        arg_periapsis = 0;
    }

    if (eccentricity != 0)
    {
        true_anomaly = rad2deg(acos(e_vec.dot(r) / (eccentricity * r_mag)));
        if (r.dot(v) < 0)
        {
            true_anomaly = 360 - true_anomaly;
        }
    }
    else
    {
        true_anomaly = rad2deg(acos(r.normalized().dot(v.normalized())));
    }

    double specific_energy = v_mag * v_mag / 2 - mu / r_mag;
    if (abs(eccentricity - 1) > 1e-8)
    {
        sma = -mu / (2 * specific_energy);
    }
    else
    {
        sma = 999999;
    }

    if (eccentricity < 1)
    {
        double E = 2 * atan(tan(deg2rad(true_anomaly) / 2) * sqrt((1 - eccentricity) / (1 + eccentricity)));
        if (E < 0)
        {
            E = E + 2 * pi_c();
        }

        mean_anomaly = rad2deg(E - eccentricity * sin(E));
    }
    else if (eccentricity > 1)
    {
        double F = 2 * atanh(tan(deg2rad(true_anomaly) / 2) * sqrt((eccentricity - 1) / (eccentricity + 1)));
        mean_anomaly = rad2deg(eccentricity * sinh(F) - F);
    }
    else
    {
        mean_anomaly = -1.0; // random val.
    }

    return std::vector<double> {sma, eccentricity, inclination, omega, arg_periapsis, true_anomaly, mean_anomaly};

}

double computeEccentricAnomaly(double M, double e, double tol=1E-10, double MAXITER=100)
{
    double ecc_anomaly = 0;
    if (e < 0.8)
    {
        ecc_anomaly = M;
    }
    else
    {
        ecc_anomaly = pi_c();
    }

    for (int i = 0; i < MAXITER; i++)
    {
        double delta = (ecc_anomaly - e * sin(ecc_anomaly) - M) / (1 - e * cos(ecc_anomaly));
        ecc_anomaly -= delta;
        if (abs(delta) < tol)
        {
            break;
        }
    }
    
    return ecc_anomaly;
}

std::pair<Vec3, Vec3> Kepler2StateVector(double a_AU, double e, double i_deg, double omega_deg, 
    double arg_peri_deg, double M_deg, double epoch_JD, double mu = 1.3271244004193938E+11)
{
    // unit conversions
    double M = deg2rad(M_deg);
    double a = a_AU * AU;
    double i = deg2rad(i_deg);
    double omega = deg2rad(omega_deg);
    double arg_peri = deg2rad(arg_peri_deg);
    SpiceDouble et = JDToEt(epoch_JD);

    // E
    double ecc_anomaly = computeEccentricAnomaly(M, e);
    double cosE = cos(ecc_anomaly);
    double sinE = sin(ecc_anomaly);

    // pos. and vel. in perifocal frame
    double r = a * (1 - e * cosE);

    double xp = a * (cosE - e);
    double yp = a * sqrt(1 - e * e) * sinE;

    double vxp = -sqrt(mu * a) / r * sinE;
    double vyp = sqrt(mu * a * (1 - e * e)) / r * cosE;

    double r_pf[3] = {xp, yp, 0.0};
    double v_pf[3] = {vxp, vyp, 0.0};

    // rot. mtx. perifocal --> inertial
    double cosO = cos(omega);
    double sinO = sin(omega);
    double cos_peri = cos(arg_peri);
    double sin_peri = sin(arg_peri);
    double cosi = cos(i);
    double sini = sin(i);

    double R[3][3] = { cosO * cos_peri - sinO * sin_peri * cosi, -cosO * sin_peri - sinO * cos_peri * cosi,  sinO * sini ,
                       sinO * cos_peri + cosO * sin_peri * cosi, -sinO * sin_peri + cosO * cos_peri * cosi, -cosO * sini ,
                       sin_peri * sini,                           cos_peri* sini,                            cosi };


    SpiceDouble r_inertial[3];
    SpiceDouble v_inertial[3];

    mxv_c(R, r_pf, r_inertial);
    mxv_c(R, v_pf, v_inertial);

    // rot. mtx. ecliptic --> equatorial
    SpiceDouble R_eq[3][3];
    pxform_c("ECLIPJ2000", "J2000", et, R_eq);

    SpiceDouble r_equat[3], v_equat[3];

    mxv_c(R_eq, r_inertial, r_equat);
    mxv_c(R_eq, v_inertial, v_equat);

    Vec3 r_res = Vec3(r_equat[0], r_equat[1], r_equat[2]);
    Vec3 v_res = Vec3(v_equat[0], v_equat[1], v_equat[2]);

    return std::pair<Vec3, Vec3> {r_res, v_res};
}

std::vector<double> cartezian2spherical(Vec3 v)
{
    double d = v.mag();
    double RA = rad2deg(atan2(v.y, v.x));
    double DEC = rad2deg(asin(v.z / d));

    if (RA < 0)
    {
        RA = RA + 360;
    }

    return std::vector<double> {RA, DEC};
}

std::string RATohms(double ra_deg) {
    double total_hours = ra_deg / 15.0;
    int hours = static_cast<int>(total_hours);
    double remainder_minutes = (total_hours - hours) * 60.0;
    int minutes = static_cast<int>(remainder_minutes);
    double seconds = (remainder_minutes - minutes) * 60.0;

    std::ostringstream oss;
    oss << std::setfill('0') << std::setw(2) << hours << ":"
        << std::setw(2) << minutes << ":"
        << std::fixed << std::setprecision(1) << std::setw(4) << seconds;
    return oss.str();
}

std::string DECToDms(double dec_deg) {
    char sign = (dec_deg >= 0) ? '+' : '-';
    double abs_deg = std::fabs(dec_deg);
    int degrees = static_cast<int>(abs_deg);
    double remainder_minutes = (abs_deg - degrees) * 60.0;
    int minutes = static_cast<int>(remainder_minutes);
    double seconds = (remainder_minutes - minutes) * 60.0;

    std::ostringstream oss;
    oss << sign << std::setfill('0') << std::setw(2) << degrees << ":"
        << std::setw(2) << minutes << ":"
        << std::fixed << std::setprecision(1) << std::setw(4) << seconds;
    return oss.str();
}

std::string doubleToResidualStr(double value) {
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(1) << std::abs(value);
    // why is the sign printed at the end of the number on MPECs, wtf
    oss << (value >= 0 ? '+' : '-');
    return oss.str();
}

// I could've made a KeplerOrbit class or whatever but this works too
void printOrbitalElements(std::vector<double> orbel)
{
    std::cout << "a: " << orbel[0] / AU << " e: " << orbel[1] << " i: " << orbel[2] << " Node: " << orbel[3] << " Peri: " << orbel[4] << " M: " << orbel[6] << "\n";
}

std::vector<double> flattenResiduals(const std::pair<std::vector<double>, std::vector<double>>& residual_pair) {
    const std::vector<double>& ra_res = residual_pair.first;
    const std::vector<double>& dec_res = residual_pair.second;
    std::vector<double> flat;
    flat.reserve(ra_res.size() + dec_res.size());

    for (std::size_t i = 0; i < ra_res.size(); ++i) {
        flat.push_back(ra_res[i]);
        flat.push_back(dec_res[i]);
    }

    return flat;
}

// --- --- --- MISC UTILS --- --- --- [END]

// --- --- --- OBSERVATIONS AND OBSERVERS HANDLING --- --- --- [START]
std::vector<Observatory> readObsCodes(const std::string& filepath = "data/ObsCodes.dat")
{
    std::cout << "Reading ObsCodes... " << std::flush;
    std::vector<Observatory> obscodes;
    int skips = 0;

    std::ifstream infile(filepath);
    if (!infile.is_open())
    {
        std::cerr << "Could not open file: " << filepath << "\n";
        return obscodes;
    }

    std::string line;
    std::getline(infile, line);

    while (std::getline(infile, line))
    {
        try
        {
            std::string code_str = line.substr(0, 3);
            std::string l_str = line.substr(6, 7);
            std::string c_str = line.substr(13, 8);
            std::string s_str = line.substr(21, 9);
            std::string name_str = trim(line.substr(30));

            obscodes.push_back(Observatory(name_str, code_str, strtod(l_str.c_str(), NULL), strtod(c_str.c_str(), NULL), strtod(s_str.c_str(), NULL)));
        }
        catch (...)
        {
            ++skips;
            continue;
        }
    }

    std::cout << "Done, skipped " << skips << " obscodes.\n";
    return obscodes;
}

std::vector<Observation> readObsFile(const std::string& filename = "primary.obs") {
    std::cout << "Reading observations file: " << filename << " ... " << std::flush;
    std::vector<Observation> observations;

    std::ifstream infile(filename);
    if (!infile.is_open()) {
        std::cerr << "Cannot open file: " << filename << "\n";
        return observations;
    }

    std::string line;
    while (std::getline(infile, line)) {
        try {
            Observation obs(line);
            observations.push_back(obs);
        }
        catch (...) {
            // Could not parse line, skip
        }
    }

    std::cout << "Done.\n\n";
    return observations;
}

std::vector<SpiceDouble> readDateFile(const std::string& filename = "dates.txt")
{
    std::cout << "Reading ephemeris dates file: " << filename << " ... " << std::flush;
    std::vector<SpiceDouble> ts;

    std::ifstream infile(filename);
    if (!infile.is_open()) {
        std::cerr << "Cannot open file: " << filename << "\n";
        return ts;
    }

    std::string line;
    while (std::getline(infile, line)) {
        try {
            SpiceDouble et;
            str2et_c(line.c_str(), &et);
            ts.push_back(et);
        }
        catch (...) {
            // Could not parse line, skip
        }
    }

    std::cout << "Done.\n\n";
    return ts;
}

// returns { permanent code, provisional code, orbital elements, p0, v0, epoch_JD }
std::tuple<std::string, std::string, std::vector<double>, Vec3, Vec3, double> readMinorPlanetFile(const std::string& filename = "mp.json")
{
    std::cout << "Reading minor planet data file: " << filename << "... " << std::flush;
    
    std::string perm;
    std::string prov;
    std::vector<double> orbital_elements; orbital_elements.resize(6);
    Vec3 p0;
    Vec3 v0;
    double epoch_JD;

    std::ifstream infile(filename);
    if (!infile) {
        std::cerr << "Could not open minor planet data file: " << filename << ".\n";
    }

    json j;
    infile >> j;

    perm = j["permanent_designation"];
    prov = j["provisional_designation"];
    orbital_elements[0] = j["orbital_elements"]["a"] * AU;
    orbital_elements[1] = j["orbital_elements"]["e"];
    orbital_elements[2] = j["orbital_elements"]["i"];
    orbital_elements[3] = j["orbital_elements"]["node"];
    orbital_elements[4] = j["orbital_elements"]["peri"];
    orbital_elements[5] = j["orbital_elements"]["M"];
    p0.x = j["state_vectors"]["x"];
    p0.y = j["state_vectors"]["y"];
    p0.z = j["state_vectors"]["z"];
    v0.x = j["state_vectors"]["vx"];
    v0.y = j["state_vectors"]["vy"];
    v0.z = j["state_vectors"]["vz"];
    epoch_JD = j["epoch"];

    std::cout << "Done.\n";

    return std::tuple<std::string, std::string, std::vector<double>, Vec3, Vec3, double> {perm, prov, orbital_elements, p0, v0, epoch_JD};
}

// lat - lon to an geocentric Earth-fixed vector (hence ECEF)
Vec3 geodeticToECEF(double lat_deg, double lon_deg)
{
    double a = 6378.137;
    double f = 1 / 298.257223563; // flattening. honestly not even needed.
    double e2 = f * (2 - f);

    double lat = deg2rad(lat_deg);
    double lon = deg2rad(lon_deg);

    double slat = sin(lat);
    double N = a / sqrt(1 - e2 * slat * slat);

    double clat = cos(lat);

    double x = N * clat * cos(lon);
    double y = N * clat * sin(lon);
    double z = N * (1 - e2) * slat;

    return Vec3(x, y, z);
}

// give obscode and time, get heliocentric position (equatorial J2000)
Vec3 getObserverPos(Observatory obsv, SpiceDouble et)
{
    SpiceDouble state[6], lt;
    spkezr_c("EARTH", et, "J2000", "NONE", "SOLAR SYSTEM BARYCENTER", state, &lt);
    Vec3 earth_pos = Vec3(state[0], state[1], state[2]);

    double lat_rad = atan2(obsv.s, obsv.c);
    double lat_deg = lat_rad * 180 / pi_c();

    Vec3 ecef = geodeticToECEF(lat_deg, obsv.l);
    SpiceDouble rot[3][3];
    pxform_c("ITRF93", "J2000", et, rot);

    SpiceDouble ecef_sd[3] = { ecef.x, ecef.y, ecef.z };
    SpiceDouble observer_j2000[3];
    mxv_c(rot, ecef_sd, observer_j2000);

    return earth_pos + Vec3(ecef_sd[0], ecef_sd[1], ecef_sd[2]);
}

// give obscode and time, get heliocentric velocity (equatorial J2000)
Vec3 getObserverVel(Observatory obsv, SpiceDouble et)
{
    SpiceDouble state[6], lt;
    spkezr_c("EARTH", et, "J2000", "NONE", "SOLAR SYSTEM BARYCENTER", state, &lt);
    Vec3 earth_vel = Vec3(state[3], state[4], state[5]);

    double lat_rad = atan2(obsv.s, obsv.c);
    double lat_deg = lat_rad * 180 / pi_c();

    Vec3 ecef = geodeticToECEF(lat_deg, obsv.l);
    SpiceDouble xform[6][6];
    sxform_c("ITRF93", "J2000", et, xform);

    SpiceDouble ecef_mxvform[6] = { ecef.x, ecef.y, ecef.z, 0, 0, 0 };
    SpiceDouble mxvg_out[6];
    mxvg_c(xform, ecef_mxvform, 6, 6, mxvg_out);

    Vec3 obsv_geo_vel = Vec3(mxvg_out[3], mxvg_out[4], mxvg_out[5]);

    Vec3 obsv_vel = earth_vel + obsv_geo_vel;
    return obsv_vel;
}

// give time, observatory and heliocentric object position, get RA and DEC
std::vector<double> getRADEC(SpiceDouble et, Observatory obsv, Vec3 p)
{
    Vec3 obsv_pos = getObserverPos(obsv, et);
    Vec3 rel_pos = p - obsv_pos;

    return cartezian2spherical(rel_pos);
}

std::vector<double> getGeocentricRADEC(SpiceDouble et, Vec3 p)
{
    SpiceDouble state[6], lt;
    spkezr_c("EARTH", et, "J2000", "NONE", "SOLAR SYSTEM BARYCENTER", state, &lt);
    Vec3 earth_pos = Vec3(state[0], state[1], state[2]);

    Vec3 rel_pos = p - earth_pos;

    return cartezian2spherical(rel_pos);
}
// --- --- --- OBSERVATIONS AND OBSERVERS HANDLING --- --- --- [END]

// --- --- --- ORBIT PROPAGATION --- --- --- [START]
// I do regret doing this the way I did it, but it doesn't bite, it can stay
StateMatrix getSolarSystemStates(SpiceDouble et)
{
    const char* bodies[9] = {
        "SUN",                 // index 0
        "MERCURY BARYCENTER",  // index 1
        "VENUS BARYCENTER",    // index 2
        "EARTH BARYCENTER",    // index 3
        "MARS BARYCENTER",     // index 4
        "JUPITER BARYCENTER",  // index 5
        "SATURN BARYCENTER",   // index 6
        "URANUS BARYCENTER",   // index 7
        "NEPTUNE BARYCENTER"   // index 8
    };

    StateMatrix states{};

    for (int i = 0; i < 9; ++i)
    {
        SpiceDouble state[6];
        SpiceDouble lt;

        spkezr_c(bodies[i], et, "J2000", "NONE", "SOLAR SYSTEM BARYCENTER", state, &lt);

        for (int j = 0; j < 3; ++j)
        {
            states[i][0][j] = state[j];     // Position (km)
            states[i][1][j] = state[j + 3]; // Velocity (km/s)
        }
    }

    return states;
}

// not running numbers on Schwardschild metrics
Vec3 gravAccel(MinorPlanet& mp, std::vector<Body>& bodies)
{
    Vec3 accel{};

    for (Body& body : bodies)
    {
        double dist = (body.pos - mp.pos).mag();
        Vec3 grav_dir = (body.pos - mp.pos) / dist;
        double grav_mag = body.GM / (dist * dist);

        accel = accel + grav_dir * grav_mag;
    }

    return accel;
}

// symplectic 8th order orbit integrator
// RKN or whatever also works probably
// I like this because I already use this a lot in other projects too, 
// so I don't have to spend brain power again implementing something else.
void stepYoshida8(MinorPlanet& mp, std::vector<Body>& bodies, SpiceDouble date0_et, double dt)
{
    constexpr double w1 = 0.311790812418427e0;
    constexpr double w2 = -0.155946803821447e1;
    constexpr double w3 = -0.167896928259640e1;
    constexpr double w4 = 0.166335809963315e1;
    constexpr double w5 = -0.106458714789183e1;
    constexpr double w6 = 0.136934946416871e1;
    constexpr double w7 = 0.629030650210433e0;
    constexpr double w0 = 1.65899088454396;

    constexpr double ds[15] = { w7, w6, w5, w4, w3, w2, w1, w0, w1, w2, w3, w4, w5, w6, w7 };
    constexpr double cs[16] = { 0.3145153251052165, 0.9991900571895715, 0.15238115813844, 0.29938547587066, -0.007805591481624963,
          -1.619218660405435, -0.6238386128980216, 0.9853908484811935, 0.9853908484811935, -0.6238386128980216,
          -1.619218660405435, -0.007805591481624963, 0.29938547587066, 0.15238115813844, 0.9991900571895715,
          0.3145153251052165 };

    SpiceDouble et = date0_et;

    for (int i = 0; i < 15; ++i)
    {
        mp.pos += mp.vel * (cs[i] * dt);

        // Advance SPICE time
        et += cs[i] * dt;

        // Update body positions
        for (Body& body : bodies)
        {
            SpiceDouble state[6], lt;
            spkezr_c(body.name.c_str(), et, "J2000", "NONE", "SOLAR SYSTEM BARYCENTER", state, &lt);
            body.pos = { state[0], state[1], state[2] };
        }

        Vec3 accel = gravAccel(mp, bodies);
        mp.vel += accel * (ds[i] * dt);
    }

    mp.pos += mp.vel * (cs[15] * dt);
}

// orbit propagator main function
std::vector<Vec3> propagate(Vec3 p0, Vec3 v0, SpiceDouble date_init, SpiceDouble date_final, double dt)
{
    MinorPlanet mp = MinorPlanet(p0, v0);

    const char* body_names[9] = {
        "SUN",
        "MERCURY BARYCENTER",
        "VENUS BARYCENTER",
        "EARTH BARYCENTER",
        "MARS BARYCENTER",
        "JUPITER BARYCENTER",
        "SATURN BARYCENTER",
        "URANUS BARYCENTER",
        "NEPTUNE BARYCENTER"
    };

    const double body_GMs[9] = {
        1.3271244004193938E+11,
        2.2031780000000021E+04,
        3.2485859200000006E+05,
        4.0350323550225981E+05,
        4.2828375214000022E+04,
        1.2671276480000021E+08,
        3.7940585200000003E+07,
        5.7945486000000080E+06,
        6.8365271005800236E+06
    };

    std::vector<Body> bodies;
    StateMatrix system_state = getSolarSystemStates(date_init);

    for (int i = 0; i < 8; i++)
    {
        bodies.push_back(Body(body_names[i], body_GMs[i]));
        bodies[i].pos = Vec3(system_state[i][0]);
    }

    double time_interval = date_final - date_init;

    if (dt <= 0)
    {
        dt = time_interval / 4;
        while (dt > 10 * 86400)
        {
            dt = dt / 2;
        }
    }

    int N_cycles = (int)(time_interval / dt);
    double dt_last = time_interval - N_cycles * dt;

    for (int cycle = 0; cycle < N_cycles; cycle++)
    {
        SpiceDouble cycle_date = date_init + cycle * dt;
        stepYoshida8(mp, bodies, cycle_date, dt);
    }

    // fractional last step if needed
    if (dt_last > 1e-6)
    {
        SpiceDouble last_step_date = date_init + N_cycles * dt;
        stepYoshida8(mp, bodies, last_step_date, dt_last);
    }

    return std::vector<Vec3> {Vec3(mp.pos), Vec3(mp.vel)};
}

// very little difference with propagate(), only in handling of dt (time step) value
std::vector<Vec3> backpropagate(Vec3 p0, Vec3 v0, SpiceDouble date_init, SpiceDouble date_final, double dt)
{
    MinorPlanet mp = MinorPlanet(p0, v0);

    const char* body_names[9] = {
        "SUN",
        "MERCURY BARYCENTER",
        "VENUS BARYCENTER",
        "EARTH BARYCENTER",
        "MARS BARYCENTER",
        "JUPITER BARYCENTER",
        "SATURN BARYCENTER",
        "URANUS BARYCENTER",
        "NEPTUNE BARYCENTER"
    };

    const double body_GMs[9] = {
        1.3271244004193938E+11,
        2.2031780000000021E+04,
        3.2485859200000006E+05,
        4.0350323550225981E+05,
        4.2828375214000022E+04,
        1.2671276480000021E+08,
        3.7940585200000003E+07,
        5.7945486000000080E+06,
        6.8365271005800236E+06
    };

    std::vector<Body> bodies;
    StateMatrix system_state = getSolarSystemStates(date_init);

    for (int i = 0; i < 8; i++)
    {
        bodies.push_back(Body(body_names[i], body_GMs[i]));
        bodies[i].pos = Vec3(system_state[i][0]);
    }

    double time_interval = date_final - date_init;

    if (dt >= 0)
    {
        dt = time_interval / 4;
        while (dt > -10 * 86400)
        {
            dt /= 2;
        }
    }

    int N_cycles = (int)(time_interval / dt);
    double dt_last = time_interval - N_cycles * dt;

    for (int cycle = 0; cycle < N_cycles; cycle++)
    {
        SpiceDouble cycle_date = date_init + cycle * dt;
        stepYoshida8(mp, bodies, cycle_date, dt);
    }

    if (abs(dt_last) > 1e-6)
    {
        SpiceDouble last_step_date = date_init + N_cycles * dt;
        stepYoshida8(mp, bodies, last_step_date, dt_last);
    }

    return std::vector<Vec3> {Vec3(mp.pos), Vec3(mp.vel)};
}
// --- --- --- ORBIT PROPAGATION --- --- --- [END]

// get true state vector at observation time given apparent position
std::vector<Vec3> computeTrueStateGeocentric(Vec3 pa, Vec3 va, SpiceDouble et_obs, double tol = 1, int MAXITER = 5)
{
    SpiceDouble lt_ignore;
    SpiceDouble obsv_sv[6];
    spkezr_c("EARTH", et_obs, "J2000", "NONE", "SOLAR SYSTEM BARYCENTER", obsv_sv, &lt_ignore);
    Vec3 obsv_pos = Vec3(obsv_sv[0], obsv_sv[1], obsv_sv[2]); // this is known
    Vec3 obsv_vel = Vec3(obsv_sv[3], obsv_sv[4], obsv_sv[5]);

    double c = 299792.458; // [km s-1] speed of light

    Vec3 p_rel = pa - obsv_pos;
    SpiceDouble pobj[3] = { p_rel.x, p_rel.y, p_rel.z };
    SpiceDouble vobs[3] = { obsv_vel.x, obsv_vel.y, obsv_vel.z };
    SpiceDouble corpos[3];
    stlabx_c(pobj, vobs, corpos);  // <--- UNDO stellar aberration

    Vec3 p_corr = obsv_pos + Vec3(corpos[0], corpos[1], corpos[2]);

    Vec3 pt = p_corr;
    Vec3 vt = va;

    double lt = (pa - obsv_pos).mag() / c; // initial light time guess
    double lt_prev = 0;
    for (int idx_iter = 0; idx_iter < MAXITER; idx_iter++)
    {

        /*
        *
        *  Do NOT propagate this with a full orbit propagator! Way too expensive!
        *
        std::vector<Vec3> sv_true = propagate(pa, va, et_obs, et_obs + lt, -1);
        pt = sv_true[0];
        vt = sv_true[1];
        *
        * Just use a linear approximation, it's good enough for nearly all cases
        *
        */
        pt = p_corr + va * lt;

        lt = (pt - obsv_pos).mag() / c;
        double lt_err = abs(lt - lt_prev);

        if (lt_err < tol)
        {
            break;
        }

        lt_prev = lt;
    }

    return { pt, vt };
}

// get apparent position at given observation time and true position
std::vector<Vec3> computeApparentStateGeocentric(Vec3 pt, Vec3 vt, SpiceDouble et_obs, double tol = 1, int MAXITER = 5)
{
    SpiceDouble lt_ignore;
    SpiceDouble obsv_sv[6];
    spkezr_c("EARTH", et_obs, "J2000", "NONE", "SOLAR SYSTEM BARYCENTER", obsv_sv, &lt_ignore);
    Vec3 obsv_pos = Vec3(obsv_sv[0], obsv_sv[1], obsv_sv[2]); // this is known
    double c = 299792.458; // [km s-1] speed of light

    Vec3 pa = pt;
    Vec3 va = vt;

    double lt = (pt - obsv_pos).mag() / c; // initial light time guess
    double lt_prev = 0;
    for (int idx_iter = 0; idx_iter < MAXITER; idx_iter++)
    {
        /*
        *
        * No no no, full propagation is too expensive for this inner iteration
        *
        std::vector<Vec3> sv_apparent = backpropagate(pt, vt, et_obs, et_obs - lt, -1);
        pa = sv_apparent[0];
        va = sv_apparent[1];
        *
        * just use linear approx., good enough for most cases
        *
        */
        pa = pt - vt * lt;

        lt = (pa - obsv_pos).mag() / c;
        double lt_err = abs(lt - lt_prev);

        if (lt_err < tol)
        {
            break;
        }

        lt_prev = lt;
    }

    // now we know the so-called apparent position (apparent position corrected for light-time)
    // now apply stellar aberration (not the same thing as so-called planetary abberation)
    Vec3 p_rel = pa - obsv_pos;
    SpiceDouble pobj[3] = { p_rel.x, p_rel.y, p_rel.z };
    SpiceDouble vobs[3] = { obsv_sv[3], obsv_sv[4], obsv_sv[5] };
    SpiceDouble apppos[3]; // stellar-aberration corrected position relative to observer
    stelab_c(pobj, vobs, apppos);

    Vec3 pa_corr = obsv_pos + Vec3(apppos[0], apppos[1], apppos[2]);

    return { pa_corr, va };
}

// get true state vector at observation time given apparent position
std::vector<Vec3> computeTrueState(Vec3 pa, Vec3 va, SpiceDouble et_obs, Observatory& obsv, double tol = 1, int MAXITER = 5)
{
    Vec3 obsv_pos = getObserverPos(obsv, et_obs); // this is known
    Vec3 obsv_vel = getObserverVel(obsv, et_obs);
    double c = 299792.458; // [km s-1] speed of light

    // undo stellar aberration
    Vec3 p_rel = pa - obsv_pos;
    SpiceDouble pobj[3] = { p_rel.x, p_rel.y, p_rel.z };
    SpiceDouble vobs[3] = { obsv_vel.x, obsv_vel.y, obsv_vel.z };
    SpiceDouble corpos[3];
    stlabx_c(pobj, vobs, corpos);  // <--- UNDO stellar aberration

    Vec3 p_corr = obsv_pos + Vec3(corpos[0], corpos[1], corpos[2]);

    Vec3 pt = p_corr;
    Vec3 vt = va;

    double lt = (pa - obsv_pos).mag() / c; // initial light time guess
    double lt_prev = 0;
    for (int idx_iter = 0; idx_iter < MAXITER; idx_iter++)
    {

        /*
        *
        *  Do NOT propagate this with a full orbit propagator! Way too expensive!
        *
        std::vector<Vec3> sv_true = propagate(pa, va, et_obs, et_obs + lt, -1);
        pt = sv_true[0];
        vt = sv_true[1];
        *
        * Just use a linear approximation, it's good enough for nearly all cases
        *
        */
        pt = p_corr + va * lt;

        lt = (pt - obsv_pos).mag() / c;
        double lt_err = abs(lt - lt_prev);

        if (lt_err < tol)
        {
            break;
        }

        lt_prev = lt;
    }

    return { pt, vt };
}

// get apparent position at given observation time and true position
std::vector<Vec3> computeApparentState(Vec3 pt, Vec3 vt, SpiceDouble et_obs, Observatory& obsv, double tol = 1, int MAXITER = 5)
{
    Vec3 obsv_pos = getObserverPos(obsv, et_obs); // this is known
    double c = 299792.458; // [km s-1] speed of light

    Vec3 pa = pt;
    Vec3 va = vt;

    double lt = (pt - obsv_pos).mag() / c; // initial light time guess
    double lt_prev = 0;
    for (int idx_iter = 0; idx_iter < MAXITER; idx_iter++)
    {
        /*
        *
        * No no no, full propagation is too expensive for this inner iteration
        *
        std::vector<Vec3> sv_apparent = backpropagate(pt, vt, et_obs, et_obs - lt, -1);
        pa = sv_apparent[0];
        va = sv_apparent[1];
        *
        * just use linear approx., good enough for most cases
        *
        */
        pa = pt - vt * lt;

        lt = (pa - obsv_pos).mag() / c;
        double lt_err = abs(lt - lt_prev);

        if (lt_err < tol)
        {
            break;
        }

        lt_prev = lt;
    }

    // now we know the so-called apparent position (apparent position corrected for light-time)
    // now apply stellar aberration (not the same thing as so-called planetary abberation)
    Vec3 obsv_vel = getObserverVel(obsv, et_obs);
    Vec3 p_rel = pa - obsv_pos;
    SpiceDouble pobj[3] = { p_rel.x, p_rel.y, p_rel.z };
    SpiceDouble vobs[3] = { obsv_vel.x, obsv_vel.y, obsv_vel.z };
    SpiceDouble apppos[3]; // stellar-aberration corrected position relative to observer
    stelab_c(pobj, vobs, apppos);

    Vec3 pa_corr = obsv_pos + Vec3(apppos[0], apppos[1], apppos[2]);

    return { pa_corr, va };
}


// feed apparent position to this function
std::pair<double, double> computePhaseAndElongationGeocentric(Vec3 pa, SpiceDouble et)
{
    SpiceDouble sv_sun[6], lt_ignore;
    spkezr_c("SUN", et, "J2000", "NONE", "SOLAR SYSTEM BARYCENTER", sv_sun, &lt_ignore);
    Vec3 p_sun = Vec3(sv_sun[0], sv_sun[1], sv_sun[2]);

    SpiceDouble sv_earth[6];
    spkezr_c("EARTH", et, "J2000", "NONE", "SOLAR SYSTEM BARYCENTER", sv_earth, &lt_ignore);
    Vec3 p_earth = Vec3(sv_earth[0], sv_earth[1], sv_earth[2]);

    Vec3 obj_to_sun = p_sun - pa;
    Vec3 obj_to_obs = p_earth - pa;

    double cos_phase = obj_to_sun.normalized().dot(obj_to_obs.normalized());
    double phase_angle_deg = rad2deg(acos(clamp(cos_phase, -1.0, 1.0)));

    Vec3 obs_to_sun = p_sun - p_earth;
    Vec3 obs_to_obj = pa - p_earth;

    double cos_elong = obs_to_sun.normalized().dot(obs_to_obj.normalized());
    double elongation_deg = rad2deg(acos(clamp(cos_elong, -1.0, 1.0)));

    return { phase_angle_deg, elongation_deg };
}

// feed apparent position to this function
std::pair<double, double> computePhaseAndElongation(Vec3 pa, SpiceDouble et, Observatory& obsv)
{
    SpiceDouble sv_sun[6], lt_ignore;
    spkezr_c("SUN", et, "J2000", "NONE", "SOLAR SYSTEM BARYCENTER", sv_sun, &lt_ignore);
    Vec3 p_sun = Vec3(sv_sun[0], sv_sun[1], sv_sun[2]);

    Vec3 p_obsv = getObserverPos(obsv, et);

    Vec3 obj_to_sun = p_sun - pa;
    Vec3 obj_to_obsv = p_obsv - pa;

    double cos_phase = obj_to_sun.normalized().dot(obj_to_obsv.normalized());
    double phase_angle_deg = rad2deg(acos(clamp(cos_phase, -1.0, 1.0)));

    Vec3 obs_to_sun = p_sun - p_obsv;
    Vec3 obsv_to_obj = pa - p_obsv;

    double cos_elong = obs_to_sun.normalized().dot(obsv_to_obj.normalized());
    double elongation_deg = rad2deg(acos(clamp(cos_elong, -1.0, 1.0)));

    return { phase_angle_deg, elongation_deg };
}

std::pair<std::vector<Vec3>, double> propagateToNextMidnight(Vec3 p0, Vec3 v0, SpiceDouble t0)
{
    double JD_0 = etToJD(t0);
    double JD_1 = std::floor(JD_0) + 0.5 + std::ceil(JD_0 - std::floor(JD_0 + 0.5)); // go to next mindight
    SpiceDouble t1 = JDToEt(JD_1);

    return { propagate(p0, v0, t0, t1, -1), JD_1 };
}

std::tuple<std::vector<std::string>, std::vector<std::string>, std::vector<std::string>,
    std::vector<double>, std::vector<double>, std::vector<double>, std::vector<std::string>, std::vector<std::string>> 
    generateEphemeris(SpiceDouble t0, Vec3 p0, Vec3 v0, std::vector<SpiceDouble> ts)
{
    std::vector<double> RA_result;
    std::vector<double> DEC_result;
    std::vector<double> range_result;
    std::vector<double> RA_rate_result;
    std::vector<double> DEC_rate_result;
    std::vector<double> phase_result;
    std::vector<double> elongation_result;

    double eps_rate = 60; // [s], time for numerical RA-DEC rate calculation

    for (int idx_tdiff = 0; idx_tdiff < ts.size(); idx_tdiff++)
    {
        if (ts[idx_tdiff] > t0)
        {
            std::vector<Vec3> sv_prop = propagate(p0, v0, t0, ts[idx_tdiff], -1);
            std::vector<Vec3> sv_apparent = computeApparentStateGeocentric(sv_prop[0], sv_prop[1], ts[idx_tdiff]);
            Vec3 p_t = sv_prop[0];
            Vec3 v_t = sv_prop[1];
            Vec3 p_a = sv_apparent[0];
            Vec3 v_a = sv_apparent[1];

            std::vector<double> RADEC = getGeocentricRADEC(ts[idx_tdiff], p_a);
            RA_result.push_back(RADEC[0]);
            DEC_result.push_back(RADEC[1]);

            // get range
            SpiceDouble state[6], lt;
            spkezr_c("EARTH", ts[idx_tdiff], "J2000", "NONE", "SOLAR SYSTEM BARYCENTER", state, &lt);
            Vec3 earth_pos = Vec3(state[0], state[1], state[2]);

            double range = (p_a - earth_pos).mag(); // [km]
            range_result.push_back(range);

            // get rates
            std::vector<Vec3> sv_prop2 = propagate(p_t, v_t, ts[idx_tdiff], ts[idx_tdiff] + eps_rate, -1);
            std::vector<Vec3> sv_apparent2 = computeApparentStateGeocentric(sv_prop2[0], sv_prop2[1], ts[idx_tdiff] + eps_rate);
            Vec3 p_a2 = sv_apparent2[0];

            std::vector<double> RADEC2 = getGeocentricRADEC(ts[idx_tdiff] + eps_rate, p_a2);
            RA_rate_result.push_back((RADEC2[0] - RADEC[0]) / eps_rate); // [deg s-1]
            DEC_rate_result.push_back((RADEC2[1] - RADEC[1]) / eps_rate); // [deg s-1]

            // get phase and elong.
            std::pair<double, double> phase_elong = computePhaseAndElongationGeocentric(p_a, ts[idx_tdiff]);
            phase_result.push_back(phase_elong.first);
            elongation_result.push_back(phase_elong.second);
        }
        else if (ts[idx_tdiff] == t0) // the date is epoch date
        {
            std::vector<Vec3> sv_apparent = computeApparentStateGeocentric(p0, v0, ts[idx_tdiff]);
            Vec3 p_a = sv_apparent[0];
            Vec3 v_a = sv_apparent[1];

            std::vector<double> RADEC = getGeocentricRADEC(ts[idx_tdiff], p_a);
            RA_result.push_back(RADEC[0]);
            DEC_result.push_back(RADEC[1]);

            // get range
            SpiceDouble state[6], lt;
            spkezr_c("EARTH", t0, "J2000", "NONE", "SOLAR SYSTEM BARYCENTER", state, &lt);
            Vec3 earth_pos = Vec3(state[0], state[1], state[2]);

            double range = (p_a - earth_pos).mag(); // [km]
            range_result.push_back(range);

            // get rates
            std::vector<Vec3> sv_prop2 = propagate(p0, v0, t0, t0 + eps_rate, -1);
            std::vector<Vec3> sv_apparent2 = computeApparentStateGeocentric(sv_prop2[0], sv_prop2[1], t0 + eps_rate);
            Vec3 p_a2 = sv_apparent2[0];

            std::vector<double> RADEC2 = getGeocentricRADEC(t0 + eps_rate, p_a2);
            RA_rate_result.push_back((RADEC2[0] - RADEC[0]) / eps_rate); // [deg s-1]
            DEC_rate_result.push_back((RADEC2[1] - RADEC[1]) / eps_rate); // [deg s-1]

            // get phase and elong.
            std::pair<double, double> phase_elong = computePhaseAndElongationGeocentric(p_a, ts[idx_tdiff]);
            phase_result.push_back(phase_elong.first);
            elongation_result.push_back(phase_elong.second);
        }
        else
        {
            std::vector<Vec3> sv_prop = backpropagate(p0, v0, t0, ts[idx_tdiff], 1);
            std::vector<Vec3> sv_apparent = computeApparentStateGeocentric(sv_prop[0], sv_prop[1], ts[idx_tdiff]);
            Vec3 p_t = sv_prop[0];
            Vec3 v_t = sv_prop[1];
            Vec3 p_a = sv_apparent[0];
            Vec3 v_a = sv_apparent[1];

            std::vector<double> RADEC = getGeocentricRADEC(ts[idx_tdiff], p_a);
            RA_result.push_back(RADEC[0]);
            DEC_result.push_back(RADEC[1]);

            // get range
            SpiceDouble state[6], lt;
            spkezr_c("EARTH", ts[idx_tdiff], "J2000", "NONE", "SOLAR SYSTEM BARYCENTER", state, &lt);
            Vec3 earth_pos = Vec3(state[0], state[1], state[2]);

            double range = (p_a - earth_pos).mag(); // [km]
            range_result.push_back(range);

            // get rates
            std::vector<Vec3> sv_prop2 = propagate(p_t, v_t, ts[idx_tdiff], ts[idx_tdiff] + eps_rate, -1);
            std::vector<Vec3> sv_apparent2 = computeApparentStateGeocentric(sv_prop2[0], sv_prop2[1], ts[idx_tdiff] + eps_rate);
            Vec3 p_a2 = sv_apparent2[0];

            std::vector<double> RADEC2 = getGeocentricRADEC(ts[idx_tdiff] + eps_rate, p_a2);
            RA_rate_result.push_back((RADEC2[0] - RADEC[0]) / eps_rate); // [deg s-1]
            DEC_rate_result.push_back((RADEC2[1] - RADEC[1]) / eps_rate); // [deg s-1]

            // get phase and elong.
            std::pair<double, double> phase_elong = computePhaseAndElongationGeocentric(p_a, ts[idx_tdiff]);
            phase_result.push_back(phase_elong.first);
            elongation_result.push_back(phase_elong.second);
        }
    }

    // convert RA and DEC from decimal degrees to hh:mm:ss and sDD:mm:ss
    std::vector<std::string> RA_hms;
    std::vector<std::string> DEC_Dms;
    std::vector<std::string> dates;
    std::vector<double> ranges;
    std::vector<double> RA_rate_arcsec_per_min;
    std::vector<double> DEC_rate_arcsec_per_min;
    std::vector<std::string> phase_p2;
    std::vector<std::string> elong_p2;

    for (int idx_result = 0; idx_result < RA_result.size(); idx_result++)
    {
        SpiceChar utcstr[20];
        et2utc_c(ts[idx_result], "ISOC", 0, 20, utcstr);
        dates.push_back(utcstr);
        RA_hms.push_back(RATohms(RA_result[idx_result]));
        DEC_Dms.push_back(DECToDms(DEC_result[idx_result]));
        ranges.push_back(range_result[idx_result] / AU);
        RA_rate_arcsec_per_min.push_back(RA_rate_result[idx_result] * 3600 * 60); // convert from [deg s-1] to [arcsec min-1]
        DEC_rate_arcsec_per_min.push_back(DEC_rate_result[idx_result] * 3600 * 60); // convert from [deg s-1] to [arcsec min-1]

        std::ostringstream oss;
        oss << std::fixed << std::setprecision(1) << phase_result[idx_result];
        if (phase_result[idx_result] < 100.0)
        {
            oss << " ";
        }
        if (phase_result[idx_result] < 10.0)
        {
            oss << " ";
        }
        std::string phase_str = oss.str();
        phase_p2.push_back(phase_str);

        std::ostringstream oss2;
        oss2 << std::fixed << std::setprecision(1) << elongation_result[idx_result];
        if (elongation_result[idx_result] < 100.0)
        {
            oss2 << " ";
        }
        if (elongation_result[idx_result] < 10.0)
        {
            oss2 << " ";
        }
        std::string elong_str = oss2.str();
        elong_p2.push_back(elong_str);
    }

    return std::tuple<std::vector<std::string>, std::vector<std::string>, std::vector<std::string>,
        std::vector<double>, std::vector<double>, std::vector<double>, std::vector<std::string>, std::vector<std::string>> 
        {dates, RA_hms, DEC_Dms, ranges, RA_rate_arcsec_per_min, DEC_rate_arcsec_per_min, elong_p2, phase_p2};
}

std::tuple<std::vector<std::string>, std::vector<std::string>, std::vector<std::string>, 
    std::vector<double>, std::vector<double>, std::vector<double>, std::vector<std::string>, std::vector<std::string>> 
    generateEphemerisObsv(SpiceDouble t0, Vec3 p0, Vec3 v0, std::vector<SpiceDouble> ts, std::string obscode, std::unordered_map<std::string, Observatory> obscode_map)
{
    std::vector<double> RA_result;
    std::vector<double> DEC_result;
    std::vector<double> range_result;
    std::vector<double> RA_rate_result;
    std::vector<double> DEC_rate_result;
    std::vector<double> phase_result;
    std::vector<double> elongation_result;

    double eps_rate = 60;

    for (int idx_tdiff = 0; idx_tdiff < ts.size(); idx_tdiff++)
    {
        if (ts[idx_tdiff] > t0)
        {
            std::vector<Vec3> sv_prop = propagate(p0, v0, t0, ts[idx_tdiff], -1);
            std::vector<Vec3> sv_apparent = computeApparentState(sv_prop[0], sv_prop[1], ts[idx_tdiff], obscode_map[obscode]);
            Vec3 p_t = sv_prop[0];
            Vec3 v_t = sv_prop[1];
            Vec3 p_a = sv_apparent[0];
            Vec3 v_a = sv_apparent[1];

            std::vector<double> RADEC = getRADEC(ts[idx_tdiff], obscode_map[obscode], p_a);
            RA_result.push_back(RADEC[0]);
            DEC_result.push_back(RADEC[1]);

            // get range
            double range = (p_a - getObserverPos(obscode_map[obscode], ts[idx_tdiff])).mag(); // [km]
            range_result.push_back(range);

            // get rates
            std::vector<Vec3> sv_prop2 = propagate(p_t, v_t, ts[idx_tdiff], ts[idx_tdiff] + eps_rate, -1);
            std::vector<Vec3> sv_apparent2 = computeApparentState(sv_prop2[0], sv_prop2[1], ts[idx_tdiff] + eps_rate, obscode_map[obscode]);
            Vec3 p_a2 = sv_apparent2[0];
            std::vector<double> RADEC2 = getRADEC(ts[idx_tdiff] + eps_rate, obscode_map[obscode], p_a2);
            RA_rate_result.push_back((RADEC2[0] - RADEC[0]) / eps_rate); // [deg s-1]
            DEC_rate_result.push_back((RADEC2[1] - RADEC[1]) / eps_rate); // [deg s-1]

            // get phase and elong.
            std::pair<double, double> phase_elong = computePhaseAndElongation(p_a, ts[idx_tdiff], obscode_map[obscode]);
            phase_result.push_back(phase_elong.first);
            elongation_result.push_back(phase_elong.second);
        }
        else if (ts[idx_tdiff] == t0) // the date is epoch date
        {
            std::vector<Vec3> sv_apparent = computeApparentState(p0, v0, t0, obscode_map[obscode]);
            Vec3 p_a = sv_apparent[0];
            Vec3 v_a = sv_apparent[1];
            std::vector<double> RADEC = getRADEC(ts[idx_tdiff], obscode_map[obscode], p_a);
            RA_result.push_back(RADEC[0]);
            DEC_result.push_back(RADEC[1]);

            // get range
            double range = (p_a - getObserverPos(obscode_map[obscode], t0)).mag(); // [km]
            range_result.push_back(range);

            // get rates
            std::vector<Vec3> sv_prop2 = propagate(p0, v0, t0, t0 + eps_rate, -1);
            std::vector<Vec3> sv_apparent2 = computeApparentState(sv_prop2[0], sv_prop2[1], t0 + eps_rate, obscode_map[obscode]);
            Vec3 p_a2 = sv_apparent2[0];
            std::vector<double> RADEC2 = getRADEC(t0, obscode_map[obscode], p_a2);
            RA_rate_result.push_back((RADEC2[0] - RADEC[0]) / eps_rate); // [deg s-1]
            DEC_rate_result.push_back((RADEC2[1] - RADEC[1]) / eps_rate); // [deg s-1]

            // get phase and elong.
            std::pair<double, double> phase_elong = computePhaseAndElongation(p_a, ts[idx_tdiff], obscode_map[obscode]);
            phase_result.push_back(phase_elong.first);
            elongation_result.push_back(phase_elong.second);
        }
        else
        {
            std::vector<Vec3> sv_prop = backpropagate(p0, v0, t0, ts[idx_tdiff], 1);
            std::vector<Vec3> sv_apparent = computeApparentState(sv_prop[0], sv_prop[1], ts[idx_tdiff], obscode_map[obscode]);
            Vec3 p_t = sv_prop[0];
            Vec3 v_t = sv_prop[0];
            Vec3 p_a = sv_apparent[0];
            Vec3 v_a = sv_apparent[1];

            std::vector<double> RADEC = getRADEC(ts[idx_tdiff], obscode_map[obscode], p_a);
            RA_result.push_back(RADEC[0]);
            DEC_result.push_back(RADEC[1]);

            // get range
            double range = (p_a - getObserverPos(obscode_map[obscode], ts[idx_tdiff])).mag(); // [km]
            range_result.push_back(range);

            // get rates
            std::vector<Vec3> sv_prop2 = propagate(p_t, v_t, ts[idx_tdiff], ts[idx_tdiff] + eps_rate, -1);
            std::vector<Vec3> sv_apparent2 = computeApparentState(sv_prop2[0], sv_prop2[1], ts[idx_tdiff] + eps_rate, obscode_map[obscode]);
            Vec3 p_a2 = sv_apparent2[0];
            std::vector<double> RADEC2 = getRADEC(ts[idx_tdiff] + eps_rate, obscode_map[obscode], p_a2);
            RA_rate_result.push_back((RADEC2[0] - RADEC[0]) / eps_rate); // [deg s-1]
            DEC_rate_result.push_back((RADEC2[1] - RADEC[1]) / eps_rate); // [deg s-1]

            // get phase and elong.
            std::pair<double, double> phase_elong = computePhaseAndElongation(p_a, ts[idx_tdiff], obscode_map[obscode]);
            phase_result.push_back(phase_elong.first);
            elongation_result.push_back(phase_elong.second);
        }
    }

    // convert RA and DEC from decimal degrees to hh:mm:ss and sDD:mm:ss
    std::vector<std::string> RA_hms;
    std::vector<std::string> DEC_Dms;
    std::vector<std::string> dates;
    std::vector<double> ranges;
    std::vector<double> RA_rate_arcsec_per_min;
    std::vector<double> DEC_rate_arcsec_per_min;
    std::vector<std::string> phase_p2;
    std::vector<std::string> elong_p2;

    for (int idx_result = 0; idx_result < RA_result.size(); idx_result++)
    {
        SpiceChar utcstr[20];
        et2utc_c(ts[idx_result], "ISOC", 0, 20, utcstr);
        dates.push_back(utcstr);
        RA_hms.push_back(RATohms(RA_result[idx_result]));
        DEC_Dms.push_back(DECToDms(DEC_result[idx_result]));
        ranges.push_back(range_result[idx_result] / AU);
        RA_rate_arcsec_per_min.push_back(RA_rate_result[idx_result] * 3600 * 60); // convert from [deg s-1] to [arcsec min-1]
        DEC_rate_arcsec_per_min.push_back(DEC_rate_result[idx_result] * 3600 * 60); // convert from [deg s-1] to [arcsec min-1]

        std::ostringstream oss;
        oss << std::fixed << std::setprecision(1) << phase_result[idx_result];
        if (phase_result[idx_result] < 100.0)
        {
            oss << " ";
        }
        if (phase_result[idx_result] < 10.0)
        {
            oss << " ";
        }
        std::string phase_str = oss.str();
        phase_p2.push_back(phase_str);

        std::ostringstream oss2;
        oss2 << std::fixed << std::setprecision(1) << elongation_result[idx_result];
        if (elongation_result[idx_result] < 100.0)
        {
            oss2 << " ";
        }
        if (elongation_result[idx_result] < 10.0)
        {
            oss2 << " ";
        }
        std::string elong_str = oss2.str();
        elong_p2.push_back(elong_str);
    }

    return std::tuple<std::vector<std::string>, std::vector<std::string>, std::vector<std::string>,
        std::vector<double>, std::vector<double>, std::vector<double>, std::vector<std::string>, std::vector<std::string>> 
        {dates, RA_hms, DEC_Dms, ranges, RA_rate_arcsec_per_min, DEC_rate_arcsec_per_min, elong_p2, phase_p2};
}

/*
std::tuple<std::vector<std::string>, std::vector<std::string>, std::vector<std::string>,
    std::vector<double>, std::vector<double>, std::vector<double>> generateEphemerisMPEC(Vec3 p0, Vec3 v0, SpiceDouble epoch_et)
{
    // get today's date with clock at midnight
    std::time_t now = std::time(nullptr);
    std::tm utc_tm{};
    errno_t err = gmtime_s(&utc_tm, &now);
    if (err) {
        throw std::runtime_error("gmtime_s failed");
    }

    utc_tm.tm_hour = 0;
    utc_tm.tm_min = 0;
    utc_tm.tm_sec = 0;

    std::time_t utc_midnight = _mkgmtime(&utc_tm);

    char utc_str[32];
    std::strftime(utc_str, sizeof(utc_str), "%Y-%m-%dT%H:%M:%S", &utc_tm);

    SpiceDouble et_today;
    utc2et_c(utc_str, &et_today);

    std::vector<SpiceDouble> target_ets;
    target_ets.push_back(et_today - 30 * 86400.0);
    target_ets.push_back(et_today - 15 * 86400.0);
    target_ets.push_back(et_today - 7 * 86400.0);
    target_ets.push_back(et_today - 1 * 86400.0);
    target_ets.push_back(et_today);
    target_ets.push_back(et_today + 1 * 86400.0);
    target_ets.push_back(et_today + 7 * 86400.0);
    target_ets.push_back(et_today + 15 * 86400.0);
    target_ets.push_back(et_today + 30 * 86400.0);

    return generateEphemeris(epoch_et, p0, v0, target_ets);
}
*/

void printEphemeris(std::string designation, std::string obscode, std::string obscode_name, 
    std::tuple<std::vector<std::string>, std::vector<std::string>, std::vector<std::string>,
    std::vector<double>, std::vector<double>, std::vector<double>, std::vector<std::string>, 
    std::vector<std::string>> ephem_data,
    std::ostream& out = std::cout)
{
    std::string date_init = std::get<0>(ephem_data)[0];
    std::string date_final = std::get<0>(ephem_data).back();
    out << "\nEphemeris for " << designation << " between " << date_init << " and " << date_final
        << " for observatory " << obscode << " (" << obscode_name << "):\n\n";
    out << "Date-Time                R.A. (J2000) Decl.       Range (AU)   RA rate (\"/min)  DEC rate (\"/min)  Elong. (deg) Phase (deg.)\n";

    for (int idx_ephem = 0; idx_ephem < std::get<0>(ephem_data).size(); idx_ephem++)
    {
        out << std::get<0>(ephem_data)[idx_ephem] << "   " << std::get<1>(ephem_data)[idx_ephem] << "   " << std::get<2>(ephem_data)[idx_ephem] << "     " 
            << std::fixed << std::setprecision(6) << std::get<3>(ephem_data)[idx_ephem] << "        "
            << std::fixed << std::setprecision(6) << std::get<4>(ephem_data)[idx_ephem] << "       "
            << std::fixed << std::setprecision(6) << std::get<5>(ephem_data)[idx_ephem] << "        " 
            << std::get<6>(ephem_data)[idx_ephem] << "        " << std::get<7>(ephem_data)[idx_ephem] << "\n";
    }

    out << "\nEnd of ephemeris output.\nMPEP v0.3.0\n";
}

void printHelpMsg()
{
    std::cout << " === MPEP HELP ===\n\n";
    std::cout << "MPEP is an ephemeris generation software for minor planets in heliocentric orbit. It computes the RA-DEC sky positions of objects "
        << "between given dates at given intervals, as well as the elongation and phase of the object.\n\n";

    std::cout << "It requires a minor planet data file in JSON format, some SPICE kernels "
        << "and an ObsCodes file which holds the observatory codes used by the MPC.\n\n";

    std::cout << "The minimal invocation is merely 'mpep' - this assumes defaults for all parameters, which are:\n\n";
    std::cout << "    Input data file: mp.json\n";
    std::cout << "    Ephemeris initial date: Today at 00:00 UTC\n";
    std::cout << "    Ephemeris final date: Tomorrow at 00:00 UTC\n";
    std::cout << "    Ephemeris step size: 1 hour\n";
    std::cout << "    Ephemeris observatory: Geocentric\n";
    std::cout << "    Date input file: None (Ignored)\n";
    std::cout << "    SPICE Folder: data/SPICE/\n";
    std::cout << "    MPC Observatory Codes file: data/ObsCodes.dat\n";
    std::cout << "    Output filename: ephemeris.txt\n\n";

    std::cout << "The parameters can be adjusted using the following arguments:\n";
    std::cout << "    -mp: Minor planet data JSON file\n";
    std::cout << "    -init_date: Ephemeris initial date (YYYY-MM-DDThh:mm:ss)\n";
    std::cout << "    -final_date: Ephemeris final date (YYYY-MM-DDThh:mm:ss)\n";
    std::cout << "    -step_day: Step size in days (cumulative with other step size parameters)\n";
    std::cout << "    -step_hour: Step size in hours (cumulative with other step size parameters)\n";
    std::cout << "    -step_minute: Step size in minutes (cumulative with other step size parameters)\n";
    std::cout << "    -step_second: Step size in seconds (cumulative with other step size parameters)\n";
    std::cout << "    -obscode: MPC observatory codes data filepath\n";
    std::cout << "    -spice: SPICE kernels directory path\n";
    std::cout << "    -inp: ephemeris date input filepath**\n";
    std::cout << "    -out: ephemeris output filepath\n\n";

    std::cout << "    ** If you use the '-inp' argument; initial date, final date and step size parameters will be ignored!\n\n";

    std::cout << "Observatory code file contents can be found in MPC website --> https://www.minorplanetcenter.net/iau/lists/ObsCodes.html\n";
    std::cout << "You can just copy the entire page into a text file.\n\n";

    std::cout << "SPICE kernels can be obtained from NAIF --> https://naif.jpl.nasa.gov/pub/naif/generic_kernels/\n";
    std::cout << "You can get whichever ones are most suitable, but please ensure you have at least the planet ephemerides kernel, a leapseconds kernel, "
        << "a generic text PCK and a binary Earth PCK.\n If you do not know what to get, you can get the following:\n";
    std::cout << "    de440.bsp (https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de440.bsp)\n";
    std::cout << "    naif0012.tls (https://naif.jpl.nasa.gov/pub/naif/generic_kernels/lsk/naif0012.tls)\n";
    std::cout << "    pck00011.tpc (https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/pck00011.tpc)\n";
    std::cout << "    earth_1962_240827_2124_combined.bpc (https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/earth_1962_240827_2124_combined.bpc)\n";
    std::cout << "        (the last one might get updated often, just pick one that looks like it)\n\n";

    std::cout << "MPEP was developed by H. A. Guler.\n";
    std::cout << "MPEP is licensed under GNU General Public License version 2.0 (GPL-2.0 License)\n\n";
}

int main(int argc, char* argv[])
{
    std::cout << "MPEP v0.3.0\n\n";

    // default parameters
    std::string mp_path = "mp.json";
    std::string date_init = "None";
    std::string date_final = "None";
    double step_size = 0;
    std::string user_obscode = "Geocentric";
    std::string in_path = "None";
    std::string obscode_path = "data/ObsCodes.dat";
    std::string spice_path = "data/SPICE/";
    int param_maxiter = 20;
    std::string out_path = "ephemeris.txt";

    // handle command line arguments
    // there is a more compact version of doing this but this is easier for my brain
    int argtype = 0;
    for (int idx_cmd = 1; idx_cmd < argc; idx_cmd++)
    {
        if (!strcmp(argv[idx_cmd], "-mp")) // also the default argument
        {
            argtype = 0;
        }
        else if (!strcmp(argv[idx_cmd], "-init_date"))
        {
            argtype = 1;
        }
        else if (!strcmp(argv[idx_cmd], "-final_date"))
        {
            argtype = 2;
        }
        else if (!strcmp(argv[idx_cmd], "-step_day"))
        {
            argtype = 3;
        }
        else if (!strcmp(argv[idx_cmd], "-step_hour"))
        {
            argtype = 4;
        }
        else if (!strcmp(argv[idx_cmd], "-step_minute"))
        {
            argtype = 5;
        }
        else if (!strcmp(argv[idx_cmd], "-step_second"))
        {
            argtype = 6;
        }
        else if (!strcmp(argv[idx_cmd], "-obscode"))
        {
            argtype = 7;
        }
        else if (!strcmp(argv[idx_cmd], "-spice"))
        {
            argtype = 8;
        }
        else if (!strcmp(argv[idx_cmd], "-inp")) // forget all arguments, read from input file instead
        {
            argtype = 9;
        }
        else if (!strcmp(argv[idx_cmd], "-out") || !strcmp(argv[idx_cmd], "-output"))
        {
            argtype = 10;
        }
        else if (!strcmp(argv[idx_cmd], "-h") || !strcmp(argv[idx_cmd], "--help")) // two dashes because people are more used to it
        {
            printHelpMsg();
            return 0;
        }
        else
        {
            switch (argtype)
            {
            case 0:
                mp_path = argv[idx_cmd];
                break;
            case 1:
                date_init = argv[idx_cmd];
                break;
            case 2:
                date_final = argv[idx_cmd];
                break;
            case 3:
                step_size += strtod(argv[idx_cmd], NULL) * 86400.0; // day
                break;
            case 4:
                step_size += strtod(argv[idx_cmd], NULL) * 60 * 60; // hour
                break;
            case 5:
                step_size += strtod(argv[idx_cmd], NULL) * 60; // minute
                break;
            case 6:
                step_size += strtod(argv[idx_cmd], NULL); // second
                break;
            case 7:
                user_obscode = argv[idx_cmd];
                break;
            case 8:
                spice_path = argv[idx_cmd];
                break;
            case 9:
                in_path = argv[idx_cmd];
                break;
            case 10:
                out_path = argv[idx_cmd];
            }
        }
    }

    std::cout << "Loading SPICE kernels... ";
    loadAllKernels(spice_path);
    std::cout << "Done\n";

    // read and parse MPC ObsCodes
    std::vector<Observatory> obscodes = readObsCodes(obscode_path);
    std::unordered_map<std::string, Observatory> obscode_map;

    // mapify obscodes by... well... code
    for (int i = 0; i < obscodes.size(); i++)
    {
        obscode_map.insert({ obscodes[i].code, obscodes[i] });
    }

    std::tuple<std::string, std::string, std::vector<double>, Vec3, Vec3, double> MPdata = readMinorPlanetFile(mp_path);
    // returns {perm, prov, orbital_elements, p0, v0, epoch_JD}

    std::string perm = std::get<0>(MPdata);
    std::string prov = std::get<1>(MPdata);
    std::vector<double> orbital_elements = std::get<2>(MPdata);
    Vec3 p0 = std::get<3>(MPdata);
    Vec3 v0 = std::get<4>(MPdata);
    double epoch_JD = std::get<5>(MPdata);

    SpiceDouble t0 = JDToEt(epoch_JD);
    std::vector<SpiceDouble> ts;

    // construct ephemeris times
    if (!strcmp(in_path.c_str(), "None"))
    {
        if (!strcmp(date_init.c_str(), "None"))
        {
            // no date given, generate ephemeris for today
            std::cout << "No ephemeris date requested, generating ephemeris for current date...\n";

            // get today's date with clock at midnight
            std::time_t now = std::time(nullptr);
            std::tm utc_tm{};
            errno_t err = gmtime_s(&utc_tm, &now);
            if (err) {
                throw std::runtime_error("gmtime_s failed");
            }

            utc_tm.tm_hour = 0;
            utc_tm.tm_min = 0;
            utc_tm.tm_sec = 0;

            std::time_t utc_midnight = _mkgmtime(&utc_tm);

            char utc_str[32];
            std::strftime(utc_str, sizeof(utc_str), "%Y-%m-%dT%H:%M:%S", &utc_tm);

            SpiceDouble et_today;
            utc2et_c(utc_str, &et_today);

            for (int idx_hour = 0; idx_hour < 25; idx_hour++)
            {
                ts.push_back(et_today + idx_hour * 60 * 60);
            }
        }
        else
        {
            // we are reading datetime parameters from arguments
            SpiceDouble et0, etf;
            str2et_c(date_init.c_str(), &et0);
            str2et_c(date_final.c_str(), &etf);

            SpiceDouble c_et = et0;
            while (c_et < etf)
            {
                ts.push_back(c_et);
                c_et += step_size;
            }
        }
    }
    else
    {
        // a file was provided, ignore arguments
        ts = readDateFile(in_path);
    }

    std::tuple<std::vector<std::string>, std::vector<std::string>, std::vector<std::string>, std::vector<double>, 
        std::vector<double>, std::vector<double>, std::vector<std::string>, std::vector<std::string>> ephem_data;
    if (!strcmp(user_obscode.c_str(), "Geocentric"))
    {
        ephem_data = generateEphemeris(t0, p0, v0, ts);
        // returns {date, RA, DEC}
    }
    else
    {
        ephem_data = generateEphemerisObsv(t0, p0, v0, ts, user_obscode, obscode_map);
        // returns {date, RA, DEC}
    }

    std::string desig = perm;
    if (!strcmp(perm.c_str(), ""))
    {
        desig = prov;
    }

    printEphemeris(desig, user_obscode, obscode_map[user_obscode].name, ephem_data);

    // write to file
    std::cout << "\nWriting output ephemeris to " << out_path << "... ";
    std::ofstream outfile(out_path);
    if (outfile.is_open())
    {
        printEphemeris(desig, user_obscode, obscode_map[user_obscode].name, ephem_data, outfile);
        outfile.close();
    }

    std::cout << "Done. Program end.\n";
}
