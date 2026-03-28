/*
 * SISI Simulator - Signal Interference from Surface Interactions
 * Copyright (c) 2023-2026, Finn Linzer, TU Wien.
 * * This file is part of the SISI Simulator, licensed under GPLv3.
 * For full license, scientific disclaimer, and third-party credits (PBRT),
 * please refer to the sisi.cpp or the LICENSE file in the project root.
 *
 */

#ifndef SISI_HPP
#define SISI_HPP

// Math-Defines for Windows (M_PI)
#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif

#include <algorithm>
#include <cmath>
#include <complex>
#include <random>
#include <string>
#include <vector>
#include <cctype>


// Add the following line to create your own distribution
// (my_custom_distribution(..)) -> recompile
// #define USE_OWN_DISTRIBUTION

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Enum, active mode
enum class SimulationMode { LASERSCANNER, PATTERSON, MATERIAL_IMPULSE };

// Parameters applicable to all modes
struct SimulationConfig {
  bool use_pseudo_random_seed;
  uint32_t n_max_detections;
};

struct LightSourceConfig {
  double wavelength_nm;
  uint32_t n_photons;
};

struct MaterialConfig {
  double n_refractive_index;
  double absorption_coefficient_per_mm;
  double scattering_coefficient_per_mm;
  double anisotropy_g;
  double specimen_thickness_m;
  double microfacet_std_deviation_deg;
};

// Modi specific structs
struct OpticConfig {
  double radius_m;
  double working_distance_m;
  double pinhole_sigma_factor;
  double rho_m;
};

struct BeamConfig {
  double sigma_x_m;
  double sigma_y_m;
  double sigma_z_m;
};

struct PhysicsOverrideConfig {
  bool force_100_percent_transmission;
};

struct ScanConfig {
  std::vector<double> angles_deg;
};

struct DiagnosticsConfig {
  bool track_photon_paths;
  bool print_results;
};

struct OutputConfig {
  bool save_results;
  bool save_path_lengths;
  bool save_registered_spot;
  bool save_hemisphere;
  bool save_photon_paths;
  std::string result_path;
  double hemisphere_radius_m;
};

// The same for all modes
struct ModeConfig {
  OpticConfig optic;
  BeamConfig beam;
  PhysicsOverrideConfig physics_override;
  ScanConfig scan;
  DiagnosticsConfig diagnostics;
  OutputConfig output;
};

// The main struct that contains all
struct SisiConfig {
  SimulationMode mode;

  SimulationConfig simulation;
  LightSourceConfig light_source;
  MaterialConfig material;

  // Which mode was chosen?
  ModeConfig active_mode;
};

// YAML END

struct Position3D {
  double x, y, z;
}; // A single photons position (no direction)

struct Direction3D {
  double x, y, z;
};

struct Direction4D {
  double x, y, z;
  int i; // The “i” stands for information
};

struct Position4D {
  double x, y, z;
  int i;
}; // i marks whether the photodiode has been hit

inline Position3D subtract(Position3D a, Position3D b) {
  Position3D result = {a.x - b.x, a.y - b.y, a.z - b.z};
  return result;
}

inline Position3D add(Position3D a, Direction3D b) {
  Position3D result = {a.x + b.x, a.y + b.y, a.z + b.z};
  return result;
}

inline Direction3D scalarMultiply(Direction3D v, double s) {
  Direction3D result = {v.x * s, v.y * s, v.z * s};
  return result;
}

inline double dot(const Direction3D &a, const Direction3D &b) {
  return a.x * b.x + a.y * b.y + a.z * b.z;
}

inline double norm(Direction3D v) { return sqrt(dot(v, v)); }

inline Direction3D normalize(Direction3D v) {
  double n = norm(v);
  Direction3D norm_vec;
  norm_vec.x = v.x / n;
  norm_vec.y = v.y / n;
  norm_vec.z = v.z / n;
  return norm_vec;
}

inline Direction3D crossProduct(Direction3D a, Direction3D b) {
  Direction3D res;
  res.x = a.y * b.z - a.z * b.y;
  res.y = a.z * b.x - a.x * b.z;
  res.z = a.x * b.y - a.y * b.x;
  return res;
} // Cross product of two 3D vectors

////// GRANITE 0 Density function
// Generate Normal Distribution (Gaussian) Density function
#include <iostream>
#include <numeric>
#include <random>
// Struct for storing the parameters of each normal distribution
struct NormalDist {
  double mu;    // mean
  double sigma; // std.dev.
};
// Function for extracting from own distribution
double sampleFromMixture(const std::vector<NormalDist> &distributions,
                         const std::vector<double> &weights,
                         std::mt19937 &rng) {

  // Create a distribution based on the weights
  std::discrete_distribution<> d(weights.begin(), weights.end());

  // Select one of the normal distributions based on the weights
  int index = d(rng);

  // Draw a sample from this normal distribution
  std::normal_distribution<> distribution(distributions[index].mu,
                                          distributions[index].sigma);
  return distribution(rng);
}
double my_custom_distribution(double mu, double sigma, std::mt19937 &rng) {
  // doing nothing with mu and sigma in this very special case of granit0
  // Density Function
  // To generate these values, use: sisi_nur1_oberfl.m
  static const std::vector<NormalDist> distributions = {
      {-58.208641, 3.0}, {-56.254428, 3.0}, {-54.300215, 3.0},
      {-52.346001, 3.0}, {-50.391788, 3.0}, {-48.437574, 3.0},
      {-46.483361, 3.0}, {-44.529148, 3.0}, {-42.574934, 3.0},
      {-40.620721, 3.0}, {-38.666507, 3.0}, {-36.712294, 3.0},
      {-34.758081, 3.0}, {-32.803867, 3.0}, {-30.849654, 3.0},
      {-28.895440, 3.0}, {-26.941227, 3.0}, {-24.987013, 3.0},
      {-23.032800, 3.0}, {-21.078587, 3.0}, {-19.124373, 3.0},
      {-17.170160, 3.0}, {-15.215946, 3.0}, {-13.261733, 3.0},
      {-11.307520, 3.0}, {-9.353306, 3.0},  {-7.399093, 3.0},
      {-5.444879, 3.0},  {-3.490666, 3.0},  {-1.536453, 3.0},
      {0.417761, 3.0},   {2.371974, 3.0},   {4.326188, 3.0},
      {6.280401, 3.0},   {8.234614, 3.0},   {10.188828, 3.0},
      {12.143041, 3.0},  {14.097255, 3.0},  {16.051468, 3.0},
      {18.005682, 3.0},  {19.959895, 3.0},  {21.914108, 3.0},
      {23.868322, 3.0},  {25.822535, 3.0},  {27.776749, 3.0},
      {29.730962, 3.0},  {31.685175, 3.0},  {33.639389, 3.0},
      {35.593602, 3.0},  {37.547816, 3.0},  {39.502029, 3.0},
      {41.456242, 3.0},  {43.410456, 3.0},  {45.364669, 3.0},
      {47.318883, 3.0},  {49.273096, 3.0},  {51.227309, 3.0},
      {53.181523, 3.0},  {55.135736, 3.0},  {57.089950, 3.0},
  };

  static const std::vector<double> weights = {
      0.000033, 0.000133, 0.001033, 0.001467, 0.001633, 0.002700, 0.002600,
      0.002567, 0.003067, 0.003233, 0.004233, 0.004100, 0.004900, 0.005267,
      0.008367, 0.007634, 0.008067, 0.013534, 0.014300, 0.014034, 0.018867,
      0.017734, 0.023834, 0.028268, 0.033901, 0.023601, 0.040235, 0.055669,
      0.065302, 0.070769, 0.076403, 0.066569, 0.049235, 0.037501, 0.042401,
      0.037301, 0.030601, 0.025868, 0.024634, 0.017501, 0.014434, 0.012934,
      0.010967, 0.012534, 0.009700, 0.006600, 0.007867, 0.005667, 0.006367,
      0.003233, 0.004733, 0.003300, 0.002200, 0.002433, 0.002967, 0.001500,
      0.002033, 0.000933, 0.000400, 0.000100,
  };

  return ((sampleFromMixture(distributions, weights, rng) / 200.) * M_PI);
}
//////////////////////////////////////////////// GRANIT 0 Density Function
/// END! //////////////

inline double sampleTilt(double mu, double sigma, std::mt19937 &rng) {
#ifdef USE_OWN_DISTRIBUTION
  // Use a custom distribution
  return my_custom_distribution(mu, sigma, rng);
#else
  // C++ Normal-Distribution
  static std::normal_distribution<double> normalDist(0.0, 1.0);
  return mu + sigma * normalDist(rng);
#endif
}

inline Direction3D positional2directional(Position3D a, bool b) {
  // If bool b is false (0):  Calculate vector direction from origin (0,0,0) to
  // Position3D If bool b is true  (1):  Vice Versa

  Direction3D result = {a.x, a.y, a.z};
  if (b)
    result = {-a.x, -a.y, -a.z};

  return result;
}

Position3D determinePhotonEmissionLocation(Position3D opticalCenter,
                                           double sigma_x_m, double sigma_y_m,
                                           double sigma_z_m, Direction3D basisX,
                                           Direction3D basisY,
                                           Direction3D basisZ,
                                           std::mt19937 &rng) {
  // You may test other distributions for research purpose
  // static std::uniform_real_distribution<double> uniform01(0.0, 1.0);
  static std::normal_distribution<double> normalDist(0.0, 1.0);
  Position3D variation;

  // For X and Y we use different standard deviations
  variation.x = sigma_x_m * normalDist(rng);
  variation.y = sigma_y_m * normalDist(rng);
  // Assuming the radiation-time of the photons is not perfect (dirac), this can
  // be represented by variation.z Assuming a pulse shape (Gaussian) along the
  // launch axis Z
  variation.z = sigma_z_m * normalDist(rng);

  Position3D positionOfEmission;
  positionOfEmission.x = opticalCenter.x + variation.x * basisX.x +
                         variation.y * basisY.x + variation.z * basisZ.x;
  positionOfEmission.y = opticalCenter.y + variation.x * basisX.y +
                         variation.y * basisY.y + variation.z * basisZ.y;
  positionOfEmission.z = opticalCenter.z + variation.x * basisX.z +
                         variation.y * basisY.z + variation.z * basisZ.z;

  return positionOfEmission;
}

Position3D determinePhotonHitSurfaceLocation(Position3D pointOfEmission,
                                             Direction3D normal) {
  Position3D hitOnTarget;

  double scale = -pointOfEmission.z / normal.z;
  hitOnTarget.x = pointOfEmission.x + scale * normal.x;
  hitOnTarget.y = pointOfEmission.y + scale * normal.y;
  hitOnTarget.z = 0; // Target always flat (currently at least)

  return hitOnTarget;
}

// Hit Photodiode?
Position4D backscatteringMayHitOptic(Position3D shooter_position,
                                     Direction3D shooter_direction,
                                     Position3D optic_center,
                                     double optic_radius, SimulationMode mode) {
  // Calculate the normal of the tilted optic to the plane
  Direction3D normal;
  if (mode == SimulationMode::PATTERSON) {
    normal = {0.0, 0.0, 1.0};
  } else {
    normal = normalize({-optic_center.x, -optic_center.y, -optic_center.z});
  }

  // Find intersection point using dot product with plane's normal
  Direction3D direction_photon2optics =
      positional2directional(subtract(optic_center, shooter_position), 0);

  double t =
      dot(direction_photon2optics, normal) / dot(shooter_direction, normal);
  double distance_to_sensor = sqrt(
      pow(optic_center.x, 2) + pow(optic_center.y, 2) + pow(optic_center.z, 2));

  if (dot(shooter_direction, normal) > 0.0 || t > 3 * distance_to_sensor) {
    // Photon trajectory points away from the target plane; intersection is
    // impossible. This can occur if the surface/microfacet is extremely tilted
    // relative to the photon's direction. Discarding ray.
    return {0, 0, 0, -1};
  }

  Position3D intersection_point =
      add(shooter_position, scalarMultiply(shooter_direction, t));

  // Check if the intersection point is within the optic
  double distance_to_center =
      sqrt(pow(intersection_point.x - optic_center.x, 2) +
           pow(intersection_point.y - optic_center.y, 2) +
           pow(intersection_point.z - optic_center.z, 2));
  return {intersection_point.x, intersection_point.y, intersection_point.z,
          (distance_to_center <= optic_radius)};
}

// Fresnel reflection for two complex refractive indices
std::pair<double, double> fresnel_reflection(double theta_i_rad, double n1,
                                             double k1, // incident medium
                                             double n2,
                                             double k2 // transmission medium
) {
  using cd = std::complex<double>;

  // Complex refractive indices (n - i k)
  const cd n1c(n1, -k1);
  const cd n2c(n2, -k2);

  const double cos_theta_i = std::cos(theta_i_rad);
  const double sin_theta_i = std::sin(theta_i_rad);

  // Snell's law (complex)
  const cd sin_theta_t = (n1c / n2c) * sin_theta_i;
  const cd cos_theta_t = std::sqrt(cd(1.0, 0.0) - sin_theta_t * sin_theta_t);

  // Fresnel coefficients
  const cd rs = (n1c * cos_theta_i - n2c * cos_theta_t) /
                (n1c * cos_theta_i + n2c * cos_theta_t);

  const cd rp = (n2c * cos_theta_i - n1c * cos_theta_t) /
                (n2c * cos_theta_i + n1c * cos_theta_t);

  // Reflectance (energy)
  double R_s = std::norm(rs);
  double R_p = std::norm(rp);

  // Clamp for numerical safety
  R_s = std::min(1.0, std::max(0.0, R_s));
  R_p = std::min(1.0, std::max(0.0, R_p));

  return std::make_pair(R_s, R_p);
}

double incident_angle(const Direction3D &incident, const Direction3D &normal) {
  Direction3D incident_norm = normalize(incident);
  Direction3D normal_norm = normalize(normal);
  double cos_theta = std::acos(dot(incident_norm, normal_norm));
  return cos_theta;
}

double normalize_angle(double angle) {
  while (angle > M_PI)
    angle -= 2 * M_PI;
  while (angle < -M_PI)
    angle += 2 * M_PI;
  return angle;
}

std::pair<double, double> relative_tilt(const Direction3D &v1,
                                        const Direction3D &v2) {
  double tilt_x = atan2(v2.y, v2.z) - atan2(v1.y, v1.z);
  double tilt_y = atan2(v2.x, v2.z) - atan2(v1.x, v1.z);
  tilt_x = normalize_angle(tilt_x);
  tilt_y = normalize_angle(tilt_y);
  return {tilt_x, tilt_y};
}

Direction3D rnd_microfacet(double std_deviation_rad, double direction,
                           std::mt19937 &rng) {

  double nx = 0, ny = 0, nz = 0;
  double sigma_rad = std_deviation_rad;

  do {
    double tiltX = sampleTilt(0.0, sigma_rad, rng); // Tilt X-Direction [rad]
    double tiltY = sampleTilt(0.0, sigma_rad, rng); // Tilt Y-Direction [rad]

    // Normal-Vector of tilted Surface
    nx = sin(tiltX);
    ny = sin(tiltY);
    nz = sqrt(1.0 - nx * nx - ny * ny); // |n| = 1
  } while ((nx * nx + ny * ny) > 1);

  double length = sqrt(nx * nx + ny * ny + nz * nz);
  nx /= length;
  ny /= length;
  nz /= length;

  // If photon direction is negative, its about to leave the specimen
  // Microfacet normal should point into the same direction
  if (direction < 0) {
    nz = -nz;
  } else {
    //    std::cout << nx << " " << ny << " " << nz << std::endl;
  }
  return {nx, ny, nz};
}

// Compute cos θ for Henyey–Greenstein sample
// Book: "Physically Based Rendering" - Third Edition, Pharr et al, Page 899
inline double computeCosThetaHG(const double &g, std::mt19937 &rng) {
  static std::uniform_real_distribution<double> uniform01(0.0, 1.0);
  const double u0 = uniform01(rng); // corresponds to u[0] in PBRT

  if (std::abs(g) < 1e-3)
    return 1.0 - 2.0 * u0;
  else {
    const double sqrTerm = (1.0 - g * g) / (1.0 - g + 2.0 * g * u0);
    return (1.0 + g * g - sqrTerm * sqrTerm) / (2.0 * g);
  }
}

// Coordinate System from a Vector, Chapter 2.2.4.
// Book: "Physically Based Rendering" - Third Edition, Pharr et al, Page 67
inline void coordinateSystem(const Direction3D &v1, Direction3D &v2,
                             Direction3D &v3) {
  if (std::abs(v1.x) > std::abs(v1.y)) {
    const double invLen = 1.0 / std::sqrt(v1.x * v1.x + v1.z * v1.z);

    v2 = {-v1.z * invLen, 0.0, v1.x * invLen};
  } else {
    const double invLen = 1.0 / std::sqrt(v1.y * v1.y + v1.z * v1.z);

    v2 = {0.0, v1.z * invLen, -v1.y * invLen};
  }

  v3 = crossProduct(v1, v2);
}

// Takes three basis vectors representing the x, y, and z axes of the
// coordinate system and returns the appropriate direction vector
// with respect to the coordinate frame defined by them
// Book: "Physically Based Rendering" - Third Edition, Pharr et al, Page 346
inline Direction3D sphericalDirection(const double &sinTheta,
                                      const double &cosTheta, const double &phi,
                                      const Direction3D &x,
                                      const Direction3D &y,
                                      const Direction3D &z) {
  return {sinTheta * std::cos(phi) * x.x + sinTheta * std::sin(phi) * y.x +
              cosTheta * z.x,
          sinTheta * std::cos(phi) * x.y + sinTheta * std::sin(phi) * y.y +
              cosTheta * z.y,
          sinTheta * std::cos(phi) * x.z + sinTheta * std::sin(phi) * y.z +
              cosTheta * z.z};
}

// The Refract-Function at MaterialTransition
Direction3D refract(const Direction3D &incident, const Direction3D &normal,
                    double n1, double n2) {
  Direction3D incident_norm = normalize(incident);
  Direction3D normal_norm = normalize(normal);
  double cos_theta_i = dot(incident_norm, normal_norm);

  // Photon direction shows outside top of material
  if (cos_theta_i < 0) {
    cos_theta_i = -cos_theta_i;
    std::swap(n1, n2);
    normal_norm.x = -normal_norm.x;
    normal_norm.y = -normal_norm.y;
    normal_norm.z = -normal_norm.z;
  }

  double sin_theta_i = std::sqrt(1.0 - cos_theta_i * cos_theta_i);
  double sin_theta_r = (n1 / n2) * sin_theta_i;

  // Total reflection
  if (sin_theta_r >= 1.0) {
    // std::cerr << "Total internal reflection occurs, no refraction" <<
    // std::endl;
    return {-9, -9, -9};
  }

  double cos_theta_r = std::sqrt(1.0 - sin_theta_r * sin_theta_r);
  Direction3D refracted;
  refracted.x = incident_norm.x * (n1 / n2) -
                normal_norm.x * ((n1 / n2) * cos_theta_i - cos_theta_r);
  refracted.y = incident_norm.y * (n1 / n2) -
                normal_norm.y * ((n1 / n2) * cos_theta_i - cos_theta_r);
  refracted.z = incident_norm.z * (n1 / n2) -
                normal_norm.z * ((n1 / n2) * cos_theta_i - cos_theta_r);

  //  Normalizing the refracted vector
  return normalize(refracted);
}

Direction3D reflect(const Direction3D &incident,
                    const Direction3D &microfacetNormal) {
  // Calculate the reflected vector
  double dot =
      2.0 * (incident.x * microfacetNormal.x + incident.y * microfacetNormal.y +
             incident.z * microfacetNormal.z);

  Direction3D reflected;
  reflected.x = incident.x - dot * microfacetNormal.x;
  reflected.y = incident.y - dot * microfacetNormal.y;
  reflected.z = incident.z - dot * microfacetNormal.z;

  return reflected;
}

Direction4D
materialTransition(const Direction3D photon_direction, const double &n1,
                   const double &n_refractive_index, const double &k,
                   const Direction3D &microfacetNormal,
                   const Direction3D &incidentPolarization,
                   bool force_100_percent_transmission, std::mt19937 &rng) {

  // Uniform Distribution between 0 and 1
  static std::uniform_real_distribution<double> uniform01(0.0, 1.0);

  //// Relative Tilt
  // if (photon_direction.z < 0) {
  // std::cout << microfacetNormal.x << " " << microfacetNormal.y << " " <<
  // microfacetNormal.z << " " << photon_direction.x << " " <<
  // photon_direction.y
  // << " " << photon_direction.z std::endl;
  //}
  std::pair<double, double> tilt =
      relative_tilt(microfacetNormal, photon_direction);
  double relative_tilt_x = tilt.first;
  double relative_tilt_y = tilt.second;

  // if (photon_direction.z < 0) {
  // std::cout << relative_tilt_x << " " << relative_tilt_y << std::endl;
  //}

  // std::cerr << "PhotonPath: " <<  photon_direction.x << " " <<
  // photon_direction.y << " " << photon_direction.z << std::endl; std::cerr <<
  // "Microfacet: " << microfacetNormal.x << " " << microfacetNormal.y << " " <<
  // microfacetNormal.z << std::endl; std::cerr << "Incidence Angle: " <<
  // incident_angle(photon_direction, microfacetNormal) << std::endl;

  // Electric-Field E is oscillating on Y axis
  // On incidence, the electric fields polarisation is "s-polarized"
  // (perpendicular, TE) against the specimens general plane, but not regarding
  // the appearing microfacet
  double incidenceAngle_rad =
      incident_angle(photon_direction, microfacetNormal);

  // if (photon_direction.z < 0) {
  // std::cout << incidenceAngle_rad << std::endl;
  //}
  // Fresnel Percentages
  std::pair<double, double> reflectionResults;
  double R_s;
  double R_p;
  double fresnel_ratio_transmission;
  if (photon_direction.z < 0) {
    // set k to 0, as we exit into air
    // Mean of fresnel, as light leaving material is not polarized anymore
    // We will ignore k in the material; doesn't really matter for the
    // simulation WRONG reflectionResults =
    // fresnel_reflection(incidenceAngle_rad, n1, n1/n_refractive_index, 0);
    reflectionResults =
        fresnel_reflection(incidenceAngle_rad, n_refractive_index, k, n1, 0.0);
    R_s = reflectionResults.first;
    R_p = reflectionResults.second;
    // fresnel_ratio_transmission = 1 - ((R_s+R_p)/2.);
    fresnel_ratio_transmission = 1.0 - 0.5 * (R_s + R_p);

  } else {
    // Fresnel on incidence
    reflectionResults =
        fresnel_reflection(incidenceAngle_rad, n1, 0.0, n_refractive_index, k);
    R_s = reflectionResults.first;
    R_p = reflectionResults.second;
    // fresnel_ratio_transeission = 1 - ((R_s+R_p)/2.);
    // Considerate percantage of parallel and vertical polarization
    // Vertikal Tilt -> fresnel_ratio_transmission = 1 -
    // ((R_s*(1-cos(relative_tilt_x))) + (R_p*cos(relative_tilt_x))); Horizontal
    // Tilt:
    // fresnel_ratio_transmission = 1 - ((R_p*(1-cos(relative_tilt_x))) +
    // (R_s*cos(relative_tilt_x))); fresnel_ratio_transmission = 1.0 - (R_s *
    // std::cos(relative_tilt_x) * std::cos(relative_tilt_x)
    //                            + R_p * std::sin(relative_tilt_x) *
    //                            std::sin(relative_tilt_x));

    // Calculate the angle of incidence (corrected using abs for high angles)
    Direction3D inc_norm = normalize(photon_direction);
    Direction3D n_norm = normalize(microfacetNormal);
    double cos_theta_val = dot(inc_norm, n_norm);
    double theta_rad = std::acos(std::abs(cos_theta_val)); // 0..90 degree

    // Retrieve basic Fresnel values (material properties)
    auto fresnel =
        fresnel_reflection(theta_rad, n1, 0.0, n_refractive_index, k);
    double R_s = fresnel.first;
    double R_p = fresnel.second;

    // POLARISATIONS-MIXING
    // Reference for S-polarization (perpendicular to the plane of incidence)
    Direction3D s_reference = crossProduct(photon_direction, microfacetNormal);

    double len_sq = dot(s_reference, s_reference);
    double s_fraction = 1.0; // Default: 100% S (safe for vertical incidence)

    // Project only if the beam is not exactly perpendicular to the facet
    if (len_sq > 1e-8) {
      // Normalize
      double inv_len = 1.0 / std::sqrt(len_sq);
      s_reference.x *= inv_len;
      s_reference.y *= inv_len;
      s_reference.z *= inv_len;

      // Projection: How much of the laser (polarization) falls on the
      // S-reference
      double projection = dot(incidentPolarization, s_reference);

      // Intensity = amplitude squared
      s_fraction = projection * projection;
    }

    // Mix effective reflection
    double R_effective = (s_fraction * R_s) + ((1.0 - s_fraction) * R_p);

    // Calculate the transmission
    fresnel_ratio_transmission = 1.0 - R_effective;
  }

  double rnd_ = uniform01(rng);

  // IF MATERIAL IMPULSE RESPONSE IS TESTED, ALWAYS ENTER MATERIAL!
  // rnd_ is set 0.0 and therefore always lower than transmission ratio
  if (force_100_percent_transmission && photon_direction.z > 0) {
    rnd_ = 0.0;
  }

  ////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////
  //// Refract or reflect, if reflected and photon inside the material, it might
  /// redo this process

  // refraction is always possible to calculate (but still with fresnel vs.
  // snellius error/mismatch)
  Direction3D refracted;
  if (photon_direction.z < 0) {
    // if photon direction negativ: direction specimen surface
    refracted =
        refract({photon_direction.x, photon_direction.y, photon_direction.z},
                microfacetNormal, n_refractive_index, n1);
    // std::cout << (incidenceAngle_rad/3.141)*180 << " " <<
    // fresnel_ratio_transmission << " " << n_refractive_index << " " << n1
  } else {
    // if photon direction positive: about to enter specimen
    refracted =
        refract({photon_direction.x, photon_direction.y, photon_direction.z},
                microfacetNormal, n1, n_refractive_index);
  }
  // Total internal reflection (.x == -9), therefore always use reflection
  if (refracted.x == -9) {
    // std::cout << "Totalreflektion" << std::endl;
    // if (fresnel_ratio_transmission > 0.0001) std::cerr << "Err: " <<
    // fresnel_ratio_transmission << std::endl;
    if (refracted.x == -9) {
      fresnel_ratio_transmission = 0.0;
      rnd_ = 1.0; // always reflect
    }
    rnd_ = 1.;
  }

  //// Refract (into the opposite Material)
  //// If light not refracted, its reflected!
  if (rnd_ < fresnel_ratio_transmission) {
    Direction4D directionWithInfo;
    directionWithInfo.x = refracted.x;
    directionWithInfo.y = refracted.y;
    directionWithInfo.z = refracted.z;
    directionWithInfo.i = 2;

    // DEPRECATED, may delete
    // TODO: Debug these returns:
    // RESOLVED (28.03.2026):
    //  Removed global Z-direction checks: The main transport loop natively 
    //  handles resulting immediate boundary collisions.
    /*if (refracted.z > 0 && photon_direction.z > 0) {
      return directionWithInfo;
    } else if (refracted.z < 0 && photon_direction.z < 0) {
      return directionWithInfo;
    } else { // std::cout << "No valid Refraction!" << std::endl);
      return {directionWithInfo.x, directionWithInfo.y, directionWithInfo.z,
              -9};
    }*/
    return directionWithInfo;
  }
  //// Reflect (at surface, normal and totalreflection alike)
  else {
    Direction3D reflected;
    reflected = reflect(photon_direction, microfacetNormal);
    Direction4D directionWithInfo;
    directionWithInfo.x = reflected.x;
    directionWithInfo.y = reflected.y;
    directionWithInfo.z = reflected.z;
    directionWithInfo.i = 1;

    // TODO: Debug these returns:
    // RESOLVED (28.03.2026):
    //  Removed global Z-direction checks: The main transport loop natively 
    //  handles resulting immediate boundary collisions.
    /*if (reflected.z > 0 && photon_direction.z < 0) {
      // Photons direction negativ: on its way out of the specimen
      return directionWithInfo;
    } else if (reflected.z < 0 && photon_direction.z > 0) {
      // Photon about to enter the specimen
      return directionWithInfo;
    } else { // std::cout << "No valid Reflection!" << std::endl);
      return {directionWithInfo.x, directionWithInfo.y, directionWithInfo.z,
              -9};
    }*/
    return directionWithInfo;
  }
}

// Sampling Phase Functions, Chapter 15.2.3
// Based on HenyeyGreenstein::Sample_p()
// Book: "Physically Based Rendering" - Third Edition, Pharr et al, Page 898
inline Direction3D scatter(const Direction3D &wo, const double &g,
                           std::mt19937 &rng) {
  static std::uniform_real_distribution<double> uniform01(0.0, 1.0);

  const Direction3D woNorm = normalize(wo);

  // Compute cosTheta for Henyey-Greenstein sample
  double cosTheta = computeCosThetaHG(g, rng);
  cosTheta = std::max(-1.0, std::min(cosTheta, 1.0));

  // Compute sinTheta for Henyey-Greenstein sample
  const double sinTheta = std::sqrt(std::max(0.0, 1.0 - cosTheta * cosTheta));

  const double phi = 2.0 * M_PI * uniform01(rng);

  // Construct local coordinate system
  Direction3D v1, v2;
  coordinateSystem(woNorm, v1, v2);

  // Compute scattered direction
  Direction3D wi = sphericalDirection(sinTheta, cosTheta, phi, v1, v2, woNorm);

  return normalize(wi);
}

// relative paths
bool isAbsolutePath(const std::string& path);
std::string getDirectoryName(const std::string& path);
std::string joinPaths(const std::string& base, const std::string& rel);
#endif // SISI_HPP
