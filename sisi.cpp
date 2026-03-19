/*
 * ==============================================================================
 * SISI Simulator - Signal Interference from Surface Interactions
 * Imperial grace, in metric space.
 * Copyright (c) 2023-2026, Finn Linzer, TU Wien, Engineering Geodesy.
 * ==============================================================================
 *
 * Hello There!
 * This software represents the physical models and Monte Carlo transport
 * architecture developed and validated within my dissertation. While I have
 * made every effort to ensure the scientific integrity of this simulation,
 * it remains an open research tool and is provided without legal warranties.
 * Users are expected to apply their own scientific rigor when configuring
 * simulations and interpreting the results. I assume no liability for
 * experimental, physical, or commercial damages resulting from the use
 * of this software.
 * -----------------------------------------------------------------------------
 * License & Commercial Use:
 * SISI is free software: you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later
 * version. Full license text: See LICENSE file.
 * -----------------------------------------------------------------------------
 * Acknowledgments:
 * - Phase Function Sampling (Henyey-Greenstein) and Coordinate System
 * routines are adapted from "Physically Based Rendering: From Theory To
 * Implementation" (Third Edition) by Matt Pharr, Wenzel Jakob, and Greg
 * Humphreys, which is licensed under the BSD License (https://pbrt.org).
 *
 * As required by their BSD License, the following notice applies exclusively to
 * those specific parts of the code:
 *
 * Copyright (c) 1998-2015, Matt Pharr, Greg Humphreys, and Wenzel Jakob.
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation
 * and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 * - The conceptual structure of the physical Monte Carlo light transport loop
 * was inspired by the educational tutorials provided by www.scratchapixel.com.
 * (https://www.scratchapixel.com/lessons/mathematics-physics-for-computer-graphics/monte-carlo-methods-in-practice/monte-carlo-simulation.html)
 */

#include "sisi.hpp"
#include "photon_log.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <yaml-cpp/yaml.h>

// Suppress Windows cmd freeze
#ifdef _WIN32
#include <windows.h>
#endif

// Enable make directory ability
#ifdef _WIN32
#include <direct.h> //  _mkdir for Windows
#else
#include <sys/stat.h> //  mkdir for Linux
#endif

// Diagnostic/Output switches from .yaml (#define look-a-like)
bool TRACK_PHOTON_PATHS = false;
bool PRINT_RESULTS = false;

bool SAVE_PATH_LENGTHS = false;
bool SAVE_REGISTERED_SPOT = false;
bool SAVE_HEMISPHERE = false;
bool SAVE_PHOTON_PATHS = false;
bool SAVE_RESULTS = false;


bool isAbsolutePath(const std::string& path) {
  if (path.empty()) return false;

#ifdef _WIN32
  // C:\...  or  \\server\share\...
  if (path.size() >= 2 &&
      std::isalpha(static_cast<unsigned char>(path[0])) &&
      path[1] == ':') {
    return true;
  }
  if (path.size() >= 2 &&
      (path[0] == '\\' || path[0] == '/') &&
      (path[1] == '\\' || path[1] == '/')) {
    return true;
  }
  return false;
#else
  return path[0] == '/';
#endif
}

std::string getDirectoryName(const std::string& path) {
  std::string normalized = path;
  std::replace(normalized.begin(), normalized.end(), '\\', '/');

  std::size_t pos = normalized.find_last_of('/');
  if (pos == std::string::npos) {
    return ".";
  }
  if (pos == 0) {
    return "/";
  }
  return normalized.substr(0, pos);
}

std::string joinPaths(const std::string& base, const std::string& rel) {
  if (base.empty() || base == ".") return rel;
  if (base.back() == '/' || base.back() == '\\') return base + rel;
  return base + "/" + rel;
}

///////////////////////////////////////
// Photon-Log (mainly for sisi_vis.py)
///////////////////////////////////////

PhotonScriptLogger g_log;
static uint64_t g_pid_base = 0; // Photon ID across all angles

static std::string make_end_label(bool via_direct, bool via_subsurface,
                                  bool had_TIR, EndReason reason) {
  if (reason == Absorbed)
    return "Absorbed";

  if (via_direct) {
    return (reason == Hit) ? "DirectReflectedAndHit" : "DirectReflectedAndGone";
  }

  if (via_subsurface) {
    if (had_TIR) {
      return (reason == Hit)
                 ? "SubsurfaceReflectedAndInnerTotalreflectionAndHit"
                 : "SubsurfaceReflectedAndInnerTotalreflectionAndGone";
    } else {
      return (reason == Hit) ? "SubsurfaceReflectedAndHit"
                             : "SubsurfaceReflectedAndGone";
    }
  }

  // Fallback
  switch (reason) {
  case Hit:
    return "Hit";
  case Gone:
    return "Gone";
  case Absorbed:
    return "Absorbed";
  default:
    return "Gone";
  }
}

static std::string make_end_note(bool via_direct, bool via_subsurface,
                                 bool had_TIR, EndReason reason, int scat_count,
                                 int refl_count) {
  std::ostringstream meta;
  meta << make_end_label(via_direct, via_subsurface, had_TIR, reason)
       << ";scat=" << scat_count << ";refl=" << refl_count
       << ";tir=" << (had_TIR ? 1 : 0);
  return meta.str();
}

///////////////////////////////////////
// YAML Parser Start
///////////////////////////////////////
SisiConfig loadConfig(const std::string &filename) {
  YAML::Node yaml = YAML::LoadFile(filename);
  SisiConfig config;

  // 1. Read Mode (laserscanner, material_impulse, patterson)
  std::string mode_str = yaml["mode"].as<std::string>();
  YAML::Node mode_node;

  if (mode_str == "laserscanner") {
    config.mode = SimulationMode::LASERSCANNER;
    mode_node = yaml["laserscanner"];
  } else if (mode_str == "material_impulse") {
    config.mode = SimulationMode::MATERIAL_IMPULSE;
    mode_node = yaml["material_impulse"];
  } else if (mode_str == "patterson") {
    config.mode = SimulationMode::PATTERSON;
    mode_node = yaml["patterson"];
  } else {
    throw std::runtime_error("ERROR: Unknown mode in YAML file: " + mode_str);
  }

  // 2. Read common parameter set
  config.simulation.use_pseudo_random_seed =
      yaml["simulation"]["use_pseudo_random_seed"].as<bool>();
  config.simulation.n_max_detections =
      yaml["simulation"]["n_max_detections"].as<uint32_t>();

  config.light_source.wavelength_nm =
      yaml["light_source"]["wavelength_nm"].as<double>();
  config.light_source.n_photons =
      yaml["light_source"]["n_photons"].as<uint32_t>();

  config.material.n_refractive_index =
      yaml["material"]["n_refractive_index"].as<double>();
  config.material.absorption_coefficient_per_mm =
      yaml["material"]["absorption_coefficient_per_mm"].as<double>();
  config.material.scattering_coefficient_per_mm =
      yaml["material"]["scattering_coefficient_per_mm"].as<double>();
  config.material.anisotropy_g = yaml["material"]["anisotropy_g"].as<double>();
  config.material.specimen_thickness_m =
      yaml["material"]["specimen_thickness_m"].as<double>();
  config.material.microfacet_std_deviation_deg =
      yaml["material"]["microfacet_std_deviation_deg"].as<double>();

  // 3. Read active mode parameter (mode_node)
  config.active_mode.optic.radius_m =
      mode_node["optic"]["radius_m"].as<double>();
  // pinhole_sigma_factor only valid for mode: laserscanner, material_impulse
  // if not available, its set to 0.0 (deactivated)
  config.active_mode.optic.pinhole_sigma_factor =
      mode_node["optic"]["pinhole_sigma_factor"]
          ? mode_node["optic"]["pinhole_sigma_factor"].as<double>()
          : 0.0;
  // Fresnel-Equation disabled if true -> Always Transmission into material
  config.active_mode.physics_override.force_100_percent_transmission =
      mode_node["physics_override"]["force_100_percent_transmission"]
          .as<bool>();

  if (config.mode == SimulationMode::PATTERSON) {
    // In Patterson mode, the detector is placed directly on the surface
    // (z = 0). A tiny offset (0.1 µm) is hard-coded to avoid numerical
    // boundary collisions.
    config.active_mode.optic.working_distance_m = 0.0000001;
    config.active_mode.optic.rho_m = mode_node["optic"]["rho_m"].as<double>();
    // In Patterson mode, all photons originate from the exact same point
    config.active_mode.beam.sigma_x_m = 0.0;
    config.active_mode.beam.sigma_y_m = 0.0;
    config.active_mode.beam.sigma_z_m = 0.0;
    // In Patterson mode, an angle of incidence is not applicable
    config.active_mode.scan.angles_deg.push_back(0.0);
  } else if (config.mode == SimulationMode::LASERSCANNER ||
             config.mode == SimulationMode::MATERIAL_IMPULSE) {
    config.active_mode.optic.working_distance_m =
        mode_node["optic"]["working_distance_m"].as<double>();
    config.active_mode.optic.rho_m = 0.0;

    config.active_mode.beam.sigma_x_m =
        mode_node["beam"]["sigma_x_m"].as<double>();
    config.active_mode.beam.sigma_y_m =
        mode_node["beam"]["sigma_y_m"].as<double>();
    config.active_mode.beam.sigma_z_m =
        mode_node["beam"]["sigma_z_m"].as<double>();

    // Read angle of incidence as list
    for (auto angle : mode_node["scan"]["angles_deg"]) {
      config.active_mode.scan.angles_deg.push_back(angle.as<double>());
    }
  } else {
    throw std::runtime_error("Fatal: Unhandled SimulationMode in config.yaml!");
  }

  // 4. Diagnostics FLAG
  config.active_mode.diagnostics.track_photon_paths =
      mode_node["diagnostics"]["track_photon_paths"].as<bool>();
  config.active_mode.diagnostics.print_results =
      mode_node["diagnostics"]["print_results"].as<bool>();

  // 5. Output FLAG
  config.active_mode.output.save_path_lengths =
      mode_node["output"]["save_path_lengths"].as<bool>();
  config.active_mode.output.save_registered_spot =
      mode_node["output"]["save_registered_spot"].as<bool>();
  config.active_mode.output.save_hemisphere =
      mode_node["output"]["save_hemisphere"].as<bool>();
  config.active_mode.output.save_photon_paths =
      mode_node["output"]["save_photon_paths"].as<bool>();
  config.active_mode.output.save_results =
      mode_node["output"]["save_results"].as<bool>();
  config.active_mode.output.result_path =
      mode_node["output"]["result_path"].as<std::string>();

  config.active_mode.output.hemisphere_radius_m =
      mode_node["output"]["hemisphere_radius_m"].as<double>();

  return config;
}
// YAML Parser END

void SISI_Simulation(double theta_in_deg, const SisiConfig &config,
                     std::mt19937 &rng) {
  ////////////////////////////////////////////////////////////////
  // Yaml/Read config start
  uint32_t nPhotons = config.light_source.n_photons;
  double n_refractive_index = config.material.n_refractive_index;
  // Laser Wavelength
  const double lambda_m = config.light_source.wavelength_nm * 1e-9;
  // Absorption and Scattering Coefficient
  double mu_a = config.material.absorption_coefficient_per_mm * 1000.0; // [1/m]
  double mu_s = config.material.scattering_coefficient_per_mm * 1000.0; // [1/m]
  // Geometry and Optics
  double working_distance_m = config.active_mode.optic.working_distance_m;
  double radius_m = config.active_mode.optic.radius_m;
  double pinhole_sigma_factor = config.active_mode.optic.pinhole_sigma_factor;
  // Laser-beam parameter, normal-distribution values
  double sigma_x_m = config.active_mode.beam.sigma_x_m;
  double sigma_y_m = config.active_mode.beam.sigma_y_m;
  double sigma_z_m = config.active_mode.beam.sigma_z_m;
  // Material
  // _deg for external readability
  double microfacet_std_deviation_deg =
      config.material.microfacet_std_deviation_deg;
  // _rad for internal/math calculations
  double microfacet_std_deviation_rad =
      (config.material.microfacet_std_deviation_deg * M_PI) / 180.0;
  double specimen_thickness_m = config.material.specimen_thickness_m;
  double hg_parameter = config.material.anisotropy_g;
  // Control
  std::string result_path = config.active_mode.output.result_path;
  uint32_t n_max_detections = config.simulation.n_max_detections;
  bool force_transmission =
      config.active_mode.physics_override.force_100_percent_transmission;
  double hemisphere_radius_m = config.active_mode.output.hemisphere_radius_m;
  // Patterson
  double rho_m = config.active_mode.optic.rho_m;
  ///////////////////////////////////////////////////////////////
  // Yaml Read End

  // Extinction coefficient
  double mu_t = mu_a + mu_s; // [1/m]
  // k (Imaginary part of the refractive index for Fresnel)
  // Derived physically from μ_a and the wavelength
  double k = (mu_a * lambda_m) / (4.0 * M_PI); // No Dimension
  // Hard-Coded refraction index of air
  double n1 = 1.00028;

  // clang-format off
  // Global photon fate counters
  uint32_t nPhotonsShot                    = 0; // Total number of photons emitted into the scene
  uint32_t nPhotonsAbsorbed                = 0; // Terminated inside the material due to absorption (mu_a)
  uint32_t nPhotonsTransmitted             = 0; // Exited through to the opposite side of the specimen
  uint32_t nPhotonsReflected               = 0; // Bounced directly off the incident surface
  uint32_t nPhotonsBackscatteredFromWithin = 0; // Scattered back to the detector side from within

  // Detector specific counters (HO = Hit Optic) 
  uint32_t nPhotonsHitOptic                = 0; // Total number of photons registered by the detector
  uint32_t nPhotonsHO_Refl                 = 0; // Sensor hits originating from direct surface reflection
  uint32_t nPhotonsHO_scat                 = 0; // Sensor hits originating from subsurface scattering
  // clang-format on

  double theta_in_rad = (theta_in_deg / 180.0) * M_PI;

  // Determine the detector's position (Optical Center)
  Position3D opticalCenter;
  if (config.mode == SimulationMode::PATTERSON) {
    // working_distance_m is set very small loading config above
    // Detector is placed at a radial distance of rho_m.
    // IMPORTANT: optic.radius_m in config.yaml must remain small (e.g., ~2 mm)!
    opticalCenter = {rho_m, 0.0, -working_distance_m};
  } else if (config.mode == SimulationMode::LASERSCANNER ||
             config.mode == SimulationMode::MATERIAL_IMPULSE) {
    // Standard setup: 2D rotation of the scanner around the Y-axis.
    // (Rotations around X or Z are currently not supported).
    opticalCenter = {-std::sin(theta_in_rad) * working_distance_m, 0.0,
                     -std::cos(theta_in_rad) * working_distance_m};
  } else {
    throw std::runtime_error(
        "Fatal: Unhandled SimulationMode for optical center calculation!");
  }

  // Uniform Distribution between 0 and 1
  std::uniform_real_distribution<double> uniform01(0.0, 1.0);
  // Exp Distribution (derived from mu_t)
  std::exponential_distribution<double> expDist(mu_t);

  ///////////////////////////////////////
  // 1x Make Directory + 4x DEBUG Data filestream provider
  ///////////////////////////////////////

  // Generate config specific sub-folder
  char folderPart[256];
  snprintf(folderPart, sizeof(folderPart), "%.2f_%.4f_%.4f_%.2f",
           n_refractive_index, mu_a / 1000., mu_s / 1000.,
           microfacet_std_deviation_deg);

  std::string full_path = result_path + "/" + folderPart;
// Debug output paths
//std::cout << "full_path   = " << full_path << std::endl;

// Create the output directory
// Cross-platform handling for Windows and Linux)
#ifdef _WIN32
  _mkdir(full_path.c_str());
#else
  mkdir(full_path.c_str(), 0777);
#endif

  // 1. SAVE_PATH_LENGTHS (Path lengths & scatter counts)
  char filenameDist[512];
  snprintf(filenameDist, sizeof(filenameDist),
           "%s/dist_%.2f_%.4f_%.4f_%.2f_%.2f.txt", full_path.c_str(),
           n_refractive_index, mu_a / 1000., mu_s / 1000.,
           microfacet_std_deviation_deg, theta_in_deg);

  std::ofstream distFile;
  if (SAVE_PATH_LENGTHS) {
    distFile.open(filenameDist, std::ios::out);
    distFile << std::fixed << std::setprecision(6);
    if (!distFile.is_open()) {
      std::cout << "Failed to open " << filenameDist << std::endl;
      std::cout
          << "Make sure result_path as stated in config.yaml is available: "
          << result_path << std::endl;
    }
  }

  // 2. SAVE_REGISTERED_SPOT (3D detection spot)
  char filenameSpot[512];
  snprintf(filenameSpot, sizeof(filenameSpot),
           "%s/spot_%.2f_%.4f_%.4f_%.2f_%.2f.txt", full_path.c_str(),
           n_refractive_index, mu_a / 1000., mu_s / 1000.,
           microfacet_std_deviation_deg, theta_in_deg);

  std::ofstream spotFile;
  if (SAVE_REGISTERED_SPOT) {
    spotFile.open(filenameSpot, std::ios::out | std::ios::app);
    spotFile << "// x_surf y_surf z_surf x_vol y_vol z_vol "
            "inner_dist scatter_count type\n";
    spotFile << std::fixed << std::setprecision(6);
    if (!spotFile.is_open()) {
      std::cout << "Failed to open " << filenameSpot << std::endl;
      std::cout
          << "Make sure result_path as stated in config.yaml is available: "
          << result_path << std::endl;
    }
  }

  // 3. SAVE_HEMISPHERE (5m-Half-Sphere for CloudCompare)
  char filenameHemi[512];
  snprintf(filenameHemi, sizeof(filenameHemi),
           "%s/hemisphere_%.2f_%.4f_%.4f_%.2f_%.2f.txt", full_path.c_str(),
           n_refractive_index, mu_a / 1000., mu_s / 1000.,
           microfacet_std_deviation_deg, theta_in_deg);

  std::ofstream hemiFile;
  if (SAVE_HEMISPHERE) {
    hemiFile.open(filenameHemi, std::ios::out);
    hemiFile << "// X Y Z R G B OriginType HitStatus InnerDist OuterDist "
                "TotalDist\n";
    hemiFile << std::fixed << std::setprecision(6);
    if (!hemiFile.is_open()) {
      std::cout << "Failed to open " << filenameHemi << std::endl;
      std::cout
          << "Make sure result_path as stated in config.yaml is available: "
          << result_path << std::endl;
    }
  }

  // 4. SAVE_PHOTON_PATHS (Subsurface Debugging / Paths)
  char filenamePaths[512];
  snprintf(filenamePaths, sizeof(filenamePaths),
           "%s/paths_%.2f_%.4f_%.4f_%.2f_%.2f.txt", full_path.c_str(),
           n_refractive_index, mu_a / 1000., mu_s / 1000.,
           microfacet_std_deviation_deg, theta_in_deg);

  std::ofstream pathsFile;
  if (SAVE_PHOTON_PATHS) {
    pathsFile.open(filenamePaths, std::ios::out);
    pathsFile << std::fixed << std::setprecision(6);
    if (!pathsFile.is_open()) {
      std::cout << "Failed to open " << filenamePaths << std::endl;
      std::cout
          << "Make sure result_path as stated in config.yaml is available: "
          << result_path << std::endl;
    }
  }

  ///////////////////////////////////////
  // Single Photon Production Centre SPPC
  ///////////////////////////////////////
  for (uint32_t n = 0; n < nPhotons; ++n) {
    const uint64_t pid = g_pid_base + static_cast<uint64_t>(n);

    bool end_via_direct = false;            // Exits via direct reflection
    bool end_via_subsurface = false;        // Exits via subsurface scattering
    bool had_inner_totalreflection = false; // Whether TIR occurred

    int scat_count = 0; // Scatter count in the material
    int refl_count = 0; // Number of reflection events (surface + TIR)

    uint16_t s_total_count = 0; // Number of distinct path segments flown
    double s_total = 0; // Cumulative distance traveled inside the material

    nPhotonsShot++;
    try {

      double photon_inner_life_distance = 0;
      double photon_outer_life_distance = 0;

      Direction4D directionAfterMaterialTransition;
      Position3D currentPos;
      Direction3D currentDir;

      ///////////////////////////////////////
      // Patterson Mode Start
      ///////////////////////////////////////
      if (config.mode == SimulationMode::PATTERSON) {
        // Calculate the starting depth z0
        // Equals z0 = 1 / mu_s' (with mu_s' = mu_s * (1 - g))
        double mu_s_prime = mu_s * (1.0 - hg_parameter);
        double z0 = 1.0 / mu_s_prime;

        currentPos = {0.0, 0.0, z0};

        // Generate isotropic vector for initial scattering
        double rnd_phi = 2.0 * M_PI * uniform01(rng);
        double rnd_cos_theta = 2.0 * uniform01(rng) - 1.0;
        double rnd_sin_theta = std::sqrt(1.0 - rnd_cos_theta * rnd_cos_theta);

        currentDir = {rnd_sin_theta * std::cos(rnd_phi),
                      rnd_sin_theta * std::sin(rnd_phi), rnd_cos_theta};

        // Initialize photon trajectory directly within the medium.
        // The state (.i = 2) bypasses the initial surface boundary check
        // and starts the subsurface transport loop.
        directionAfterMaterialTransition = {currentDir.x, currentDir.y,
                                            currentDir.z, 2};

        g_log.log(pid, InitPhotonPose,
                  Pose{currentPos.x, currentPos.y, currentPos.z, currentDir.x,
                       currentDir.y, currentDir.z},
                  "PattersonZ0");
      } else if (config.mode == SimulationMode::LASERSCANNER ||
                 config.mode == SimulationMode::MATERIAL_IMPULSE) {
        ///////////////////////////////////////
        // laserscanner or material_impulse Start
        ///////////////////////////////////////

        // Initialize a new 3D coordinate system by defining three orthogonal
        // unit vectors (basisX, basisY, basisZ). This setup is used for
        // positioning a new photon with Gaussian distribution around its
        // origin. Gram–Schmidt orthogonalization:
        Direction3D basisZ =
            normalize(positional2directional(opticalCenter, 1));
        Direction3D basisY;
        // Gimbal Lock Check
        if (std::abs(basisZ.y) > 0.99999) {
          basisY = {1.0, 0.0, 0.0};
          // Jump in direction is negligible, as when IA near 0, R_p = R_s
        } else {
          basisY = {0.0, 1.0, 0.0};
        }
        Direction3D basisX = crossProduct(basisY, basisZ);
        basisX = normalize(basisX);
        basisY = crossProduct(basisZ, basisX);
        basisY = normalize(basisY);

        Position3D newPhotonEmissionLocation = determinePhotonEmissionLocation(
            opticalCenter, sigma_x_m, sigma_y_m, sigma_z_m, basisX, basisY,
            basisZ, rng);

        g_log.log(pid, InitPhotonPose,
                  Pose{newPhotonEmissionLocation.x, newPhotonEmissionLocation.y,
                       newPhotonEmissionLocation.z, basisZ.x, basisZ.y,
                       basisZ.z},
                  "0");

        Position3D newPhotonHitSurfaceLocation =
            determinePhotonHitSurfaceLocation(newPhotonEmissionLocation,
                                              basisZ);

        photon_outer_life_distance += sqrt(
            pow(newPhotonEmissionLocation.x - newPhotonHitSurfaceLocation.x,
                2) +
            pow(newPhotonEmissionLocation.y - newPhotonHitSurfaceLocation.y,
                2) +
            pow(newPhotonEmissionLocation.z - newPhotonHitSurfaceLocation.z,
                2));

        // -------------------------------------------------------------------
        // Material Boundary Transition
        // Return states (.i):
        //   1 : Reflected back at surface
        //   2 : Refracted into the other medium
        //  -9 : Physically invalid geometry (e.g. into shadow), retry
        // -------------------------------------------------------------------
        do {
          Direction3D microfacetNormal;
          // Generate a random microfacet normal
          // Reject if normal shows away from the incoming photon
          do {
            microfacetNormal =
                rnd_microfacet(microfacet_std_deviation_rad, basisZ.z, rng);
          } while (dot(basisZ, microfacetNormal) < 0);
          // Calculate reflection or refraction based on Fresnel equations
          directionAfterMaterialTransition = materialTransition(
              basisZ, n1, n_refractive_index, k, microfacetNormal, basisY,
              force_transmission, rng);
        } while (directionAfterMaterialTransition.i == -9);

        // directionAfterMaterialTransition could be reflection or refraction
        currentPos = newPhotonHitSurfaceLocation;
        currentDir = {directionAfterMaterialTransition.x,
                      directionAfterMaterialTransition.y,
                      directionAfterMaterialTransition.z};
      }

      Direction3D backscatteredPhoton_direction;
      Position3D backscatteredPhoton_position;

      ///////////////////////////////////////
      // Photon inside Material
      ///////////////////////////////////////
      if (directionAfterMaterialTransition.i == 2) {

        g_log.log(pid, RefractedIntoMaterial,
                  Pose{currentPos.x, currentPos.y, currentPos.z, currentDir.x,
                       currentDir.y, currentDir.z},
                  "i1");

        while (1) {
          // Travelled distance to new event (mu_t = mu_a + mu_s)
          double s = expDist(rng);

          s_total += s;
          s_total_count += 1;
          double distToBoundary = 0;

          // When photon leaves specimen in the back (which is unlikely in case
          // of a high thickness)
          if (currentDir.z > 0) {
            distToBoundary = (specimen_thickness_m - currentPos.z) /
                             currentDir.z; // Points towards the bottom surface
            if (s > distToBoundary) {
              nPhotonsTransmitted++;
              break;
            }
          } else if (currentDir.z <= 0) {
            // Points towards the top surface
            distToBoundary = -currentPos.z / currentDir.z;
          }
          // Photon may leave specimen on top surface
          if (s > distToBoundary && (currentDir.z < 0)) {
            do {
              Direction3D microfacetNormal;
              // Generate microfacet normal for the "inner" surface
              // Reject impossible microfacets
              do {
                microfacetNormal = rnd_microfacet(microfacet_std_deviation_rad,
                                                  currentDir.z, rng);
              } while (dot(currentDir, microfacetNormal) < 0);

              directionAfterMaterialTransition = materialTransition(
                  currentDir, n1, n_refractive_index, k, microfacetNormal,
                  {-1.0, -10.00, -100.000}, force_transmission, rng);
              // Dummy vector: Polarization is lost after internal scattering
            } while (directionAfterMaterialTransition.i == -9);

            // Photon is leaving the specimen and is refracted at currentPos
            if (directionAfterMaterialTransition.i == 2) {
              // Photon leaves specimen, volumetric position
	      double lastX = currentPos.x;
	      double lastY = currentPos.y;
              double lastZ = currentPos.z;
              // Photon leaves specimen, surface position (genuine laserspot)
              currentPos =
                  add(currentPos, scalarMultiply(currentDir, distToBoundary));

              // PINHOLE LOGIC
              // ---------------------------------------------------------
              bool passed_spatial_filter = true;

              // The check is performed only if the factor is greater than 0
              if (pinhole_sigma_factor > 0.0001) {
                // Check the elliptical distance to the center (0,0)
                // Formel: (x/Rx)^2 + (y/Ry)^2 <= 1

                double limit_x = sigma_x_m * pinhole_sigma_factor;
                double limit_y = sigma_y_m * pinhole_sigma_factor;

                // Avoiding division by zero
                if (limit_x > 0 && limit_y > 0) {
                  double dist_sq =
                      (currentPos.x * currentPos.x) / (limit_x * limit_x) +
                      (currentPos.y * currentPos.y) / (limit_y * limit_y);

                  if (dist_sq > 1.0) {
                    passed_spatial_filter = false;
                  }
                }
              }

              backscatteredPhoton_position = currentPos;
              backscatteredPhoton_direction = {
                  directionAfterMaterialTransition.x,
                  directionAfterMaterialTransition.y,
                  directionAfterMaterialTransition.z};

              photon_inner_life_distance += (s_total - s + distToBoundary);
              g_log.log(pid, PhotonLeavesSpecimen,
                        Pose{currentPos.x, currentPos.y, currentPos.z,
                             backscatteredPhoton_direction.x,
                             backscatteredPhoton_direction.y,
                             backscatteredPhoton_direction.z},
                        "gs" + std::to_string(s_total_count));

              /////////////////////////////////////////////////////////////////////////////
              // Surface Exit Logging
              // Track all backscattered photons exiting the material
              // An infinite radius (999999.9) ensures every exiting photon is
              // projected onto the surface for the pathsFile log.
              Position4D backscatteredPhoton_all = backscatteringMayHitOptic(
                  backscatteredPhoton_position, backscatteredPhoton_direction,
                  {0, 0, -working_distance_m}, 999999.9, config.mode);
              // Berechne die Distanz von der Oberfläche bis zu dieser Ebene
              double dist_surface_to_plane =
                  sqrt(pow(backscatteredPhoton_position.x -
                               backscatteredPhoton_all.x,
                           2) +
                       pow(backscatteredPhoton_position.y -
                               backscatteredPhoton_all.y,
                           2) +
                       pow(backscatteredPhoton_position.z -
                               backscatteredPhoton_all.z,
                           2));
              if (SAVE_PHOTON_PATHS) {
                pathsFile << "SubSurfaceScatteredALL "
                          << backscatteredPhoton_all.x << " "
                          << backscatteredPhoton_all.y << " "
                          << backscatteredPhoton_all.z << " "
                          << "0 0 " //  Placeholder (formerly diodeX, diodeY)
                          << currentDir.x << " " << currentDir.y << " "
                          << currentDir.z << " " << photon_inner_life_distance
                          << " " << s_total_count << "\n";
              }
              /////////////////////////////////////////////////////////////////////////////

              // Position4D backscatteredPhoton =
              // backscatteringMayHitOptic(backscatteredPhoton_position,
              // backscatteredPhoton_direction, opticalCenter,radius_m,
              // config.mode);
              Position4D backscatteredPhoton;

              // -------------- SPATIAL FILTER
              if (passed_spatial_filter) {
                //  Calculate only when the photon passes through the pinhole
                backscatteredPhoton = backscatteringMayHitOptic(
                    backscatteredPhoton_position, backscatteredPhoton_direction,
                    opticalCenter, radius_m, config.mode);
              } else {
                // Spatial filter blocked -> Photon is ignored (no hit)
                backscatteredPhoton = {0, 0, 0, -1};
              }
              if (SAVE_HEMISPHERE) {
                const double R_hemi = hemisphere_radius_m;
                Position3D P = backscatteredPhoton_position;
                Direction3D D = backscatteredPhoton_direction;

                double b = 2.0 * (P.x * D.x + P.y * D.y + P.z * D.z);
                double c =
                    (P.x * P.x + P.y * P.y + P.z * P.z) - R_hemi * R_hemi;
                double discriminant = b * b - 4.0 * c;
                double d_hemi = (discriminant >= 0)
                                    ? (-b + std::sqrt(discriminant)) / 2.0
                                    : 0;

                Position3D shell_hit = {P.x + d_hemi * D.x, P.y + d_hemi * D.y,
                                        P.z + d_hemi * D.z};

                int hit_status = (backscatteredPhoton.i == 1) ? 1 : 0;
                double total_outer = photon_outer_life_distance + d_hemi;
                double total_path = photon_inner_life_distance + total_outer;

                hemiFile << shell_hit.x << " " << shell_hit.y << " "
                         << shell_hit.z << " "
                         << (hit_status ? "0 255 0 "
                                        : "255 0 0 ") // green Hit!, Red Miss X
                         << "2 " << hit_status
                         << " " // OriginType 2 = Subsurface
                         << photon_inner_life_distance << " " << total_outer
                         << " " << total_path << "\n";
              }

              // HIT!
              if (backscatteredPhoton.i == 1) {
                nPhotonsHitOptic++;
                nPhotonsHO_scat++;
                double photon_from_scatter_to_optics =
                    sqrt(pow(currentPos.x - backscatteredPhoton.x, 2) +
                         pow(currentPos.y - backscatteredPhoton.y, 2) +
                         pow(currentPos.z - backscatteredPhoton.z, 2));
                photon_outer_life_distance += photon_from_scatter_to_optics;

                if (SAVE_PATH_LENGTHS) {
                  distFile << photon_outer_life_distance << " "
                           << photon_inner_life_distance << " " << s_total_count
                           << std::endl;
                }

                if (SAVE_REGISTERED_SPOT) {
		  spotFile << currentPos.x << " " << currentPos.y << " " << currentPos.z << " "
         << lastX << " " << lastY << " " << lastZ << " " << photon_inner_life_distance << " " << s_total_count << " 2 " << std::endl;
                }

                // DiodeHitPlane?:
                g_log.log(pid, SubSurfaceScattered,
                          Pose{backscatteredPhoton.x, backscatteredPhoton.y,
                               backscatteredPhoton.z, currentDir.x,
                               currentDir.y, currentDir.z},
                          "hs" + std::to_string(s_total_count));
                // g_log.logEnd(pid, Hit);

                end_via_subsurface = true;

                g_log.logEnd(pid, Hit,
                             make_end_note(end_via_direct, end_via_subsurface,
                                           had_inner_totalreflection, Hit,
                                           scat_count, refl_count));

                if (SAVE_PHOTON_PATHS) {
                  pathsFile
                      << "HIT " << currentPos.x << " " << currentPos.y << " "
                      << currentPos.z << " " << backscatteredPhoton.x << " "
                      << backscatteredPhoton.y << " " << backscatteredPhoton.z
                      << " " << photon_outer_life_distance << " "
                      << photon_inner_life_distance << " "
                      << (photon_outer_life_distance +
                          photon_inner_life_distance)
                      << " " << s_total_count << "\n";
                }
              } else { // NO HIT!
                // DiodeHitPlane?:
                g_log.log(pid, SubSurfaceScattered,
                          Pose{currentPos.x + backscatteredPhoton_direction.x *
                                                  working_distance_m,
                               currentPos.y + backscatteredPhoton_direction.y *
                                                  working_distance_m,
                               currentPos.z + backscatteredPhoton_direction.z *
                                                  working_distance_m,
                               currentDir.x, currentDir.y, currentDir.z},
                          "gs" + std::to_string(s_total_count));
                end_via_subsurface = true;

                g_log.logEnd(pid, Gone,
                             make_end_note(end_via_direct, end_via_subsurface,
                                           had_inner_totalreflection, Gone,
                                           scat_count, refl_count));
              }

              nPhotonsBackscatteredFromWithin++;
              break;
            } else if (directionAfterMaterialTransition.i == 1) {

              // Photon got reflected from the inner surface
              // Reflect: direction * s;
              currentPos =
                  add(currentPos, scalarMultiply(currentDir, distToBoundary));
              g_log.log(pid, FresnelReflPointIN,
                        Pose{currentPos.x, currentPos.y, currentPos.z,
                             currentDir.x, currentDir.y, currentDir.z},
                        "f" + std::to_string(s_total_count));
              had_inner_totalreflection = true;
              refl_count++;

              currentDir = {directionAfterMaterialTransition.x,
                            directionAfterMaterialTransition.y,
                            directionAfterMaterialTransition.z};
              currentPos = add(
                  currentPos, scalarMultiply(currentDir, (s - distToBoundary)));

              g_log.log(pid, FresnelReflPointENDIN,
                        Pose{currentPos.x, currentPos.y, currentPos.z,
                             currentDir.x, currentDir.y, currentDir.z},
                        "o" + std::to_string(s_total_count));

            } else {
              /*not_possible*/
            }
          } // if (s > distToBoundary && (muz < 0))
          // Photon does not leave specimen on top surface, its still inside
          else {
            currentPos = add(currentPos, scalarMultiply(currentDir, s));
            // First direction, after going into the material
            // if ( s_total_count == 1) { std::cout << currentPos.x << " " <<
            // currentPos.y << " " << currentPos.z << std::endl; }
          }

          if (uniform01(rng) < mu_a / mu_t) {
            nPhotonsAbsorbed++;

            if (SAVE_PHOTON_PATHS) {
              pathsFile
                  << "ABSORBED " << currentPos.x << " " << currentPos.y << " "
                  << currentPos.z << " "
                  << "0 0 0 0 0 0 " // Set the remaining values to 0 (as before)
                  << s_total_count << "\n";
            }

            g_log.log(pid, AbsorbedInMaterial,
                      Pose{currentPos.x, currentPos.y, currentPos.z,
                           currentDir.x, currentDir.y, currentDir.z},
                      "a" + std::to_string(s_total_count));
            // g_log.logEnd(pid, Absorbed);
            g_log.logEnd(pid, Absorbed,
                         make_end_note(end_via_direct, end_via_subsurface,
                                       had_inner_totalreflection, Absorbed,
                                       scat_count, refl_count));

            break;
          }

          currentDir = scatter(currentDir, hg_parameter, rng);
          g_log.log(pid, ScatteredInMaterial,
                    Pose{currentPos.x, currentPos.y, currentPos.z, currentDir.x,
                         currentDir.y, currentDir.z},
                    "s" + std::to_string(s_total_count));
          scat_count++;

          if (SAVE_PHOTON_PATHS) {
            pathsFile << "SCATTERED " << currentPos.x << " " << currentPos.y
                      << " " << currentPos.z << " "
                      << "0 0 0 0 0 0 " << s_total_count << "\n";
          }
        }
      }
      // Reflection from Surface
      // ///////////////////////////////////////////////////////
      else if (directionAfterMaterialTransition.i == 1) {

        g_log.log(pid, ReflectedOnMaterial,
                  Pose{currentPos.x, currentPos.y, currentPos.z, currentDir.x,
                       currentDir.y, currentDir.z},
                  "r");
        refl_count++;

        nPhotonsReflected++;
        // Photon got backscattered right from the surface
        Position4D backscatteredPhoton = backscatteringMayHitOptic(
            currentPos, currentDir, opticalCenter, radius_m, config.mode);
        if (SAVE_HEMISPHERE) {
          const double R_hemi = hemisphere_radius_m;
          Position3D P = currentPos;
          Direction3D D = currentDir;

          double b = 2.0 * (P.x * D.x + P.y * D.y + P.z * D.z);
          double c = (P.x * P.x + P.y * P.y + P.z * P.z) - R_hemi * R_hemi;
          double discriminant = b * b - 4.0 * c;
          double d_hemi =
              (discriminant >= 0) ? (-b + std::sqrt(discriminant)) / 2.0 : 0;

          Position3D shell_hit = {P.x + d_hemi * D.x, P.y + d_hemi * D.y,
                                  P.z + d_hemi * D.z};

          int hit_status = (backscatteredPhoton.i == 1) ? 1 : 0;
          double total_outer = photon_outer_life_distance + d_hemi;

          hemiFile << shell_hit.x << " " << shell_hit.y << " " << shell_hit.z
                   << " " << (hit_status ? "0 255 0 " : "255 0 0 ") << "1 "
                   << hit_status << " "            // OriginType 1 = Direct
                   << "0.0 " << total_outer << " " // Inner is always 0
                   << total_outer << "\n";         // Total = Outer
        }

        backscatteredPhoton_position = currentPos;
        backscatteredPhoton_direction = {directionAfterMaterialTransition.x,
                                         directionAfterMaterialTransition.y,
                                         directionAfterMaterialTransition.z};

        if (backscatteredPhoton.i == 1) {

          nPhotonsHO_Refl++;
          nPhotonsHitOptic++;
          double photon_from_scatter_to_optics =
              sqrt(pow(currentPos.x - backscatteredPhoton.x, 2) +
                   pow(currentPos.y - backscatteredPhoton.y, 2) +
                   pow(currentPos.z - backscatteredPhoton.z, 2));
          photon_outer_life_distance += photon_from_scatter_to_optics;

          if (SAVE_PATH_LENGTHS) {
            distFile << photon_outer_life_distance << " "
                     << photon_inner_life_distance << " " << s_total_count
                     << std::endl;
          }

          if (SAVE_REGISTERED_SPOT) {
	    spotFile << currentPos.x << " " << currentPos.y << " " << currentPos.z << " "
         <<  currentPos.x<< " " << currentPos.y << " " << currentPos.z << " " << " " << photon_inner_life_distance << " " << s_total_count << " 1 " << std::endl;
          }

          g_log.log(pid, ReflectedDirect,
                    Pose{currentPos.x + currentDir.x * working_distance_m,
                         currentPos.y + currentDir.y * working_distance_m,
                         currentPos.z + currentDir.z * working_distance_m,
                         currentDir.x, currentDir.y, currentDir.z},
                    "hr");

          // g_log.logEnd(pid, Hit);
          end_via_direct = true;

          g_log.logEnd(pid, Hit,
                       make_end_note(end_via_direct, end_via_subsurface,
                                     had_inner_totalreflection, Hit, scat_count,
                                     refl_count));
        } else { // NO HIT!
          g_log.log(pid, ReflectedDirect,
                    Pose{currentPos.x + currentDir.x * working_distance_m,
                         currentPos.y + currentDir.y * working_distance_m,
                         currentPos.z + currentDir.z * working_distance_m,
                         currentDir.x, currentDir.y, currentDir.z},
                    "gr");
          // g_log.logEnd(pid, Gone);
          end_via_direct = true;

          g_log.logEnd(pid, Gone,
                       make_end_note(end_via_direct, end_via_subsurface,
                                     had_inner_totalreflection, Gone,
                                     scat_count, refl_count));
        }
        backscatteredPhoton_position = currentPos;
        backscatteredPhoton_direction = {directionAfterMaterialTransition.x,
                                         directionAfterMaterialTransition.y,
                                         directionAfterMaterialTransition.z};
      }

      // throw exception;
    } // try
    catch (...) {
      std::cerr << "Error " << n << std::endl;
    }
    if (nPhotonsHitOptic >= n_max_detections)
      break;
  } // for (int n = 0; n < nPhotons; ++n)

  // clang-format off
  if (PRINT_RESULTS) {
    static bool header_printed = false;
    if (!header_printed) {
        std::cout << "\n" << std::left 
                  << std::setw(8)  << "Angle" 
                  << std::setw(10) << "nPhot" 
                  << std::setw(10) << "nHIT" 
                  << std::setw(10) << "nBS" 
                  << std::setw(10) << "nRefl" 
                  << std::setw(10) << "nAbs" 
                  << std::setw(10) << "nHORefl" 
                  << std::setw(10) << "nHOScat" 
                  << std::setw(10) << "nShot" 
                  << std::setw(10) << "nTrans" << std::endl;
        std::cout << std::string(98, '-') << std::endl;
        header_printed = true;
    }

    std::cout << std::left 
              << std::setw(8)  << theta_in_deg 
              << std::setw(10) << nPhotons 
              << std::setw(10) << nPhotonsHitOptic
              << std::setw(10) << nPhotonsBackscatteredFromWithin 
              << std::setw(10) << nPhotonsReflected 
              << std::setw(10) << nPhotonsAbsorbed 
              << std::setw(10) << nPhotonsHO_Refl 
              << std::setw(10) << nPhotonsHO_scat 
              << std::setw(10) << nPhotonsShot 
              << std::setw(10) << nPhotonsTransmitted << std::endl;
  }
  // clang-format on

  // SAVE_RESULTS
  char filename[512];
  snprintf(filename, sizeof(filename), "%s/%.2f_%.4f_%.4f_%.2f.txt",
           full_path.c_str(), n_refractive_index, mu_a / 1000., mu_s / 1000.,
           microfacet_std_deviation_deg);

  std::ofstream sFile;
  if (SAVE_RESULTS) {
    sFile.open(filename, std::ios::out | std::ios::app);
    sFile << std::fixed << std::setprecision(6);

    if (!sFile.is_open()) {
      std::cout << "Failed to open " << filename << std::endl;
      std::cout
          << "Make sure result_path as stated in config.yaml is available: "
          << result_path << std::endl;
    } else {
      // Write the exact values from PRINT_RESULTS to the file
      sFile << theta_in_deg << " " << nPhotons << " " << nPhotonsHitOptic << " "
            << nPhotonsBackscatteredFromWithin << " " << nPhotonsReflected
            << " " << nPhotonsAbsorbed << " " << nPhotonsHO_Refl << " "
            << nPhotonsHO_scat << " " << nPhotonsShot << " "
            << nPhotonsTransmitted << "\n";

      sFile.close();
    }
  }
  // Close all open files
  if (SAVE_PATH_LENGTHS && distFile.is_open()) {
    distFile.close();
  }
  if (SAVE_REGISTERED_SPOT && spotFile.is_open()) {
    spotFile.close();
  }
  if (SAVE_HEMISPHERE && hemiFile.is_open()) {
    hemiFile.close();
  }
  if (SAVE_PHOTON_PATHS && pathsFile.is_open()) {
    pathsFile.close();
  }
} // SISI_Simulation(...)

int main(int argc, char **argv) {
  if (argc != 2) {
    std::cerr << "=========================================================="
              << std::endl;
    std::cerr << "ERROR: No yaml-configuration file specified!" << std::endl;
    std::cerr << "Usage: sisi.exe <path_to_config.yaml>" << std::endl;
    std::cerr << "\nI am looking for the file in this directory:" << std::endl;
    std::cerr << argv[0] << std::endl;
    std::cerr << "=========================================================="
              << std::endl;

#ifdef _WIN32
    std::cout << "\nPress Enter to exit..." << std::endl;
    std::cin.get();
#endif
    return 1;
  }

#ifdef _WIN32
  // Dont freeze Console Window on Windows
  HANDLE hInput = GetStdHandle(STD_INPUT_HANDLE);
  DWORD prev_mode;
  GetConsoleMode(hInput, &prev_mode);
  // Disable QuickEdit-Mode
  SetConsoleMode(hInput,
                 ENABLE_EXTENDED_FLAGS | (prev_mode & ~ENABLE_QUICK_EDIT_MODE));
#endif

std::string config_file = argv[1];
SisiConfig config = loadConfig(config_file);

std::string config_dir = getDirectoryName(config_file);
std::string& result_path = config.active_mode.output.result_path;

if (!isAbsolutePath(result_path)) {
  result_path = joinPaths(config_dir, result_path);
}


// Debug output paths
//std::cout << "config_file = " << config_file << std::endl;
//std::cout << "config_dir  = " << config_dir << std::endl;
//std::cout << "result_path = " << result_path << std::endl;

  // Global Switch, Debugging
  TRACK_PHOTON_PATHS = config.active_mode.diagnostics.track_photon_paths;
  PRINT_RESULTS = config.active_mode.diagnostics.print_results;

  SAVE_PATH_LENGTHS = config.active_mode.output.save_path_lengths;
  SAVE_REGISTERED_SPOT = config.active_mode.output.save_registered_spot;
  SAVE_HEMISPHERE = config.active_mode.output.save_hemisphere;
  SAVE_RESULTS = config.active_mode.output.save_results;
  SAVE_PHOTON_PATHS = config.active_mode.output.save_photon_paths;

  // C++ style Random Number Generator (rng)
  std::mt19937 rng;
  if (config.simulation.use_pseudo_random_seed) {
    rng.seed(1337);
  } else {
    rng.seed(std::random_device{}());
  }

  if (TRACK_PHOTON_PATHS)
    g_log.setEnabled(true);
  else
    g_log.setEnabled(false);

  for (size_t i = 0; i < config.active_mode.scan.angles_deg.size(); i++) {
    double angle = config.active_mode.scan.angles_deg[i];
    SISI_Simulation(angle, config, rng);
    g_pid_base += config.light_source.n_photons;
  }

  if (TRACK_PHOTON_PATHS) {
    g_log.writeJSONL(config.active_mode.output.result_path +
                     "/photon_paths.jsonl");
  }
#ifdef _WIN32
  std::cout << "\nSimulation finished. Press Enter to exit..." << std::endl;
  std::cin.get();
#endif

  return 0;
}
