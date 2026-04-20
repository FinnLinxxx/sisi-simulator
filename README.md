<div align="center">

# SISI Simulator 

✧ **S**ignal **I**nterference from **S**urface **I**nteractions ✧

<br>

*A C++ Monte Carlo photon transport simulator for modeling light propagation in various materials and spatial setups* <br> *— in an SI world.*

For the blog entry and web version of SISI (2D) [click here](https://finnlinxxx.github.io/sisi/index.html)

<br>

[![C++11](https://img.shields.io/badge/C++-11-00599C?style=flat-square&logo=c%2B%2B&logoColor=white)](https://isocpp.org/)
[![CMake](https://img.shields.io/badge/CMake-Build-064F8C?style=flat-square&logo=cmake&logoColor=white)](https://cmake.org/)

</div>

---

## Description
The **SISI Simulator** is a C++ Monte Carlo photon transport simulator, testing propagation in various materials and setups. It supports three operational modes: classic laserscanner , material impulse response testing, and the [Patterson measurement setup](https://opg.optica.org/ao/fulltext.cfm?uri=ao-28-12-2331). The simulation relies on a configuration file (`config.yaml`) to dictate physical parameters, optics, beam profiles, and output formatting.

## Installation

<details open>
<summary><b> Linux </b></summary>
<br>

To install and build the project on Linux, use a terminal to clone the repository and build it via CMake:

```bash
# Clone repository
$ cd ~
$ git clone https://git.geo.tuwien.ac.at/finnlinxxx/sisi-simulator.git
$ cd sisi-simulator

# Copy your own configuration file from the example
$ cp config.example.yaml config.yaml
(you may adapt the config.yaml to your setup, see below)

# Build directory and compile
$ mkdir build
$ cd build 
$ cmake ..
$ make

# Run the simulator using the config file as an argument
$ ./sisi ../config.yaml

# Alternatively, copy the executable to the main folder and run it there
$ cd ..
$ cp build/sisi .
$ ./sisi config.yaml
```

</details>

<details>
<summary><b> Windows </b></summary>
<br>

On Windows, build and run the simulator was tested using Visual Studio 2017 with its built-in CMake support.

**1. Clone the repository**

Open PowerShell and clone the project to your desired directory:
```powershell
git clone https://git.geo.tuwien.ac.at/finnlinxxx/sisi-simulator.git
cd sisi-simulator
```

**2. Setup Configuration**
Create your own configuration file by copying the example:
Copy `config.example.yaml` and rename it to `config.yaml` in the `sisi-simulator` folder.

**3. Build in Visual Studio**
* Open Visual Studio.
* Go to **File > Open > Folder...** and select the `sisi-simulator` folder.
* Wait a moment for Visual Studio to parse the CMake environment.
* At the top toolbar, set the Target to **`x64-Release`**.
* In the top menu, click **CMake > Build All** (Building takes ~30 seconds). You should see "Build succeeded".

**4. Setup Run/Debug Configuration**
* On the right side in the **Solution Explorer** (Folder View), find `CMakeLists.txt`.
* Right-click it > **Debug and Launch Settings** > select `sisi.exe`.
* This opens the `launch.vs.json` file. Replace its entire content with the text below (which is also available in `sisi-simulator/tools/visual_studio/launch.vs.example.json`):

```json
{
  "version": "0.2.1",
  "configurations": [
    {
      "type": "default",
      "project": "CMakeLists.txt",
      "projectTarget": "sisi.exe",
      "name": "Run sisi (workspace config)",
      "args": [ "${workspaceRoot}\\config.yaml" ]
    }
  ]
}
```
* Save (Ctrl+S) and close the file.

**5. Run the Simulator**
* At the top toolbar, next to the green "Play" button, click the small dropdown arrow ("Select Startup Item").
* Choose **"Run sisi (workspace config)"**.
* Press **Play**.
* A console will start and process your `config.yaml`. Results will be saved relative to your workspace in `sisi_data/...` (Note: the actual `.exe` is stored in a hidden `CMakeBuild` AppData folder, but the launch configuration ensures it runs in the correct workspace directory).
* Press Enter or X to close the console when finished. 
* **Warning:** Running the simulation again will overwrite existing files in the output directory!

</details>


The simulator is controlled via a `config.yaml` file. You can choose between three main modes: `laserscanner`, `material_impulse`, and `patterson`. 


### Configuration Parameters Explained

A brief overview of what each parameter in the `config.yaml` file does:

* **`mode`**: Defines which simulation type runs (`laserscanner`, `material_impulse`, or `patterson`).
* **`simulation`**:
    * `use_pseudo_random_seed`: If `true`, the simulation uses a fixed seed, ensuring identical results on repeat runs.
    * `n_max_detections`: The maximum number of photon detections allowed before the simulation forcibly stops.
* **`light_source`**:
    * `wavelength_nm`: The wavelength of the simulated light (in nanometers).
    * `n_photons`: The total number of photons to simulate.
* **`material`**: Defines the optical properties of the target.
    * `n_refractive_index`: Index of refraction.
    * `absorption_coefficient_per_mm`: Likelihood of a photon being absorbed per mm.
    * `scattering_coefficient_per_mm`: Likelihood of a photon scattering per mm.
    * `anisotropy_g`: The scattering directionality (0.0 = isotropic, 1.0 = mainly forward).
    * `specimen_thickness_m`: Thickness of the material (in meter).
    * `microfacet_std_deviation_deg`: Surface roughness standard deviation in degrees.
* **`optic` / `beam`**:
    * `radius_m`: Radius of the receiving optic.
    * `working_distance_m`: Distance from the optic to the specimen.
    * `pinhole_sigma_factor`: Filters stray photons (Focus); `0.0` disables it.
    * `rho_m`: Specific parameter for Patterson mode defining receiver distance offset.
    * `sigma_x_m / y_m / z_m`: Standard deviations defining the laser beam's gaussian profile.
* **`physics_override`**:
    * `force_100_percent_transmission`: If `true`, forces photons to enter the material, ignoring specular reflection.
* **`scan`**:
    * `angles_deg`: A list of incidence angles to simulate (e.g., `[0, 1, 2]`).
* **`diagnostics` & `output`**:
    * Toggles (`true`/`false`) to print output to the console or save results: track photon paths, save registered spots, log paths or impact on a theoretical hemisphere.
    * `result_path`: Directory where output files will be saved (e.g., `sisi_data/laserscanner`). Warning: files here get overwritten on subsequent runs.
