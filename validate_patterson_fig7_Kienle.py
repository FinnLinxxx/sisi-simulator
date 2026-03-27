import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import c as c_vacuum

# ==========================================
# 1. KONFIGURATION
# ==========================================
# Dateipfad zur Simulation
data_file = "/home/finn/sisi-simulator/sisi_data/patterson/1.42_0.0045_0.6000_0.10/dist_1.42_0.0045_0.6000_0.10_0.00.txt"

# ANPASSUNG AN FIGURE 7
# ---------------------
# 1. Y-Achse: Skalierung auf ca. 4500 Counts (visueller Peak im Paper)
TARGET_PEAK_COUNTS = 4500.0 

# 2. Form: Faltung mit der Unschärfe des Lasers (180 ps FWHM)
FWHM_IRF_ps = 180.0
sigma_IRF_ns = (FWHM_IRF_ps / 1000.0) / 2.355 

# 3. Position: Zeit-Versatz (Offset)
# Der "System"-Puls (blaue Linie in Fig 7) hat seinen Peak bei ca. 0.25 ns.
# Um diesen Wert schieben wir unsere t=0 Simulation nach rechts.
TIME_OFFSET_NS = 0.25 

# Physikalische Parameter (Wet Wood, parallel)
mu_a_mm = 0.0045
sigma_s_mm = 6.0    # Eingabe Simulation (mu_s)
g = 0.9             # Anisotropie
n = 1.41            # Brechungsindex
rho = 25.0          # Abstand Quelle-Detektor

# ==========================================
# 2. PHYSIK BERECHNEN (Patterson)
# ==========================================
c = (c_vacuum / n) * 1e-9 
mu_s_prime = sigma_s_mm * (1 - g) # Sollte 0.6 ergeben
D = 1.0 / (3.0 * (mu_a_mm + mu_s_prime))
z0 = 1.0 / mu_s_prime

# Theorie-Funktion (Impulsantwort)
def patterson_R_rho_t(t_ps):
    t = np.maximum(t_ps, 1e-5)
    term1 = (4 * np.pi * D * c)**(-1.5)
    term2 = z0 * (t**(-2.5))
    term3 = np.exp(-mu_a_mm * c * t)
    dist_sq = rho**2 + z0**2
    term4 = np.exp(-dist_sq / (4 * D * c * t))
    return term1 * term2 * term3 * term4

# Zeitachse berechnen (Fein aufgelöst für saubere Faltung)
# Wir berechnen von 0 bis 6 ns, Plotten später aber verschoben
t_theory_ns = np.linspace(0.0, 6.0, 2000) 
t_theory_ps = t_theory_ns * 1000.0
R_theory_ideal = patterson_R_rho_t(t_theory_ps)

# ==========================================
# 3. FALTUNG (CONVOLUTION)
# ==========================================
# Gauß-Kern erstellen
dt_ns = t_theory_ns[1] - t_theory_ns[0]
kernel_radius = int(5 * sigma_IRF_ns / dt_ns)
x_kernel = np.arange(-kernel_radius, kernel_radius+1) * dt_ns
kernel = np.exp(-x_kernel**2 / (2 * sigma_IRF_ns**2))
kernel /= np.sum(kernel)

# Theorie falten (Glätten)
R_theory_convolved = np.convolve(R_theory_ideal, kernel, mode='same')

# Skalieren auf Zielhöhe (Counts)
R_theory_final = (R_theory_convolved / np.max(R_theory_convolved)) * TARGET_PEAK_COUNTS

# ==========================================
# 4. SIMULATION DATEN
# ==========================================
print("Lade Simulationsdaten...")
try:
    raw_data = np.loadtxt(data_file)
    if raw_data.ndim == 1: raw_data = raw_data.reshape(1, -1)
    dist_meters = raw_data[:, 1]
except Exception as e:
    print(f"WARNUNG: Nutze Dummy-Daten ({e})")
    # Fallback damit Code läuft
    dist_meters = np.random.lognormal(np.log(0.5), 0.6, 5000) * 0.3 

# Umrechnung in Zeit (ns)
times_ns = (dist_meters * 1000.0 / c) / 1000.0

# Histogramm erstellen
bins_count = 120
hist_counts, bin_edges = np.histogram(times_ns, bins=bins_count, range=(0, 6.0))
bin_centers = 0.5 * (bin_edges[1:] + bin_edges[:-1])

# Simulation falten (Gleiche Physik wie Theorie)
bin_width_ns = bin_edges[1] - bin_edges[0]
k_bins_radius = int(4 * sigma_IRF_ns / bin_width_ns)
if k_bins_radius < 1: k_bins_radius = 1
k_bins = np.arange(-k_bins_radius, k_bins_radius+1) * bin_width_ns
kernel_hist = np.exp(-k_bins**2 / (2 * sigma_IRF_ns**2))
kernel_hist /= np.sum(kernel_hist)

hist_convolved = np.convolve(hist_counts, kernel_hist, mode='same')

# Skalieren
if np.max(hist_convolved) > 0:
    hist_final = (hist_convolved / np.max(hist_convolved)) * TARGET_PEAK_COUNTS
else:
    hist_final = hist_convolved

# ==========================================
# 5. PLOTTEN MIT OFFSET
# ==========================================
plt.figure(figsize=(10, 7))

# Wir addieren hier TIME_OFFSET_NS auf die X-Achse!

# A) Simulation (Rot)
mask = hist_final > 1e-1 
plt.semilogy(bin_centers[mask] + TIME_OFFSET_NS,  # <--- HIER IST DER SHIFT
             hist_final[mask], 
             'rd', markersize=4, label='MC Simulation (Convolved + Shifted)')

# B) Theorie (Schwarz)
valid = t_theory_ns > 0.05
plt.semilogy(t_theory_ns[valid] + TIME_OFFSET_NS, # <--- HIER IST DER SHIFT
             R_theory_final[valid], 
             'k-', linewidth=2.5, label='Patterson Model (Convolved + Shifted)')

# Layout wie im Paper
plt.title(f'Validierung Wet Wood ($\\rho={rho}$ mm)\n'
          f'Offset={TIME_OFFSET_NS}ns | IRF={int(FWHM_IRF_ps)}ps')
plt.xlabel('time (ns)')
plt.ylabel('relative reflectance [counts]')
plt.xlim(0, 5)
plt.ylim(1e0, 1e4)
plt.legend()
plt.grid(True, which="both", linestyle='--', alpha=0.5)
plt.tight_layout()
plt.show()
