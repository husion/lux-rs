import numpy as np
import matplotlib.pyplot as plt
import os

# Load original data
data = np.loadtxt('../playground/data/Native-R.csv', delimiter=',')
wl_orig = data[:, 0]
spd_orig = data[:, 1]
max_orig = spd_orig.max()
spd_orig_norm = spd_orig / max_orig

# Define wavelength grid for smooth plot
wl = np.linspace(380, 780, 1000)

def gauss(wl, peak, fwhm, height):
    sig = fwhm / (2.0 * np.sqrt(2.0 * np.log(2.0)))
    return height * np.exp(-0.5 * ((wl - peak) / sig) ** 2)

# Fitted peaks
peaks = [453.962, 540.000, 613.197, 631.868, 648.087]
fwhms = [21.195, 60.000, 13.718, 9.938, 13.113]
heights = [0.196576, 0.356536, 5.37827, 15.90185, 2.25871]

# Compute model spectra
model_spd = np.zeros_like(wl)
individual_peaks = []
for p, f, h in zip(peaks, fwhms, heights):
    pk = gauss(wl, p, f, h)
    model_spd += pk
    individual_peaks.append(pk)

max_model = model_spd.max()
model_spd_norm = model_spd / max_model
individual_peaks_norm = [pk / max_model for pk in individual_peaks]

# Plot styling
fig, ax = plt.subplots(figsize=(10, 6), dpi=150)
ax.grid(True, linestyle=':', alpha=0.6)

# Plot original
ax.plot(wl_orig, spd_orig_norm, label='Original Native-R (Measured)', color='#7f8c8d', linestyle='--', alpha=0.7, linewidth=1.5)

# Plot components
colors = ['#3498db', '#2ecc71', '#e67e22', '#e74c3c', '#9b59b6']
labels = ['Blue Leakage (~454nm)', 'Green Leakage (~540nm)', 'KSF Red 1 (~613nm)', 'KSF Red 2 (~632nm)', 'KSF Red 3 (~648nm)']
for pk, col, lbl in zip(individual_peaks_norm, colors, labels):
    ax.plot(wl, pk, color=col, alpha=0.6, linewidth=1.2, label=lbl)
    ax.fill_between(wl, pk, 0, color=col, alpha=0.06)

# Plot combined model
ax.plot(wl, model_spd_norm, label='Fitted Multi-peak Model', color='#c0392b', linewidth=2.5)

# Labels & title
ax.set_title('Red Primary Spectrum: Native-R vs. Multi-peak Fitted Model', fontsize=13, fontweight='bold', pad=15)
ax.set_xlabel('Wavelength (nm)', fontsize=11)
ax.set_ylabel('Normalized Spectral Power Distribution (a.u.)', fontsize=11)
ax.set_xlim(380, 780)
ax.set_ylim(0, 1.05)
ax.legend(frameon=True, facecolor='white', framealpha=0.9, fontsize=9, loc='upper right')

plt.tight_layout()
os.makedirs('/Users/daniel/.gemini/antigravity-cli/brain/03a82375-9361-41db-8bdc-a4d36f442cd9', exist_ok=True)
plt.savefig('/Users/daniel/.gemini/antigravity-cli/brain/03a82375-9361-41db-8bdc-a4d36f442cd9/fitted_red_primary.png')
print('Plot saved successfully.')
