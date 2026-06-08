import numpy as np
import matplotlib.pyplot as plt
import os

wl = np.linspace(380, 780, 1000)

def gauss(wl, peak, fwhm, height=1.0):
    sig = fwhm / (2.0 * np.sqrt(2.0 * np.log(2.0)))
    return height * np.exp(-0.5 * ((wl - peak) / sig) ** 2)

# Isolated Red (KSF phosphor)
red = 0.340923 * gauss(wl, 613.197, 13.718) + 1.008005 * gauss(wl, 631.868, 9.938) + 0.143178 * gauss(wl, 648.087, 13.113)
red /= red.max()

# Isolated Green
green = gauss(wl, 531.4658, 45.1469)
green /= green.max()

# Isolated Blue
blue = gauss(wl, 454.9756, 26.1472)
blue /= blue.max()

# Plot
fig, ax = plt.subplots(figsize=(10, 6), dpi=150)
ax.grid(True, linestyle=':', alpha=0.6)

# Plot curves
ax.plot(wl, blue, color='#3498db', linewidth=2.5, label='Isolated Blue Primary (peak=455.0nm, fwhm=26.1nm)')
ax.fill_between(wl, blue, 0, color='#3498db', alpha=0.1)

ax.plot(wl, green, color='#2ecc71', linewidth=2.5, label='Isolated Green Primary (peak=531.5nm, fwhm=45.1nm)')
ax.fill_between(wl, green, 0, color='#2ecc71', alpha=0.1)

ax.plot(wl, red, color='#e74c3c', linewidth=2.5, label='Isolated Red Primary (KSF, peak=631.9nm)')
ax.fill_between(wl, red, 0, color='#e74c3c', alpha=0.1)

# Labels & title
ax.set_title('RGB Isolated Primary Spectra (Leakage-Free)', fontsize=13, fontweight='bold', pad=15)
ax.set_xlabel('Wavelength (nm)', fontsize=11)
ax.set_ylabel('Normalized Spectral Power Distribution (a.u.)', fontsize=11)
ax.set_xlim(380, 780)
ax.set_ylim(0, 1.05)
ax.legend(frameon=True, facecolor='white', framealpha=0.9, fontsize=9.5)

plt.tight_layout()
os.makedirs('/Users/daniel/.gemini/antigravity-cli/brain/03a82375-9361-41db-8bdc-a4d36f442cd9', exist_ok=True)
plt.savefig('/Users/daniel/.gemini/antigravity-cli/brain/03a82375-9361-41db-8bdc-a4d36f442cd9/rgb_isolated_primaries.png')
print('RGB plot saved successfully.')
