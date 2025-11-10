import pandas as pd
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter, butter, filtfilt

# Load your data file (adjust path if needed)
df = pd.read_csv("fit.data", comment="#", delim_whitespace=True, header=None, names=["t", "d", "d_fit"])

# --- 1. Moving average smoothing ---
df["d_movavg"] = df["d"].rolling(window=200, center=True).mean()

# --- 2. Savitzky-Golay filter ---
df["d_savgol"] = savgol_filter(df["d"], window_length=11, polyorder=4)

# --- 3. Butterworth low-pass filter ---
def butter_lowpass_filter(data, cutoff, fs, order=4):
    nyquist = 0.5 * fs
    normal_cutoff = cutoff / nyquist
    b, a = butter(order, normal_cutoff, btype='low', analog=False)
    return filtfilt(b, a, data)

fs = 1 / (df["t"].iloc[1] - df["t"].iloc[0])  # Sampling frequency
cutoff = 0.05  # Low cutoff frequency for smoothing
df["d_butter"] = butter_lowpass_filter(df["d"], cutoff, fs)

# --- Plot all results ---
plt.figure(figsize=(10, 6))
plt.plot(df["t"], df["d"], alpha=0.4, label="Original (noisy)")
plt.plot(df["t"], df["d_movavg"], label="Moving Average",linestyle='-',linewidth=0.7)
plt.plot(df["t"], df["d_savgol"], label="Savitzky-Golay",linestyle=':' ,linewidth=0.7)
plt.plot(df["t"], df["d_butter"], label="Butterworth Low-Pass",linestyle=':',linewidth=0.7)
plt.xlabel("Time (t)")
plt.ylabel("d")
plt.title("Comparison of Smoothing Methods for Noisy Data")
plt.legend()
plt.grid(True)
plt.show()

