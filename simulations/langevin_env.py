import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

def get_rates_from_sad(k_G, zeta_G, alpha_G, k_T, b1=1.0, r=1000.0):
    """
    Calcola i rate del sistema invertendo i parametri della SAD.
    Parametrizzazione basata su: sigma_d^2 = b1 + d1
    """
    # 1. Calcolo d1 dal rapporto (alpha_G - 1) / k_G
    # Con alpha_G = 1 + 2(d1-b1)/sigma_e^2 e sigma_e^2=(b1+d1)/k_G
    R = (alpha_G - 1.0) / k_G
    d1 = b1 * (2.0 + R) / (2.0 - R)
    
    # 2. Calcolo dei rumori
    sig_d_sq = b1 + d1 # Rumore demografico
    sig_e_sq = sig_d_sq / k_G # Rumore ambientale (parametro libero)
    
    # 3. Calcolo b0 (fissando d0=1)
    d0 = 1.0 
    b0 = 0.5 * zeta_G * sig_d_sq + d0
    
    # 4. Calcolo gamma per mantenere k_T = k_G * r / gamma
    gamma = r * (k_G / k_T)
    
    return {'b0': b0, 'b1': b1, 'd0': d0, 'd1': d1, 
            'sig_e_sq': sig_e_sq, 'sig_d_sq': sig_d_sq, 
            'r': r, 'gamma': gamma}

def run_simulation(p, dt=0.001, total_steps=1000000, skip=1000):
    """
    Simulazione Eulero-Maruyama con campionamento decorrelato (skip).
    """
    xs = np.zeros(total_steps // skip)
    ys = np.zeros(total_steps // skip)
    
    x, y = 1.0, 1.0 # Condizioni iniziali
    
    # Pre-generazione rumore per velocità
    noise = np.random.normal(0, np.sqrt(dt), (total_steps, 2))
    
    idx = 0
    for t in range(total_steps):
        # Drift e Diffusione DNA (x)
        drift_x = (p['b0'] - p['d0']) + (p['b1'] - p['d1']) * x
        diff_x = np.sqrt(max(0, p['sig_d_sq'] * x + p['sig_e_sq'] * x**2))
        
        # Drift e Diffusione mRNA (y)
        drift_y = p['r'] * x - p['gamma'] * y
        diff_y = np.sqrt(max(0, p['r'] * x + p['gamma'] * y))
        
        # Update
        x += drift_x * dt + diff_x * noise[t, 0]
        y += drift_y * dt + diff_y * noise[t, 1]
        
        # Boundary (densità positive)
        x = max(x, 1e-9); y = max(y, 1e-9)
        
        if t % skip == 0 and idx < len(xs):
            xs[idx], ys[idx] = x, y
            idx += 1
            
    return xs, ys

# --- Esempio di utilizzo ---
# Input SAD
params_sad = {'k_G': 1.2, 'zeta_G': 0.4, 'alpha_G': 1.15, 'k_T': 12.0}

# Ottenimento Rate
rates = get_rates_from_sad(**params_sad)

# Simulazione
x_sim, y_sim = run_simulation(rates)

# Salvataggio dati decorrelati
pd.DataFrame({'DNA': x_sim, 'mRNA': y_sim}).to_csv('langevin_decorrelated_data.csv', index=False)

# Plotting
fig, ax = plt.subplots(1, 2, figsize=(14, 5))

# SAD Plot
sns.kdeplot(x_sim, ax=ax[0], log_scale=True, label='DNA', fill=True)
sns.kdeplot(y_sim, ax=ax[0], log_scale=True, label='mRNA', fill=True)
ax[0].set_title("Species Abundance Distributions")
ax[0].legend()

# Heatmap Correlazione
hb = ax[1].hexbin(x_sim, y_sim, gridsize=30, cmap='magma', bins='log', xscale='log', yscale='log')
fig.colorbar(hb, ax=ax[1], label='log10(count)')
ax[1].set_xlabel("DNA concentration (x)")
ax[1].set_ylabel("mRNA concentration (y)")
ax[1].set_title(f"Correlation Heatmap (r={rates['r']}, γ={rates['gamma']:.1f})")

plt.tight_layout()
plt.show()
