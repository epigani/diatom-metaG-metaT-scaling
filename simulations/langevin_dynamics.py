import numpy as np
import numba as nb
import os
from concurrent.futures import ProcessPoolExecutor
import sys
import logging
import scipy.optimize as opt
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.special import hyp2f1

@nb.njit
def stochastic_simulation(x0, y0, r0, r1, gamma1, b0, b1, d0, d1, Delta, T, dt):
    def W10(x, y):
        if y < Delta:
            return b0 + b1 * x
        elif (y >= Delta) & (x <= 0):
            return b0
        else:
            return b0 + b1 * x - d0 - d1 * x

    def W01(x, y):
        if y < Delta:
            return r1 * y - gamma1 * y + (b0 + b1 * x) * Delta + r0 * x
        elif (y >= Delta) & (x <= 0):
            return r1 * y - gamma1 * y + b0 * Delta
        else:
            return r1 * y - gamma1 * y + (b0 + b1 * x - d0 - d1 * x) * Delta + r0 * x

    def W20(x, y):
        if y < Delta:
            return b0 + b1 * x
        elif (y >= Delta) & (x <= 0):
            return b0
        else:
            return b0 + b1 * x + d0 + d1 * x

    def W02(x, y):
        if y < Delta:
            return r1 * y  + gamma1 * y + (b0 + b1 * x) * (Delta ** 2) + r0 * x
        elif (y >= Delta) & (x <= 0):
            return r1 * y + gamma1 * y + b0 * (Delta ** 2)
        else:
            return r1 * y + gamma1 * y + (b0 + b1 * x + d0 + d1 * x) * (Delta ** 2) + r0 * x

    def W11(x, y):
        if y < Delta:
            return (b0 + b1 * x) * Delta
        elif (y >= Delta) & (x <= 0):
            return (b0) * Delta
        else:
            return (b0 + b1 * x + d0 + d1 * x) * Delta

    x = np.zeros(T)
    y = np.zeros(T)
    x[0] = x0
    y[0] = y0

    IdMatrix = np.eye(2)
    noise = np.zeros(2)
    G = np.zeros((2, 2))

    for t in range(1, T):
        dx = W10(x[t-1], y[t-1]) * dt
        dy = W01(x[t-1], y[t-1]) * dt

        W20_element = W20(x[t-1], y[t-1])
        W02_element = W02(x[t-1], y[t-1])
        W11_element = W11(x[t-1], y[t-1])
        
        if t == 1:
            W20_el = W20_element
            W02_el = W02_element
            W11_el = W11_element

        if Delta == 0:
            G[0, 0] = np.sqrt(W20_element)
            G[1, 1] = np.sqrt(W02_element)
            G[0, 1] = 0.0
            G[1, 0] = 0.0
        else:
            trace_W2 = W20_element + W02_element
            det_W2 = W20_element * W02_element - W11_element ** 2

            sqrt_detW2 = np.sqrt(det_W2)
            factor = 1.0 / np.sqrt(trace_W2 + 2 * sqrt_detW2)

            G[0, 0] = factor * (W20_element + sqrt_detW2)
            G[1, 1] = factor * (W02_element + sqrt_detW2)
            G[0, 1] = factor * W11_element
            G[1, 0] = factor * W11_element

        noise[0] = np.random.normal(0, np.sqrt(dt))
        noise[1] = np.random.normal(0, np.sqrt(dt))

        dx += (G[0, 0] * noise[0] + G[0, 1] * noise[1])
        dy += (G[1, 0] * noise[0] + G[1, 1] * noise[1])

        x[t] = max(x[t-1] + dx, 0)
        y[t] = max(y[t-1] + dy, 0)

    return x, y, W20_el, W02_el, W11_el

def find_uncorrelated_interval(data, dt, threshold=1e-2):
    def autocorrelation_fft(x):
        n = len(x)
        x = x - np.mean(x)
        f = np.fft.fft(x, n=2*n)
        acf = np.fft.ifft(f * np.conjugate(f))[:n].real
        acf /= acf[0]
        return acf

    subsampled_data = data[::int(1 / dt)]
    autocorr = autocorrelation_fft(subsampled_data)
    try:
        tau_time = np.argmax(autocorr < threshold)
    except ValueError:
        tau_time = len(autocorr)
    tau_steps = int(tau_time / dt)
    return tau_steps, tau_time

def log_pdf(data, nbins=20):
    if len(data) < 3:
        return np.array([]), np.array([]), np.array([]), np.array([])
    else:
        positive_data = data[data > 0]
        bins = np.logspace(np.log10(np.min(positive_data)), np.log10(np.max(positive_data)), nbins)
        hist, _ = np.histogram(positive_data, bins=bins, density=False)
        prob = hist / np.sum(hist)
        pdf = prob / np.diff(bins)
        x_plot = np.sqrt(bins[1:] * bins[:-1])
        y_plot = pdf
        x_plot = x_plot[y_plot > 0]
        y_plot = y_plot[y_plot > 0]
        return x_plot, y_plot, pdf, bins 

def fit_log_generalized_gamma(x, y):
    def log_generalized_gamma(X, a, b, C):
        return np.log(C) + (-1+a)*X - b*np.exp(X)
    if len(x) < 3:
        return np.array([np.nan, np.nan, np.nan])
    else:
        X = np.log(x)
        Y = np.log(y)
        popt, _ = opt.curve_fit(log_generalized_gamma, X, Y, p0=[1, 1, 1], bounds=([-5, 0, 1e-9], [5, np.inf, np.inf]))
        return popt

def plotSAD(r0, r1, gamma1, b0, b1, d0, d1, Delta, n_x, PDF_x, n_y, PDF_y, popt_x, popt_y, chain=0):
    def generalized_gamma(x, a, b, C):
        return C * x**(-1+a) * np.exp(-b*x)
    fig, axes = plt.subplots(1, 2, figsize=(figsizes["2 columns"][0], 0.5*figsizes["2 columns"][1]), constrained_layout=True)
    
    title = f'r0={r0:.1e}, r1={r1:.1e}, gamma1={gamma1:.1e}, b0={b0:.1e}, b1={b1:.1e}, d0={d0:.1e}, d1={d1:.1e}, Delta={Delta:.1e}'
    fig.suptitle(title, fontsize=7)

    axes[0].set_xscale('log')
    axes[0].set_yscale('log')
    text = f'$a={popt_x[0]:.2f}$\n$b={popt_x[1]:.2e}$\n$C={popt_x[2]:.2e}$'
    axes[0].text(0.1, 0.1, text, transform=axes[0].transAxes, fontsize=8, verticalalignment='bottom')
    axes[0].scatter(n_x, PDF_x, color=light_palette[0], alpha=1, edgecolor='k', lw=0.5, s=15)  
    if not np.isnan(popt_x).any(): 
        axes[0].plot(n_x, generalized_gamma(n_x, *popt_x), color=light_palette[0], alpha=0.5)  
        axes[0].set_xlabel("DNA")
        axes[0].set_ylabel("PDF")

    axes[1].set_xscale('log')
    axes[1].set_yscale('log')
    text = f'$a={popt_y[0]:.2f}$\n$b={popt_y[1]:.2e}$\n$C={popt_y[2]:.2e}$'
    axes[1].text(0.1, 0.1, text, transform=axes[1].transAxes, fontsize=8, verticalalignment='bottom')
    axes[1].scatter(n_y, PDF_y, color=light_palette[1], alpha=1, edgecolor='k', lw=0.5, s=15)    
    if not np.isnan(popt_y).any():
        axes[1].plot(n_y, generalized_gamma(n_y, *popt_y), color=light_palette[1], alpha=0.5)  
        axes[1].set_xlabel("RNA")
        axes[1].set_ylabel("PDF")

    figname = f'figures/r0_{r0:.4f}_r1_{r1:.4f}_gamma1_{gamma1:.4f}_b0_{b0:.4f}_b1_{b1:.4f}_d0_{d0:.4f}_d1_{d1:.4f}_Delta_{Delta:.4f}_chain_{chain}.png'
    fig.savefig(figname, dpi=200)
    plt.close()
    
def plotCORR(r0, r1, gamma1, b0, b1, d0, d1, Delta, DNA_xplot, RNA_yplot, chain=0):
    fig, ax = plt.subplots(1, figsize=(figsizes["1 column"][0], figsizes["1 column"][0]))
    mask = (DNA_xplot > 0) & (RNA_yplot > 0)
    DNA_xplot = DNA_xplot[mask]
    RNA_yplot = RNA_yplot[mask]
    ax.scatter(DNA_xplot,RNA_yplot, color='forestgreen', alpha=0.95, s=5, lw=0.01, edgecolor='k')
    r2 = np.corrcoef(DNA_xplot, RNA_yplot)[0, 1]**2
    r2_log = np.corrcoef(np.log(DNA_xplot), np.log(RNA_yplot))[0, 1]**2
    text = f'$R^2={r2:.2f}$\n$R^2_{{\log}}={r2_log:.2f}$'
    ax.text(0.93, 0.07, text, transform=ax.transAxes, fontsize=8, verticalalignment='bottom', horizontalalignment='right', bbox=dict(facecolor='white', edgecolor='black', boxstyle='square, pad=1.'))

    ax.set_xlabel("DNA")
    ax.set_ylabel("RNA")
    try:
        max_val = max(DNA_xplot.max(), RNA_yplot.max())
        ax.set_xlim(1e-2, 1.05*max_val)
        ax.set_ylim(1e-2, 1.05*max_val)
    except ValueError:
        print('No data to plot')
        pass
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.grid(True)

    figname = f'figures/scatter_r0_{r0:.4f}_r1_{r1:.4f}_gamma1_{gamma1:.4f}_b0_{b0:.4f}_b1_{b1:.4f}_d0_{d0:.4f}_d1_{d1:.4f}_Delta_{Delta:.4f}_chain_{chain}.png'
    fig.savefig(figname, dpi=300)
    plt.close()
    
def expected_mean_without_r0(b0, b1, d0, d1):
    P0 = 1 / hyp2f1(1, b0 / b1, 1 + d0 / d1, b1 / d1)
    return (b0-d0*(1-P0))/(d1-b1)

def run_simulation(x0, y0, r0, r1, gamma1, b0, b1, d0, d1, Delta, T, dt, chain=0):
    logging.info(f'Running simulations with parameters: {x0:.4f}, {y0:.4f}, {r0:.4f}, {r1:.4f}, {gamma1:.4f}, {b0:.4f}, {b1:.4f}, {d0:.4f}, {d1:.4f}, {Delta:.4f}, {T}, {dt}')
    x, y, W20el, W02el, W11el = stochastic_simulation(x0, y0, r0, r1, gamma1, b0, b1, d0, d1, Delta, T, dt)
    logging.info(f'Simulation completed with parameters: {x0:.4f}, {y0:.4f}, {r0:.4f}, {r1:.4f}, {gamma1:.4f}, {b0:.4f}, {b1:.4f}, {d0:.4f}, {d1:.4f}, {Delta:.4f}, {T}, {dt}')
    
    store = False
    if store:
        filename = f'W_elements/{r0:1.4f}_{r1:1.4f}_{gamma1:1.4f}_{b0:1.4f}_{b1:1.4f}_{d0:1.4f}_{d1:1.4f}_{Delta:1.4f}_{T}_{dt:1.4f}_{chain}.txt'
        os.makedirs('W_elements', exist_ok=True)
        with open(filename, 'w') as f:
            f.write(f'W20el: {W20el}\n')
            f.write(f'W02el: {W02el}\n')
            f.write(f'W11el: {W11el}\n')
        
    y_max_nonnan = np.max(y[~np.isnan(y)])
    y_last = y[np.where(y == y_max_nonnan)[0][0]-50:np.where(y == y_max_nonnan)[0][0]]
    logging.info(f'Last 50 values of y before reaching maximum: {y_last}')   
     
    x_tau_steps, _ = find_uncorrelated_interval(x, dt)
    y_tau_steps, _ = find_uncorrelated_interval(y, dt)
    
    logging.info(f'Autocorrelation intervals determined: {x_tau_steps}, {y_tau_steps}')
    
    try:
        max_tau_steps = max(x_tau_steps, y_tau_steps)
    except ValueError:
        max_tau_steps = x_tau_steps
        
    logging.info(f'Uncorrelated interval determined: {max_tau_steps}')
    
    x_uncorrelated = x[::max_tau_steps]
    y_uncorrelated = y[::max_tau_steps]
    
    store = True
    if store:
        filename = f'mean_var/{r0:1.4f}_{r1:1.4f}_{gamma1:1.4f}_{b0:1.4f}_{b1:1.4f}_{d0:1.4f}_{d1:1.4f}_{Delta:1.4f}_{T}_{dt:1.4f}_{chain}.txt'
        os.makedirs('mean_var', exist_ok=True)
        with open(filename, 'w') as f:
            f.write(f'mean_x: {np.mean(x_uncorrelated):1.4f}\n')
            f.write(f'var_x: {np.var(x_uncorrelated):1.4f}\n')
            f.write(f'mean_y: {np.mean(y_uncorrelated):1.4f}\n')
            f.write(f'var_y: {np.var(y_uncorrelated):1.4f}\n')
            f.write(f'expected_mean_x: {expected_mean_without_r0(b0, b1, d0, d1):1.4f}\n')
        
    n_x, PDF_x, _, _ = log_pdf(x_uncorrelated, nbins=25)
    popt_x = fit_log_generalized_gamma(n_x, PDF_x)
    
    n_y, PDF_y, _, _ = log_pdf(y_uncorrelated, nbins=25)
    popt_y = fit_log_generalized_gamma(n_y, PDF_y)
    
    plotSAD(r0, r1, gamma1, b0, b1, d0, d1, Delta, n_x, PDF_x, n_y, PDF_y, popt_x, popt_y, chain)
    plotCORR(r0, r1, gamma1, b0, b1, d0, d1, Delta, x_uncorrelated, y_uncorrelated, chain)
    
    new_data = {
        'r0': r0, 'r1': r1, 'gamma1': gamma1, 
        'b0': b0, 'b1': b1, 'd0': d0, 
        'd1': d1, 'Delta': Delta, 'T': T, 'dt': dt, 'chain': chain, 'max_tau_steps': max_tau_steps,
        'popt_x_1': popt_x[0], 'popt_x_2': popt_x[1], 'popt_x_3': popt_x[2],
        'popt_y_1': popt_y[0], 'popt_y_2': popt_y[1], 'popt_y_3': popt_y[2]
    }
    
    filename = 'fits/gamma.csv'
    if os.path.exists(filename):
        df = pd.read_csv(filename)
    else:
        columns = ['r0', 'r1', 'gamma1', 'b0', 'b1', 'd0', 'd1', 'Delta', 'T', 'dt', 'chain', 'max_tau_steps',
                   'popt_x_1', 'popt_x_2', 'popt_x_3', 
                   'popt_y_1', 'popt_y_2', 'popt_y_3']
        df = pd.DataFrame(columns=columns)
        
    new_data_df = pd.DataFrame([new_data])
    df = pd.concat([df, new_data_df], ignore_index=True)
    df.to_csv(filename, index=False)

    filename = f'data/{r0:1.4f}_{r1:1.4f}_{gamma1:1.4f}_{b0:1.4f}_{b1:1.4f}_{d0:1.4f}_{d1:1.4f}_{Delta:1.4f}_{T}_{dt:1.4f}_{chain}.npz'
    format_scientific = np.vectorize(lambda x: float(np.format_float_scientific(x, precision=3, exp_digits=2)))
    x_uncorrelated_sc = format_scientific(x_uncorrelated)
    y_uncorrelated_sc = format_scientific(y_uncorrelated)
    np.savez(filename, x_uncorrelated=x_uncorrelated_sc, y_uncorrelated=y_uncorrelated_sc, max_tau_steps=max_tau_steps, x_tau_steps=x_tau_steps, y_tau_steps=y_tau_steps)
    
    return x_uncorrelated, y_uncorrelated, max_tau_steps, x_tau_steps, y_tau_steps

def setup_logging(stdout_file='stdout.log', stderr_file='stderr.log'):
    sys.stdout = open(stdout_file, 'w')
    sys.stderr = open(stderr_file, 'w')
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler("execution.log"),
            logging.StreamHandler(sys.stdout)
        ]
    )

def main():
    x0, y0 = 1., 1.
    Delta = 0.
    T = int(2e7)
    dt = 0.01
    number_of_chains = 200
    
    r0s = [500]
    r1s = [0.0]
    gamma1s = [50.01]
    b0s = [4.5]
    b1s = [1.]
    d0s = [5.5]
    d1s = [1.0001]
    print(Delta)
    
    parameter_sets = []
    for r0 in r0s:
        for r1 in r1s:
            for gamma1 in gamma1s:
                for b0 in b0s:
                    for b1 in b1s:
                        for d0 in d0s:
                            for d1 in d1s:
                                for chain in range(number_of_chains):
                                    parameter_sets.append((x0, y0, r0, r1, gamma1, b0, b1, d0, d1, Delta, T, dt, chain))
    print(parameter_sets)

    num_workers = 24
    logging.info(f'Starting simulations with {num_workers} workers.')
    with ProcessPoolExecutor(max_workers=num_workers) as executor:
        futures = [executor.submit(run_simulation, *params) for params in parameter_sets[:]]
        results = [future.result() for future in futures]
    logging.info('Simulations completed.')

if __name__ == "__main__":
    os.makedirs('logs', exist_ok=True)
    os.makedirs('figures', exist_ok=True)
    os.makedirs('data', exist_ok=True)
    os.makedirs('fits', exist_ok=True)
    setup_logging(stdout_file='logs/stdout.log', stderr_file='logs/stderr.log')
    logging.info('='*50)
    logging.info('Starting script execution.')

    import time
    start_time = time.time()
    
    light_palette = ['#4479C4', '#BD4135', '#948d99']
    dark_palette = ['#285DB1', '#AC3127', '#c1bbb0']
    color_clusters=['#648fff', '#ffb000', '#948d99']

    fig_formats = ['.pdf', '.eps', '.tiff']
    golden = (1 + 5 ** 0.5) / 2
    cms = 0.393701
    figsize = (golden*8*cms,8*cms)

    figsizes = {
        "1 column": (3.43, 3.43/golden),
        "1.5 columns": (4.49, 4.49/golden),
        "2 columns": (7.01, 7.01/golden)
    }

    sns.set_theme(
        rc={
            'figure.figsize':figsizes["2 columns"],
            'figure.dpi': 200,
            'savefig.dpi': 300},
        font="Helvetica Neue",
        font_scale=1.3,
        style="ticks")

    plt.rcParams['legend.edgecolor'] = 'k'
    plt.rcParams['legend.facecolor'] = 'w'
    plt.rcParams['legend.frameon'] = True
    plt.rcParams['legend.framealpha'] = 1
    plt.rcParams['legend.fancybox'] = False
    plt.rcParams['legend.fontsize'] = 10
    plt.rcParams['axes.linewidth'] = 1.5  
    plt.rcParams['axes.edgecolor'] = 'k'
    plt.rcParams['xtick.labelsize'] = 12
    plt.rcParams['ytick.labelsize'] = 12
    plt.rcParams['axes.labelsize'] = 12

    try:
        main()
    except Exception as e:
        logging.error(f'An error occurred: {str(e)}')
    finally:
        end_time = time.time()
        elapsed_time = end_time - start_time
        logging.info(f"Total execution time: {elapsed_time:.2f} seconds\n\n")
        sys.stdout.close()
        sys.stderr.close()
