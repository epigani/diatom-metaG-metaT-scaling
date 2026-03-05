import numpy as np
import numba as nb
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.optimize as opt
from scipy.signal import fftconvolve
from pathlib import Path
from tqdm import tqdm


##############################################################
# PARAMETERS FROM SAD → RATES
##############################################################

def sad_to_rates(kG, zetaG, alphaG, sigma_e2):

    sigma_d2 = kG * sigma_e2
    mu = 0.5 * zetaG * sigma_d2
    d = 0.5 * (alphaG - 1) * sigma_e2

    d1 = (sigma_d2 + d)/2
    b1 = (sigma_d2 - d)/2

    b0 = max(mu,0)
    d0 = max(-mu,0)

    return dict(
        b0=b0,
        b1=b1,
        d0=d0,
        d1=d1,
        sigma_d2=sigma_d2,
        sigma_e2=sigma_e2
    )


##############################################################
# RNA PARAMETERS
##############################################################

def set_rna_rates(kG, kT, r=1000):

    gamma = kG * r / kT

    return dict(r=r, gamma=gamma)


##############################################################
# LANGEVIN SIMULATION
##############################################################

@nb.njit
def simulate_xy(x0,y0,
                b0,b1,d0,d1,
                sigma_d2,sigma_e2,
                r,gamma,
                dt,T):
    
    # print(f"Running simulation... Parameters: b0={b0:.4f}, b1={b1:.4f}, d0={d0:.4f}, d1={d1:.4f}, sigma_d2={sigma_d2:.4f}, sigma_e2={sigma_e2:.4f}, r={r:.4f}, gamma={gamma:.4f}")
    # print("This may take a while...")
    
    counter=0
    
    x=np.zeros(T)
    y=np.zeros(T)

    x[0]=x0
    y[0]=y0

    for t in range(1,T):

        mu_x = (b1-d1)*x[t-1] + (b0-d0)
        mu_y = r*x[t-1] - gamma*y[t-1]

        noise_x = np.sqrt(sigma_d2*x[t-1] + sigma_e2*x[t-1]**2)
        noise_y = np.sqrt(r*x[t-1] + gamma*y[t-1])

        dx = mu_x*dt + noise_x*np.random.normal()*np.sqrt(dt)
        dy = mu_y*dt + noise_y*np.random.normal()*np.sqrt(dt)
        
        new_x = x[t-1] + dx
        new_y = y[t-1] + dy
        
        if new_x < 0:
            # new_x = x[t-1]
            new_x = 0.0
            counter += 1
        if new_y < 0:
            # new_y = y[t-1]
            new_y = 0.0
            counter += 1

        x[t] = new_x
        y[t] = new_y

    return x,y,counter


##############################################################
# AUTOCORRELATION
##############################################################

def autocorr_fft(x):

    x=x-np.mean(x)

    if np.allclose(x, 0.0):
        return np.ones(1)

    result=fftconvolve(x,x[::-1],mode='full')
    result=result[result.size//2:]

    return result/result[0]


def decorrelation_time(x,threshold=0.01):

    ac=autocorr_fft(x)

    idx=np.where(ac<threshold)[0]

    if len(idx)==0:
        return len(ac)

    return idx[0]


def decorrelation_time_chunked(
    x,
    threshold=0.01,
    chunk_size=int(2e6),
    min_chunk_size=4096,
    aggregate_quantile=0.75,
    max_chunks=20,
    saturation_threshold=0.25,
    max_chunk_growth_steps=3,
    max_chunk_size=None,
    return_diagnostics=False
):

    n=len(x)
    chunk_size=max(int(chunk_size), 1)
    max_chunks=max(int(max_chunks), 1)
    min_chunk_size=max(int(min_chunk_size), 1)
    max_chunk_growth_steps=max(int(max_chunk_growth_steps), 0)
    if max_chunk_size is None:
        max_chunk_size=n
    else:
        max_chunk_size=max(int(max_chunk_size), 1)
    max_chunk_size=min(max_chunk_size, n)
    effective_chunk_size=min(chunk_size, max_chunk_size)

    if n <= effective_chunk_size:
        tau=max(decorrelation_time(x, threshold=threshold), 1)
        if not return_diagnostics:
            return tau
        diagnostics={
            "chunk_size_initial": chunk_size,
            "chunk_size_used": n,
            "growth_steps": 0,
            "n_chunks_used": 1,
            "saturation_fraction": float(tau >= n),
            "adaptation_exhausted": False
        }
        return tau, diagnostics

    growth_steps=0
    tau_estimate=1
    saturation_fraction=0.0
    n_chunks_used=0
    adaptation_exhausted=False

    while True:
        starts=np.arange(0, n, effective_chunk_size, dtype=np.int64)
        if len(starts) > max_chunks:
            select_idx=np.linspace(0, len(starts) - 1, max_chunks, dtype=np.int64)
            starts=starts[select_idx]

        taus=[]
        saturated_chunks=0
        for start in starts:
            stop=min(start + effective_chunk_size, n)
            chunk=x[start:stop]
            if len(chunk) < min_chunk_size:
                continue
            tau_chunk=max(decorrelation_time(chunk, threshold=threshold), 1)
            taus.append(tau_chunk)
            if tau_chunk >= len(chunk):
                saturated_chunks += 1

        n_chunks_used=len(taus)
        if n_chunks_used == 0:
            tau_estimate=1
            saturation_fraction=0.0
        else:
            tau_estimate=max(int(np.ceil(np.quantile(taus, aggregate_quantile))), 1)
            saturation_fraction=saturated_chunks / n_chunks_used

        can_grow=(
            growth_steps < max_chunk_growth_steps
            and effective_chunk_size < max_chunk_size
            and effective_chunk_size < n
        )
        if saturation_fraction > saturation_threshold and can_grow:
            next_chunk_size=min(effective_chunk_size * 2, max_chunk_size, n)
            if next_chunk_size == effective_chunk_size:
                break
            effective_chunk_size=next_chunk_size
            growth_steps += 1
            continue

        adaptation_exhausted=(saturation_fraction > saturation_threshold) and (not can_grow)
        break

    if not return_diagnostics:
        return tau_estimate

    diagnostics={
        "chunk_size_initial": chunk_size,
        "chunk_size_used": effective_chunk_size,
        "growth_steps": growth_steps,
        "n_chunks_used": n_chunks_used,
        "saturation_fraction": saturation_fraction,
        "adaptation_exhausted": adaptation_exhausted
    }
    return tau_estimate, diagnostics


##############################################################
# LOG BINNED PDF
##############################################################

def log_pdf(data,nbins=30):

    data=data[data>0]

    bins=np.logspace(np.log10(data.min()),np.log10(data.max()),nbins)

    hist,_=np.histogram(data,bins)

    prob=hist/np.sum(hist)

    pdf=prob/np.diff(bins)

    x=np.sqrt(bins[1:]*bins[:-1])

    mask=pdf>0

    return x[mask],pdf[mask]


##############################################################
# HEATMAP
##############################################################

def heatmap_xy(x,y,bins=60):

    mask=(x>0)&(y>0)

    X=np.log10(x[mask])
    Y=np.log10(y[mask])

    H,xedges,yedges=np.histogram2d(X,Y,bins=bins,density=True)

    return H.T,xedges,yedges


##############################################################
# OUTPUT
##############################################################

def format_scientific_array(values, precision=3):

    formatter=np.vectorize(
        lambda x: float(np.format_float_scientific(x, precision=precision, exp_digits=2))
    )

    return formatter(values)


def build_output_filename(params, chain=0, outdir=None):

    if outdir is None:
        outdir = Path(__file__).resolve().parent / "data" / "langevin_environmental"
    else:
        outdir = Path(outdir)

    outdir.mkdir(parents=True, exist_ok=True)

    ordered_keys=[
        "kG", "zetaG", "alphaG", "kT",
        "sigma_e2", "b1_target", "r",
        "x0", "y0", "T", "dt", "decor_thresh"
    ]

    pieces=[]
    for key in ordered_keys:
        value=params[key]
        if key == "T":
            pieces.append(f"{key}={int(value)}")
        elif key == "decor_thresh":
            pieces.append(f"{key}={value:.2e}")
        else:
            pieces.append(f"{key}={value:.4f}")

    pieces.append(f"chain={int(chain)}")

    return str(outdir / ("__".join(pieces) + ".npz"))


##############################################################
# MAIN
##############################################################

def main():

    ##########################################################
    # TARGET PARAMETERS
    ##########################################################

    x0=1.0
    y0=1.0
    kG=1.0
    zetaG=1.0
    alphaG=1.5
    kT=20.0
    sigma_e2=1.0
    b1_target=1.0
    r=500.0
    dt=1e-5
    T=int(1e9)
    decor_thresh=0.01
    decor_chunk_size=int(2e6)
    decor_aggregate_quantile=0.75
    decor_max_chunks=20
    decor_saturation_threshold=0.25
    decor_max_chunk_growth_steps=3
    chain=0

    ##########################################################
    # BUILD RATES
    ##########################################################

    rates=sad_to_rates(kG, zetaG, alphaG, sigma_e2)

    rna=set_rna_rates(kG, kT, r=r)

    print("Rates:")
    print(rates)
    print(rna)

    ##########################################################
    # SIMULATION
    ##########################################################

    x,y, counts_boundary=simulate_xy(
        x0,y0,
        rates["b0"],rates["b1"],
        rates["d0"],rates["d1"],
        rates["sigma_d2"],rates["sigma_e2"],
        rna["r"],rna["gamma"],
        dt=dt,
        T=T
    )
    print(f"Simulation completed. Boundary hits: {counts_boundary}")

    ##########################################################
    # DECORRELATION
    ##########################################################

    x_tau_steps,x_decor_diag=decorrelation_time_chunked(
        x,
        threshold=decor_thresh,
        chunk_size=decor_chunk_size,
        aggregate_quantile=decor_aggregate_quantile,
        max_chunks=decor_max_chunks,
        saturation_threshold=decor_saturation_threshold,
        max_chunk_growth_steps=decor_max_chunk_growth_steps,
        return_diagnostics=True
    )
    y_tau_steps,y_decor_diag=decorrelation_time_chunked(
        y,
        threshold=decor_thresh,
        chunk_size=decor_chunk_size,
        aggregate_quantile=decor_aggregate_quantile,
        max_chunks=decor_max_chunks,
        saturation_threshold=decor_saturation_threshold,
        max_chunk_growth_steps=decor_max_chunk_growth_steps,
        return_diagnostics=True
    )
    max_tau_steps=max(max(x_tau_steps, y_tau_steps), 1)

    x=x[::max_tau_steps]
    y=y[::max_tau_steps]

    print("decorrelation step:", max_tau_steps)
    print("decorrelation chunk size:", decor_chunk_size)
    print("decorrelation aggregation quantile:", decor_aggregate_quantile)
    print("decorrelation max chunks:", decor_max_chunks)
    print("decorrelation saturation threshold:", decor_saturation_threshold)
    print("decorrelation max chunk growth steps:", decor_max_chunk_growth_steps)
    print(
        "x decor diagnostics:",
        f"chunk_size_used={x_decor_diag['chunk_size_used']},",
        f"growth_steps={x_decor_diag['growth_steps']},",
        f"saturation_fraction={x_decor_diag['saturation_fraction']:.2f},",
        f"n_chunks_used={x_decor_diag['n_chunks_used']}"
    )
    print(
        "y decor diagnostics:",
        f"chunk_size_used={y_decor_diag['chunk_size_used']},",
        f"growth_steps={y_decor_diag['growth_steps']},",
        f"saturation_fraction={y_decor_diag['saturation_fraction']:.2f},",
        f"n_chunks_used={y_decor_diag['n_chunks_used']}"
    )
    if x_decor_diag["adaptation_exhausted"]:
        print(
            "WARNING: x decorrelation estimate may be censored."
            " Increase decor_chunk_size or decor_max_chunk_growth_steps."
        )
    if y_decor_diag["adaptation_exhausted"]:
        print(
            "WARNING: y decorrelation estimate may be censored."
            " Increase decor_chunk_size or decor_max_chunk_growth_steps."
        )

    params = {
        "x0": x0,
        "y0": y0,
        "kG": kG,
        "zetaG": zetaG,
        "alphaG": alphaG,
        "kT": kT,
        "sigma_e2": sigma_e2,
        "b1_target": b1_target,
        "r": r,
        "dt": dt,
        "T": T,
        "decor_thresh": decor_thresh,
        "decor_chunk_size": decor_chunk_size,
        "decor_aggregate_quantile": decor_aggregate_quantile,
        "decor_max_chunks": decor_max_chunks,
        "decor_saturation_threshold": decor_saturation_threshold,
        "decor_max_chunk_growth_steps": decor_max_chunk_growth_steps
    }

    filename = build_output_filename(params=params, chain=chain)
    x_uncorrelated_sc=format_scientific_array(x)
    y_uncorrelated_sc=format_scientific_array(y)

    np.savez(
        filename,
        x_uncorrelated=x_uncorrelated_sc,
        y_uncorrelated=y_uncorrelated_sc,
        max_tau_steps=max_tau_steps,
        x_tau_steps=x_tau_steps,
        y_tau_steps=y_tau_steps,
        x_decor_chunk_size_used=x_decor_diag["chunk_size_used"],
        y_decor_chunk_size_used=y_decor_diag["chunk_size_used"],
        x_decor_growth_steps=x_decor_diag["growth_steps"],
        y_decor_growth_steps=y_decor_diag["growth_steps"],
        x_decor_saturation_fraction=x_decor_diag["saturation_fraction"],
        y_decor_saturation_fraction=y_decor_diag["saturation_fraction"],
        x_decor_n_chunks_used=x_decor_diag["n_chunks_used"],
        y_decor_n_chunks_used=y_decor_diag["n_chunks_used"],
        x_decor_adaptation_exhausted=x_decor_diag["adaptation_exhausted"],
        y_decor_adaptation_exhausted=y_decor_diag["adaptation_exhausted"],
        chain=chain,
        gamma=rna["gamma"],
        b0=rates["b0"],
        b1=rates["b1"],
        d0=rates["d0"],
        d1=rates["d1"],
        sigma_d2=rates["sigma_d2"],
        sigma_e2_eff=rates["sigma_e2"],
        **params
    )
    print(f"Saved uncorrelated data: {filename}")

    ##########################################################
    # SAD
    ##########################################################

    nx,px=log_pdf(x)
    ny,py=log_pdf(y)

    ##########################################################
    # PLOTS
    ##########################################################

    fig,ax=plt.subplots(1,2,figsize=(10,4))

    ax[0].loglog(nx,px,'o')
    ax[0].set_title("DNA SAD")

    ax[1].loglog(ny,py,'o')
    ax[1].set_title("RNA SAD")

    plt.show()

    ##########################################################
    # CORRELATION
    ##########################################################

    fig,ax=plt.subplots()

    ax.scatter(x,y,s=5,alpha=0.3)

    ax.set_xscale("log")
    ax.set_yscale("log")

    ax.set_xlabel("DNA")
    ax.set_ylabel("RNA")

    plt.show()

    ##########################################################
    # HEATMAP
    ##########################################################

    H,xedges,yedges=heatmap_xy(x,y)

    fig,ax=plt.subplots()

    im=ax.imshow(
        H,
        origin="lower",
        extent=[xedges[0],xedges[-1],yedges[0],yedges[-1]],
        aspect='auto'
    )

    ax.set_xlabel("log DNA")
    ax.set_ylabel("log RNA")

    plt.colorbar(im)

    plt.show()


if __name__=="__main__":
    main()
