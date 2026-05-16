"""Step 1 (FINAL): Data, SIR demo, prior predictive
Model: cases ~ NegBin(rho * I(t), phi)  with N=11750, I0=1
Params: R0, gamma, rho, phi
"""
import numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt, os

# Get the directory of this script for relative paths
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_DIR = os.path.dirname(SCRIPT_DIR)

np.random.seed(42)
PLOT_DIR = os.path.join(PROJECT_DIR, "outputs", "plots")
DATA_PATH = os.path.join(PROJECT_DIR, "data", "clean_epidemic_dataset_2025.csv")
os.makedirs(PLOT_DIR, exist_ok=True)
plt.rcParams.update({'font.family':'serif','font.size':11,'axes.grid':True,'grid.alpha':0.3,'grid.color':'gray','figure.facecolor':'white'})

df = pd.read_csv(DATA_PATH)
weeks = df['week_index'].values; cases = df['cases'].values; n_weeks = len(cases)
N = 11750; I0 = 1

def solve_sir(beta, gamma, N, I0, n_weeks, dt=0.1):
    spw = int(1.0/dt)
    S,I = float(N-I0),float(I0)
    S_o,I_o,R_o = np.zeros(n_weeks+1),np.zeros(n_weeks+1),np.zeros(n_weeks+1)
    S_o[0],I_o[0] = S,I
    for w in range(n_weeks):
        for _ in range(spw):
            dS=-beta*S*I/N; dI=beta*S*I/N-gamma*I
            S2,I2=S+dt/2*dS,I+dt/2*dI; dS2=-beta*S2*I2/N; dI2=beta*S2*I2/N-gamma*I2
            S3,I3=S+dt/2*dS2,I+dt/2*dI2; dS3=-beta*S3*I3/N; dI3=beta*S3*I3/N-gamma*I3
            S4,I4=S+dt*dS3,I+dt*dI3; dS4=-beta*S4*I4/N; dI4=beta*S4*I4/N-gamma*I4
            S+=dt/6*(dS+2*dS2+2*dS3+dS4); I+=dt/6*(dI+2*dI2+2*dI3+dI4)
            S,I=max(S,0),max(I,0)
        S_o[w+1],I_o[w+1]=S,I; R_o[w+1]=N-S-I
    return S_o,I_o,R_o

# Plot 1: Raw data
fig,ax = plt.subplots(figsize=(10,5))
ax.scatter(weeks,cases,color='black',s=30,zorder=5)
ax.set_xlabel("Week Index"); ax.set_ylabel("Number of Cases")
ax.set_title("Weekly Reported Cases Over Time")
plt.tight_layout(); fig.savefig(f"{PLOT_DIR}/01_raw_data.png",dpi=200,bbox_inches='tight'); plt.close()

# Plot 2: SIR compartments with best-fit example
beta_ex,gamma_ex,rho_ex = 0.748, 0.44, 0.049
S_e,I_e,R_e = solve_sir(beta_ex,gamma_ex,N,I0,n_weeks)
fig,axes = plt.subplots(1,2,figsize=(12,5))
axes[0].plot(np.arange(n_weeks+1),S_e,color='black',lw=2,label='S')
axes[0].plot(np.arange(n_weeks+1),I_e,color='black',lw=2,ls='--',label='I')
axes[0].plot(np.arange(n_weeks+1),R_e,color='black',lw=2,ls=':',label='R')
axes[0].set_xlabel("Week"); axes[0].set_ylabel("Individuals"); axes[0].set_title("SIR Compartments (full scale)")
axes[0].legend(frameon=True,edgecolor='black')
axes[1].plot(weeks,rho_ex*I_e[1:],color='black',lw=2,label=f'rho*I(t), rho={rho_ex}')
axes[1].scatter(weeks,cases,color='black',s=30,zorder=5,label='Data')
axes[1].set_xlabel("Week"); axes[1].set_ylabel("Cases"); axes[1].set_title("Observed vs Predicted (rho*I(t))")
axes[1].legend(frameon=True,edgecolor='black')
plt.tight_layout(); fig.savefig(f"{PLOT_DIR}/02_sir_compartments.png",dpi=200,bbox_inches='tight'); plt.close()

# Prior predictive: R0~LN(log(1.7),0.25), gamma~LN(log(0.44),0.3), rho~LN(log(0.05),0.5)
prior_curves = []
for _ in range(200):
    r0 = np.random.lognormal(np.log(1.7), 0.25)
    g = np.random.lognormal(np.log(0.44), 0.3)
    rho = np.random.lognormal(np.log(0.05), 0.5)
    b = r0*g
    if b>5 or g>3 or rho>1: continue
    _,Ip,_ = solve_sir(b,g,N,I0,n_weeks)
    prior_curves.append(rho*Ip[1:])

fig,ax = plt.subplots(figsize=(10,5))
for c in prior_curves: ax.plot(weeks,c,color='gray',alpha=0.12,lw=0.8)
ax.scatter(weeks,cases,color='black',s=30,zorder=5,label='Observed data')
ax.set_xlabel("Week Index"); ax.set_ylabel("Cases")
ax.set_title("Prior Predictive Check"); ax.legend(frameon=True,edgecolor='black'); ax.set_ylim(bottom=0)
plt.tight_layout(); fig.savefig(f"{PLOT_DIR}/03_prior_predictive.png",dpi=200,bbox_inches='tight'); plt.close()

# Check prior predictive coverage
pa = np.array(prior_curves)
cov = sum(1 for t in range(n_weeks) if np.percentile(pa[:,t],5)<=cases[t]<=np.percentile(pa[:,t],95))
print(f'Prior predictive coverage: {cov/n_weeks*100:.0f}% of data in 90% interval => {"PASS" if cov/n_weeks>0.4 else "FAIL"}')
np.savez(f"{PLOT_DIR}/../step1_data.npz", weeks=weeks, cases=cases, n_weeks=n_weeks, N=N, I0=I0)
print("Step 1 done.")
