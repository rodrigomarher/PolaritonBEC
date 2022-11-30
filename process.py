import numpy as np
import matplotlib.pyplot as plt
import imageio as iio
import glob
import gc
from joblib import Parallel, delayed
import sys

njobs=16
path = "results/"
wfc_names = []
wfx_names = []
pump_names = []

for i,name in enumerate(glob.iglob(path+"pump*")):
    pump_names.append(name)
for i,name in enumerate(glob.iglob(path+"wfc*")):
    wfc_names.append(name)
for i,name in enumerate(glob.iglob(path+"wfx*")):
    wfx_names.append(name)

pump_names.sort()
wfc_names.sort()
wfx_names.sort()
nimages = len(wfx_names)

x = np.linspace(-60,60,500)
y = np.linspace(-60,60,500)
t = np.linspace(0,20,500)

pump = np.zeros(nimages)
photon = np.zeros(nimages)
exciton = np.zeros(nimages)
for i in range(nimages):
    field_pump = np.loadtxt(pump_names[i],dtype=np.complex128)
    field_wfc = np.loadtxt(wfc_names[i],dtype=np.complex128)
    field_wfx = np.loadtxt(wfx_names[i],dtype=np.complex128)
    pump[i] = np.max(np.abs(field_pump))
    photon[i] = np.max(np.abs(field_wfc))
    exciton[i] = np.max(np.abs(field_wfx))

max_p = pump.max()
max_c = photon.max()
max_x = exciton.max()

def plot(ti):
    wfx = np.genfromtxt(wfx_names[ti],dtype=np.complex128)
    wfc = np.genfromtxt(wfc_names[ti],dtype=np.complex128)
    pump = np.genfromtxt(pump_names[ti],dtype=np.complex128)

    fig = plt.figure(figsize=(18,8))
    spec = fig.add_gridspec(nrows=2, ncols=9,width_ratios=[1,0.05,0.5,1,0.05,0.5,1,0.05,0.5],wspace=0.2,hspace=0.3)
    ax_fig = []
    ax_cbar = []
    ax_empty = []
    cbar = []
    for i in range(2):
        for j in range(3):
            ax_fig.append(fig.add_subplot(spec[i,3*j]))
            ax_cbar.append(fig.add_subplot(spec[i,3*j+1]))
            ax_empty.append(fig.add_subplot(spec[i,3*j+2]))

    for i in range(6):
        ax_empty[i].axis('off')

    ax_fig[0].set_title('Pump field',fontsize=17)
    ax_fig[1].set_title("Photonic field",fontsize=17)
    ax_fig[2].set_title("Excitonic field",fontsize=17)
    for i in range(6):
        ax_fig[i].set_xlabel('x (μm)', fontsize=16)
        ax_fig[i].set_ylabel('y (μm)', fontsize=16)
        ax_fig[i].tick_params(axis='both', which='major', labelsize=16)

    im = ax_fig[0].imshow(np.abs(pump)**2,vmax=max_p**2,vmin=0,extent=[x[0],x[-1],y[0],y[-1]],cmap='turbo',aspect='auto',origin='lower')
    cbar.append(fig.colorbar(im,cax=ax_cbar[0]))
    im = ax_fig[1].imshow(np.abs(wfc)**2*1E2,vmax=max_c**2*1E2,vmin=0,extent=[x[0],x[-1],y[0],y[-1]],cmap='turbo',aspect='auto',origin='lower')
    cbar.append(fig.colorbar(im,cax=ax_cbar[1]))
    im = ax_fig[2].imshow(np.abs(wfx)**2*1E2,vmax=max_x**2*1E2,vmin=0,extent=[x[0],x[-1],y[0],y[-1]],cmap='turbo',aspect='auto',origin='lower')
    cbar.append(fig.colorbar(im,cax=ax_cbar[2]))
    im = ax_fig[3].imshow(np.angle(pump),vmax=np.pi,vmin=-np.pi,cmap='twilight',aspect='auto',origin='lower',extent=[x[0],x[-1],y[0],y[-1]])
    cbar.append(fig.colorbar(im,cax=ax_cbar[3]))
    im = ax_fig[4].imshow(np.angle(wfc),vmax=np.pi,vmin=-np.pi,cmap='twilight',aspect='auto',origin='lower',extent=[x[0],x[-1],y[0],y[-1]])
    cbar.append(fig.colorbar(im,cax=ax_cbar[4]))
    im = ax_fig[5].imshow(np.angle(wfx),vmax=np.pi,vmin=-np.pi,cmap='twilight',aspect='auto',origin='lower',extent=[x[0],x[-1],y[0],y[-1]])
    cbar.append(fig.colorbar(im,cax=ax_cbar[5]))
    fig.tight_layout()
    fig.suptitle(f"t = {round(t[ti],3)} ps", fontsize=17)
    plt.savefig(f"tmp/snapshot_{str(ti).zfill(8)}.png")
    plt.cla()
    plt.clf()
    plt.close('all')
    gc.collect()

results = Parallel(n_jobs=6)(delayed(plot)(ti) for ti in range(nimages))

filenames = [i for i in glob.iglob('tmp/*')]
filenames.sort()

if len(sys.argv)>1: 
    gifname = sys.argv[1]
else:
    gifname = "condensate.gif"

with iio.get_writer(gifname, mode='I') as writer:
    for filename in filenames:
        image = iio.imread(filename)
        writer.append_data(image)

