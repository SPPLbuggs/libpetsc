import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.gridspec as gridspec

size = 12
med_size = 13
big_size = 13

plt.rc('font', size=size)
plt.rc('axes', titlesize=size)
plt.rc('axes', labelsize=med_size)
plt.rc('xtick', labelsize=size)
plt.rc('ytick', labelsize=size)
plt.rc('legend', fontsize=size)
plt.rc('figure', titlesize=big_size)
plt.rcParams['figure.figsize'] = (4.5, 3)

# Uncomment for LaTex fonts:
#plt.rcParams['font.family'] = 'serif'
#plt.rcParams['font.serif'] = ['Computer Modern']
#mpl.rc('text', usetex=True)

cm_subsection = np.linspace(0.0, 1.0, 4) 
colors = [ mpl.cm.viridis(x) for x in cm_subsection ]

path = 'Output'

x = np.fromfile(path + '/meshx.dat',dtype=float)
y = np.fromfile(path + '/meshy.dat',dtype=float)
t = np.fromfile(path + '/time.dat', dtype=float)

nx = len(x)
ny = max(len(y),1)
nt = len(t)

temp = np.fromfile(path + '/f1.dat',dtype=float)
f1 = temp.reshape([nt, ny, nx])

temp = np.fromfile(path + '/f2.dat',dtype=float)
f2 = temp.reshape([nt, ny, nx])

tt, xx = np.meshgrid(t,x)
yloc = 0

tloc = [nt/10, nt/2, -1]
tloc[0] = np.argmin(np.abs(t-0.2))


# Figure 1
fig = plt.figure(figsize=(5.25,5.25))
gs = gridspec.GridSpec(2,2)
gs.set_width_ratios([0.97,0.03])
ax0 = fig.add_subplot(gs[0,0])
ax1 = fig.add_subplot(gs[1,0], sharex=ax0)
ax2 = fig.add_subplot(gs[1,1])

plt.setp(ax0.get_xticklabels(), visible=False)
ax0.spines['right'].set_visible(False)
ax0.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
ax1.spines['top'].set_visible(False)

ax0.set_ylabel('Amplitude')
ax1.set_xlabel('Position')
ax1.set_ylabel('Time')

ax0.plot(x, f1[tloc[0],yloc,:], label='{:.1f}'.format(t[tloc[0]]), color=colors[0])
ax0.plot(x, f1[tloc[1],yloc,:], label='{:.1f}'.format(t[tloc[1]]), color=colors[1])
ax0.plot(x, f1[tloc[2],yloc,:], label='{:.1f}'.format(t[tloc[2]]), color=colors[2])
ax0.legend(bbox_to_anchor = (1.325, 0.55), loc = 5, frameon=False, title='Time')

im = ax1.contourf(xx, tt, f1[:,yloc,:].T, 30, cmap='viridis')
fig.colorbar(im, cax = ax2)

ax0.set_ylim([-0.05, 1.05])
ax1.set_xlim([-0.05, 1.05])
ax1.set_ylim([-0.5, 10.5])

gs.tight_layout(fig, rect=[0, 0, 1, 1])

# Figure 2
fig2 = plt.figure(figsize=(5.25,5.25))
gs2 = gridspec.GridSpec(2,2)
gs2.set_width_ratios([0.97,0.03])
ax0 = fig2.add_subplot(gs2[0,0])
ax1 = fig2.add_subplot(gs2[1,0], sharex=ax0)
ax2 = fig2.add_subplot(gs2[1,1])

plt.setp(ax0.get_xticklabels(), visible=False)
ax0.spines['right'].set_visible(False)
ax0.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
ax1.spines['top'].set_visible(False)

ax0.set_ylabel('Amplitude')
ax1.set_xlabel('Position')
ax1.set_ylabel('Time')

ax0.plot(x, f2[tloc[0],yloc,:], label='{:.1f}'.format(t[tloc[0]]), color=colors[0])
ax0.plot(x, f2[tloc[1],yloc,:], label='{:.1f}'.format(t[tloc[1]]), color=colors[1])
ax0.plot(x, f2[tloc[2],yloc,:], label='{:.1f}'.format(t[tloc[2]]), color=colors[2])
ax0.legend(bbox_to_anchor = (1.325, 0.55), loc = 5, frameon=False, title='Time')

im = ax1.contourf(xx, tt, f2[:,yloc,:].T, 30, cmap='viridis')
fig2.colorbar(im, cax = ax2)

ax0.set_ylim([-1.05, 1.05])
ax1.set_xlim([-0.05, 1.05])
ax1.set_ylim([-0.5, 10.5])

gs2.tight_layout(fig2, rect=[0, 0, 1, 1])

plt.show()
