import matplotlib.pyplot as plt
import numpy as np
import flopy
ws = './ex1'
name = 'mymodel'

def get_z(z0, a, b, alpha, x):
    return z0 + x * np.tan(alpha) + a * np.sin(b * x / np.cos(alpha)) / np.cos(alpha)

Lx = 20000
Lz = 10000
nlay = 400
nrow = 1
ncol = 200
top = Lz
bot = 0
dx = Lx / ncol
dz = (top - bot) / nlay
botm = [top - (b + 1) * dz for b in range(nlay)]
a = 200.
alpha = np.arctan2(1000, Lx)
period = 5000.
b = 2 * np.pi / period
x = np.arange(dx / 2, Lx + dx / 2, dx)
z = get_z(Lz, a, b, alpha, x)
chdspd = [[(0, 0, j), z[j]] for j in range(ncol)]

sim = flopy.mf6.MFSimulation(sim_name=name, sim_ws=ws, exe_name='mf6')
tdis = flopy.mf6.ModflowTdis(sim)
ims = flopy.mf6.ModflowIms(sim, print_option="all", inner_maximum=100)
gwf = flopy.mf6.ModflowGwf(sim, modelname=name, save_flows=True)
dis = flopy.mf6.ModflowGwfdis(
    gwf, nlay=nlay, nrow=nrow, ncol=ncol, top=top, botm=botm, delr=dx)
ic = flopy.mf6.ModflowGwfic(gwf, strt=top)
npf = flopy.mf6.ModflowGwfnpf(gwf, save_specific_discharge=True)
chd = flopy.mf6.ModflowGwfchd(gwf, stress_period_data=chdspd)
oc = flopy.mf6.ModflowGwfoc(gwf,
                            budget_filerecord=f"{name}.bud",
                            head_filerecord=f"{name}.hds",
                            printrecord=[('HEAD', 'ALL'), ('BUDGET', 'ALL')],
                            saverecord=[('HEAD', 'ALL'), ('BUDGET', 'ALL')])
sim.write_simulation()
sim.run_simulation()

head = gwf.output.head().get_data()
bud = gwf.output.budget()
spdis = bud.get_data(text='DATA-SPDIS')[0]
qx, qy, qz = flopy.utils.postprocessing.get_specific_discharge(spdis, gwf)
u = qx.reshape((nlay, ncol))
phi = u[-1::-1].cumsum(axis=0)
phi = np.flipud(phi)

fig, ax = plt.subplots(figsize=(10, 10))
ax.set_aspect(1)
pxs = flopy.plot.PlotCrossSection(gwf, ax=ax, line={"row":0})
pxs.contour_array(
    head, 
    levels=np.arange(Lz, z.max(), 25), 
    linewidths=1., 
    colors="k",
    linestyles="dashed",
)
pxs.contour_array(phi, levels=np.linspace(phi.min(), phi.max(), 10))
#pxs.plot_vector(qx, qy, qz, normalize=True, color="black")
ax.set_xlim(0, 20000)
ax.set_ylim(0, 11000)
ax.plot(x, z)