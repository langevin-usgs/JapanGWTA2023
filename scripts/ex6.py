import pathlib as pl
import matplotlib.pyplot as plt
import matplotlib.animation
import numpy as np
from shapely.geometry import Polygon, LineString, Point
import flopy
import flopy.utils.triangle
import flopy.utils.voronoi

ws = './ex6'
name = 'mymodel'
tempdir = pl.Path(ws, "temp")
tempdir.mkdir(parents=True, exist_ok=True)

domain = Polygon([
    [1831.38, 6335.54],
    [4337.73, 6851.13],
    [6428.74, 6707.91],
    [8662.98, 6493.08],
    [9350.43, 5891.56],
    [9235.86, 4717.15],
    [8963.74, 3685.97],
    [8691.62, 2783.68],
    [8047.13, 2038.94],
    [7416.96, 578.09],
    [6414.42, 105.46],
    [5354.59, 205.72],
    [4624.17, 363.26],
    [3363.83, 563.77],
    [1330.11, 1809.78],
    [399.18, 2998.51],
    [914.77, 5132.49],
])
area_max = 1 * 100.0**2
tri = flopy.utils.triangle.Triangle(maximum_area=area_max, angle=30, model_ws=tempdir)
tri.add_polygon(domain)
tri.build(verbose=False)
vor = flopy.utils.voronoi.VoronoiGrid(tri)
disv_gridprops = vor.get_disv_gridprops()

top = 50.
botm = [-100.]
nlay = 1

# create and run steady-state flow simulation
ws_flow = pl.Path(ws, "flow")
sim = flopy.mf6.MFSimulation(sim_name=name, sim_ws=ws_flow, exe_name='mf6')
tdis = flopy.mf6.ModflowTdis(sim)
ims = flopy.mf6.ModflowIms(sim, print_option="all", inner_maximum=100)
gwf = flopy.mf6.ModflowGwf(sim, modelname=name, save_flows=True)
disv = flopy.mf6.ModflowGwfdisv(
    gwf, nlay=nlay, **vor.get_disv_gridprops(), top=top, botm=botm)
ic = flopy.mf6.ModflowGwfic(gwf, strt=top)
npf = flopy.mf6.ModflowGwfnpf(
    gwf, save_specific_discharge=True, save_saturation=True, icelltype=1,
)

gi = flopy.utils.GridIntersect(gwf.modelgrid)
ls = LineString([p for p in domain.exterior.coords])
chdspd = [[(0, j), 1.] for j in np.array(gi.intersects(ls)["cellids"], dtype=int)]
chd = flopy.mf6.ModflowGwfchd(gwf, stress_period_data=chdspd)
rch = flopy.mf6.ModflowGwfrcha(gwf, recharge=0.001)
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

fig, axes = plt.subplots(2, 1, figsize=(8, 11))
ax = axes[0]
ax.set_title("Map View")
ax.set_aspect(1)
ax.set_xlabel("x")
ax.set_ylabel("y")
pmv = flopy.plot.PlotMapView(gwf, ax=ax)
pmv.plot_bc(ftype="CHD")
pmv.plot_grid(linewidth=0.5)
#pmv.plot_array(head)
pmv.plot_vector(qx, qy, normalize=True, color="black", istep=1, jstep=1)
pmv.contour_array(head)
ax.plot(*domain.exterior.xy, color="black", linewidth=3.)

ax = axes[1]
ax.set_title("Cross Section")
ax.set_aspect(10.)
ax.set_xlabel("x")
ax.set_ylabel("z")
pxs = flopy.plot.PlotCrossSection(gwf, ax=ax, line={"line": [(0, 4500), (9500, 4500)]})
pxs.plot_bc(ftype="CHD")
pxs.plot_array(head, head=head)
pxs.plot_grid(color="black", linewidth=0.5)
plt.tight_layout()

# create and run transient transport simulation using flow results
ws_transport = pl.Path(ws, "transport")
sim = flopy.mf6.MFSimulation(sim_name=name, sim_ws=ws_transport, exe_name='mf6')
tdis = flopy.mf6.ModflowTdis(sim, nper=1, perioddata=[(36500., 100, 1.0)])
ims = flopy.mf6.ModflowIms(sim, print_option="all", linear_acceleration="BICGSTAB", inner_maximum=100)
gwt = flopy.mf6.ModflowGwt(sim, modelname=name, save_flows=True)
disv = flopy.mf6.ModflowGwtdisv(
    gwt, nlay=nlay, **vor.get_disv_gridprops(), top=top, botm=botm)
ic = flopy.mf6.ModflowGwtic(gwt, strt=0.)
mst = flopy.mf6.ModflowGwtmst(gwt, porosity=0.25)
adv = flopy.mf6.ModflowGwtadv(gwt, scheme="upstream")
dsp = flopy.mf6.ModflowGwtdsp(gwt, alh=1., ath1=0.1)
sourcerecarray = [()]
ssm = flopy.mf6.ModflowGwtssm(gwt, sources=sourcerecarray)
pd = [
    ("GWFHEAD", f"../flow/{name}.hds"),
    ("GWFBUDGET", f"../flow/{name}.bud"),
]
fmi = flopy.mf6.ModflowGwtfmi(gwt, packagedata=pd)
gi = flopy.utils.GridIntersect(gwt.modelgrid)
ix = gi.intersect(Point(6000., 4000.))
cnc_node = ix["cellids"][0]
cnc = flopy.mf6.ModflowGwtcnc(gwt, stress_period_data=[((0, cnc_node), 100.),])
oc = flopy.mf6.ModflowGwtoc(
    gwt,
    budget_filerecord=f"{name}.cbc",
    concentration_filerecord=f"{name}.ucn",
    saverecord=[("CONCENTRATION", "ALL"), ("BUDGET", "LAST")],
    printrecord=[("CONCENTRATION", "LAST"), ("BUDGET", "LAST")],
)
sim.write_simulation()
sim.run_simulation()

cobj = gwt.output.concentration()
conc = cobj.get_data()
conc_alldata = cobj.get_alldata()
times = cobj.times
budgwt = gwt.output.budget()

levels = [0.001, 0.01, 0.1, 1., 10., 100.]
cmap = plt.colormaps["jet"].with_extremes(under="white", over="white")
ca_dict = {
    "cmap": cmap, 
    "extend": "both", 
    "levels":levels, 
    "filled":True,
    "norm": "log",
}

fig, ax = plt.subplots(figsize=(8, 8))
ax.set_aspect(1)
pmv = flopy.plot.PlotMapView(gwt, ax=ax)
pmv.plot_grid(linewidth=0.5)
pmv.contour_array(head)
pmv.contour_array(conc, **ca_dict)
ax.plot(*domain.exterior.xy, color="black", linewidth=3.)

fig, ax = plt.subplots(figsize=(8, 6))
ax.set_aspect(1)
ax.set_xlabel(r'x')
ax.set_ylabel(r'y')
title = ax.set_title(f"Time = {times[0]} days")

# plot persistent items
pmv = flopy.plot.PlotMapView(gwt, ax=ax)
pmv.plot_grid(linewidth=0.1)
pmv.contour_array(head)
levels = [0.001, 0.01, 0.1, 1., 10., 100.]
cmap = plt.colormaps["jet"].with_extremes(under="white", over="white")
ca_dict = {
    "cmap": cmap, 
    "extend": "both", 
    "levels":levels, 
    "filled":True,
    "norm": "log",
}
cont = pmv.contour_array(conc_alldata[0], **ca_dict)
fig.colorbar(cont, shrink=0.5)
ax.plot(*domain.exterior.xy, color="black", linewidth=3.)

def animate(i):
    global cont
    global title
    for c in cont.collections:
        c.remove()  # removes only the contours, leaves the rest intact
    cont = pmv.contour_array(conc_alldata[i], **ca_dict)
    title = ax.set_title(f"Time = {times[i]} days")
    return cont

ani = matplotlib.animation.FuncAnimation(fig, animate, frames=conc_alldata.shape[0])
plt.close()

# in notebook use this method
#from IPython.display import HTML
#HTML(ani.to_jshtml())

# can use this command to write animation to file (use .avi on windows)
print ("saving animation...")
ani.save("ex6-animation.mp4")
print ("done...")
