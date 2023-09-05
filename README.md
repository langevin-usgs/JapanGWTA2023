![alt](images/header.png)

# JapanGWTA2023

This repository contains information for a lecture to the Japan Ground Water Technology Association.  The purpose of the lecture is to present recent advances for the MODFLOW 6 groundwater simulator.

![alt](images/grid.png)

## Software Installation

See software installation instructions [here](./software.md).

## Examples

Example 1 -- Toth Flow System [script](./scripts/ex1.py) and [notebook](./notebooks/ex1.ipynb)

## Examples using FloPy and MODFLOW 6
### Example 1 -- Toth Flow System

```
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
```

### Example 2 -- Regular Grid

```
import pathlib as pl
import matplotlib.pyplot as plt
import numpy as np
from shapely.geometry import Polygon, LineString
import flopy

ws = './ex2'
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

Lx = 10000
Ly = 8000
nlay = 1
nrow = 32
ncol = 40
delr = Lx / ncol * np.ones(ncol, dtype=float)
delc = Ly / nrow * np.ones(nrow, dtype=float)
top = 50 * np.ones((nrow, ncol), dtype=float)
botm = -100 * np.ones((nlay, nrow, ncol), dtype=float)
sg = flopy.discretization.StructuredGrid(
    nlay=nlay,
    nrow=nrow,
    ncol=ncol,
    delr=delr,
    delc=delc,
    top=top,
    botm=botm
)
idomain = np.zeros((nlay, nrow, ncol), dtype=int)
gi = flopy.utils.GridIntersect(sg)

# inside domain polygon
ixp = gi.intersect(domain)
for i, j in ixp["cellids"]:
    idomain[:, i, j] = 1

# identify cells that touch domain polygon
ls = LineString([p for p in domain.exterior.coords])    
ixl = gi.intersect(ls)
for i, j in ixl["cellids"]:
    idomain[:, i, j] = 2

sim = flopy.mf6.MFSimulation(sim_name=name, sim_ws=ws, exe_name='mf6')
tdis = flopy.mf6.ModflowTdis(sim)
ims = flopy.mf6.ModflowIms(sim, print_option="all", inner_maximum=100)
gwf = flopy.mf6.ModflowGwf(sim, modelname=name, save_flows=True)
dis = flopy.mf6.ModflowGwfdis(
    gwf, nlay=nlay, nrow=nrow, ncol=ncol, delr=delr, delc=delc, 
    top=top, botm=botm, idomain=idomain)
ic = flopy.mf6.ModflowGwfic(gwf, strt=top)
npf = flopy.mf6.ModflowGwfnpf(gwf, save_specific_discharge=True, icelltype=1)
chdspd = [[(0, i, j), 1.] for i, j in ixl["cellids"]]
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
pmv.plot_inactive(color_noflow="gray")
pmv.plot_grid(color="black", linewidth=0.5)
pmv.plot_vector(qx, qy, normalize=True, color="black", istep=1, jstep=1)
pmv.contour_array(head)

ax = axes[1]
ax.set_title("Cross Section")
ax.set_aspect(10.)
ax.set_xlabel("x")
ax.set_ylabel("z")
pxs = flopy.plot.PlotCrossSection(gwf, ax=ax, line={"row": int(nrow/2)})
pxs.plot_inactive(color_noflow="gray")
pxs.plot_bc(ftype="CHD")
pxs.plot_array(head, head=head)
pxs.plot_grid(color="black", linewidth=0.5)
plt.tight_layout()
```

### Example 3 -- Quadtree Grid

```
import pathlib as pl
import matplotlib.pyplot as plt
import numpy as np
from shapely.geometry import Polygon, LineString
import flopy
from flopy.utils.gridgen import Gridgen

ws = './ex3'
name = 'mymodel'
gridgen_ws = pl.Path(ws, "gridgen")
gridgen_ws.mkdir(parents=True, exist_ok=True)

domain = [[
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
    [1831.38, 6335.54],
]]

Lx = 10000
Ly = 8000
nlay = 1
nrow = 32
ncol = 40
delr = Lx / ncol * np.ones(ncol, dtype=float)
delc = Ly / nrow * np.ones(nrow, dtype=float)
top = 50 * np.ones((nrow, ncol), dtype=float)
botm = -100 * np.ones((nlay, nrow, ncol), dtype=float)
sg = flopy.discretization.StructuredGrid(
    nlay=nlay,
    nrow=nrow,
    ncol=ncol,
    delr=delr,
    delc=delc,
    top=top,
    botm=botm
)

g = Gridgen(sg, model_ws=gridgen_ws, surface_interpolation="interpolate")
g.add_active_domain([domain], range(nlay))
g.add_refinement_features(domain, "line", level=1, layers=range(nlay))
poly = [[[(6000, 4000), (7000, 4000), (7000, 5000), (6000, 5000), (6000, 4000)]]]
g.add_refinement_features(poly, "polygon", level=2, layers=range(nlay))
g.build()

sim = flopy.mf6.MFSimulation(sim_name=name, sim_ws=ws, exe_name='mf6')
tdis = flopy.mf6.ModflowTdis(sim)
ims = flopy.mf6.ModflowIms(sim, print_option="all", inner_maximum=100)
gwf = flopy.mf6.ModflowGwf(sim, modelname=name, save_flows=True)
gridprops_disv = g.get_gridprops_disv()
disv = flopy.mf6.ModflowGwfdisv(
    gwf, **gridprops_disv)
ic = flopy.mf6.ModflowGwfic(gwf, strt=gridprops_disv["top"])
npf = flopy.mf6.ModflowGwfnpf(gwf, save_specific_discharge=True, icelltype=1)

gi = flopy.utils.GridIntersect(gwf.modelgrid)
ls = LineString(domain[0])
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
pmv.plot_grid(color="black", linewidth=0.5)
#pmv.plot_array(head)
pmv.plot_vector(qx, qy, normalize=True, color="black", istep=1, jstep=1)
pmv.contour_array(head)

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
```

### Example 4 -- Voronoi Grid

```
import pathlib as pl
import matplotlib.pyplot as plt
import numpy as np
from shapely.geometry import Polygon, LineString
import flopy
import flopy.utils.triangle
import flopy.utils.voronoi

ws = './ex4'
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
area_max = 1. * 100.0**2
tri = flopy.utils.triangle.Triangle(maximum_area=area_max, angle=30, model_ws=tempdir)
tri.add_polygon(domain)
tri.build(verbose=False)
vor = flopy.utils.voronoi.VoronoiGrid(tri)
disv_gridprops = vor.get_disv_gridprops()

top = 50.
botm = [-100.]
nlay = 1
sim = flopy.mf6.MFSimulation(sim_name=name, sim_ws=ws, exe_name='mf6')
tdis = flopy.mf6.ModflowTdis(sim)
ims = flopy.mf6.ModflowIms(sim, print_option="all", inner_maximum=100)
gwf = flopy.mf6.ModflowGwf(sim, modelname=name, save_flows=True)
disv = flopy.mf6.ModflowGwfdisv(
    gwf, nlay=nlay, **vor.get_disv_gridprops(), top=top, botm=botm)
ic = flopy.mf6.ModflowGwfic(gwf, strt=top)
npf = flopy.mf6.ModflowGwfnpf(gwf, save_specific_discharge=True, icelltype=1)

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
pmv.plot_grid(color="black", linewidth=0.5)
#pmv.plot_array(head)
pmv.plot_vector(qx, qy, normalize=True, color="black", istep=1, jstep=1)
pmv.contour_array(head)

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
```

### Example 5 -- Local Grid Refinement

```
import pathlib as pl
import matplotlib.pyplot as plt
import numpy as np
from shapely.geometry import Polygon, LineString
import flopy
from flopy.utils.lgrutil import Lgr

ws = './ex5'
name = 'mymodel'

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

Lx = 10000
Ly = 8000
nlayp = 1
nrowp = 32
ncolp = 40
delrp = Lx / ncolp * np.ones(ncolp, dtype=float)
delcp = Ly / nrowp * np.ones(nrowp, dtype=float)
topp = 50 * np.ones((nrowp, ncolp), dtype=float)
botmp = -100 * np.ones((nlayp, nrowp, ncolp), dtype=float)
idomainp = np.ones((nlayp, nrowp, ncolp), dtype=int)
idomainp[:, 11:15, 23:27] = 0

lgr = Lgr(
    nlayp,
    nrowp,
    ncolp,
    delrp,
    delcp,
    topp,
    botmp,
    idomainp,
    ncpp=3,
    ncppl=[3],
    xllp=0.0,
    yllp=0.0,
)

sg = flopy.discretization.StructuredGrid(
    nlay=nlayp,
    nrow=nrowp,
    ncol=ncolp,
    delr=delrp,
    delc=delcp,
    top=topp,
    botm=botmp
)

idomain = np.zeros((nlayp, nrowp, ncolp), dtype=int)
gi = flopy.utils.GridIntersect(sg)

# inside domain polygon
ixp = gi.intersect(domain)
for i, j in ixp["cellids"]:
    idomain[:, i, j] = 1

# touching domain polygon
ls = LineString([p for p in domain.exterior.coords])    
ixl = gi.intersect(ls)
for i, j in ixl["cellids"]:
    idomain[:, i, j] = 2
    
idomain[np.where(idomainp == 0)] = 0

sim = flopy.mf6.MFSimulation(sim_name=name, sim_ws=ws, exe_name='mf6')
tdis = flopy.mf6.ModflowTdis(sim)
ims = flopy.mf6.ModflowIms(sim)

# parent model
name = "parent"
gwf = flopy.mf6.ModflowGwf(sim, modelname=name, save_flows=True)
dis = flopy.mf6.ModflowGwfdis(
    gwf, 
    nlay=nlayp, 
    nrow=nrowp, 
    ncol=ncolp, 
    delr=delrp, 
    delc=delcp,
    top=topp,
    botm=botmp,
    idomain=idomain,
)
ic = flopy.mf6.ModflowGwfic(gwf)
npf = flopy.mf6.ModflowGwfnpf(gwf, save_specific_discharge=True, icelltype=1)
chdspd = [[(0, i, j), 1.] for i, j in ixl["cellids"]]
chd = flopy.mf6.ModflowGwfchd(gwf, stress_period_data=chdspd)
rch = flopy.mf6.ModflowGwfrcha(gwf, recharge=0.001)
budget_file = name + '.bud'
head_file = name + '.hds'
oc = flopy.mf6.ModflowGwfoc(gwf,
                            budget_filerecord=budget_file,
                            head_filerecord=head_file,
                            saverecord=[('HEAD', 'ALL'), ('BUDGET', 'ALL')])
gwf_parent = gwf

# child model
name = "child"
gwf = flopy.mf6.ModflowGwf(sim, modelname=name, 
                           save_flows=True)
dis = flopy.mf6.ModflowGwfdis(
    gwf, 
    nlay=lgr.nlay, 
    nrow=lgr.nrow, 
    ncol=lgr.ncol, 
    delr=lgr.delr, 
    delc=lgr.delc,
    top=lgr.top,
    botm=lgr.botm,
    xorigin=lgr.xll,
    yorigin=lgr.yll,
)
ic = flopy.mf6.ModflowGwfic(gwf)
npf = flopy.mf6.ModflowGwfnpf(gwf, save_specific_discharge=True, icelltype=1)
rch = flopy.mf6.ModflowGwfrcha(gwf, recharge=0.001)
budget_file = name + '.bud'
head_file = name + '.hds'
oc = flopy.mf6.ModflowGwfoc(gwf,
                            budget_filerecord=budget_file,
                            head_filerecord=head_file,
                            printrecord=[('HEAD', 'ALL'), ('BUDGET', 'ALL')],
                            saverecord=[('HEAD', 'ALL'), ('BUDGET', 'ALL')],
)
gwf_child = gwf

# setup the exchange
exchangedata = lgr.get_exchange_data(
    angldegx=True, 
    cdist=True
)
flopy.mf6.ModflowGwfgwf(sim, 
                        exgmnamea="parent",
                        exgmnameb="child",
                        xt3d=True,
                        nexg=len(exchangedata),
                        auxiliary=["ANGLDEGX", "CDIST"],
                        exchangedata=exchangedata)

sim.write_simulation()
sim.run_simulation()

fig, axes = plt.subplots(2, 1, figsize=(8, 11))
ax = axes[0]
ax.set_title("Map View")
ax.set_aspect(1)
ax.set_xlabel("x")
ax.set_ylabel("y")

for igrid, gwf in enumerate([gwf_parent, gwf_child]):
    head = gwf.output.head().get_data()
    bud = gwf.output.budget()
    spdis = bud.get_data(text='DATA-SPDIS')[0]
    qx, qy, qz = flopy.utils.postprocessing.get_specific_discharge(spdis, gwf)

    pmv = flopy.plot.PlotMapView(gwf, ax=ax)
    if igrid == 0:
        # parent
        pmv.plot_inactive(color_noflow="gray")
        pmv.plot_bc(ftype="CHD")
        pmv.plot_vector(qx, qy, normalize=True, color="black")
    if igrid == 1:
        # child
        c = pmv.plot_array(head, masked_values=[1.e30])
    pmv.plot_grid(colors='black', linewidths=0.5)
    pmv.contour_array(head, colors="blue")
ax.set_xlim(0, 10000)
ax.set_ylim(0, 8000)

ax = axes[1]
ax.set_title("Cross Section")
ax.set_aspect(10.)
ax.set_xlabel("x")
ax.set_ylabel("z")
for igrid, gwf in enumerate([gwf_parent, gwf_child]):
    head = gwf.output.head().get_data()
    bud = gwf.output.budget()
    spdis = bud.get_data(text='DATA-SPDIS')[0]
    qx, qy, qz = flopy.utils.postprocessing.get_specific_discharge(spdis, gwf)
    pxs = flopy.plot.PlotCrossSection(gwf, ax=ax, line={"line": [(0, 5000), (10000, 5000)]}, geographic_coords=True)
    if igrid == 0:
        pxs.plot_inactive(color_noflow="gray")
    c = pxs.plot_array(head, masked_values=[1.e30], head=head)
    pxs.plot_grid(colors='black', linewidths=0.5)
ax.set_xlim(0, 10000)
ax.set_ylim(-100, 50)
plt.tight_layout()
```

### Example 6 -- Voronoi Flow and Transport

This problem is divided into several parts, which can go into separate jupyter notebook cells.

#### Create and Run the Flow Model
```
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
```

#### Create and Run the Transport Model

```
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
```

#### Post-Process the Results into an Animation

```
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

from IPython.display import HTML
HTML(ani.to_jshtml())

# can use this command to write animation to file
#ani.save("ex6-animation.avi")
```

### Example 7 -- Saltwater Intrusion Model


#### Create and Run the Saltwater Intrusion Model
```
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation
import flopy

ws = './ex7'
name = 'henry'
gwfname = 'gwf_' + name
gwtname = 'gwt_' + name

nlay = 20
ncol = 40
fx = 0.5
fz = 0.5
lx = 2
lz = 1
sealevel = 0.85
amplitude = 0.14
frequency = 4
wellfact = 0.25

def get_idomain(nlay, nrow, ncol, lx, lz, fx, fz):
    idomain = np.ones((nlay, nrow, ncol), dtype=int)
    x1 = fx * lx
    y1 = lz
    x2 = lx
    y2 = fz * lz
    slope = (y2 - y1) / (x2 - x1)
    b = y1 - slope * x1

    delr = lx / ncol
    delv = lz / nlay
    xcenters = np.linspace(delr / 2, lx - delr / 2, ncol)
    zcenters = np.linspace(lz - delv / 2, delv / 2, nlay)

    for k in range(nlay):
        zc = zcenters[k]
        for j in range(ncol):
            xc = xcenters[j]
            zedge = slope * xc + b
            if zc > zedge:
                idomain[k, 0, j] = 0

    kidm0, iidmn0, jidmn0 = np.where(idomain == 0)
    for k, j in zip(kidm0, jidmn0):
        if idomain[k, 0, j] == 0 and idomain[k, 0, j - 1] == 1:
            idomain[k, 0, j - 1] = 2

    for k, j in zip(kidm0, jidmn0):
        if idomain[k, 0, j] == 0 and idomain[k + 1, 0, j] == 1:
            idomain[k + 1, 0, j] = 3

    return idomain

def sinfunc(a, b, c, d, x):
    return a * np.sin(b * (x - c)) + d

nrow = 1
delr = lx / ncol
delc = 1.0
top = lz
delz = lz / nlay
botm = list(top - np.arange(delz, nlay * delz + delz, delz))

perlen = [0.25] + 1000 * [0.001]
nper = len(perlen)
nstp = [250] + 1000 * [1]
tsmult = 1.0
tdis_rc = []
for i in range(nper):
    tdis_rc.append((perlen[i], nstp[i], tsmult))

nouter, ninner = 200, 50
hclose, rclose, relax = 1e-7, 1e-5, 0.97

sim = flopy.mf6.MFSimulation(sim_name=name, sim_ws=ws, exe_name='mf6')
tdis = flopy.mf6.ModflowTdis(
    sim, time_units='DAYS', nper=nper, perioddata=tdis_rc
)

gwf = flopy.mf6.ModflowGwf(
    sim, 
    modelname=gwfname, 
    newtonoptions="NEWTON UNDER_RELAXATION",
)
idomain = get_idomain(nlay, nrow, ncol, lx, lz, fx=fx, fz=fz)
dis = flopy.mf6.ModflowGwfdis(
    gwf, 
    nlay=nlay, 
    nrow=nrow, 
    ncol=ncol,
    delr=delr, 
    delc=delc,
    top=top, 
    botm=botm,
    idomain=idomain,
)
ic = flopy.mf6.ModflowGwfic(gwf, strt=lz)
npf = flopy.mf6.ModflowGwfnpf(
    gwf, 
    xt3doptions=False,
    save_flows=True,
    save_specific_discharge=True,
    icelltype=1,
    k=864.
)
sto = flopy.mf6.ModflowGwfsto(
    gwf, 
    iconvert=1,
    ss=1.e-3, 
    sy=0.35,
    steady_state=[False],
    transient=[True]
)
buy = flopy.mf6.ModflowGwfbuy(
    gwf, 
    packagedata=[(0, 0.7, 0.0, f"{gwtname}", "none")], 
    denseref=1000.0, 
)

# drn and ghb
kidx, iidx, jidx = np.where(idomain > 1)
xcellcenters = gwf.modelgrid.xcellcenters
zcellcenters = gwf.modelgrid.zcellcenters
botm = dis.botm.get_data()
dt = 0.001
times = np.arange(dt, 1.0 + dt, dt)
sealevelts = [sealevel] + list(
    sinfunc(amplitude, frequency * 2 * np.pi, 0, sealevel, times)
)
ghbspd = {}
drnspd = {}
for kper in range(nper):
    if kper == 0:
        sl = sealevel
    else:
        sl = sealevelts[kper]
    ghblist = []
    drnlist = []
    for k, i, j in zip(kidx, iidx, jidx):
        zcell = zcellcenters[k, i, j]
        cond = 864.0 * (delz * delc) / (0.5 * delr)
        if zcell > sl:
            drnlist.append([(k, i, j), zcell, 864.0, 0.0])
        else:
            ghblist.append([(k, i, j), sl, 864.0, 35.0, 1024.5])
    if len(ghblist) > 0:
        ghbspd[kper] = ghblist
    if len(drnlist) > 0:
        drnspd[kper] = drnlist

# drn
drn1 = flopy.mf6.ModflowGwfdrn(
    gwf,
    stress_period_data=drnspd,
    print_input=True,
    print_flows=True,
    save_flows=False,
    pname="DRN-1",
    auxiliary="CONCENTRATION",
)

# ghb
ghb1 = flopy.mf6.ModflowGwfghb(
    gwf,
    stress_period_data=ghbspd,
    print_input=True,
    print_flows=True,
    save_flows=False,
    pname="GHB-1",
    auxiliary=["CONCENTRATION", "DENSITY"],
)

wellist1 = []
qwell = 5.7024 / nlay
for k in range(nlay):
    wellist1.append([(k, 0, 0), qwell, 0.])
wel = flopy.mf6.ModflowGwfwel(
    gwf,
    stress_period_data=wellist1,
    save_flows=False,
    auxiliary='CONCENTRATION',
    pname="WEL-1",
)
oc = flopy.mf6.ModflowGwfoc(
    gwf,
    budget_filerecord=f"{gwfname}.cbc",
    head_filerecord=f"{gwfname}.hds",
    headprintrecord=[
        ('COLUMNS', 10, 'WIDTH', 15,
         'DIGITS', 6, 'GENERAL')],
    saverecord=[('HEAD', 'ALL'),
                ('BUDGET', 'ALL')],
    printrecord=[('HEAD', 'LAST'),
                 ('BUDGET', 'LAST')])

gwt = flopy.mf6.ModflowGwt(sim, modelname=gwtname)
dis = flopy.mf6.ModflowGwtdis(
    gwt, 
    nlay=nlay, 
    nrow=nrow, 
    ncol=ncol,
    delr=delr, 
    delc=delc,
    top=top, 
    botm=botm,
    idomain=idomain,
)
ic = flopy.mf6.ModflowGwtic(gwt, strt=35.)
adv = flopy.mf6.ModflowGwtadv(gwt, scheme='UPSTREAM')

diffusion_only = True
if diffusion_only:
    diffc = 0.57024
    alh = 0.0
    ath = 0.0
    xt3d = False
else:
    diffc = 0.0
    alh = 0.1
    ath = 0.01
    xt3d = True
xt3d_off = not xt3d
dsp = flopy.mf6.ModflowGwtdsp(
    gwt, xt3d_off=xt3d_off, diffc=diffc, alh=alh, ath1=ath
)

mst = flopy.mf6.ModflowGwtmst(gwt, porosity=0.35)
sourcerecarray = [
    ("WEL-1", "AUX", "CONCENTRATION"),
    ("GHB-1", "AUX", "CONCENTRATION"),
]
ssm = flopy.mf6.ModflowGwtssm(gwt, sources=sourcerecarray)
fmi = flopy.mf6.ModflowGwtfmi(gwt, flow_imbalance_correction=True)
oc = flopy.mf6.ModflowGwtoc(
    gwt,
    budget_filerecord=f"{gwtname}.cbc",
    concentration_filerecord=f"{gwtname}.ucn",
    concentrationprintrecord=[
        ('COLUMNS', 10, 'WIDTH', 15,
         'DIGITS', 6, 'GENERAL')],
    saverecord=[('CONCENTRATION', 'ALL')],
    printrecord=[('CONCENTRATION', 'LAST'),
                 ('BUDGET', 'LAST')])
gwfgwt = flopy.mf6.ModflowGwfgwt(
    sim, 
    exgtype='GWF6-GWT6',
    exgmnamea=gwfname, 
    exgmnameb=gwtname,
)

imsgwf = flopy.mf6.ModflowIms(
    sim, 
    print_option='ALL',
    outer_dvclose=hclose,
    outer_maximum=nouter,
    under_relaxation='NONE',
    inner_maximum=ninner,
    inner_dvclose=hclose, 
    rcloserecord=rclose,
    linear_acceleration='BICGSTAB',
    scaling_method='NONE',
    reordering_method='NONE',
    relaxation_factor=relax,
    filename=f"{gwf.name}.ims",
)
imsgwt = flopy.mf6.ModflowIms(
    sim, 
    print_option='ALL',
    outer_dvclose=hclose,
    outer_maximum=nouter,
    under_relaxation='NONE',
    inner_maximum=ninner,
    inner_dvclose=hclose, 
    rcloserecord=rclose,
    linear_acceleration='BICGSTAB',
    scaling_method='NONE',
    reordering_method='NONE',
    relaxation_factor=relax,
    filename=f"{gwt.name}.ims",
)
sim.register_ims_package(imsgwf, [gwf.name])
sim.register_ims_package(imsgwt, [gwt.name])

sim.write_simulation()
sim.run_simulation()
```

#### Create Animation 1
```
cobj = gwt.output.concentration()
conc = cobj.get_data()
conc_alldata = cobj.get_alldata()
times = cobj.times
budgwf = gwf.output.budget()

def get_sealevel(totim):
    if totim <= 0.25:
        return sealevel
    else:
        return sinfunc(amplitude, frequency * 2 * np.pi, 0, sealevel, totim - 0.25)

fig, ax = plt.subplots(figsize=(8, 6))
ax.set_aspect(1)
ax.set_xlabel(r'x')
ax.set_ylabel(r'z')
title = ax.set_title(f"Time = {times[0]:.3f} days")

# plot persistent items
pxs = flopy.plot.PlotCrossSection(gwt, ax=ax, line={"row": 0})

sl = get_sealevel(0.)
seapoly = np.array([[lx * fx, sl], [lx, sl], [lx, 0]])
seapatch = matplotlib.patches.Polygon(
    seapoly, closed=True, facecolor="darkred", zorder=0
)
ax.add_patch(seapatch)

aqpoly = np.array(
    [[0, 0], [lx, 0], [lx, fz * lz], [lx * fx, lz], [0, lz]]
)
aqpatch = matplotlib.patches.Polygon(
    aqpoly, closed=True, facecolor=".7", zorder=1
)
ax.add_patch(aqpatch)
pxs.plot_grid(linewidth=0.1)
levels = 35 * np.array([0.01, 0.1, 0.5, 0.9, 0.99])
cmap = plt.colormaps["jet"].with_extremes(under="white", over="white")
ca_dict = {
    "cmap": cmap, 
    "extend": "both", 
    "levels":levels, 
    "filled":False,
    "norm": "linear",
}

cont = pxs.contour_array(conc_alldata[0], **ca_dict)
fig.colorbar(cont, shrink=0.5)

def animate(i):
    global cont
    global title
    for c in cont.collections:
        c.remove()  # removes only the contours, leaves the rest intact
    cont = pxs.contour_array(conc_alldata[i], **ca_dict)
    title = ax.set_title(f"Time = {times[i]:.3f} days")
    
    totim = times[i]
    sl = get_sealevel(totim)
    seapoly = np.array([[lx * fx, sl], [lx, sl], [lx, 0]])
    seapatch.set_xy(seapoly)
    
    return

frames = range(250, 1250, 10)
ani = matplotlib.animation.FuncAnimation(fig, animate, frames=frames)
plt.close()

from IPython.display import HTML
HTML(ani.to_jshtml())
```

#### Create Animation 2
```
def get_patch_collection(modelgrid, head, conc, cmap="jet", zorder=None):
    # create patches for each cell
    import matplotlib.collections
    import matplotlib.patches

    xv, yv, zv = modelgrid.xyzvertices
    botm = modelgrid.botm
    patches = []
    for k in range(modelgrid.nlay):
        for j in range(modelgrid.ncol):
            x0 = xv[0, j]
            x1 = xv[0, j + 1]
            z0 = zv[k, 0, j]
            z0 = min(z0, head[k, 0, j])
            z0 = max(z0, botm[k, 0, j])
            z1 = zv[k + 1, 0, j]
            poly = [[x0, z0], [x1, z0], [x1, z1], [x0, z1], [x0, z0]]
            patch = matplotlib.patches.Polygon(
                poly, closed=True, edgecolor="k", facecolor="red"
            )
            patches.append(patch)
    pc = matplotlib.collections.PatchCollection(
        patches, cmap=cmap, zorder=zorder
    )
    pc.set_array(conc.flatten())
    return pc

hobj = gwf.output.head()
head_alldata = hobj.get_alldata()
cobj = gwt.output.concentration()
conc = cobj.get_data()
conc_alldata = cobj.get_alldata()

times = cobj.times
budgwf = gwf.output.budget()
spdis = budgwf.get_data(text='DATA-SPDIS')

def get_sealevel(totim):
    if totim <= 0.25:
        return sealevel
    else:
        return sinfunc(amplitude, frequency * 2 * np.pi, 0, sealevel, totim - 0.25)

fig, ax = plt.subplots(figsize=(8, 6))
ax.set_aspect(1)
ax.set_xlabel(r'x')
ax.set_ylabel(r'z')
title = ax.set_title(f"Time = {times[0]:.3f} days")

pxs = flopy.plot.PlotCrossSection(gwt, ax=ax, line={"row": 0})

sl = get_sealevel(0.)
seapoly = np.array([[lx * fx, sl], [lx, sl], [lx, 0]])
seapatch = matplotlib.patches.Polygon(
    seapoly, closed=True, facecolor="darkred", zorder=0
)
ax.add_patch(seapatch)

h = head_alldata[0]
c = conc_alldata[0]
conc_pc = get_patch_collection(gwf.modelgrid, h, c, cmap="jet", zorder=None)
conc_pc.set_clim(0, 35.0)
ax.add_collection(conc_pc)

pxs.plot_grid(linewidth=0.1)
fig.colorbar(conc_pc, shrink=0.5)
botm = gwf.modelgrid.botm
qx, qy, qz = flopy.utils.postprocessing.get_specific_discharge(spdis[0], gwf)
quiver = pxs.plot_vector(qx, qy, qz, normalize=False, color="white")

def animate(i):
    global conc_pc
    global title
    global quiver
    
    conc_pc.remove()
    h = head_alldata[i]
    c = conc_alldata[i]
    c = np.ma.masked_greater(c, 1e20)
    c = np.ma.masked_where(h < botm, c)
    
    conc_pc = get_patch_collection(gwf.modelgrid, h, c, cmap="jet", zorder=None)
    conc_pc.set_clim(0, 35.0)
    ax.add_collection(conc_pc)

    qx, qy, qz = flopy.utils.postprocessing.get_specific_discharge(spdis[i], gwf)
    quiver.set_UVC(qx, qz)
    
    title = ax.set_title(f"Time = {times[i]:.3f} days")
    
    totim = times[i]
    sl = get_sealevel(totim)
    seapoly = np.array([[lx * fx, sl], [lx, sl], [lx, 0]])
    seapatch.set_xy(seapoly)
    
    return

frames = range(250, 1250, 10)
ani = matplotlib.animation.FuncAnimation(fig, animate, frames=frames)
plt.close()

from IPython.display import HTML
HTML(ani.to_jshtml())
#ani.save("ex7-animation.mp4")
```

#### Plot Sealevel Fluctuation
```
s = []
for t in times:
    s.append(get_sealevel(t))
plt.plot(times, s)
plt.xlabel("time, in days")
plt.ylabel("sealevel elevation, in meters")
```
