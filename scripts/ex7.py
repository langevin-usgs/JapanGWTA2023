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

# create animation 1
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

#from IPython.display import HTML
#HTML(ani.to_jshtml())

print ("saving animation...")
ani.save("ex7-animation1.mp4")

# create animation 2
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

# use this in notebook
#from IPython.display import HTML
#HTML(ani.to_jshtml())

print ("saving animation...")
ani.save("ex7-animation2.mp4")
print ("done...")
