# Software Installation

Lecture attendees may wish to follow along by running the FloPy/MODFLOW 6 examples.  This will require that software be installed prior to the lecture.  Required software includes:

1.  Python (preferably 3.10 or newer).  [Miniconda](https://docs.conda.io/en/latest/miniconda.html) is a popular Python distribution that runs on Windows, Mac, and Linux.

2.  Use conda to install a customized environment for FloPy and MODFLOW. The first step is create a file called `flopy_environment.yml`.  This file should contain the following contents.

```
name: flopy

channels:
  - conda-forge

dependencies:
  - python>=3.8
  - numpy>=1.15.0
  - matplotlib>=1.4.0
  - flopy
  - python-dateutil>=2.4.0
  - affine
  - scipy
  - pandas
  - pyshp
  - rasterio
  - fiona
  - descartes
  - pyproj
  - shapely>=1.8
  - geos
  - geojson
  - vtk
  - rasterstats
  - pyvista
  - imageio
  - pymetis
  - trame
  - jupyter

  # pip installations
  - pip
  - pip:
      - git+https://github.com/MODFLOW-USGS/modflowapi.git
```

After this file is created, run the following command from a terminal.

```
$ conda env create -f flopy_environment.yml
```

Whenever you want to use this FloPy environment, you will need to activate it with the following command.

```
conda activate flopy
```

3.  Install MODFLOW executables and related programs.  After FloPy has been installed, MODFLOW executables can be downloaded to your computer using the [get-modflow](https://github.com/modflowpy/flopy/blob/develop/docs/get_modflow.md) utility, which is installed with FloPy. The following command is one way to use `get-modflow`, which will download the executables and make them available for use with FloPy scripts.  Make sure that the flopy environment is active before you run `get-modflow`.

```
$ get-modflow :flopy 
```
