<a name="readme-top"></a>

<!-- PROJECT SHIELDS -->

[![License](https://img.shields.io/github/license/Sad-Abd/pyPolyMesher.svg?style=for-the-badge)](https://github.com/Sad-Abd/pyPolyMesher/blob/main/LICENSE)
[![Made with love in SUT (Iran)](https://img.shields.io/badge/Made%20with%20%E2%9D%A4%EF%B8%8F%20in-SUT%20(Iran)-0c674a?style=for-the-badge)](https://sut.ac.ir)
[![GitHub Stars](https://img.shields.io/github/stars/Sad-Abd/pyPolyMesher.svg?style=for-the-badge)](https://github.com/Sad-Abd/pyPolyMesher/stargazers)

<!-- PROJECT LOGO -->

<br />
<div align="center">
  <a href="https://github.com/Sad-Abd/pyPolyMehser">
    <img src="images/Logo.png" alt="Logo" width="426" height="272">
  </a>

<h3 align="center">pyPolyMesher</h3>

  <p align="center">
    Generation of polygonal Mesh
    <br />
    <a href="https://github.com/Sad-Abd/pyPolyMehser"><strong>Explore the docs »</strong></a>
    <br />
    <br />
    <a href="https://github.com/Sad-Abd/pyPolyMehser">View Demo</a>
    ·
    <a href="https://github.com/Sad-Abd/pyPolyMehser/issues">Report Bug</a>
    ·
    <a href="https://github.com/Sad-Abd/pyPolyMehser/issues">Request Feature</a>
  </p>
</div>



<!-- TABLE OF CONTENTS -->
<details>
  <summary>Table of Contents</summary>
  <ol>
    <li>
      <a href="#about-the-project">About The Project</a>
    </li>
    <li>
      <a href="#getting-started">Getting Started</a>
      <ul>
        <li><a href="#installation">Installation</a></li>
      </ul>
    </li>
    <li><a href="#usage">Usage</a>
      <ol>
        <li><a href="#1.-read-image">Read Image</a></li>
      </ol>
    </li>
    <li><a href="#roadmap">Roadmap</a></li>
    <li><a href="#contributing">Contributing</a></li>
    <li><a href="#license">License</a></li>
    <li><a href="#contact">Contact</a></li>
    <li><a href="#acknowledgments">Acknowledgments</a></li>
  </ol>
</details>

<!-- ABOUT THE PROJECT -->
## About The Project

**pyPolyMesher** is a python package for generating unstructured polygonal meshes in arbitrarily defined 2D domains. It allows users to mathematically specify domains using signed distance functions (SDFs) and generates high-quality meshes adapted to the geometry and features of the domain. **pyPolyMesher** was initially created as a Python version of the [MATLAB PolyMesher program](http://paulino.princeton.edu/software.html) but has since been enriched with additional features.


Key capabilities:

- Define 2D domains mathematically using signed distance functions
- Built-in library of SDF primitives (circles, rectangles etc.) and operations to construct complex domains
- Ability to define custom SDFs for new domain geometries
- Generate unstructured triangular/polygonal meshes adapted to domains
- Apply boundary conditions and mark fixed points
- Assess mesh quality metrics like element aspect ratio
- Animate mesh generation process
- Import and mesh polygons from DXF files

By leveraging SDFs to represent domains, pyPolyMesher can capture intricate geometries and generate optimized meshes tailored to them, making it useful for simulations and analysis.

The package provides Lloyd's algorithm for efficient and robust meshing of arbitrary SDF-based domains. Researchers can conveniently translate geometric constructs and concepts into code using the SDF formalism.

Overall, pyPolyMesher simplifies the entire workflow - from domain specification to quality mesh generation to numerical analysis.


<p align="right">(<a href="#readme-top">back to top</a>)</p>



<!-- GETTING STARTED -->
## Getting Started

This part explains how to install and use this package.

### Installation

Since pyPolyMesher is not published on PyPI [_yet_], you need to clone the repository:

```
git clone https://github.com/Sad-Abd/pyPolyMesher.git
```

Then install it using pip:

```
pip install ./pyPolyMesher
```

<p align="right">(<a href="#readme-top">back to top</a>)</p>

<!-- USAGE EXAMPLES -->
## Usage

### Basic Usage:

`pyPolyMesher.PolyMesherPolyMesher(Domain, NElem, MaxIter, P=None, anim=False)`: Generate polygon mesh on `Domain` with `NElem` number of elements. Improve mesh for `MaxIter` iterations. Can be given an initial point set `P`. 

```python
import pyPolyMesher
from pyPolyMesher.exampleDomains import MichellDomain
MichellDomain.Plot()
Node, Element, Supp, Load, P = pyPolyMesher.PolyMesher(MichellDomain, 50, 100)
```

### Internal SDFs:
```python
from pyPolyMesher import dFunctions as DF
```
1. `DF.dLine(P, x1, y1, x2, y2)`: Calculate the signed distance from points P to a line segment defined by two endpoints (x1, y1) and (x2, y2).
2. `DF.dCircle(P, xc, yc, r)`: Calculate the signed distance from points P to a circle defined by its center (xc, yc) and radius (r).
3. `DF.dRectangle(P, x1, x2, y1, y2)`: Calculate the signed distance from points P to a rectangle defined by its bottom-left (x1, y1) and top-right (x2, y2) coordinates.
4. `DF.dPolygon(P, vertices)`: Calculate the signed distance from points P to a polygon defined by its vertices.
5. `DF.dUnion(d1, d2)`: Calculate the signed distance field resulting from the union of two distance fields (d1 and d2).
6. `DF.dIntersect(d1, d2)`: Calculate the signed distance field resulting from the intersection of two distance fields (d1 and d2).
7. `DF.dDiff(d1, d2)`: Calculate the signed distance field resulting from the difference of two distance fields (d1 and d2).

### Example Domains:

1. `pyPolyMesher.exampleDomains.MbbDomain`

![MbbDomain](images/MBB_random.png)

2. `pyPolyMesher.exampleDomains.HornDomain`

![HornDomain](images/Horn.png)

3. `pyPolyMesher.exampleDomains.WrenchDomain`

![WrenchDomain](images/Wrench.png)

4. `pyPolyMesher.exampleDomains.MichellDomain`

![MichellDomain](images/Michell.png)

5. `pyPolyMesher.exampleDomains.SuspensionDomain`

![SuspensionDomain](images/suspension.png)

6. `pyPolyMesher.exampleDomains.CookDomain`

![CooksMembrane](images/cook_mesh.png)

### Import Polygon Domain from DXF:

```python
from from pyPolyMesher import PolyMesher, Domain, mesh_assessment
from pyPolyMesher.dxfImporter import dxf_polygon
from pyPolyMesher.dFunctions import dPolygon

dxf_file_path = 'examples/polygon1.dxf'
v = dxf_polygon(dxf_file_path)

SDF = lambda P: dPolygon(P, v)
dxfDomain = Domain("DXF Polygon Domain", [0,100,0,100], SDF)
dxfDomain.Plot()
```

![polygon_dxf](images/polygon_dxf.png)

```python
Node, Element, Supp, Load, P = PolyMesher(dxfDomain, 50, 100)
mesh_assessment(Node, Element)
```
![polygon_dxf_mesh](images/polygon_dxf_mesh.png)


See [Examples.py](examples/Examples.py) for more examples.

<p align="right">(<a href="#readme-top">back to top</a>)</p>

<!-- ROADMAP -->
## Roadmap

### Section 1 - Current Focus and Issue Resolution

1. ~~Translate other domain codes.~~
2. ~~Add docstrings and comments to the following files:~~
   - ~~pyPolyMesher~~
   - ~~pydFunction~~
   - ~~Domains~~
3. ~~Make the resequence function work properly ([Issue #3](https://github.com/Sad-Abd/pyPolyMesher/issues/3)).~~
4. ~~Transform Domain definitions into classes using Object-Oriented Programming (OOP) principles.~~
5. ~~Rethink the file hierarchy and user experience.~~
6. Use Jupyter notebook for `Examples.py`.
7. Use Jupyter notebook to illustrate `Domain` creation.
8. Define `Domain` from `dxf` files
    * ~~Polygon importer~~
    * Circle importer
    * Spline importer
    * Automatic SDF for geometries
9. Add mesh quality assessments
    * ~~Aspect Ratio~~
    * ~~Standard Deviation of Elements Areas~~
10. Add some example meshes to the **README**.

### Section 2 - Upcoming Priorities
1. Enhance the **README** with more detailed information.
2. Publish the package on *PYPI* and *Zenodo* for wider distribution.
3. Add some tests.

### Section 3 - Vision and Future Prospects

1. Develop a GUI for domain definition to improve user interaction.
2. Plugin for CAD programs.
3. Explore and brainstorm alternative options for domain definition and future possible expansions.

See the [open issues](https://github.com/Sad-Abd/pyPolyMesher/issues) for a full list of proposed features (and known issues).

<p align="right">(<a href="#readme-top">back to top</a>)</p>

<!-- CONTRIBUTING -->
## Contributing

Contributions are what make the open source community such an amazing place to learn, inspire, and create. Any contributions you make are **greatly appreciated**.

If you have a suggestion that would make this better, please fork the repo and create a pull request. You can also simply open an issue with the tag "enhancement".
We appreciate your interest in pyPolyMesher!. Don't forget to give the project a star! Thanks again!


<p align="right">(<a href="#readme-top">back to top</a>)</p>

<!-- LICENSE -->
## License

This project is licensed under the GPLv3 License - see the LICENSE file for details.
Contact

<p align="right">(<a href="#readme-top">back to top</a>)</p>

## Contact

If you have any questions or feedback, feel free to reach out:

Email: AbediSadjad@gmail.com

GitHub: [Sad-Abd](https://github.com/Sad-Abd)
