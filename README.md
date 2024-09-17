# Diatomic Molecule Dissociation by Gravitational Gradient
This project was undertaken in the summer of 2024. It comprises an investigation into what conditions are required for a molecule's atoms to dissociate when close to a Kerr-Newman-Taub-NUT black hole. This was performed by first obtaining the analytical equations of motion using the Euler-Lagrange equations for a particle in the aforementioned metric and a Lennard-Jones potential in Mathematica. The expressions obtained were then compressed and manipulated using Python into a form which can be evaluated using the Euler method in C++. The results were then analysed and visualised using Python again.

## Theory
The total Lagrangian of the system is
```math
L = L_g + L_Q + L_m \,\,,
```
where $`L_g`$ is the Lagrangian for the metric in Boyer–Lindquist coordinates ( $`x_{\mu} = (t, r, \phi, \theta)`$ ), $`L_Q`$ is the Lagrangian for the influence of the electric charge of the black hole and $`L_m`$ is the Lagrangian for the Lennard-Jones potential between the atoms in the diatomic molecule.

The metric Lagrangian is given by
```math
L_g = -\frac{\Delta}{\Sigma} \big(t' - \chi\phi'\big)^2 + \Sigma \Big(\frac{r'^2}{\Delta} + \theta'^2\Big) + \frac{\sin^2(\theta)}{\Sigma} \big(a t' - (r^2 + l^2 + a^2)\phi'\big)^2 \,\,,
```
where
```math
\Sigma = r^2 + \big(l + a\cos(\theta)\big)^2 \,\,,
```
```math
\Delta = r^2 - 2Mr + a^2 - l^2 + Q^2 \,\,,
```
```math
\chi = a\sin^2(\theta) - 2l\cos(\theta) \,\,,
```
$`a`$ is the black hole angular momentum to mass ratio, $`l`$ is the black hole gravitomagnetic monopole moment and $`Q`$ is the black hole electric charge.


The black hole electric charge Lagrangian is
```math
L_Q = \frac{q}{m} \frac{Qr}{\Sigma} (-t' + \chi\phi') \,\,,
```
where $`q`$ is the particle's charge and $`m`$ is the particle's mass.


The Lagrangian for the Lennard-Jones potential is given by
```math
L_m = \frac{1}{m} V_{LJ} \beta\cdot x \,\,,
```
where $`\beta`$ is the 4-velocity of the source, $`\cdot`$ is the Minkowski product and
```math
V_{LJ} = 4\epsilon \Bigg(\Big(\frac{\sigma}{|\vec{n}|}\Big)^{12} - \Big(\frac{\sigma}{|\vec{n}|}\Big)^6\Bigg)
```
is the Lennard-Jones potential, where $`\epsilon`$ is the depth of the potential well, $`\sigma`$ is the distance at which the particle-particle potential energy $`V_{LJ}`$ is zero and $`\vec{n}`$ is a 3-vector from the source particle to the test particle.

## Requirements
### C++
* boost/multiprecision/cpp_dec_float.hpp
### Python
* Numba
* tqdm
* pathlib
* mpmath
* Matplotlib
* NumPy

## Running
Firstly, run the `.nb` Mathematica file to generate the expressions in the `.cpp` file. Then run the `Euler-Lagrange preconfigure.ipynb` file to compress the expressions and write the file in a form that can be evaluated. Next, run the `.cpp` main file given the same name as this project. Lastly, analyse your results using the `Data Analysis.ipynb` file.

Attributes of the system can be changed using the `properties.txt` file in `/data/`.

## Next Steps
Modify the project so that a variety of variables can be sampled in one run to determine the structure of the parameter space that allows molecules to dissociate.


By finding the velocities such that the molecules dissociate one should be able to combine this with data for the expected distribution of velocities $`p(\vec{v})`$ to obtain the probability of one molecule dissociating given a position in space $`\vec{r}`$ and orientation $`\hat{\omega}`$, which is the unit vector from one atom to the other. This probability would be expressed as
```math
P(\vec{r}, \hat{\omega}) = \int\limits_V d^3\vec{v} \, p(\vec{v}) \,\,,
```
with
```math
V = \Big\{\vec{v} \in  \mathbb{R}^3 \Big{|} \, |\vec{v}| < \delta_i(\hat{\omega}, \vec{r}) \cap |\vec{v}| > \delta_j(\hat{\omega}, \vec{r}) \,\, ; \, i \neq j \in \\{1,2\\} \Big\} \,\,,
```
where $`\delta_i(\hat{\omega}, \vec{r})`$ is the velocity required for particle $i$ to escape. This can next be integrated over all orientations and all radii greater than the Schwarzschild radius while also including the number density profile with distance to the black hole. By differentiating with respect to the solid angle from the origin of the system one can obtain the differential number of dissociation events with respect to the solid angle
```math
\frac{dN_{\text{diss}}}{d\Omega} = \int\limits_{r_s}^{\infty} dr \, r^2 n(r) \int\limits_{S^2} d\omega \, P(\vec{r}, \hat{\omega}) \,\,.
```
In principle, this could be used to obtain the number of expected dissociation events observed in any direction.
