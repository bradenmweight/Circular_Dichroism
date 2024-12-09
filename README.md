# Linear Absorption and Circular Dichroism Spectroscopy
    ~~~ Braden M. Weight -- May 9, 2024 ~~~

## Theory
Electric Dipole Moment:
```math
\hat{\mu}^\mathrm{el} = -\sum_{j}^{N_\mathrm{el}} {\bf \hat{r}}_j + \sum_{J}^{N_\mathrm{ions}} Z_J {\bf R}_J \otimes \mathbb{I}_\mathrm{el}
```

Magnetic Dipole Moment:
```math
\hat{\mu}^\mathrm{mag} = -\sum_{j}^{N_\mathrm{el}} {\bf \hat{r}}_j \times \frac{i}{\hbar}{\bf \hat{p}}_j
```

Note: $\hat{p}_x = -i\hbar \hat{\nabla}_x$

Oscillator Strength:
```math
f_{0\alpha} = 2 E_{0\alpha}~\hat{\mu}_{0\alpha}^\mathrm{el} \cdot \hat{\mu}_{\alpha 0}^\mathrm{el}
```

Rotary Strength: 
```math
R_{0\alpha} = 2 E_{0\alpha}~\hat{\mu}^\mathrm{el} \cdot \hat{\mu}^\mathrm{mag}, where $E_{0\alpha} = E_\alpha - E_0$
```

Absorption Spectroscopy:
```math
\begin{align}
\mathrm{ABS}(E) &= \sum_\alpha f_{0\alpha}~\delta(E - E_{0\alpha})\\
&\approx \sum_\alpha f_{0\alpha} \mathrm{e}^{-\frac{(E - E_{0\alpha})^2}{2\sigma^2}}\\
&\approx \sum_\alpha f_{0\alpha}~\frac{\frac{\sigma^2}{4}}{(E - E_{0\alpha})^2 + (\frac{\sigma}{2})^2}
\end{align}
```

Circular Dichroism Spectroscopy:
```math
\begin{align}
\mathrm{CD}(E) &= \sum_\alpha R_{0\alpha}~\delta(E - E_{0\alpha})\\
&\approx \sum_\alpha R_{0\alpha} \mathrm{e}^{-\frac{(E - E_{0\alpha})^2}{2\sigma^2}}\\
&\approx \sum_\alpha R_{0\alpha}~\frac{\frac{\sigma^2}{4}}{(E - E_{0\alpha})^2 + (\frac{\sigma}{2})^2},
\end{align}
```
where $\sigma$ is a finite broadening parameter

## The Fermionic RPA Solution
Practically speaking, 
```math
\boldsymbol{\mu}_{0\alpha}^\mathrm{el} = \mathrm{Tr}[\hat{\mu}^\mathrm{el}~\hat{\xi}_{0\alpha}] = \int dr \mu^\mathrm{el}(r)~\xi_{0\alpha}(r) = \int d\boldsymbol{r}~\psi^*_\alpha(\boldsymbol{r})~\boldsymbol{r}~\psi_0(\boldsymbol{r}')
```
and 
```math
\mu_{0\alpha}^\mathrm{mag} = \mathrm{Tr}[\hat{\mu}^\mathrm{mag}~\hat{\xi}_{0\alpha}] = \int d\boldsymbol{r} \mu^\mathrm{mag}(\boldsymbol{r})~\xi_{0\alpha}(\boldsymbol{r}) = \int d\boldsymbol{r}~\psi^*_\alpha(\boldsymbol{r})~(\boldsymbol{r}\times \boldsymbol{\hat{\nabla}})~\psi_0(\boldsymbol{r})
```
where $\xi_{0\alpha}(r,r') = \psi_\alpha^*(r)\psi_0(r')$ is the one-particle transition density matrix between the ground $0$ and excited state $\alpha$. In linear-response TDDFT, the transition density $\xi_{0\alpha} = X + Y$ where $| X,Y\rangle$ is the eigenvector of the RPA operator which satisfy,
```math
\begin{bmatrix}
A & B\\
-A^* & -B^*
\end{bmatrix}
\begin{bmatrix}
X\\
Y
\end{bmatrix} =
\omega_{0\alpha}\begin{bmatrix}
1 & 0\\
0 & -1
\end{bmatrix}
\begin{bmatrix}
X\\
Y
\end{bmatrix}
```
Here, $\omega_{0\alpha}$ is the ground-to-excited transition energy and $|\xi^{0 \alpha}\rangle \equiv | X^{0 \alpha}, Y^{0 \alpha} \rangle$ can be interpreted as excitation and de-excitation vectors of the $0\rightarrow\alpha_\mathrm{th}$ electronic transition. The normalization can be written as $\langle\xi^{0 \alpha}|\xi^{0 \beta}\rangle = \langle X^{0 \alpha} | X^{0 \beta}\rangle - \langle Y^{0 \alpha} | Y^{0 \beta}\rangle= \delta_{\alpha\beta}$.

Importantly, the RPA eigenvalue equation can be re-cast into a different form as,
```math
(A + B) (A - B) (X^{0\alpha} - Y^{0\alpha}) = \omega_{0\alpha} (X^{0\alpha} - Y^{0\alpha})
```
whose normalization is then defined as $\langle z^{0\alpha} | z^{0\beta} \rangle= \langle X^{0\alpha} - Y^{0\alpha} | X^{0\beta} - Y^{0\beta} \rangle = \delta_{\alpha\beta}$, where $|z^{0\alpha}\rangle = | X^{0\alpha} - Y^{0\alpha} \rangle$.

## Matrix Elements of an Observable
To find expectation values (needed for transition properties), one must perform the trace of the operator with the transition densities as,
```math
\begin{align}
\langle \psi_\alpha | \hat{O} | \psi_0 \rangle &= \mathrm{Tr}\big[ \hat{O}~\xi^{0\alpha} \big]\\
&= \sum_i^{N_\mathrm{occ}}\sum_a^{N_\mathrm{vir}} O_{ia}~\xi_{ia}^{0\alpha}\nonumber\\
&= \sum_{\mu\nu}^{N_\mathrm{AO}} O_{\mu\nu}~\xi_{\mu\nu}^{0\alpha}\nonumber
\end{align}
```
where $\xi_{\mu\nu}^{0\alpha} = \sum_ i^{N_\mathrm{occ}} \sum_ a^{N_\mathrm{vir }}C_{\mu i}~\xi_{ia}^{0\alpha}~C_{\nu a}$ is the atomic orbital representation of the transition density matrix and $\{C_{\mu i}\}$ are the molecular orbitals in the basis of the atomic orbitals.




## Calculating Electronic and Magnetic Transition Dipoles
```math
\begin{align}
\langle \psi_\alpha | \vec{\hat{\mu}}^\mathrm{el} | \psi_0 \rangle &= \sum_{\mu\nu}^{N_\mathrm{AO}} (\vec{\hat{\mu}}^\mathrm{el})_{\mu\nu}~\xi_{\mu\nu}^{0\alpha}\\ &= \sum_{\mu\nu}^{N_\mathrm{AO}} (\vec{\hat{\mu}}^\mathrm{el})_{\mu\nu}~(X_{\mu\nu}^{0\alpha} + Y_{\mu\nu}^{0\alpha})\nonumber\\\langle \psi_\alpha | \vec{\hat{\mu}}^\mathrm{mag} | \psi_0 \rangle &= \sum_{\mu\nu}^{N_\mathrm{AO}} (\vec{\hat{\mu}}^\mathrm{mag})_{\mu\nu}~\xi_{\mu\nu}^{0\alpha}\\&= \sum_{\mu\nu}^{N_\mathrm{AO}} (\vec{\hat{\mu}}^\mathrm{mag})_{\mu\nu}~(X_{\mu\nu}^{0\alpha} - Y_{\mu\nu}^{0\alpha})\nonumber
\end{align}
```
Since the electronic dipole operator is Hermitian, i.e., $\mu_{\mu\nu}^\mathrm{el} = \mu_{\nu\mu}^\mathrm{el}$, the density required is $X + Y$ while for the magnetic dipole operator, which is anti-Hermitian, i.e., $\mu_{\mu\nu}^\mathrm{el} = -\mu_{\nu\mu}^\mathrm{el}~$, $X-Y$ is required. (BMW ~ Are there any tricks to get around having to know whether the operator is hermi or not ?)

## Real-space Projected Transition and Dipole Density
The elements of the transition dipole matrix in the real-space representation can be written as a 6D-object,
```math
\xi(x,y,z,x',y',z') = \sum_{\mu\nu}^{N_\mathrm{AO}}\langle x y z | \mu \rangle~\xi_{\mu\nu}~\langle \nu | x' y' z' \rangle.
```
Taking only the diagonal elements, such that $|x'y'z'\rangle = |xyz\rangle$, we get $\xi(x,y,z)$, a 3D object. This 3D object is natural to cast as an isosurface (volumetric data), in a similar way as molecular orbitals, total density, or electrostatic potentials are visualized. Furthermore, one can also generate the volumetric data for any observable by simply not performing the trace and taking only the diagonal elements of the Eq. 7.
```math
\begin{align}
\vec{\mu}_{0\alpha}^\mathrm{el}(x,y,z) &= \sum_{\mu\nu\gamma}^{N_\mathrm{AO}}\langle x y z | \mu \rangle~\vec{\mu}_{\mu\nu}^\mathrm{el}~\xi^{0\alpha}_{\nu\gamma}~\langle \gamma | xyz \rangle\\
\vec{\mu}_{0\alpha}^\mathrm{mag}(x,y,z) &= \sum_{\mu\nu\gamma}^{N_\mathrm{AO}}\langle x y z | \mu \rangle~\vec{\mu}_{\mu\nu}^\mathrm{mag}~\xi^{0\alpha}_{\nu\gamma}~\langle \gamma | xyz \rangle
\end{align}
```

## Natural transition orbitals (NTOs) 
NTOs can be easily constructed by diagonalizing (i.e., in general SVD) the transition density matrix in the MO representation to obtain 
```math
\begin{align}
\xi^{0\alpha}_{ia} &= \sum_s^{N_\mathrm{occ}\times N_\mathrm{vir}} U^{0\alpha}_{ia,s}~\lambda_s^{0\alpha}~V^{0\alpha}_{ia,s}\\
\mathrm{NTO}^{0\alpha}_{s} &= \sum_{i}^{N_\mathrm{occ}}\sum_a^{N_\mathrm{vir}} U_{ia,s}^{0\alpha} \ \ \ (\mathrm{Occupied})\\
\mathrm{NTO}^{0\alpha}_{s} &= \sum_{i}^{N_\mathrm{occ}}\sum_a^{N_\mathrm{vir}} V_{ia,s}^{0\alpha} \ \ \ (\mathrm{Virtual})
\end{align}
```

## Instructions for Gaussian16 + MultiWfn
```
1. Run Gaussian16 TDDFT with "IOp(6/8=3) IOp(9/40=5)"
2. Plot ABS and CD spectra: 
  --> python3 make_spectroscopy.py
3. Generate densities:
  --> ./make_transition_density.slurm
  --> ./make_dipole_density.slurm
  --> ./make_magnetic_dipole_density.slurm
4. Generate various plots based on densities:
  --> compute_plots_projections.py
```
Note: This script computes the rotary strength density from the transition dipole and transition magnetic dipole densities as output from MultiWfn (V3.7 -- https://sobereva.com/multiwfn/) using only the TDA excitation coefficients, even if the de-excitation coefficients exist (i.e., full RPA).

## Instructions for PySCF (in development)
```
1. Run the example file at 'pyscf/run_example.py'

!!!! Warning -- Results are not yet consistent with Gaussian16/MultiWfn implementation
TODO -- Cube file generator not yet working. Atom/density not coinciding...
TODO -- Implement TDA version to check with MultiWfn (modify CI normalization in this case)
```

