# Finite volume schemes for conservation laws

## Plan for coding

The supervisors suggests to begin experimenting with implementing some finite volume schemes in order to get familiar with the numerics and Julia.
The following equations gives invcreasing difficoulty.

- [ ] Advection equation
- [ ] Burgers' equation
- [ ] Shallow-water equations

## Methods

* Lax-Friedrich
* Lax-Wendrof
* Central-upwind (se also Section 4.5 in gpu-conslaws (2))

It will also be beneficial to read about automatic differentiation (AD). The supervisors suggests reading the relevant chapters of Torjei's master thesis and play around with the Julia package ForwardDiff.jl.

## Reading material

### Textbooks and lecture notes

1. Lecture notes from Sid Mishra for [Numerical methods for conservation laws
and related equations](https://www.uio.no/studier/emner/matnat/math/MAT-IN9240/h17/pensumliste/numcl_notes.pdf), mainly chapters 1 - 5
2. "How to Solve Systems of Conservation Laws Numerically Using the Graphics Processor as a High-Performance Computational Engine" by Knut-Andreas and previous colleagues at SINTEF. Good introduction to numerical solutions of conservation laws.
3. "Finite-Volume Methods for Hyperbolic Problems" by LeVeque - textbook that can be used to look up theory or as a third alternative to the above.


### Relevant articles for motivation or for later

4. Kurganov, Noelle and Petrova - central-upwind scheme, relevant for the flood simulation (see also ch 4.5 in (2) above).
5. Kurganov, Petrova - variant of central-upwind for shallow-water equations
6. Brodtkorb, Sætra, Altinakar - Written by those from SINTEF that were in the US to simulate dam breach and flood. Uses the numerical method from Kurganov, Petrova
7. Fernandez-Pato et al. - Rainfall/runoff simulation with 2D full shallow water equations

