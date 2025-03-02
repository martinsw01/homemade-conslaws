# Homemade Conslaws

This is a package in Julia written for my specialization project during the last year of my master's degree in applied physics and mathematics at NTNU. The project aims to get an intuition about conservation laws and their numerical finite volume (FVM) solutions. Furthermore, a FVM solver inspired by SINTEF Digital's solver, [`SinFVM.jl`](https://github.com/sintefmath/SinFVM.jl) (previously `SinSWE.jl`) is implemented in order to investigate computational difficulties with shallow water simulations.

One of these problems can be seen in the following animation:

<img class="ignore-theme" src="example\_computational\_difficulties.gif" alt="Animation demonstrating computational difficulties">

The CFL-condition may in some situations severilly restrict the time step size, potentially making the simulation very slow. This and other problems are describe din detail in the final report, which can be found [here](https://martinsw01.github.io/homemade-conslaws/dev/final_report).