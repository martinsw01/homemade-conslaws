var documenterSearchIndex = {"docs":
[{"location":"generated/theory/linear_transport_eqs/#Linear-transport-equations","page":"Linear transport equations","title":"Linear transport equations","text":"","category":"section"},{"location":"generated/theory/linear_transport_eqs/","page":"Linear transport equations","title":"Linear transport equations","text":"The simplest case of the one-dimensional version of the transport equation,","category":"page"},{"location":"generated/theory/linear_transport_eqs/","page":"Linear transport equations","title":"Linear transport equations","text":"beginequation\n    U_t + a(x t) U_x = 0 quad text in  R times R_+ labeleqtransport\nendequation","category":"page"},{"location":"generated/theory/linear_transport_eqs/","page":"Linear transport equations","title":"Linear transport equations","text":"is when the velocity field a is constant. The Cauchy problem is solved using the method of characteristics, where the solution is constant along the characteristic curves x = x_0 + a t. The solution is then","category":"page"},{"location":"generated/theory/linear_transport_eqs/","page":"Linear transport equations","title":"Linear transport equations","text":"beginequation\n    U(x t) = U_0(x - a t) labeleqtransport_const\nendequation","category":"page"},{"location":"generated/theory/linear_transport_eqs/","page":"Linear transport equations","title":"Linear transport equations","text":"for any (xt) in R times R_+. We see that the initial data is transported with the velocity a. For more general cases, the characteristics may not be possible to solve explicitly. However, we can obtain some information of the sulutions with the following a priori energy estimate:","category":"page"},{"location":"generated/theory/linear_transport_eqs/","page":"Linear transport equations","title":"Linear transport equations","text":"lemma: Lemma\nAssume U is a smooth solution of the transport equation decaying to zero at infinity for all t in R_+ and a in C^1(R R_+). Then, U satisfies the energy bound    int_R U^2(x t) dd x leq e^norma_C^1t int_R U_0^2(x) dd xproof: Proof\nFollows by multiplying the transport equation by U and integrating over space. Then, use that U decays to zero at infinity and apply Grönwall's inequality.","category":"page"},{"location":"generated/theory/linear_transport_eqs/","page":"Linear transport equations","title":"Linear transport equations","text":"The lemma shows that the energy is bounded. Using another functional, the assumptions on U can be relaxed:","category":"page"},{"location":"generated/theory/linear_transport_eqs/","page":"Linear transport equations","title":"Linear transport equations","text":"lemma: Lemma\nAssume U is a smooth bounded solution of the transport equation. Then, we have    sup_x in R absU(x t) leq normU_0_L^infty(R)for any t  0. proof: Proof\nFor any (x t) in R times R_+, there exists xi in R such that U(x t) = U_0(xi)","category":"page"},{"location":"generated/theory/linear_transport_eqs/#Finite-difference-schemes-for-the-transport-equation","page":"Linear transport equations","title":"Finite difference schemes for the transport equation","text":"","category":"section"},{"location":"generated/theory/linear_transport_eqs/","page":"Linear transport equations","title":"Linear transport equations","text":"For some velocity fields a(xt), it may not be possible to derive an explicit formula for the characteristic equation. We thus use numerical methods to approximate the solutions of the transport equation.","category":"page"},{"location":"generated/theory/linear_transport_eqs/#Discretization","page":"Linear transport equations","title":"Discretization","text":"","category":"section"},{"location":"generated/theory/linear_transport_eqs/","page":"Linear transport equations","title":"Linear transport equations","text":"For simplicity, assume that the velocity field is positive. AS R is unbounded, we truncate the domain into Omega = x_L X_R. Thus, we must impose boundary conditions, which will be discussed below. For simplicity, we use a uniform mesh of mesh size Delta x with N+1 points x_j:","category":"page"},{"location":"generated/theory/linear_transport_eqs/","page":"Linear transport equations","title":"Linear transport equations","text":"    x_L = x_0  x_1  cdots  x_N = x_R","category":"page"},{"location":"generated/theory/linear_transport_eqs/","page":"Linear transport equations","title":"Linear transport equations","text":"We further choose some terminal time T and divide into M+1 points t^n = n Delta t. We set the initial approximaition U_j^0 = U_0(x_j) and update the next approximation U_j^n+1 using a finite difference scheme.","category":"page"},{"location":"generated/theory/linear_transport_eqs/#Centered-finite-difference-scheme","page":"Linear transport equations","title":"Centered finite difference scheme","text":"","category":"section"},{"location":"generated/theory/linear_transport_eqs/","page":"Linear transport equations","title":"Linear transport equations","text":"One such scheme is the forward difference in time and cetral difference in space approximating (refeqtransport_const):","category":"page"},{"location":"generated/theory/linear_transport_eqs/","page":"Linear transport equations","title":"Linear transport equations","text":"    fracU_j^n+1 - U_j^nDelta t + fraca(U_j+1^n - U_j-1^n)2 Delta x = 0\n    quad 0  j  N","category":"page"},{"location":"generated/theory/linear_transport_eqs/#Example","page":"Linear transport equations","title":"Example","text":"","category":"section"},{"location":"generated/theory/linear_transport_eqs/","page":"Linear transport equations","title":"Linear transport equations","text":"We consider the domain 0 1 with initial data","category":"page"},{"location":"generated/theory/linear_transport_eqs/","page":"Linear transport equations","title":"Linear transport equations","text":"    U_0(x) = sin(2pi x)","category":"page"},{"location":"generated/theory/linear_transport_eqs/","page":"Linear transport equations","title":"Linear transport equations","text":"u0(x) = sin(2*pi*x)","category":"page"},{"location":"generated/theory/linear_transport_eqs/","page":"Linear transport equations","title":"Linear transport equations","text":"u0 (generic function with 1 method)","category":"page"},{"location":"generated/theory/linear_transport_eqs/","page":"Linear transport equations","title":"Linear transport equations","text":"and a = 1. Since the data is periodic, we impose periodic boundary conditions. Numerically, we implement this by setting","category":"page"},{"location":"generated/theory/linear_transport_eqs/","page":"Linear transport equations","title":"Linear transport equations","text":"    U_0^n = U_N^n quad U_N+1^n = U_1^n","category":"page"},{"location":"generated/theory/linear_transport_eqs/","page":"Linear transport equations","title":"Linear transport equations","text":"Thus, on the boundary, j = 0 N, we have","category":"page"},{"location":"generated/theory/linear_transport_eqs/","page":"Linear transport equations","title":"Linear transport equations","text":"beginalign*\n    fracU_1^n+1 - U_1^nDelta t + fraca(U_2^n - U_N^n)2 Delta x = 0 \n    fracU_N^n+1 - U_N^nDelta t + fraca(U_1^n - U_N-1^n)2 Delta x = 0\nendalign*","category":"page"},{"location":"generated/theory/linear_transport_eqs/","page":"Linear transport equations","title":"Linear transport equations","text":"For the first time step n = 1, we have","category":"page"},{"location":"generated/theory/linear_transport_eqs/","page":"Linear transport equations","title":"Linear transport equations","text":"beginequation\n    U_j^1 = U_0(x_j) - fracDelta t2 Delta x (U_0(x_j+1) - U_0(x_j-1)) labeleqtransport_first\nendequation","category":"page"},{"location":"generated/theory/linear_transport_eqs/","page":"Linear transport equations","title":"Linear transport equations","text":"In addition, reshape U_j^n into a column vector U_j^n = U_j+nN. Then, then above equations can be formulated as A bm U = bm f. Then the first N entries of bm f are the right hand side of (refeqtransport_first) and the rest are zero. The matrix A will then be a block tridiagonal looking like this:","category":"page"},{"location":"generated/theory/linear_transport_eqs/","page":"Linear transport equations","title":"Linear transport equations","text":"(Image: )","category":"page"},{"location":"generated/theory/linear_transport_eqs/","page":"Linear transport equations","title":"Linear transport equations","text":"Using a grid of 50 mesh points, simulating to time T = 3, we get the following result:","category":"page"},{"location":"generated/theory/linear_transport_eqs/","page":"Linear transport equations","title":"Linear transport equations","text":"x_L, x_R = 0, 1\nT = 2\nN = 50\ndx = (x_R - x_L)/N\ndt = 1 * dx\n\nx = x_L+dx:dx:x_R\nt = dt:dt:T\nM = length(t)\n\nA, f = central_difference.central_difference_scheme(N, M, dx, dt, u0)\nU = A \\ f\n\nViz.animate_solution(reshape(U, N, M),\n                     (x, t) -> u0(x-t),\n                     x, t)","category":"page"},{"location":"generated/theory/linear_transport_eqs/","page":"Linear transport equations","title":"Linear transport equations","text":"(Image: )","category":"page"},{"location":"generated/theory/linear_transport_eqs/","page":"Linear transport equations","title":"Linear transport equations","text":"Intuitively, the information should propagate from left to right. However, the scheme uses information from both sides. This can be explained rigorously using the the discrete energy","category":"page"},{"location":"generated/theory/linear_transport_eqs/","page":"Linear transport equations","title":"Linear transport equations","text":"    E^n = frac12 Delta x sum_j (U_j^n)^2","category":"page"},{"location":"generated/theory/linear_transport_eqs/","page":"Linear transport equations","title":"Linear transport equations","text":"We know that the exact solution have bounded energy. We say that a scheme is energy stable if E^n leq E^0 for all n.","category":"page"},{"location":"generated/theory/linear_transport_eqs/","page":"Linear transport equations","title":"Linear transport equations","text":"lemma: Lemma\nLet U_j^n be the solutions computed with the central difference scheme. Then, we have the  following dicsrete energy estimate:    E^n+1 = E^n + fracDelta x2 sum_j (U_j^n+1 - U_j^n)^2Thus, the energy grows at each time step for any choice of Delta x and Delta t.proof: Proof\nSimilar to the proof of the continuous energy estimate and using the identity    d_2(d_1 - d_2) = frac12qty(d_1^2 - d_2^2)","category":"page"},{"location":"generated/theory/linear_transport_eqs/#Upwind-scheme","page":"Linear transport equations","title":"Upwind scheme","text":"","category":"section"},{"location":"generated/theory/linear_transport_eqs/","page":"Linear transport equations","title":"Linear transport equations","text":"To respect the flow of information, we can use forward and backward differences in space depending on the direction of propagation of information, i.e.","category":"page"},{"location":"generated/theory/linear_transport_eqs/","page":"Linear transport equations","title":"Linear transport equations","text":"    fracU_j^n+1 - U_j^nDelta t\n    + fraca^+ (U_j^n - U_j-1^n)Delta x + fraca^- (U_j+1^n - U_j^n)Delta x = 0","category":"page"},{"location":"generated/theory/linear_transport_eqs/","page":"Linear transport equations","title":"Linear transport equations","text":"Information is \"carried by the wind\", hence the name upwind. The above equations can be written as","category":"page"},{"location":"generated/theory/linear_transport_eqs/","page":"Linear transport equations","title":"Linear transport equations","text":"    fracU_j^n+1 - U_j^nDelta t + fraca (U_j+1^n - U_j-1^n)2 Delta x\n    = fracabsa2Delta x (U_j+1^n - 2 U_j^n + U_j-1^n)","category":"page"},{"location":"generated/theory/linear_transport_eqs/","page":"Linear transport equations","title":"Linear transport equations","text":"We see we have the central difference scheme with a diffusion term; the right hand side approximates fracDelta x absa2 U_xx. Thus, it adds numerical viscosity to the unstable central difference scheme, which will play a crucial role later.","category":"page"},{"location":"generated/theory/linear_transport_eqs/","page":"Linear transport equations","title":"Linear transport equations","text":"Now, we do the same numerical experiment as previously with a = 1:","category":"page"},{"location":"generated/theory/linear_transport_eqs/","page":"Linear transport equations","title":"Linear transport equations","text":"using homemade_conslaws.upwind\n\nx_L, x_R = 0, 1\nT = 1\nN = 50\ndx = (x_R - x_L)/N\ndt = 1.3 * dx\nx = x_L + dx:dx:x_R\nt = dt:dt:T\nM = length(t)\n\nA, f = upwind.upwind_scheme(N, M, dx, dt, u0)\nU = A \\ f\n\nViz.animate_solution(reshape(U, N, M),\n                    (x, t) -> u0(x-t),\n                    x, t)","category":"page"},{"location":"generated/theory/linear_transport_eqs/","page":"Linear transport equations","title":"Linear transport equations","text":"(Image: )","category":"page"},{"location":"generated/theory/linear_transport_eqs/","page":"Linear transport equations","title":"Linear transport equations","text":"x_L, x_R = 0, 1\nT = 1\nN = 100\ndx = (x_R - x_L)/N\ndt = 0.7 * dx\nx = x_L+dx:dx:x_R\nt = dt:dt:T\nM = length(t)\n\nA, f = upwind.upwind_scheme(N, M, dx, dt, u0)\nU = A \\ f\n\nViz.animate_solution(reshape(U, N, M),\n                     (x, t) -> u0(x-t),\n                     x, t)","category":"page"},{"location":"generated/theory/linear_transport_eqs/","page":"Linear transport equations","title":"Linear transport equations","text":"(Image: )","category":"page"},{"location":"generated/theory/linear_transport_eqs/","page":"Linear transport equations","title":"Linear transport equations","text":"We see that stability depends on the relation fracDelta tDelta x. ","category":"page"},{"location":"generated/theory/linear_transport_eqs/","page":"Linear transport equations","title":"Linear transport equations","text":"lemma: Lemma\nIf the mesh parameters satisfy the condition    fracabsa Delta tDelta x leq 1then the upwind solution satisfies the estimate    E^n+1 leq E^nso the scheme is conditionally stable.proof: Proof\nStart similarly to the proof of the unconditional unstability of the central difference scheme and use the mentioned identity several times.","category":"page"},{"location":"generated/theory/linear_transport_eqs/","page":"Linear transport equations","title":"Linear transport equations","text":"The above condition is called A CFL condition. We also have L^1 and L^infty stability:","category":"page"},{"location":"generated/theory/linear_transport_eqs/","page":"Linear transport equations","title":"Linear transport equations","text":"lemma: Lemma\nAssume the above CFL condition holds. Then, the solutions of the upwind scheme satisfy    normU^n+1_L^1 le normU^n_L^1 quad\n    normU^n+1_L^infty le normU^n_L^inftyproof: Proof\nThe first inequality follows from U_j^n+1being a convex combination of U_j-1^n, U_j^n and U_j+1^n.","category":"page"},{"location":"generated/theory/linear_transport_eqs/#Discontinuous-initial-data","page":"Linear transport equations","title":"Discontinuous initial data","text":"","category":"section"},{"location":"generated/theory/linear_transport_eqs/","page":"Linear transport equations","title":"Linear transport equations","text":"We now consider the transport equation with a = 1 in the domain 0 1 with initial data","category":"page"},{"location":"generated/theory/linear_transport_eqs/","page":"Linear transport equations","title":"Linear transport equations","text":"    U_0(x) = begincases\n        2  x  05 \n        1  x geq 05\n    endcases","category":"page"},{"location":"generated/theory/linear_transport_eqs/","page":"Linear transport equations","title":"Linear transport equations","text":"This yields the discontinuous solution U(x t) = U_0(x - t). We use the upwind scheme with N = 50 and N = 200 mesh points.","category":"page"},{"location":"generated/theory/linear_transport_eqs/","page":"Linear transport equations","title":"Linear transport equations","text":"x_L, x_R = 0, 1\nT = 0.25\nN = 50\ndx = (x_R - x_L)/N\ndt = 0.9 * dx\n\nx = x_L+dx:dx:x_R\nt = dt:dt:T\nM = length(t)\n\nu0(x) = x .< 0.5 ? 2 : 1\n\nA, f = upwind.upwind_scheme(N, M, dx, dt, u0)\nU = A \\ f\n\nViz.animate_solution(reshape(U, N, M),\n                     (x, t) -> u0(x-t),\n                     x, t)","category":"page"},{"location":"generated/theory/linear_transport_eqs/","page":"Linear transport equations","title":"Linear transport equations","text":"(Image: )","category":"page"},{"location":"generated/theory/linear_transport_eqs/","page":"Linear transport equations","title":"Linear transport equations","text":"x_L, x_R = 0, 1\nT = 0.25\nN = 200\ndx = (x_R - x_L)/N\ndt = 0.9 * dx\n\nx = x_L+dx:dx:x_R\nt = dt:dt:T\nM = length(t)\n\nA, f = upwind.upwind_scheme(N, M, dx, dt, u0)\nU = A \\ f\n\nViz.animate_solution(reshape(U, N, M),\n                     (x, t) -> u0(x-t),\n                     x, t)","category":"page"},{"location":"generated/theory/linear_transport_eqs/","page":"Linear transport equations","title":"Linear transport equations","text":"(Image: )","category":"page"},{"location":"#Homemade-Conslaws","page":"Conservation Laws","title":"Homemade Conslaws","text":"","category":"section"},{"location":"#Conservation-laws","page":"Conservation Laws","title":"Conservation laws","text":"","category":"section"},{"location":"","page":"Conservation Laws","title":"Conservation Laws","text":"Let boldsymbol U be a quantity defined on a domain Omega subset R^n. For any subdomain omega subset Omega, the temporal rate of change of bm U is equal to the amount of bm U created or destroyed and the flux going through the boundary partial omega. It can be described mathematically by","category":"page"},{"location":"","page":"Conservation Laws","title":"Conservation Laws","text":"    dvt int_omega bm U dd x\n    = - int_partial omega bm F cdot bm nu dd s\n    + int_omega bm S dd x","category":"page"},{"location":"","page":"Conservation Laws","title":"Conservation Laws","text":"where bm F is the flux and bm S is the source term. By the Gauss divergence theorem, we can rewrite this as","category":"page"},{"location":"","page":"Conservation Laws","title":"Conservation Laws","text":"    dvt int_omega bm U dd x\n    + int_omega div bm F dd x\n    = int_omega bm S dd x","category":"page"},{"location":"","page":"Conservation Laws","title":"Conservation Laws","text":"Since this equation holds for all subdomains omega subseteq Omega, we can write","category":"page"},{"location":"","page":"Conservation Laws","title":"Conservation Laws","text":"    bm U_t + div bm F = bm S quad text in  Omega times R_+","category":"page"},{"location":"","page":"Conservation Laws","title":"Conservation Laws","text":"We call this equation a conservation law.","category":"page"}]
}