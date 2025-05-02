using Gridap
using GridapEmbedded
using STLCutters
using LinearAlgebra
using Plots


# Manufactured solution
u_ex(x) = VectorValue(2*x[1] + cos(2*π*x[2]), -2 * x[2] + sin(2*π*x[1] ))
p_ex(x) = sin(2*π*x[1])

# compute pde data
# viscosity ν
# change so that ν is a function of u ? or a function of x? I think a function of u...
ν(u) = u^2      #for eksempel... Kan være man må ha noe annet her?
f(x) = -divergence(ν ⋅ ε(u_ex))(x) + ∇(p_ex)(x)
g(x) = tr(∇(u_ex)(x))

# Dirichlet data
ud(x) = u_ex(x)