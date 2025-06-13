#TODO: Debug this code
#TODO: 1. Look at pure elasticity problem/vector Laplace

# %% 
using Gridap
# import Gridap: ∇
using GridapEmbedded
# using STLCutters

# Manufactured solution
u_ex(x) = VectorValue(2*x[1]+cos(2*π*x[2]),-2*x[2]+sin(2*π*x[1]))
p_ex(x) = sin(2*π*x[1])

# %% Compute PDE data
# Viscosity
μ = 1.0 
# TODO: This alternative does not work
# f(x) = -∇⋅(strain_u)(x) + ∇(p_ex)(x)
f(x) = -divergence(ε(u_ex))(x) + ∇(p_ex)(x)
g(x) = tr(ε(u_ex)(x))

# Can use Laplace if div(u) = 0 or Stokes with Laplace operator is used
# f(x) = -0.5*Δ(u_ex)(x) + ∇(p_ex)(x)
# g(x) = tr(∇(u_ex)(x))

# Dirichlet data
u_D(x) = u_ex(x)

# %% Define geometry
R = 0.512

# Disc center
p0 = Point(0.0, 0.0)
# Take complement of disk
# geo = !disk(R, x0=p0)
geo = disk(R, x0=p0)    #kan invertere hva som er bgmodel og hva som er physical mesh ved å fjerne utropstegn

# Define background mesh
n = 60
partition = (n, n)
dim = length(partition)
pmin = p0+Point(-1.0, 0.0)
pmax = p0+Point(1.0, 1.0)
bgmodel = CartesianDiscreteModel(pmin,pmax,partition)

# Define active and physical meshes
cutgeo = cut(bgmodel, geo)

# Needed to extract (potentially cut!) boundary facets from active mesh
cutgeo_facets = cut_facets(bgmodel,geo)

Ω_bg = Triangulation(bgmodel)
Ω_act = Triangulation(cutgeo, ACTIVE)
Ω = Triangulation(cutgeo, PHYSICAL)

# Embedded boundary
# Has Dirichlet for u and Neumann for pF
Γd = EmbeddedBoundary(cutgeo)
n_Γd = get_normal_vector(Γd)

# Outer mesh boundary 
# Has Dirichlet for pF and no-stress for u
Γs = BoundaryTriangulation(cutgeo_facets, PHYSICAL_IN)
n_Γs = get_normal_vector(Γs)

# Get ghost penalty facets
Fg = GhostSkeleton(cutgeo)
n_Fg = get_normal_vector(Fg)

# All interior facets of active mesh
Fi = SkeletonTriangulation(Ω_act)
n_Fi = get_normal_vector(Fi)

# Function spaces orders
# order_u  = 1
# order_pT = order_u
order_u  = 2
order_p = order_u - 1

# Define measures
degree = 2*order_u
dΩ = Measure(Ω,degree)
dΓd = Measure(Γd, degree)
dΓs = Measure(Γs, degree)
dFg = Measure(Fg, degree)
dFi = Measure(Fi, degree)

# Write out domains 
writevtk(Ω_act, "mesh_act")
writevtk(Ω_bg, "mesh_bg")
writevtk(Ω, "mesh")
writevtk(Γd, "surface_gamma_d")
writevtk(Γs, "surface_gamma_s")
writevtk(Fg, "ghost_facets")

# %% Define weak formulation
# Define function spaces
reffe_u  = ReferenceFE(lagrangian,VectorValue{dim, Float64},order_u)
reffe_p = ReferenceFE(lagrangian,Float64, order_p)

V = TestFESpace(Ω_act, reffe_u,  conformity=:H1)
Q = TestFESpace(Ω_act, reffe_p, conformity=:H1, constraint=:zeromean)

U = TrialFESpace(V)
P = TrialFESpace(Q)

X = MultiFieldFESpace([U, P])
Y = MultiFieldFESpace([V, Q])

# Physical parameters (potentially rescaled after time-step discretization)

# Nitsche and Ghost penalty stabilization parameter
βu = 10.0*order_u^2
γu1 = 1.0
γu2 = 1.0
γp = 0.1
βp = 0.1

# Weak formulation

# mesh size
h = (pmax-pmin)[1]/partition[1]

# TODO: For now we impose Dirichlet boundary conditions on the whole boundary, so need to change b.c. on Γs
# Bilinear forms
# a(u,v) =  ( ∫( ∇(u)⊙∇(v))dΩ 
        #  + ∫( -(n_Γd⋅∇(u))⋅v - u⋅(n_Γd⋅∇(v)) + (βu/h)*u⋅v)dΓd
a(u,v) =  ( ∫( ε(u)⊙ε(v))dΩ 
         + ∫( -(n_Γd⋅ε(u))⋅v -(n_Γd⋅ε(v))⋅u + βu/h*u⋅v)dΓd
         + ∫( -(n_Γs⋅ε(u))⋅v -(n_Γs⋅ε(v))⋅u + βu/h*u⋅v)dΓs
)

b(v, p) = ( ∫(-1*(∇⋅v*p))dΩ 
            + ∫((n_Γd⋅v)*p)dΓd
            + ∫(v⋅n_Γs*p)dΓs
)

# Define ghost penalty
function jump_nn(u,n)
  return ( n.plus⋅ (n.plus⋅∇∇(u).plus) -  n.minus ⋅(n.minus ⋅ ∇∇(u).minus) )
end

g_u(u,v) = ( ∫( (γu1*h)*jump(n_Fg⋅∇(u))⋅jump(n_Fg⋅∇(v)) )dFg 
            +  
               ∫( (γu2*h^3)*jump_nn(u,n_Fg)⋅jump_nn(v,n_Fg) )dFg
)

g_p(p,q) = ∫( (γp*h^3)*jump(n_Fg⋅∇(p))*jump(n_Fg⋅∇(q)) )dFg
s_p(p,q) = ∫( (βp*h^3)*jump(n_Fi⋅∇(p))*jump(n_Fi⋅∇(q)) )dFi

# Totatl bilinear form
Ah((u,p),(v,q)) = a(u,v) + g_u(u,v) + b(v,p) + b(u,q) - g_p(p,q)
# Ah((u,p),(v,q)) = a(u,v) + g_u(u,v) + b(v,p) + b(u,q) - g_p(p,q) - s_p(p,q)

# Linear forms
lh((v,q)) = ( ∫(f⋅v - g⋅q)dΩ
# lh((v,q)) = ( ∫(fdiv⋅v - g⋅q)dΩ
            #  + ∫( -1*(n_Γd⋅∇(v))⋅u_D + βu/h*(u_D⋅v) + u_D⋅n_Γd*q)dΓd 
             + ∫( -1*(n_Γd⋅ε(v))⋅u_D + βu/h*(u_D⋅v) + u_D⋅n_Γd*q)dΓd 
             + ∫( -1*(n_Γs⋅ε(v))⋅u_D + βu/h*(u_D⋅v) + u_D⋅n_Γs*q)dΓs
)


# %% Solve and write out results
# Solve
op = AffineFEOperator(Ah,lh,X,Y)
uh, ph = solve(op)

outputfile = "stokes"
if outputfile !== nothing
  writevtk(Ω,"$(outputfile)",cellfields=["uh"=>uh,
                                         "ph"=>ph, 
                                         "u_ex"=> u_ex, 
                                         "p_ex"=> p_ex,
                                         "f"=> f,
                                         "g"=> g
                                         ])
end