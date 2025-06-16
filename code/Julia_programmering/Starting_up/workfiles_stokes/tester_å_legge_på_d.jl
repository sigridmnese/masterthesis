# %% Load required modules
using Gridap
# import Gridap: ∇
using GridapEmbedded
using STLCutters
using Test

# Manufactured solution
u_ex(x) = VectorValue(2*x[1]+cos(2*π*x[2]),-2*x[2]+sin(2*π*x[1]), 0)
p_ex(x) = sin(2*π*x[1])

# Compute PDE data
# Viscosity
μ = 1.0 
f(x) = -divergence(ε(u_ex))(x) + ∇(p_ex)(x)
g(x) = tr(ε(u_ex)(x))

# Can use Laplace if div(u) = 0 or Stokes with Laplace operator is used
# f(x) = -0.5*Δ(u_ex)(x) + ∇(p_ex)(x)
# g(x) = tr(∇(u_ex)(x))

# Dirichlet data
u_D(x) = u_ex(x)

# %% Define background mesh and geometry
# Disc center and radius
# R = 0.753
p0 = Point(0.0, 0.0, 0.0)
# geo = sphere(R, x0=p0)

# Take complement of the ball
# geo = !sphere(R, x0=p0)

# Define background mesh
n = 15
partition = (n, n, n)
dim = length(partition)
pmin = p0+Point(-1.0, -1.0, 0.0)
# pmin = p0+Point(-1.0, 0.0, 0.0)
pmax = p0+Point(1.0, 1.0, 1.0)
bgmodel = CartesianDiscreteModel(pmin, pmax, partition)
Ω_bg = Triangulation(bgmodel)
writevtk(Ω_bg, "mesh_bg")
h = 1/n
γn = 0.05  # hentet fra Josefin sin artikkel.
γt = 0.05  # hentet fra Josefin sin artikkel.
e = 0        # e er det som er epsilon i artikkelen til Josefin. Denne er 0 for unit square. 

# %% Test STL geo
geo_stl = STLGeometry("sphere.stl")
@test check_requisites(geo_stl,bgmodel)

# %% Define actual computational domain
cutgeo = cut(bgmodel, geo_stl)

# Needed to extract (potentially cut!) boundary facets from active mesh
cutgeo_facets = cut_facets(bgmodel,geo_stl)

Ω_act = Triangulation(cutgeo, ACTIVE)
Ω = Triangulation(cutgeo, PHYSICAL)

# Embedded boundary
# Has Dirichlet for u and Neumann for pF
Γd = EmbeddedBoundary(cutgeo)
n_Γd = get_normal_vector(Γd)
# lag tangent‐felt på Γd:
t_Γd = CellField(x -> VectorValue(-n_Γd(x)[2],
                                   n_Γd(x)[1]),Γd)

I₂ = one(TensorValue{3,3,Float64})

# Projeksjons‐operatorene:
# defining the Navier operators
Pn(n) = n ⊗ n        # normal projection operator
Pt(n) = I₂ - Pn(n)     # tangential projection operator

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

# mesh size
h = (pmax-pmin)[1]/partition[1]

# TODO: For now we impose Dirichlet boundary conditions on the whole boundary, so need to change b.c. on Γs
# Bilinear forms
# a(u,v) =  ( ∫( ∇(u)⊙∇(v))dΩ 
        #  + ∫( -(n_Γd⋅∇(u))⋅v - u⋅(n_Γd⋅∇(v)) + (βu/h)*u⋅v)dΓd
a(u,v) =  ( ∫( ε(u)⊙ε(v))dΩ 
         #+ ∫( -(n_Γd⋅ε(u))⋅v -(n_Γd⋅ε(v))⋅u + βu/h*u⋅v)dΓd
         + ∫( -(n_Γs⋅ε(u))⋅v -(n_Γs⋅ε(v))⋅u + βu/h*u⋅v)dΓs
)

a0(u, v) = ∫( ε(v)⊙(flux(ε(u))))dΩ              # sjekket med Hanna
a1(u, v) = ∫(((n_Γd ⋅ (flux(ε(u)))) ⋅ v))dΓd     # ok
a2(u, v) = ∫((Pn(n_Γd)⋅ u) ⋅ (n_Γd ⋅ (flux(ε(v)))))dΓd      # OK
a3(u, v) = ∫((2* nu0/(γn*h)*(n_Γd ⋅ u)) ⋅ (n_Γd ⋅ v)  )dΓd # linje 2, ledd 1                            OK
a4(u, v) = ∫((e/(e + γt *h) * ((Pt(n_Γd) ⋅  (n_Γd ⋅ (flux(ε(u)))))) ⋅ v))dΓd # linje 2, ledd 2            # OK
a5(u, v) = ∫((nu0/(e + γt*h)*(t_Γd ⋅ u)) ⋅ ( t_Γd ⋅ v))dΓd #linje 2, ledd 3                               # OK     
a6(u, v) = ∫((e*γt * h/(e + γt*h) * (Pt(n_Γd) ⋅  (n_Γd ⋅ (flux∘ε(u)))) ) ⋅ (2*n_Γd ⋅ ε(v)))dΓd            # OK
a7(u, v) = ∫((nu0 * γt * h / (e + γt *h) * (Pt(n_Γd)⋅ u)) ⋅ (2* n_Γd ⋅ ε(v)))dΓd                          # OK
ad(u, v) = a0(u, v) - a1(u, v) - a2(u, v) + a3(u, v) + a4(u, v) + a5(u, v) - a6(u, v) - a7(u, v)

ld(v, q) = (∫(-1*(n_Γd ⋅ ud) ⋅ ( n_Γd ⋅ (n_Γd ⋅ flux(ε(v)))))dΓd #OK
    + ∫(((2*nu0)/(γn *h) * (n_Γd ⋅ ud)) ⋅ (n_Γd ⋅ v) )dΓd #OK
    - ∫((n_Γd ⋅ ud) ⋅ q)dΓd  #OK
    + ∫((nu0/(e + γt * h) ⋅ (t_Γd ⋅ ud) )⋅ (t_Γd ⋅ v) )dΓd #OK
    - ∫((nu0 * γt * h/(e + γt * h) * Pt(n_Γd) ⋅ ud) ⋅ (2*n_Γd ⋅ ε(v)))dΓd # OK
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
Ah((u,p),(v,q)) = a(u,v) +ad(u, v) + g_u(u,v) + b(v,p) + b(u,q) - g_p(p,q)
# Ah((u,p),(v,q)) = a(u,v) + g_u(u,v) + b(v,p) + b(u,q) - g_p(p,q) - s_p(p,q)

# Linear forms
lh((v,q)) = ( ∫(f⋅v - g⋅q)dΩ
# lh((v,q)) = ( ∫(fdiv⋅v - g⋅q)dΩ
            #  + ∫( -1*(n_Γd⋅∇(v))⋅u_D + βu/h*(u_D⋅v) + u_D⋅n_Γd*q)dΓd 
             #+ ∫( -1*(n_Γd⋅ε(v))⋅u_D + βu/h*(u_D⋅v) + u_D⋅n_Γd*q)dΓd 
             + ∫( -1*(n_Γs⋅ε(v))⋅u_D + βu/h*(u_D⋅v) + u_D⋅n_Γs*q)dΓs
             +ld(v, q)
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