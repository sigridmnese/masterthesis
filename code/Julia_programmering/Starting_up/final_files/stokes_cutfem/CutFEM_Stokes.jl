using Gridap
using GridapEmbedded
using STLCutters
using LinearAlgebra
using Plots
using LineSearches: BackTracking
import Random
using Logging
using LoggingExtras

include("C:\\Users\\Sigri\\Documents\\Master\\report\\code\\Julia_programmering\\Starting_up\\utils\\utils.jl")
include("C:\\Users\\Sigri\\Documents\\Master\\report\\code\\Julia_programmering\\Starting_up\\workfiles_stokes\\testing_non-linear_stokes.jl")

"""CutFEM with weak boundary imposition using Nitsches method. Symmetric stokes equations with Navier boundary conditions on entire boundary. """
# Denne er med konstant viskositet, og med nitsche implementasjon av stokes med navier BC.
# Defining constants
nu0 =  100   # klarer ikke mindre enn 0.01 og klarer heller ikke større enn 100 (størrelsesordener)
r = 4/3
A = nu0
ϵ_0 = 1e-7       

# Defining manufactured solutions
u_exact(x) =  VectorValue(2*x[1] + cos(2*π*x[2]), -2*x[2] + sin(2*π*x[1]))#VectorValue(2*x[1] + exp(x[1]/2) * cos(2*π*x[2]), -2*x[2] + exp(x[2]/2) * sin(2*π*x[1])) #bytte til sin/cos-uttrykk  VectorValue(-x[2], x[1])
p_exact(x) = sin(2*π*x[1])*cos(2*π*x[2])


viskositet(εu) = nu0#^(1-r)*(ϵ_0^2 + norm(εu)^2)^((r-2)/2)
dviskositet(εu, εdu) = 0#nu0^(1-r) *(ϵ_0^2 + εu ⋅ εu)^((r-4)/2)*(εu⊙εdu) *(r-2)

flux(εu) = 2 * nu0 * εu         # det kan være jeg ikke skal gange med 2 her...?
dflux(εu, εdu) = flux∘(εdu) # eller her for øvrig

f(x) =  -divergence(flux∘ε(u_exact))(x) + ∇(p_exact)(x)
ud(x) = u_exact(x)
g(x) = tr(ε(u_exact)(x))

# solver parameters
n = 32
stabilize = true
δ = 0 
save = true
calc_condition = false
order = 2
geometry = "heart"
βu0 = 0.001
βp0 = 1  # får ok resultater når denne er 1 også??
β_1 = 1
β_2 = 1
β_3 = 0.1
γ=1
nu = 1  
γu1 = 1
γu2 = 1.0
γp = 0.1

function Stokes_CutFEM_solver(;n, u_exact, p_exact, f, g, ud, order, geometry, βu0, γu1, γu2, γp, βp0, nu, stabilize, δ, save = false, calc_condition = false)
      """
      Using a stabilized Nitsche ficticious domain method as decribed by Massing and Larson, Logg and Rognes. Using P2-P1 Taylor-Hood elements.  
      n: number of grid elements. Powers of 2 for simplicity and convergence estimates.
      u_exact: exact solution for method of manufactured solutions
      order: order of polynomial degree. 
      f: lhs for first term, -Δ u_ex + ∇p = f
      g: lhs for second term u = g
      geometry: optional between "Circle", "Flower", "Heart", "Glacier".
      stabilize: wheather to add the stabilization term or not
      δ: perturbation of cut
      """
      # Define background mesh
      partition = (n, n)
      dim = length(partition)
      a = 1.2
      pmin = Point(-a + δ, -a + δ)
      pmax = Point(a + δ, a + δ)
      bgmodel = CartesianDiscreteModel(pmin,pmax,partition)
      # mesh size
      h = (pmax-pmin)[1]/partition[1]
  
      # defining ghost penalty constants
      βu = βu0 *nu/(h^2)                # disse er de samme som Sigmund bruker...
      βp = βp0/h
  
      geo = create_geometry(geometry, n)
      # Define active and physical mesh
      cutgeo = cut(bgmodel,geo)
      cutgeo_facets = cut_facets(bgmodel,geo)
      Ω_bg = Triangulation(bgmodel)
      Ω_act = Triangulation(cutgeo, ACTIVE)
      Ω = Triangulation(cutgeo, PHYSICAL)
  
      # Embedded boundary
      # Dirichlet conditions on u
      Γd = EmbeddedBoundary(cutgeo)
      n_Γd = get_normal_vector(Γd)
  
      # Get ghost penalty facets
      Fg = GhostSkeleton(cutgeo)
      n_Fg = get_normal_vector(Fg)
  
      # Define measures
      degree = 2*order
      dΩ = Measure(Ω,degree)
      dΓd = Measure(Γd, degree)
      dFg = Measure(Fg, degree)
  
      # Define function spaces 
      reffe_u  = ReferenceFE(lagrangian,VectorValue{dim, Float64},order)
      reffe_p = ReferenceFE(lagrangian,Float64, order - 1)
  
      V = TestFESpace(Ω_act, reffe_u,  conformity=:H1)
      Q = TestFESpace(Ω_act, reffe_p, conformity=:H1, constraint=:zeromean)
  
      U = TrialFESpace(V)
      P = TrialFESpace(Q)
  
      X = MultiFieldFESpace([U, P])
      Y = MultiFieldFESpace([V, Q])
  
        γ = 10*2*2
        mu = nu0*γ/h
        #βu = 10.0*order^2
        γu1 = 1.0
        γu2 = 1.0
        γp = 0.1
        #βp = 0.1
    | # Weak formulation components
      a(u, v) = (∫(2* nu0 * ε(u) ⊙ ε(v))dΩ + ∫(-2*nu0*((n_Γd ⋅ ε(u)) ⋅ v) + (-2*nu0*(n_Γd ⋅ ε(v)) ⋅ u)+(mu * u ⋅ v))dΓd )
      b(v, p) = (∫(-1*(∇ ⋅ v*p))dΩ
                  + ∫((n_Γd ⋅ v) * p)dΓd)#
                  
    gu(u,v) = ( ∫( (γu1 * h * nu0)*jump(n_Fg⋅∇(u))⋅jump(n_Fg⋅∇(v)) )dFg 
            +  
               ∫( (γu2 * h^3 * nu0)*jump_nn(u,n_Fg)⋅jump_nn(v,n_Fg) )dFg
)

    gp(p,q) = ∫( (γp*h^3/nu0)*jump(n_Fg⋅∇(p))*jump(n_Fg⋅∇(q)) )dFg
      
  
      l2_norm(u) = (sum( ∫( u ⋅ u )*dΩ ))^(1/2)
      h1_semi(u) = sum(∫(∇(u) ⊙ ∇(u))*dΩ)^(1/2)
      
      l1((v, q)) = ∫(f ⋅ v)dΩ
      l2(v) = ∫(-2* nu0*(n_Γd ⋅ ε(v)) ⋅ ud + mu ⋅ v ⋅ ud )dΓd
      l3(q) = ∫((n_Γd ⋅ ud) *q)dΓd
  
      if stabilize
          A((u,p),(v,q)) =(a(u, v) + b(v, p) - b(u, q) 
          + gu(u,v)
          + gp(p, q)
          )
          L((v, q)) = l1((v, q)) + l2(v) - l3(q)
          op = AffineFEOperator(A,L,X,Y)
          uh, ph = solve(op)
      else
          B((u,p),(v,q)) = a(u,v) + b(v, p) - b(u, q) 
          M((v,q)) = l1((v, q)) + l2(v) - l3(q)
       # Linear forms
          op = AffineFEOperator(B,M,X,Y)
          uh, ph = solve(op)
      end
  
    errp = p_exact - ph
    erru = u_exact - uh
    
    # condition number
    if calc_condition
      condition_numb= cond(Array(get_matrix(op)),2)   # kanskje bruke infinitynormen istedenfor
    else
      condition_numb = 1
    end
  
    if save
      writevtk(Ω, "C:\\Users\\Sigri\\Documents\\Master\\report\\results\\stokes\\$n $geometry $order.vtu", cellfields=["u_ex" => u_exact, "uh"=>uh, "erru"=> erru, "p_ex" => p_exact, "ph"=>ph, "errp"=> errp, "nablau" => ∇(u_exact)]) #, "erru" => erru]) 
    end
    return uh, u_exact, erru, l2_norm(uh - u_exact), h1_semi(uh - u_exact), ph, p_exact, errp, l2_norm(ph - p_exact), h1_semi(ph - p_exact), condition_numb, Ω_act
end


solver = Stokes_CutFEM_solver
# uh, u_exact, erru, ul2_norm, uh1_semi, ph, p_exact, errp, pl2_norm, ph1_semi, condition_numb, Ω_act = solver(;n, u_exact, p_exact, f, g, ud, order, geometry, βu0, γu1, γu2, γp, βp0, nu, stabilize, δ, save, calc_condition)
# order = 2
# numb_it = 5
# uarr_l2_stab, uarr_h1_stab, parr_l2_stab, parr_h1_stab, h = convergence_stokes_weird_domain(;numb_it, u_exact, p_exact, f, g, ud, order, geometry, solver, δ, βu0, γu1, γu2, γp, βp0, nu, stabilize, save)

# stabilize = false
# uarr_l2_nostab, uarr_h1_nostab, parr_l2_nostab, parr_h1_nostab, h = convergence_stokes(;numb_it, u_exact, p_exact, f, g, ud, order, geometry, solver, δ, βu0, γu1, γu2, γp, βp0, nu, stabilize, save)

# # # # ########## velocity convergence plot:
# plot_convergence_u(
#     uarr_l2_stab[2:end], uarr_h1_stab[2:end], h[2:end];
#     uarr_l2_nostab = uarr_l2_nostab[2:end], uarr_h1_nostab = uarr_h1_nostab[2:end],
#     title_str="Convergence Stokes with Dirichlet BC"
# )

# plot_convergence_p(
#     parr_l2_stab[2:end], parr_h1_stab[2:end], h[2:end];
#     parr_l2_nostab=parr_l2_nostab[2:end],
#     parr_h1_nostab=parr_h1_nostab[2:end],
#     title_str="Convergence Stokes with Dirichlet BC"
# )

# #print eoc-verdier:
# print_eoc_latex_combined_only_eoc(h;
#     uarr_l2_stab = uarr_l2_stab,
#     uarr_h1_stab = uarr_h1_stab,
#     uarr_l2_nostab = uarr_l2_nostab,
#     uarr_h1_nostab = uarr_h1_nostab,
#     start = 1
# )

# print_eoc_latex_combined_pressure_only_eoc(h;
#     parr_l2_stab = parr_l2_stab,
#     parr_h1_stab = parr_h1_stab,
#     parr_l2_nostab = parr_l2_nostab,
#     parr_h1_nostab = parr_h1_nostab,
#     start = 1
# )

n = 16                # øke denne?
M = 2000               #full kjøring med M = 2000 med 2000 så kjører det nok i 2 timer. 
# βu0 = 1
# γu1 = 1
# γu2 = 1
# γp = 0.1
# βp0 = 0.1
stabilize = true
save = false
calc_condition = false
geometry = "circle"

arr_δ, arr_l2u, arr_h1u, arr_l2p, arr_h1p, arr_cond = sensitivity_stokes(;n, M, u_exact, p_exact, f, g, ud, order, geometry, solver, βu0, γu1, γu2, γp, βp0, nu, stabilize, save)
stabilize = false
start = 1
save = true
arr_δ_nostab, arr_l2u_nostab, arr_h1u_nostab, arr_l2p_nostab, arr_h1p_nostab, arr_cond_nostab = sensitivity_stokes(;n, M, u_exact, p_exact, f, g, ud, order, geometry, solver, βu0, γu1, γu2, γp, βp0, nu, stabilize, save)

plot_sensitivity_velocity(
    arr_δ, arr_l2u, arr_h1u,
    arr_δ_nostab, arr_l2u_nostab, arr_h1u_nostab;
    start = start,
    title_str = "Sensitivity Stokes with Dirichlet BC"
)

plot_sensitivity_pressure(
    arr_δ, arr_l2p, arr_h1p,
    arr_δ_nostab, arr_l2p_nostab, arr_h1p_nostab;
    start = start,
    title_str = "Sensitivity Stokes with Dirichlet BC"
)

plot_sensitivity_poisson(
    arr_δ, arr_l2u, arr_h1u, arr_δ_nostab, arr_l2u_nostab, arr_h1u_nostab;
    marker_l2_stab = :circle,
    marker_l2_nostab = :star5,
    marker_h1_stab = :circle,
    marker_h1_nostab = :star5,
    markstep = 50,
    title_str = "Sensitivity Stokes with Dirichlet BC",
)

plot_sensitivity_poisson(
    arr_δ, arr_l2p, arr_h1p, arr_δ_nostab, arr_l2p_nostab, arr_h1p_nostab;
    marker_l2_stab = :circle,
    marker_l2_nostab = :star5,
    marker_h1_stab = :circle,
    marker_h1_nostab = :star5,
    markstep = 50,
    title_str = "Sensitivity p-Stokes, pressure",
)