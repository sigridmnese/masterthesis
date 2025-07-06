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

"""Fitted FEM with weak boundary imposition using Nitsches method. Symmetric p-stokes equations, using Newtons method to linearize
denne funker for nu0 = 100, 1 og 0.01!!"""
#""" Lagt inn konvergensplott for denne i teams"""

# Defining constants
   # test ulike verdier for nu0!!! Da synes jeg du skal kjøre konvergenstester men også bare generelt se på hva som skjer med løsningen når du endrer nu0.
r = 4/3
A = 2^(1/(1-r))
ϵ_0 = 1e-6
nu0 =  A

# Defining manufactured solutions
u_exact(x) =  VectorValue(2*x[1] + cos(2*π*x[2]), -2*x[2] + sin(2*π*x[1]))#VectorValue(2*x[1] + exp(x[1]/2) * cos(2*π*x[2]), -2*x[2] + exp(x[2]/2) * sin(2*π*x[1])) #bytte til sin/cos-uttrykk  VectorValue(-x[2], x[1])
p_exact(x) = sin(2*π*x[1])*cos(2*π*x[2])

# Defining the problem 
flux(εu) = nu0^(1-r)*(ϵ_0^2 + tr(εu' ⋅ εu))^((r-2)/2)⋅εu 
dflux(εu, εdu)=(r-2)*nu0^(1-r)*(ϵ_0^2 + tr(εu' ⋅ εu))^((r-4)/2)*(εu⊙εdu) ⋅ εu + nu0^(1-r)*(ϵ_0^2 + tr(εu' ⋅ εu))^((r-2)/2)*εdu
viskositet(εu) = nu0^(1-r)*(ϵ_0^2 + tr(εu' ⋅ εu))^((r-2)/2)
dviskositet(εu, εdu) = nu0^(1-r) *(ϵ_0^2 + tr(εu' ⋅ εu))^((r-4)/2)*(εu⊙εdu) *(r-2)


# viskositet(εu) = 1/2*A^(1-r) * (ϵ_0^2 + 1/2 * norm(εu)^2)^((r-2)/2)
# dviskositet(εu, εdu) =  (r-2)/4*A^(1-r) *(ϵ_0^2 + 1/2* norm(εu)^2)^((r-4)/2) * (εu⊙εdu)  # eller delt på 4?
# flux(εu) = A^(1-r) * (ϵ_0^2 + 1/2 *norm(εu)^2)^((r-2)/2) ⋅ εu
# dflux(εu, εdu)= (r-2)*A^(1-r)*(ϵ_0^2 + 1/2 * norm(εu)^2)^((r-4)/2)*(εu⊙εdu) ⋅ εu + A^(1-r)*(ϵ_0^2 + 1/2 * norm(εu)^2)^((r-2)/2)*εdu


f(x) =  -divergence(flux∘ε(u_exact))(x) + ∇(p_exact)(x)      # prøver å endre f her...
ud(x) = u_exact(x)
g(x) = -tr(ε(u_exact)(x))


# solver:
#fitted FEM stokes solver:
function pstokes_FEM_newton(;n, u_exact, p_exact, f, g, ud, order, geometry, βu0, γu1, γu2, γp, βp0, nu, stabilize, δ, save = false, calc_condition = false)
    """    
    Case A er mye bedre enn Case B
    Fitted FEM with weak boundary imposition using Nitsches method. Symmetric p-stokes equations, using Newtons method to linearize.
    Should test for different values of nu0.
    Symmetric case.
    n: number of grid elements. Powers of 2 for simplicity and convergence estimates.
    u_exact: exact solution for method of manufactured solutions
    order: order of polynomial degree. 
    f: lhs for first term, -Δ u_ex + ∇p = f
    g: lhs for second term u = g
    geometry: optional between "Circle", "Flower", "Heart", "Glacier". Will not affect anything - just added so that the function takes similar function arguments as the stokes_solver
    βu0, γu1, γu2, γp, βp0: parameters for the stokes solver. Will not affect anything - just added so that the function takes similar function arguments as the stokes_solver
    stabilize: will not affect anything - just added so that the function takes similar function arguments as the stokes_solver
    δ: perturbation of cut, will not affect anything - just added so that the function takes similar function arguments as the stokes_solver
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
    
    order = 2
    # Define function spaces 
    reffe_u  = ReferenceFE(lagrangian,VectorValue{2, Float64}, order)
    reffe_p = ReferenceFE(lagrangian,Float64, order - 1)
  
    V = TestFESpace(Ω_act, reffe_u,  conformity=:H1)
    Q = TestFESpace(Ω_act, reffe_p, conformity=:H1, constraint=:zeromean)
    Y = MultiFieldFESpace([V,Q])
  
    U = TrialFESpace(V)         # har fjernet begrensningen til grensebetingelsen ud her..
    P = TrialFESpace(Q)
    X = MultiFieldFESpace([U,P])

    ######################## kvadratish domain - men bruker nitsche implementasjon. Hvilke spaces skal jeg bruke da? #########################
    # weak formulation components   

    a(u, v) = ∫(ε(v)⊙(flux∘ε(u)))dΩ  + ∫(-((n_Γd ⋅ (flux∘ε(u))) ⋅ v) + (-(n_Γd ⋅ (flux∘ε(v))) ⋅ u) + 2*γ/h * (viskositet∘ε(u) ⋅ (u ⋅ v)))dΓd      # denne må ha et ekstra boundary term. Finn ut hvordan det ser ut. 
    b(v, p) = (∫(-1*(∇ ⋅ v*p))dΩ + ∫((n_Γd ⋅ v) * p)dΓd)   # b er den samme som før. 
    l1((u, p), (v, q)) = ∫(f ⋅ v - g ⋅ q)dΩ + ∫(-(n_Γd ⋅ (flux∘ε(v))) ⋅ ud)dΓd + ∫( 2*γ/h* (viskositet∘ε(u) ⋅ (ud ⋅ v)))dΓd
    l2(q) = ∫((n_Γd ⋅ ud) * q)dΓd

    dl1((u, du), (v, q)) = ∫( 2*γ/h*(dviskositet∘(ε(u), ε(du))⋅ (ud ⋅ v)))dΓd 
  
    # a(u, v) = ∫( ε(v)⊙(flux∘ε(u)))dΩ  + ∫(-((n_Γd ⋅ (flux∘ε(u))) ⋅ v) + (-(n_Γd ⋅ (flux∘ε(v))) ⋅ u) + 2*nu0*γ/h * (u ⋅ v))dΓd      # denne må ha et ekstra boundary term. Finn ut hvordan det ser ut. 
    # b(v, p) = (∫(-1*(∇ ⋅ v*p))dΩ + ∫((n_Γd ⋅ v) * p)dΓd)   # b er den samme som før. 
    # l1((v, q)) = ∫(f ⋅ v + g ⋅ q)dΩ
    # l2(v) = ∫(-(n_Γd ⋅ (flux∘ε(v))) ⋅ ud)dΓd    
    # l3((v, q)) = ∫( 2*nu0*γ/h* (ud ⋅ v))dΓd
    # l4(q) = ∫(-1*(n_Γd ⋅ ud) * q)dΓd
    
    da(u, du, v) = ∫( ε(v)⊙(dflux∘(ε(u), ε(du))))dΩ  + ∫(-((n_Γd ⋅ (dflux∘(ε(u), ε(du)))) ⋅ v) + (-(n_Γd ⋅ (flux∘ε(v))) ⋅ du) + (2*γ/h*(dviskositet∘(ε(u), ε(du))⋅ (u ⋅ v)))   +   (2*γ/h*(viskositet∘ε(u)⋅ (du ⋅ v))))dΓd    

    g_u(u,v) = ( ∫( (γu1*h)*jump(n_Fg⋅∇(u))⋅jump(n_Fg⋅∇(v)) )dFg 
            +  
      ∫( (γu2*h^3)*jump_nn(u,n_Fg)⋅jump_nn(v,n_Fg) )dFg
)
    g_p(p,q) = ∫( (γp*h^3)*jump(n_Fg⋅∇(p))*jump(n_Fg⋅∇(q)) )dFg

    res((u,p),(v,q)) = a(u, v) + b(v, p) - b(u, q) - l1((u,p), (v, q)) + l2(q) + g_u(u, v) + g_p(p, q)
    jac((u, p), (du, dp), (v, q)) = b(v, dp) - b(du, q) + da(u, du, v) + g_u(du, v) + g_p(dp, q)  - dl1((u, du), (v, q))
      
    op = FEOperator(res, jac, X, Y)

    # non-linear phase
    nls = NLSolver(
    show_trace=true, method=:newton, linesearch=BackTracking(), iterations=50)      #prøver å legge inn et max antall iterasjoner og en lav toleranse      
    solver = FESolver(nls)

    (uh, ph) = solve(solver, op)

    errp = p_exact - ph
    erru = u_exact - uh

    l2_norm(u) = (sum( ∫( u ⋅ u )*dΩ ))
    h1_semi(u) = sum(∫(∇(u) ⊙ ∇(u))*dΩ)
    
    # condition number
    if calc_condition
      condition_numb= cond(Array(get_matrix(op)),2)   # kanskje bruke infinitynormen istedenfor
    else
      condition_numb = 1
    end
  
    if save
        writevtk(Ω, "C:\\Users\\Sigri\\Documents\\Master\\report\\results\\stokes\\$n $geometry $order.vtu", cellfields=["u_ex" => u_exact, "uh"=>uh, "erru"=> erru, "p_ex" => p_exact, "ph"=>ph, "errp"=> errp, "nablau" => ∇(u_exact)]) #, "erru" => erru]) 
    end
    return uh, u_exact, erru, l2_norm(uh - u_exact), h1_semi(uh - u_exact), ph, p_exact, errp, l2_norm(ph - p_exact), h1_semi(ph - p_exact), condition_numb, Ω
end

# solver parameters
n = 32
stabilize = true
solver = pstokes_FEM_newton
δ = 0 
save = true
calc_condition = false
order = 2
geometry = "heart"
βu0 = 1
γu1 = 0.1
γu2 = 0.1
γp = 0.1
βp0 = 0.1
β_1 = 1
β_2 = 1
β_3 = 0.1
γ=10*2*2
nu = 1      # denne brukes ikke i p_stokes_cutfem, men sendes kun inn for at fuksjonskallet skal være likt i konvergens-funksjonen. 
uh, u_exact, erru, ul2_norm, uh1_semi, ph, p_exact, errp, pl2_norm, ph1_semi, condition_numb, Ω_act = solver(;n, u_exact, p_exact, f, g, ud, order, geometry, βu0, γu1, γu2, γp, βp0, nu, stabilize, δ, save, calc_condition)

# ################################################## convergence ##########################################################
# numb_it = 6
# uarr_l2, uarr_h1, parr_l2, parr_h1, h = convergence_stokes(;numb_it, u_exact, p_exact, f, g, ud, order, geometry, solver, δ, βu0, γu1, γu2, γp, βp0, nu, stabilize, save)

# ########## velocity convergence plot:
# plot(
#     0,
#     titlefont = 16,
#     guidefont = 14,
#     tickfont = 12
# )

# plot!(h, uarr_l2, xaxis=:log, yaxis=:log, marker=:o, lw=2, label="L2")
# plot!(h, uarr_h1, marker=:o, lw=2, label="H1")
# xlabel!("Mesh size h")
# ylabel!("Velocity error")
# title!("Convergence of p-stokes fitted FEM with Nitsche's method")

# # # ########### pressure convergence plot:
# plot(
#     0,
#     titlefont = 16,
#     guidefont = 14,
#     tickfont = 12
# )

# plot!(h, parr_l2, xaxis=:log, yaxis=:log, marker=:o, lw=2, label="L2")
# plot!(h, parr_h1, marker=:o, lw=2, label="H1")
# xlabel!("Mesh size h")
# ylabel!("Pressure error")
# title!("Convergence of p-stokes fitted FEM with Nitsche's method")
