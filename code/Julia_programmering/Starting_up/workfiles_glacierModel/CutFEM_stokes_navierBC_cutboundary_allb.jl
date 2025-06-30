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
nu0 =  1   # klarer ikke mindre enn 0.01 og klarer heller ikke større enn 100 (størrelsesordener)
r = 4/3
A = nu0
ϵ_0 = 1e-7       

# Defining manufactured solutions
u_exact(x) =  VectorValue(2*x[1] + cos(2*π*x[2]), -2*x[2] + sin(2*π*x[1]))#VectorValue(2*x[1] + exp(x[1]/2) * cos(2*π*x[2]), -2*x[2] + exp(x[2]/2) * sin(2*π*x[1])) #bytte til sin/cos-uttrykk  VectorValue(-x[2], x[1])
p_exact(x) = sin(2*π*x[1])*cos(2*π*x[2])

# Defining the problem 
#flux(εu) = nu0*(ϵ_0 + norm(εu)^2)^((r-2)/2)⋅εu 
#dflux(εdu,εu)=(r-2)*nu0*(ϵ_0 + norm(εu)^2)^((r-4)/2)*(εu⊙εdu) ⋅ εu + nu0*(ϵ_0 + norm(εu)^2)^((r-2)/2)*εdu

viskositet(εu) = nu0#^(1-r)*(ϵ_0^2 + norm(εu)^2)^((r-2)/2)
dviskositet(εu, εdu) = 0#nu0^(1-r) *(ϵ_0^2 + εu ⋅ εu)^((r-4)/2)*(εu⊙εdu) *(r-2)

# dette er hvis viskositeten er en konstant:
flux(εu) = 2 * nu0 * εu         # det kan være jeg ikke skal gange med 2 her...?
dflux(εu, εdu) = flux∘(εdu) # eller her for øvrig

# dette er hvis viskositeten ikke er en konstant:
# flux(εu) = 2 * viskositet(εu) ⋅ εu         # det kan være jeg ikke skal gange med 2 her...?
# dflux(εu, εdu) = 2 * dviskositet(εu, εdu)⋅ εu + 2 * viskositet(εu) ⋅ εdu # eller her for øvrig

f(x) =  -divergence(flux∘ε(u_exact))(x) + ∇(p_exact)(x)      # prøver å endre f her...
ud(x) = u_exact(x)
g(x) = tr(ε(u_exact)(x))
# prøver også å bytte ud til g...

# solver parameters
n = 128
stabilize = true

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
γ=1
nu = 1  

# p stokes cut FEM fungerer veldig bra, men må få lagt inn rikig viskositetsavhengighet.
function stokes_navierBC_CutFEM_cutmeshboundary(;n, u_exact, p_exact, f, g, ud, order, geometry, βu0, γu1, γu2, γp, βp0, nu, stabilize, δ, save = false, calc_condition = false)
    """
    cutFEM solver for the non-linear stokes (p-stokes) equations. Using P2-P1 Taylor-Hood elements.  
    n: number of grid elements. Powers of 2 for simplicity and convergence estimates.
    u_exact: exact solution for method of manufactured solutions
    order: order of polynomial degree. 
    f: lhs for first term, -∇ (ν ∇(u_ex) + ∇p = f
    g: lhs for second term u = g
    geometry: optional between "Circle", "Flower", "Heart", "Glacier".
    stabilize: wheather to add the stabilization term or not
    δ: perturbation of cut
    """

    γn = 0.05  # hentet fra Josefin sin artikkel.
    γt = 0.05  # hentet fra Josefin sin artikkel.
    e = 0 
    # %% Define geometry
    R = 0.512

    # Disc center
    p0 = Point(0.0, 0.0)
    # Take complement of disk
    # geo = !disk(R, x0=p0)
    geo = disk(R, x0=p0)    #kan invertere hva som er bgmodel og hva som er physical mesh ved å fjerne utropstegn

    # Define background mesh
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
    Γs = EmbeddedBoundary(cutgeo)
    n_Γs = get_normal_vector(Γs)

    # Outer mesh boundary 
    # Has Dirichlet for pF and no-stress for u
    Γb = BoundaryTriangulation(cutgeo_facets, PHYSICAL_IN)
    n_Γb = get_normal_vector(Γb)

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
    dΓb = Measure(Γb, degree)
    dΓs = Measure(Γs, degree)
    dFg = Measure(Fg, degree)
    dFi = Measure(Fi, degree)

    # %% Define weak formulation
    # Define function spaces
    reffe_u  = ReferenceFE(lagrangian,VectorValue{dim, Float64},order)
    reffe_p = ReferenceFE(lagrangian,Float64, order-1)

    V = TestFESpace(Ω_act, reffe_u,  conformity=:H1)
    Q = TestFESpace(Ω_act, reffe_p, conformity=:H1, constraint=:zeromean)

    U = TrialFESpace(V)
    P = TrialFESpace(Q)

    X = MultiFieldFESpace([U, P])
    Y = MultiFieldFESpace([V, Q])

    # Physical parameters (potentially rescaled after time-step discretization)

    # Nitsche and Ghost penalty stabilization parameter
    βu = 10.0*order_u^2
    γu1 = 0.1
    γu2 = 0.1
    γp = 0.1
    βp = 0.1

    # Weak formulation

    # mesh size
    h = (pmax-pmin)[1]/partition[1]
    
    #println(sum( ∫( p_exact) * dΩ ))
    l2_norm(u) = (sum( ∫( u ⋅ u )*dΩ ))^(1/2)
    h1_semi(u) = sum(∫(∇(u) ⊙ ∇(u))*dΩ)^(1/2)

    I₂ = one(TensorValue{2,2,Float64})

    # Projeksjons‐operatorene:
    # defining the Navier operators
    Pn(n) = n ⊗ n        # normal projection operator
    Pt(n) = I₂ - Pn(n)     # tangential projection operator
    

    # Legger på base boundary condition på hele
    a0(u, v) = ∫( ε(v)⊙(flux(ε(u))))dΩ              # sjekket med Hanna
    a1(u, v) = (∫(-((n_Γb ⋅ (flux(ε(u)))) ⋅ v))dΓb
    - ∫(((n_Γs ⋅ (flux(ε(u)))) ⋅ v))dΓs 
    )
    a2(u, v) = (∫(-(Pn(n_Γb)⋅ u) ⋅ (n_Γb ⋅ (flux(ε(v)))))dΓb
    - ∫((Pn(n_Γs)⋅ u) ⋅ (n_Γs ⋅ (flux(ε(v)))))dΓs
    )
    a3(u, v) = (∫((2* nu0/(γn*h)*(n_Γb ⋅ u)) ⋅ (n_Γb ⋅ v)  )dΓb
    + ∫((2* nu0/(γn*h)*(n_Γs ⋅ u)) ⋅ (n_Γs ⋅ v)  )dΓs
    ) 
    a4(u, v) = (∫((e/(e + γt *h) * ((Pt(n_Γb) ⋅  (n_Γb ⋅ (flux(ε(u)))))) ⋅ v))dΓb
    + ∫((e/(e + γt *h) * ((Pt(n_Γs) ⋅  (n_Γs ⋅ (flux(ε(u)))))) ⋅ v))dΓs
    )
    a5(u, v) = (∫((nu0/(e + γt*h)*(Pt(n_Γb) ⋅ u)) ⋅ ( Pt(n_Γb) ⋅ v))dΓb 
    + ∫((nu0/(e + γt*h)*(Pt(n_Γs) ⋅ u)) ⋅ ( Pt(n_Γs) ⋅ v))dΓs
    ) 
    a6(u, v) = (∫(-(e*γt * h/(e + γt*h) * (Pt(n_Γb) ⋅  (n_Γb ⋅ (flux∘ε(u)))) ) ⋅ (2*n_Γb ⋅ ε(v)))dΓb         
    - ∫((e*γt * h/(e + γt*h) * (Pt(n_Γs) ⋅  (n_Γs ⋅ (flux∘ε(u)))) ) ⋅ (2*n_Γs ⋅ ε(v)))dΓs
    )
    a7(u, v) = (∫(-(nu0 * γt * h / (e + γt *h) * (Pt(n_Γb)⋅ u)) ⋅ (2* n_Γb ⋅ ε(v)))dΓb  
    - ∫((nu0 * γt * h / (e + γt *h) * (Pt(n_Γs)⋅ u)) ⋅ (2* n_Γs ⋅ ε(v)))dΓs
    )
    
    b(p, v) = (
        ∫(-1*(∇ ⋅ v*p))dΩ 
    + ∫((n_Γb ⋅ v) * p)dΓb 
    + ∫((n_Γs ⋅ v) * p)dΓs
    )

    l(v, q) = (∫(f ⋅ v - g ⋅ q)dΩ)     # OK

    lb(v, q) =(∫(-(n_Γb ⋅ ud) ⋅ ( n_Γb ⋅ (n_Γb ⋅ flux(ε(v)))))dΓb #OK
    + ∫(((2*nu0)/(γn *h) * (n_Γb ⋅ ud)) ⋅ (n_Γb ⋅ v) )dΓb #OKb
    - ∫((n_Γb ⋅ ud) ⋅ q)dΓb  #OK
    + ∫((nu0/(e + γt * h) ⋅ (Pt(n_Γb) ⋅ ud) )⋅ (Pt(n_Γb) ⋅ v) )dΓb #OK    # endrer her, og endrer her... skal dette være odot?
    - ∫((nu0 * γt * h/(e + γt * h) * Pt(n_Γb) ⋅ ud) ⋅ (2*n_Γb ⋅ ε(v)))dΓb # OK
    )

    ls(v, q) =(
    ∫(-(n_Γs⋅ ud) ⋅ ( n_Γs ⋅ (n_Γs ⋅ flux(ε(v)))))dΓs #OK
    + ∫(((2*nu0)/(γn *h) * (n_Γs ⋅ ud)) ⋅ (n_Γs ⋅ v) )dΓs #OK
    - ∫((n_Γs ⋅ ud) ⋅ q)dΓs  #OK
    + ∫((nu0/(e + γt * h) ⋅ (Pt(n_Γs) ⋅ ud) )⋅ (Pt(n_Γs) ⋅ v) )dΓs #OK    # endrer her, og endrer her... skal dette være odot?
    - ∫((nu0 * γt * h/(e + γt * h) * Pt(n_Γs) ⋅ ud) ⋅ (2*n_Γs ⋅ ε(v)))dΓs # OK
    )
    
    # stabiliserer alle facets
    g_u(u,v) = ( ∫( (γu1*h*nu0)*jump(n_Fg⋅∇(u))⋅jump(n_Fg⋅∇(v)) )dFg 
            +  
               ∫( (γu2*h^3*nu0)*jump_nn(u,n_Fg)⋅jump_nn(v,n_Fg) )dFg
)
    g_p(p,q) = ∫( (γp*h^3/nu0)*jump(n_Fg⋅∇(p))*jump(n_Fg⋅∇(q)) )dFg

    if stabilize # have tested with different combinations of calling the gp in the residual and jacobian, but this one seems to work best.
       # hvis jeg vil teste uten newton-løser:
        A((u,p),(v,q)) = a0(u, v) + a1(u, v) + a2(u, v) + a3(u, v) + a4(u, v) + a5(u, v) + a6(u, v) + a7(u, v) + b(p, v) - b(q, u) + g_u(u,v) + g_p(p,q)
        L((v, q)) = l(v, q) + lb(v, q) + ls(v, q)
        op = AffineFEOperator(A, L, X, Y)
        (uh, ph) = solve(op)

    else
        B((u,p),(v,q)) = a0(u, v) - a1(u, v) - a2(u, v) + a3(u, v) + a4(u, v) + a5(u, v) - a6(u, v) - a7(u, v) + b(p, v) - b(q, u) 
        M((v, q)) = l(v, q) + lb(v, q) + ls(v, q)
        op = AffineFEOperator(B, M, X, Y)
        (uh, ph) = solve(op)
    end

    errp = p_exact - ph
    erru = u_exact - uh
    
    condition_numb = 1
    if save
        # writevtk(bgmodel, "C:\\Users\\Sigri\\Documents\\Master\\report\\figures\\Domenefigurer\\mesh_bg$geometry $δ.vtu")
        # writevtk(Fg, "C:\\Users\\Sigri\\Documents\\Master\\report\\figures\\Domenefigurer\\Fg$geometry $δ.vtu")
        # writevtk(Γb, "C:\\Users\\Sigri\\Documents\\Master\\report\\figures\\Domenefigurer\\surface_gamma_b_$geometry $δ.vtu")
        # writevtk(Γs, "C:\\Users\\Sigri\\Documents\\Master\\report\\figures\\Domenefigurer\\surface_gamma_s_$geometry $δ.vtu")
        # writevtk(Ω_act, "C:\\Users\\Sigri\\Documents\\Master\\report\\figures\\Domenefigurer\\Ω_act$geometry $δ.vtu")
        writevtk(Ω, "C:\\Users\\Sigri\\Documents\\Master\\report\\results\\stokes\\$n $geometry $order $δ allb.vtu", cellfields=["u_ex" => u_exact, "uh"=>uh, "erru"=> erru, "p_ex" => p_exact, "ph"=>ph, "errp"=> errp, "nablau" => ∇(u_exact), "flux" => flux∘ε(u_exact), "viskositet" => viskositet∘ε(u_exact)]) #, "erru" => erru]) 
    end 
    return uh, u_exact, erru, l2_norm(uh - u_exact), h1_semi(uh - u_exact), ph, p_exact, errp, l2_norm(ph - p_exact), h1_semi(ph - p_exact), condition_numb, Ω
end
solver = stokes_navierBC_CutFEM_cutmeshboundary
uh, u_exact, erru, ul2_norm, uh1_semi, ph, p_exact, errp, pl2_norm, ph1_semi, condition_numb, Ω_act = solver(;n, u_exact, p_exact, f, g, ud, order, geometry, βu0, γu1, γu2, γp, βp0, nu, stabilize, δ, save, calc_condition)
# order = 2
# numb_it = 4
# uarr_l2_stab, uarr_h1_stab, parr_l2_stab, parr_h1_stab, h_arr = convergence_stokes_weird_domain(;numb_it, u_exact, p_exact, f, g, ud, order, geometry, solver, δ, βu0, γu1, γu2, γp, βp0, nu, stabilize, save)

# stabilize = false
# uarr_l2_nostab, uarr_h1_nostab, parr_l2_nostab, parr_h1_nostab, h = convergence_stokes(;numb_it, u_exact, p_exact, f, g, ud, order, geometry, solver, δ, βu0, γu1, γu2, γp, βp0, nu, stabilize, save)

# ########## velocity convergence plot:
# plot_convergence_u(
#     uarr_l2_stab, uarr_h1_stab, h_arr;
#     #uarr_l2_nostab = uarr_l2_nostab, uarr_h1_nostab = uarr_h1_nostab,
#     title_str="Convergence of CutFEM Stokes Navier BC"
# )

# plot_convergence_p(
#     parr_l2_stab, parr_h1_stab, h_arr;
#     #parr_l2_nostab=parr_l2_nostab,
#     #parr_h1_nostab=parr_h1_nostab,
#     title_str="Convergence of CutFEM Stokes Navier BC"
# )

# # #print eoc-verdier:
# print_eoc_latex_combined(h;
#     uarr_l2_stab = uarr_l2_stab,
#     uarr_h1_stab = uarr_h1_stab,
#     uarr_l2_nostab = uarr_l2_nostab,
#     uarr_h1_nostab = uarr_h1_nostab,
#     start = 1
# )

# print_eoc_latex_combined_pressure(h;
#     parr_l2_stab = parr_l2_stab,
#     parr_h1_stab = parr_h1_stab,
#     parr_l2_nostab = parr_l2_nostab,
#     parr_h1_nostab = parr_h1_nostab,
#     start = 1
# )

# n = 16                # øke denne?
# M = 100               #full kjøring med M = 2000 med 2000 så kjører det nok i 2 timer. 
# # βu0 = 1
# # γu1 = 1
# # γu2 = 1
# # γp = 0.1
# # βp0 = 0.1
# stabilize = true
# save = false
# calc_condition = false

# arr_δ, arr_l2u, arr_h1u, arr_l2p, arr_h1p, arr_cond = sensitivity_stokes(;n, M, u_exact, p_exact, f, g, ud, order, geometry, solver, βu0, γu1, γu2, γp, βp0, nu, stabilize, save)
# stabilize = false
# start = 1
# save = true
# arr_δ_nostab, arr_l2u_nostab, arr_h1u_nostab, arr_l2p_nostab, arr_h1p_nostab, arr_cond_nostab = sensitivity_stokes(;n, M, u_exact, p_exact, f, g, ud, order, geometry, solver, βu0, γu1, γu2, γp, βp0, nu, stabilize, save)

# plot_sensitivity_velocity(
#     arr_δ, arr_l2u, arr_h1u,
#     arr_δ_nostab, arr_l2u_nostab, arr_h1u_nostab;
#     start = start,
#     title_str = "Sensitivity of CutFEM Stokes Navier BC"
# )

# plot_sensitivity_pressure(
#     arr_δ, arr_l2p, arr_h1p,
#     arr_δ_nostab, arr_l2p_nostab, arr_h1p_nostab;
#     start = start,
#     title_str = "Sensitivity of CutFEM Stokes Navier BC"
# )
