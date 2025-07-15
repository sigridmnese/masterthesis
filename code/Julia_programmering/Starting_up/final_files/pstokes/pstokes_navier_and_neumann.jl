using Gridap
using GridapEmbedded
using STLCutters
using LinearAlgebra
using Plots
using LineSearches: BackTracking
import Random
using Logging
using LoggingExtras
using PGFPlotsX


include("C:\\Users\\Sigri\\Documents\\Master\\report\\code\\Julia_programmering\\Starting_up\\utils\\utils.jl")
include("C:\\Users\\Sigri\\Documents\\Master\\report\\code\\Julia_programmering\\Starting_up\\workfiles_stokes\\testing_non-linear_stokes.jl")

# Defining constants
r = 4/3
A = 2^(1/(1-r))
ϵ_0 = 1e-6
nu0 =  A

# Defining manufactured solutions
u_exact(x) =  VectorValue(2*x[1] + cos(2*π*x[2]), -2*x[2] + sin(2*π*x[1]))#VectorValue(2*x[1] + exp(x[1]/2) * cos(2*π*x[2]), -2*x[2] + exp(x[2]/2) * sin(2*π*x[1])) #bytte til sin/cos-uttrykk  VectorValue(-x[2], x[1])
p_exact(x) = sin(2*π*x[1])*cos(2*π*x[2])

# Defining the problem 
flux(εu) = nu0^(1-r)*(ϵ_0^2 + 1/2* tr(εu' ⋅ εu))^((r-2)/2) ⋅ εu         #nu0^(1-r)*(ϵ_0^2 + tr(εu' ⋅ εu))^((r-2)/2)⋅εu 
dflux(εu, εdu)= (r-2)/2 * nu0^(1-r) *(ϵ_0^2 + 1/2 * tr(εu' ⋅ εu))^((r-4)/2)*(εu⊙εdu) ⋅ εu + nu0^(1-r)*(ϵ_0^2 + 1/2* tr(εu' ⋅ εu))^((r-2)/2) ⋅ εdu #(r-2)*nu0^(1-r)*(ϵ_0^2 + tr(εu' ⋅ εu))^((r-4)/2)*(εu⊙εdu) ⋅ εu + nu0^(1-r)*(ϵ_0^2 + tr(εu' ⋅ εu))^((r-2)/2)*εdu
viskositet(εu) = 1/2 * nu0^(1-r)*(ϵ_0^2 + 1/2* tr(εu' ⋅ εu))^((r-2)/2)
dviskositet(εu, εdu) = (r-2)/4 * nu0^(1-r) *(ϵ_0^2 + 1/2 * tr(εu' ⋅ εu))^((r-4)/2)*(εu⊙εdu) 

f(x) =  -divergence(flux∘ε(u_exact))(x) + ∇(p_exact)(x)      # prøver å endre f her...
ud(x) = u_exact(x)
g(x) = -tr(ε(u_exact)(x))

# solver:
function pstokes_CutFEM_navierandneumann(;n, u_exact, p_exact, f, g, ud, order, geometry, βu0, γu1, γu2, γp, βp0, nu, stabilize, δ, save = false, calc_condition = false)
    """    
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
    geo = !disk(R, x0=p0)    #kan invertere hva som er bgmodel og hva som er physical mesh ved å fjerne utropstegn

    # Define background mesh
    partition = (n, n)
    dim = length(partition)
    a = 1
    pmin = p0 + Point(-a + δ, -a + δ)
    pmax = p0 + Point(a + δ, a + δ)
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
    Γd = BoundaryTriangulation(cutgeo_facets, PHYSICAL_IN)
    n_Γd = get_normal_vector(Γd)

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

    reffe_u  = ReferenceFE(lagrangian,VectorValue{dim, Float64},order)
    reffe_p = ReferenceFE(lagrangian,Float64, order-1)

    V = TestFESpace(Ω_act, reffe_u,  conformity=:H1)
    Q = TestFESpace(Ω_act, reffe_p, conformity=:H1) #, constraint=:zeromean)

    U = TrialFESpace(V)
    P = TrialFESpace(Q)

    X = MultiFieldFESpace([U, P])
    Y = MultiFieldFESpace([V, Q])
 
    # mesh size
    h = (pmax-pmin)[1]/partition[1]
    
    #println(sum( ∫( p_exact) * dΩ ))
    l2_norm(u) = (sum( ∫( u ⋅ u )*dΩ ))^(1/2)
    h1_semi(u) = sum(∫(∇(u) ⊙ ∇(u))*dΩ)^(1/2)

    I₂ = one(TensorValue{2,2,Float64})

    # Projeksjons‐operatorene:
    Pn(n) = n ⊗ n        # normal projection operator
    Pt(n) = I₂ - Pn(n)     # tangential projection operator
    
    # Defining the weak form
    a0(u, v) = ∫(ε(v)⊙(flux∘ε(u)))dΩ              # sjekket med Hanna
    a1(u, v) = (∫(-((n_Γd ⋅ (flux∘(ε(u)))) ⋅ v))dΓd
    )
    a2(u, v) = (∫(-(Pn(n_Γd)⋅ u) ⋅ (n_Γd ⋅ (flux∘ε(v))))dΓd
    )
    a3(u, v) = (∫((2/(γn*h)*(viskositet∘ε(u) ⋅ ((n_Γd ⋅ u)) ⋅ (n_Γd ⋅ v))) )dΓd
    ) 
    a4(u, v) = (∫((e/(e + γt *h) * ((Pt(n_Γd) ⋅  (n_Γd ⋅ (flux∘(ε(u)))))) ⋅ v))dΓd
    )
    a5(u, v) = (∫(1/(e + γt*h)* ( viskositet∘ε(u) ⋅ (Pt(n_Γd) ⋅ u) ⋅ ( Pt(n_Γd) ⋅ v)))dΓd 
    ) 
    a6(u, v) = (∫(-(e*γt * h/(e + γt*h) * (Pt(n_Γd) ⋅  (n_Γd ⋅ (flux∘ε(u)))) ) ⋅ (2*n_Γd ⋅ ε(v)))dΓd         
    )
    a7(u, v) = (∫(-(γt * h / (e + γt *h) *( viskositet∘ε(u) ⋅ ((Pt(n_Γd)⋅ u) ⋅ (2* n_Γd ⋅ ε(v))))))dΓd  
    )

    da0(u, du, v) = ∫( ε(v)⊙(dflux∘(ε(u), ε(du))))dΩ              # sjekket med Hanna
    da1(u, du, v) = (∫(-((n_Γd ⋅ ((dflux∘(ε(u), ε(du))))) ⋅ v))dΓd
    )
    da2(u, du, v) = (∫(-(Pn(n_Γd)⋅ du) ⋅ (n_Γd ⋅ ((flux∘(ε(v))))))dΓd
    )
    da3(u, du, v) = (∫((2/(γn*h)*(dviskositet∘(ε(u), ε(du)) ⋅ ((n_Γd ⋅ u)) ⋅ (n_Γd ⋅ v))) )dΓd
                    + 
                    ∫((2/(γn*h)*(viskositet∘ε(u) ⋅ ((n_Γd ⋅ du)) ⋅ (n_Γd ⋅ v))) )dΓd
    ) 
    da4(u, du, v) = (∫((e/(e + γt *h) * ((Pt(n_Γd) ⋅  (n_Γd ⋅ ((dflux∘(ε(u), ε(du))))))) ⋅ v))dΓd
    )
    da5(u, du, v) = (∫(1/(e + γt*h)* ( dviskositet∘(ε(u), ε(du)) ⋅ (Pt(n_Γd) ⋅ u) ⋅ ( Pt(n_Γd) ⋅ v)))dΓd
    +
    ∫(1/(e + γt*h)* ( viskositet∘ε(u) ⋅ (Pt(n_Γd) ⋅ du) ⋅ ( Pt(n_Γd) ⋅ v)))dΓd 
    ) 
    da6(u, du, v) = (∫(-(e*γt * h/(e + γt*h) * (Pt(n_Γd) ⋅  (n_Γd ⋅ (dflux∘(ε(u), ε(du))))) ) ⋅ (2*n_Γd ⋅ ε(v)))dΓd         
    )

    da7(u, du, v) = (∫(-(γt * h / (e + γt *h) *( dviskositet∘(ε(u), ε(du)) ⋅ ((Pt(n_Γd)⋅ u) ⋅ (2* n_Γd ⋅ ε(v))))))dΓd
    +
    ∫(-(γt * h / (e + γt *h) *( viskositet∘ε(u) ⋅ ((Pt(n_Γd)⋅ du) ⋅ (2* n_Γd ⋅ ε(v))))))dΓd  
    )

    ####################### differentiating a##################################
    
    b(p, v) = (
        ∫(-1*(∇ ⋅ v*p))dΩ 
    + ∫((n_Γd ⋅ v) * p)dΓd 
    )

    bu(q, u) = (
        ∫(-(∇ ⋅ u*q))dΩ             #+ ∫((n_Γb ⋅ u) * q)dΓb   
        )

    l(v, q) = (∫(f ⋅ v - g ⋅ q)dΩ)     # OK

    lb(u, v, q) =(∫(-(n_Γd ⋅ ud) ⋅ ( n_Γd ⋅ (n_Γd ⋅ (flux∘ε(v)))))dΓd #OK
    + ∫((2/(γn *h) * (viskositet∘ε(u) ⋅ (((n_Γd ⋅ ud)) ⋅ (n_Γd ⋅ v)))))dΓd #OKb
    - ∫((n_Γd ⋅ ud) ⋅ q)dΓd  #OK
    + ∫(1/(e + γt * h) ⋅ (viskositet∘ε(u) ⋅ ((Pt(n_Γd) ⋅ ud) ⋅ (Pt(n_Γd) ⋅ v)) ))dΓd #OK    # endrer her, og endrer her... skal dette være odot?
    - ∫((γt * h/(e + γt * h) * (viskositet∘ε(u) ⋅ (( Pt(n_Γd) ⋅ ud) ⋅ (2*n_Γd ⋅ ε(v))))))dΓd # OK
    )

    dlb(u, du,  v, q) =(#∫(-(n_Γd ⋅ ud) ⋅ ( n_Γd ⋅ (n_Γd ⋅ (flux∘ε(v)))))dΓd #OK
        ∫(2/(γn *h) * (dviskositet∘(ε(u), ε(du)) ⋅ ((n_Γd ⋅ ud) ⋅ (n_Γd ⋅ v))))dΓd #OKb
        #- ∫((n_Γd ⋅ ud) ⋅ q)dΓd  #OK
        + ∫((1/(e + γt * h) * (dviskositet∘(ε(u), ε(du)) ⋅ ((Pt(n_Γd) ⋅ ud) ⋅ (Pt(n_Γd) ⋅ v)))))dΓd #OK
        - ∫((γt * h/(e + γt * h) * (dviskositet∘(ε(u), ε(du)) ⋅ (( Pt(n_Γd) ⋅ ud) ⋅ (2*(n_Γd ⋅ ε(v)))))))dΓd # OK
    )

    ls(v, q) =(
    ∫(((  (n_Γs ⋅(flux∘ε(u_exact))) - n_Γs ⋅ I₂ ⋅ p_exact) ⋅  v))dΓs)     #-1*q * (n_Γs ⋅ ud) +
    
    
    # stabiliserer alle facets
    g_u(u,v) = ( ∫( (γu1*h)*jump(n_Fg⋅∇(u))⋅jump(n_Fg⋅∇(v)) )dFg 
            +  
               ∫( (γu2*h^3)*jump_nn(u,n_Fg)⋅jump_nn(v,n_Fg) )dFg
)
    g_p(p,q) = ∫( (γp*h^3)*jump(n_Fg⋅∇(p))*jump(n_Fg⋅∇(q)) )dFg

    
    if stabilize
      res((u,p),(v,q)) = (a0(u, v) 
      + a1(u, v) 
      + a2(u, v) 
      + a3(u, v) 
      + a4(u, v) 
      + a5(u, v) 
      + a6(u, v) 
      + a7(u, v) 
      + b(p, v) 
      - b(q, u)
      + g_u(u,v) 
      + g_p(p,q) 
      - l(v, q) 
      - lb(u, v, q)
      - ls(v, q)
      )
      jac((u, p), (du, dp), (v, q)) = (da0(u, du, v) 
      + da1(u, du, v) 
      + da2(u, du, v) 
      + da3(u, du, v) 
      + da4(u, du, v) 
      + da5(u, du, v) 
      + da6(u, du, v) 
      + da7(u, du, v) 
      + b(dp, v)
      - b(q, du)
      + g_u(du, v) 
      + g_p(dp, q) 
      - dlb(u, du, v, q)
      )

      op = FEOperator(res, jac, X, Y)

      # non-linear phase
      nls = NLSolver(
      show_trace=true, method=:newton, linesearch=BackTracking(), iterations=100)      #prøver å legge inn et max antall iterasjoner og en lav toleranse      
      solver = FESolver(nls)

      (uh, ph) = solve(solver, op)

    else
      res2((u,p),(v,q)) = a0(u, v) + a1(u, v) + a2(u, v) + a3(u, v) + a4(u, v) + a5(u, v) + a6(u, v) + a7(u, v) + b(p, v) - b(q, u) - l(v, q) - lb(u, v, q)
      jac2((u, p), (du, dp), (v, q)) = da0(u, du, v) + da1(u, du, v) + da2(u, du, v) + da3(u, du, v) + da4(u, du, v) + da5(u, du, v) + da6(u, du, v) + da7(u, du, v) + b(dp, v) - b(q, du) - dlb(u, du, v, q)
    
      op = FEOperator(res2, jac2, X, Y)

      # non-linear phase
      nls = NLSolver(
      show_trace=true, method=:newton, linesearch=BackTracking(), iterations=100)      #prøver å legge inn et max antall iterasjoner og en lav toleranse      
      solver = FESolver(nls)

      (uh, ph) = solve(solver, op)
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
        writevtk(Ω, "C:\\Users\\Sigri\\Documents\\Master\\report\\results\\stokes\\$n $geometry $order.vtu", cellfields=["u_ex" => u_exact, "uh"=>uh, "erru"=> erru, "p_ex" => p_exact, "ph"=>ph, "errp"=> errp, "nablau" => ∇(u_exact), "viskositet" => viskositet∘ε(u_exact)]) #, "erru" => erru]) 
    end
    return uh, u_exact, erru, l2_norm(uh - u_exact), h1_semi(uh - u_exact), ph, p_exact, errp, l2_norm(ph - p_exact), h1_semi(ph - p_exact), condition_numb, Ω
end

# solver parameters
n = 64
stabilize = true
solver = pstokes_CutFEM_navierandneumann
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

# # ################################################## convergence ##########################################################
# numb_it = 7
# stabilize = true
# uarr_l2, uarr_h1, parr_l2, parr_h1, harr = convergence_stokes(;numb_it, u_exact, p_exact, f, g, ud, order, geometry, solver, δ, βu0, γu1, γu2, γp, βp0, nu, stabilize, save)

# # stabilize = false
# # uarr_l2_nostab, uarr_h1_nostab, parr_l2_nostab, parr_h1_nostab, harr = convergence_stokes(;numb_it, u_exact, p_exact, f, g, ud, order, geometry, solver, δ, βu0, γu1, γu2, γp, βp0, nu, stabilize, save)

# df = DataFrame(harr = harr, uarr_l2 = uarr_l2, uarr_h1 = uarr_h1, parr_l2 = parr_l2, parr_h1 = parr_h1)

# # # ########## velocity convergence plot:
# plot_convergence_u(
#     uarr_l2, uarr_h1, harr;
#     # uarr_l2_nostab=uarr_l2_nostab,
#     # uarr_h1_nostab=uarr_h1_nostab,
#     title_str="Convergence p-Stokes"
# )

# plot_convergence_p(
#     parr_l2, parr_h1, harr;
#     # parr_l2_nostab=parr_l2_nostab,
#     # parr_h1_nostab=parr_h1_nostab,
#     title_str="Convergence p-Stokes"
# )

# # # # # #print eoc-verdier:
# print_eoc_latex_combined_only_eoc(harr;
#     uarr_l2_stab = uarr_l2,
#     uarr_h1_stab = uarr_h1,
#     # uarr_l2_nostab = uarr_l2_nostab,
#     # uarr_h1_nostab = uarr_h1_nostab,
#     start = 1
# )

# print_eoc_latex_combined_pressure_only_eoc(harr;
#     parr_l2_stab = parr_l2,
#     parr_h1_stab = parr_h1,
#     # parr_l2_nostab = parr_l2_nostab,
#     # parr_h1_nostab = parr_h1_nostab,
#     start = 1
# )

# # n = 16                # øke denne?
# # M = 10               #full kjøring med M = 2000 med 2000 så kjører det nok i 2 timer. 
# # # βu0 = 1
# # # γu1 = 1
# # # γu2 = 1
# # # γp = 0.1
# # # βp0 = 0.1
# # stabilize = true
# # save = false
# # calc_condition = false

# # arr_δ, arr_l2u, arr_h1u, arr_l2p, arr_h1p, arr_cond = sensitivity_stokes(;n, M, u_exact, p_exact, f, g, ud, order, geometry, solver, βu0, γu1, γu2, γp, βp0, nu, stabilize, save)
# # stabilize = false
# # start = 1
# # save = false
# # arr_δ_nostab, arr_l2u_nostab, arr_h1u_nostab, arr_l2p_nostab, arr_h1p_nostab, arr_cond_nostab = sensitivity_stokes(;n, M, u_exact, p_exact, f, g, ud, order, geometry, solver, βu0, γu1, γu2, γp, βp0, nu, stabilize, save)


# # # df = DataFrame(arr_δ = arr_δ,
# # # arr_l2u = arr_l2u,
# # # arr_h1u = arr_h1u,
# # # arr_l2p = arr_l2p,
# # # arr_h1p = arr_h1p,
# # # arr_δ_nostab = arr_δ_nostab,
# # # arr_l2u_nostab = arr_l2u_nostab,
# # # arr_h1u_nostab = arr_h1u_nostab,
# # # arr_l2p_nostab = arr_l2p_nostab,
# # # arr_h1p_nostab = arr_h1p_nostab
# # # )

# # plot_sensitivity_velocity(
# #     arr_δ, arr_l2u, arr_h1u,
# #     arr_δ_nostab, arr_l2u_nostab, arr_h1u_nostab;
# #     start = start,
# #     title_str = "Sensitivity p-Stokes"
# # )

# # plot_sensitivity_pressure(
# #     arr_δ, arr_l2p, arr_h1p,
# #     arr_δ_nostab, arr_l2p_nostab, arr_h1p_nostab;
# #     start = start,
# #     title_str = "Sensitivity p-Stokes"
# # )

# # plot_sensitivity_poisson(
# #     arr_δ, arr_l2u, arr_h1u, arr_δ_nostab, arr_l2u_nostab, arr_h1u_nostab;
# #     marker_l2_stab = :circle,
# #     marker_l2_nostab = :star5,
# #     marker_h1_stab = :circle,
# #     marker_h1_nostab = :star5,
# #     markstep = 50,
# #     title_str = "Sensitivity p-Stokes, velocity",
# # )

# # plot_sensitivity_poisson(
# #     arr_δ, arr_l2p, arr_h1p, arr_δ_nostab, arr_l2p_nostab, arr_h1p_nostab;
# #     marker_l2_stab = :circle,
# #     marker_l2_nostab = :star5,
# #     marker_h1_stab = :circle,
# #     marker_h1_nostab = :star5,
# #     markstep = 50,
# #     title_str = "Sensitivity p-Stokes, pressure",
# # )