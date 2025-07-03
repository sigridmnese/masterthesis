using Gridap
using GridapEmbedded
using STLCutters
using LinearAlgebra
using Plots
using LineSearches: BackTracking
import Random
using Logging
using Test
using DataFrames, CSV

const CONV_FILE = "convergence.csv"

function append_convergence_row(n::Int; stabilize::Bool = true)
    # --- kall solver ---
    uh, u_ex, erru, ul2_norm, uh1_semi,
    ph, p_ex, errp, pl2_norm, ph1_semi,
    condition_numb, Ω_act =
      solver(; n,
               u_exact, p_exact, f, g, ud,
               order, geometry,
               βu0, γu1, γu2, γp, βp0,
               nu,
               stabilize,    # <-- toggles stab/no-stab
               δ,
               save=false,
               calc_condition=false)

    # Beregn h
    h = 1/n

    # --- les inn eksisterende eller lag nytt DataFrame ---
    df = isfile(CONV_FILE) ? CSV.read(CONV_FILE, DataFrame) : DataFrame(
        n = Int[],
        stabilize = Bool[],
        h = Float64[],
        u_l2 = Float64[],
        u_h1 = Float64[],
        p_l2 = Float64[],
        p_h1 = Float64[],
    )

    # --- legg til ny rad og skriv ut ---
    push!(df, (
        n = n,
        stabilize = stabilize,
        h = h,
        u_l2 = ul2_norm,
        u_h1 = uh1_semi,
        p_l2 = pl2_norm,
        p_h1 = ph1_semi,
    ))
    CSV.write(CONV_FILE, df)
    println("Appended: n=$n, stabilize=$(stabilize), h=$(h), u_l2=$(ul2_norm), u_h1=$(uh1_semi), p_l2=$(pl2_norm), p_h1=$(ph1_semi)")
end

function jump_nn(u,n)
    return ( n.plus ⋅ (n.plus⋅∇∇(u).plus) - n.minus ⋅ (n.minus ⋅ ∇∇(u).minus) )       # andre ordens hopp... Forklare dette skikkelig. 
end

# Denne er med konstant viskositet, og med nitsche implementasjon av stokes med navier BC.
# Defining constants
nu0 =  1   # klarer ikke mindre enn 0.01 og klarer heller ikke større enn 100 (størrelsesordener)
r = 4/3
A = nu0
ϵ_0 = 1e-7       

# Defining manufactured solutions
u_exact(x) =  VectorValue(2*x[1] + cos(2*π*x[2]), -2*x[2] + sin(2*π*x[1]), 0)#VectorValue(2*x[1] + exp(x[1]/2) * cos(2*π*x[2]), -2*x[2] + exp(x[2]/2) * sin(2*π*x[1])) #bytte til sin/cos-uttrykk  VectorValue(-x[2], x[1])
p_exact(x) = sin(2*π*x[1])*cos(2*π*x[2])

viskositet(εu) = nu0#^(1-r)*(ϵ_0^2 + norm(εu)^2)^((r-2)/2)
dviskositet(εu, εdu) = 0#nu0^(1-r) *(ϵ_0^2 + εu ⋅ εu)^((r-4)/2)*(εu⊙εdu) *(r-2)

# dette er hvis viskositeten er en konstant:
flux(εu) = 2 * nu0 * εu         # det kan være jeg ikke skal gange med 2 her...?
dflux(εu, εdu) = flux∘(εdu) # eller her for øvrig

f(x) =  -divergence(flux∘ε(u_exact))(x) + ∇(p_exact)(x)      # prøver å endre f her...
ud(x) = u_exact(x)
g(x) = tr(ε(u_exact)(x))
# prøver også å bytte ud til g...

# solver parameters
n = 32
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
function stokes_navierBC_CutFEM_3d(;n, u_exact, p_exact, f, g, ud, order, geometry, βu0, γu1, γu2, γp, βp0, nu, stabilize, δ, save = false, calc_condition = false)
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
    println(n)
    γn = 0.05  # hentet fra Josefin sin artikkel.
    γt = 0.05  # hentet fra Josefin sin artikkel.
    e = 0 

    # %% Define background mesh and geometry
    # Disc center and radius
    # R = 0.753
    p0 = Point(0.0, 0.0, 0.0)
    # geo = sphere(R, x0=p0)

    # Take complement of the ball
    # geo = !sphere(R, x0=p0)

    # Define background mesh
    partition = (n, n, n)
    dim = length(partition)
    a = 1
    pmin = p0+Point(-a + δ, -a + δ, δ)
    pmax = p0+Point(a + δ, a + δ, a + δ)
    bgmodel = CartesianDiscreteModel(pmin, pmax, partition)
    Ω_bg = Triangulation(bgmodel)
    writevtk(Ω_bg, "mesh_bg")

    # mesh size
    h = (pmax-pmin)[1]/partition[1]

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
    Γb = EmbeddedBoundary(cutgeo)
    n_Γb = get_normal_vector(Γb)

    # Outer mesh boundary 
    # Has Dirichlet for pF and no-stress for u
    Γs = BoundaryTriangulation(cutgeo_facets, PHYSICAL_IN)
    n_Γs = get_normal_vector(Γs)

    # Get ghost penalty facets
    Fg = GhostSkeleton(cutgeo)
    n_Fg = get_normal_vector(Fg)

    # # All interior facets of active mesh
    # Fi = SkeletonTriangulation(Ω_act)
    # n_Fi = get_normal_vector(Fi)

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
    # dFi = Measure(Fi, degree)

    # Write out domains 
    writevtk(Ω_act, "mesh_act")
    writevtk(Ω, "mesh")
    writevtk(Γb, "surface_gamma_d")
    writevtk(Γs, "surface_gamma_s")
    writevtk(Fg, "ghost_facets")

    # %% Define weak formulation
    # Define function spaces
    reffe_u  = ReferenceFE(lagrangian,VectorValue{dim, Float64},order_u)
    reffe_p = ReferenceFE(lagrangian,Float64, order_p)

    V = TestFESpace(Ω_act, reffe_u,  conformity=:H1)
    Q = TestFESpace(Ω_act, reffe_p, conformity=:H1)     

    U = TrialFESpace(V)
    P = TrialFESpace(Q)

    X = MultiFieldFESpace([U, P])
    Y = MultiFieldFESpace([V, Q])
    
    #println(sum( ∫( p_exact) * dΩ ))
    l2_norm(u) = (sum( ∫( u ⋅ u )*dΩ ))^(1/2)
    h1_semi(u) = sum(∫(∇(u) ⊙ ∇(u))*dΩ)^(1/2)
    
    ϵ_0 = 1e-6

    #defining viskosity parameter:
    # You have to chage the call for γ in three functions: a - last term, l3 - first term, da - last term
    γ = nu0*order*order
    #γ(∇u) = viskositet(∇u)
    #γ(u) = viskositet∘ε(u)

    # weak formulation components    
    # a(u, v) = ∫( ε(v)⊙(flux∘ε(u)))dΩ  + ∫(-((n_Γd ⋅ (flux∘ε(u))) ⋅ v) + (-(n_Γd ⋅ (flux∘ε(v))) ⋅ u) + viskositet∘ε(u)/h * (u ⋅ v))dΓd      # denne må ha et ekstra boundary term. Finn ut hvordan det ser ut. 
    # b(v, p) = (∫(-1*(∇ ⋅ v*p))dΩ + ∫((n_Γd ⋅ v) * p)dΓd)   # b er den samme som før. 
    # l1((v, q)) = ∫(f ⋅ v)dΩ
    # l2(v) = ∫(-(n_Γd ⋅ (flux∘ε(v))) ⋅ ud)dΓd    # har brukt ud som dirichlet grense
    # l3((u, p), (v, q)) = ∫( viskositet∘ε(u)/h* (ud ⋅ v))dΓd
    # l4(q) = ∫((n_Γd ⋅ ud) * q)dΓd               # kan det være at denne skal være annerledes? mtp flux * ... og ikke n_Γd ⋅ ud?

    I₂ = one(TensorValue{3,3,Float64})

    # Projeksjons‐operatorene:
    # defining the Navier operators
    Pn(n) = n ⊗ n        # normal projection operator
    Pt(n) = I₂ - Pn(n)     # tangential projection operator
    
    # Legger på base boundary condition på hele
    a0(u, v) = ∫( ε(v)⊙(flux(ε(u))))dΩ              
    a1(u, v) = (∫(-((n_Γb ⋅ (flux(ε(u)))) ⋅ v))dΓb
    )
    a2(u, v) = (∫(-(Pn(n_Γb)⋅ u) ⋅ (n_Γb ⋅ (flux(ε(v)))))dΓb
    )
    a3(u, v) = (∫((2* nu0/(γn*h)*(n_Γb ⋅ u)) ⋅ (n_Γb ⋅ v)  )dΓb
    ) 
    a4(u, v) = (∫((e/(e + γt *h) * ((Pt(n_Γb) ⋅  (n_Γb ⋅ (flux(ε(u)))))) ⋅ v))dΓb
    )
    a5(u, v) = (∫((nu0/(e + γt*h)*(Pt(n_Γb) ⋅ u)) ⋅ ( Pt(n_Γb) ⋅ v))dΓb 
    ) 
    a6(u, v) = (∫(-(e*γt * h/(e + γt*h) * (Pt(n_Γb) ⋅  (n_Γb ⋅ (flux∘ε(u)))) ) ⋅ (2*n_Γb ⋅ ε(v)))dΓb         
    )
    a7(u, v) = (∫(-(nu0 * γt * h / (e + γt *h) * (Pt(n_Γb)⋅ u)) ⋅ (2* n_Γb ⋅ ε(v)))dΓb  
    )
    
    b(p, v) = (
        ∫(-1*(∇ ⋅ v*p))dΩ 
    + ∫((n_Γb ⋅ v) * p)dΓb 
    )

    bu(q, u) = (
        ∫(-(∇ ⋅ u*q))dΩ             #+ ∫((n_Γb ⋅ u) * q)dΓb   
        )

    l(v, q) = (∫(f ⋅ v - g ⋅ q)dΩ)     # OK        

    lb(v, q) =(∫(-(n_Γb ⋅ ud) ⋅ ( n_Γb ⋅ (n_Γb ⋅ flux(ε(v)))))dΓb #OK
    + ∫(((2*nu0)/(γn *h) * (n_Γb ⋅ ud)) ⋅ (n_Γb ⋅ v) )dΓb #OKb
    - ∫((n_Γb ⋅ ud) ⋅ q)dΓb  #OK
    + ∫((nu0/(e + γt * h) ⋅ (Pt(n_Γb) ⋅ ud) )⋅ (Pt(n_Γb) ⋅ v) )dΓb #OK    # endrer her, og endrer her... skal dette være odot?
    - ∫((nu0 * γt * h/(e + γt * h) * Pt(n_Γb) ⋅ ud) ⋅ (2*n_Γb ⋅ ε(v)))dΓb # OK
    )

    ls(v, q) =(
    ∫(((  2*nu0 *(n_Γs ⋅(ε(u_exact))) - n_Γs ⋅ I₂ ⋅ p_exact) ⋅  v))dΓs)     #-1*q * (n_Γs ⋅ ud) +
    
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
        B((u,p),(v,q)) = a0(u, v) + a1(u, v) + a2(u, v) + a3(u, v) + a4(u, v) + a5(u, v) + a6(u, v) + a7(u, v) + b(p, v) - b(q, u) 
        M((v, q)) = l(v, q) + lb(v, q) + ls(v, q)
        op = AffineFEOperator(B, M, X, Y)
        (uh, ph) = solve(op)
    end


    errp = p_exact - ph
    erru = u_exact - uh
    
    condition_numb = 1
    if save
        #writevtk(bgmodel, "C:\\Users\\Sigri\\Documents\\Master\\report\\figures\\Domenefigurer\\mesh_bg$geometry $δ.vtu")
        # writevtk(Fg, "C:\\Users\\Sigri\\Documents\\Master\\report\\figures\\Domenefigurer\\Fg$geometry $δ.vtu")
        # writevtk(Γd, "C:\\Users\\Sigri\\Documents\\Master\\report\\figures\\Domenefigurer\\surface_gamma_d_$geometry $δ.vtu")
        # writevtk(Ω_act, "C:\\Users\\Sigri\\Documents\\Master\\report\\figures\\Domenefigurer\\Ω_act$geometry $δ.vtu")
        # writevtk(Ω_act, "C:\\Users\\Sigri\\Documents\\Master\\report\\figures\\Domenefigurer\\Ω_act$geometry $δ.vtu")
        writevtk(Ω, "$n.vtu", cellfields=["u_ex" => u_exact, "uh"=>uh, "erru"=> erru, "p_ex" => p_exact, "ph"=>ph, "errp"=> errp, "nablau" => ∇(u_exact), "flux" => flux∘ε(u_exact), "viskositet" => viskositet∘ε(u_exact)]) #, "erru" => erru]) 
    end 
    return uh, u_exact, erru, l2_norm(uh - u_exact), h1_semi(uh - u_exact), ph, p_exact, errp, l2_norm(ph - p_exact), h1_semi(ph - p_exact), condition_numb, Ω
end
solver = stokes_navierBC_CutFEM_3d
# n = 16
# uh, u_exact, erru, ul2_norm, uh1_semi, ph, p_exact, errp, pl2_norm, ph1_semi, condition_numb, Ω_act = solver(;n, u_exact, p_exact, f, g, ud, order, geometry, βu0, γu1, γu2, γp, βp0, nu, stabilize, δ, save, calc_condition)


append_convergence_row(16; stabilize=true)



# ###################################  sensitivity_stokes test ########################################## 
# n = 16                # øke denne?
# M = 200               #full kjøring med M = 2000 med 2000 så kjører det nok i 2 timer. 
# βu0 = 1
# γu1 = 1
# γu2 = 1
# γp = 0.1
# βp0 = 0.1
# stabilize = true
# save = false
# calc_condition = false

# arr_δ, arr_l2u, arr_h1u, arr_l2p, arr_h1p, arr_cond = sensitivity_stokes(;n, M, u_exact, p_exact, f, g, ud, order, geometry, solver, βu0, γu1, γu2, γp, βp0, nu, stabilize, save)
# stabilize = false
# start = 1
# save = false
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
