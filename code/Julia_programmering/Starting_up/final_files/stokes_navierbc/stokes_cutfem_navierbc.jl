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

"""CutFEM with weak boundary imposition using Nitsches method. Symmetric stokes equations with Navier boundary conditions. Lagt inn plott i teams?"""
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
g(x) = VectorValue(0, 0)#tr(ε(u_exact)(x))
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
function stokes_navierBC_CutFEM(;n, u_exact, p_exact, f, g, ud, order, geometry, βu0, γu1, γu2, γp, βp0, nu, stabilize, δ, save = false, calc_condition = false)
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
  
    # Define function spaces 
    reffe_u  = ReferenceFE(lagrangian,VectorValue{dim, Float64},order)
    reffe_p = ReferenceFE(lagrangian,Float64, order - 1)
  
    V = TestFESpace(Ω_act, reffe_u,  conformity=:H1)
    Q = TestFESpace(Ω_act, reffe_p, conformity=:H1, constraint=:zeromean)
  
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
    

    I₂ = one(TensorValue{2,2,Float64})

    # Projeksjons‐operatorene:
    # defining the Navier operators
    Pn(n) = n ⊗ n        # normal projection operator
    Pt(n) = I₂ - Pn(n)     # tangential projection operator
    
    # da(u, du, v) = ∫( ε(v)⊙(dflux∘(ε(du), ε(u))))dΩ  + ∫(-((n_Γd ⋅ (dflux∘(ε(du), ε(u)))) ⋅ v) + (-(n_Γd ⋅ (flux∘ε(v))) ⋅ du)+(viskositet∘ε(u)/h*(du ⋅ v)))dΓd       
    a0(u, v) = ∫( ε(v)⊙(flux(ε(u))))dΩ              # sjekket med Hanna
    a1(u, v) = ∫(((n_Γd ⋅ (flux(ε(u)))) ⋅ v))dΓd     # ok
    a2(u, v) = ∫((Pn(n_Γd)⋅ u) ⋅ (n_Γd ⋅ (flux(ε(v)))))dΓd      # OK
    a3(u, v) = ∫((2* nu0/(γn*h)*(n_Γd ⋅ u)) ⋅ (n_Γd ⋅ v)  )dΓd # linje 2, ledd 1                            OK
    a4(u, v) = ∫((e/(e + γt *h) * ((Pt(n_Γd) ⋅  (n_Γd ⋅ (flux(ε(u)))))) ⋅ v))dΓd # linje 2, ledd 2            # OK
    a5(u, v) = ∫((nu0/(e + γt*h)*(Pt(n_Γd) ⋅ u)) ⋅ ( Pt(n_Γd) ⋅ v))dΓd #linje 2, ledd 3     må jeg endre her til odot?                              # OK     
    a6(u, v) = ∫((e*γt * h/(e + γt*h) * (Pt(n_Γd) ⋅  (n_Γd ⋅ (flux∘ε(u)))) ) ⋅ (2*n_Γd ⋅ ε(v)))dΓd            # OK
    a7(u, v) = ∫((nu0 * γt * h / (e + γt *h) * (Pt(n_Γd)⋅ u)) ⋅ (2* n_Γd ⋅ ε(v)))dΓd                          # OK   
    
    b(p, v) = (∫(-1*(∇ ⋅ v*p))dΩ + ∫((n_Γd ⋅ v) * p)dΓd)

    l(v, q) = (∫(f ⋅ v)dΩ     # OK
    - ∫((n_Γd ⋅ ud) ⋅ ( n_Γd ⋅ (n_Γd ⋅ flux(ε(v)))))dΓd #OK
    + ∫(((2*nu0)/(γn *h) * (n_Γd ⋅ ud)) ⋅ (n_Γd ⋅ v) )dΓd #OK
    - ∫((n_Γd ⋅ ud) ⋅ q)dΓd  #OK
    + ∫((nu0/(e + γt * h) ⋅ (Pt(n_Γd) ⋅ ud) )⋅ (Pt(n_Γd) ⋅ v) )dΓd #OK    # endrer her, og endrer her... skal dette være odot?
    - ∫((nu0 * γt * h/(e + γt * h) * Pt(n_Γd) ⋅ ud) ⋅ (2*n_Γd ⋅ ε(v)))dΓd # OK
    )
    g_u(u,v) = ( ∫( (γu1*h)*jump(n_Fg⋅∇(u))⋅jump(n_Fg⋅∇(v)) )dFg 
            +  
               ∫( (γu2*h^3)*jump_nn(u,n_Fg)⋅jump_nn(v,n_Fg) )dFg
)
    g_p(p,q) = ∫( (γp*h^3)*jump(n_Fg⋅∇(p))*jump(n_Fg⋅∇(q)) )dFg

    if stabilize # have tested with different combinations of calling the gp in the residual and jacobian, but this one seems to work best.
       # hvis jeg vil teste uten newton-løser:
        A((u,p),(v,q)) = a0(u, v) - a1(u, v) - a2(u, v) + a3(u, v) + a4(u, v) + a5(u, v) - a6(u, v) - a7(u, v) + b(p, v) - b(q, u) + g_u(u,v) + g_p(p,q)
        L((v, q)) = l(v, q)
        op = AffineFEOperator(A, L, X, Y)
        (uh, ph) = solve(op)

    else
        B((u,p),(v,q)) = a0(u, v) - a1(u, v) - a2(u, v) + a3(u, v) + a4(u, v) + a5(u, v) - a6(u, v) - a7(u, v) + b(p, v) - b(q, u) 
        M((v, q)) = l(v, q)
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
        writevtk(Ω, "C:\\Users\\Sigri\\Documents\\Master\\report\\results\\stokes\\$n $geometry $order $δ.vtu", cellfields=["u_ex" => u_exact, "uh"=>uh, "erru"=> erru, "p_ex" => p_exact, "ph"=>ph, "errp"=> errp, "nablau" => ∇(u_exact), "flux" => flux∘ε(u_exact), "viskositet" => viskositet∘ε(u_exact)]) #, "erru" => erru]) 
    end 
    return uh, u_exact, erru, l2_norm(uh - u_exact), h1_semi(uh - u_exact), ph, p_exact, errp, l2_norm(ph - p_exact), h1_semi(ph - p_exact), condition_numb, Ω
end
solver = stokes_navierBC_CutFEM
uh, u_exact, erru, ul2_norm, uh1_semi, ph, p_exact, errp, pl2_norm, ph1_semi, condition_numb, Ω_act = solver(;n, u_exact, p_exact, f, g, ud, order, geometry, βu0, γu1, γu2, γp, βp0, nu, stabilize, δ, save, calc_condition)



# ################################################################ Convergence test ##########################################################
# numb_it = 6
# uarr_l2_stab, uarr_h1_stab, parr_l2_stab, parr_h1_stab, h = convergence_stokes(;numb_it, u_exact, p_exact, f, g, ud, order, geometry, solver, δ, βu0, γu1, γu2, γp, βp0, nu, stabilize, save)


# stabilize = false
# uarr_l2_nostab, uarr_h1_nostab, parr_l2_nostab, parr_h1_nostab, h = convergence_stokes(;numb_it, u_exact, p_exact, f, g, ud, order, geometry, solver, δ, βu0, γu1, γu2, γp, βp0, nu, stabilize, save)

# # ########## velocity convergence plot:
# plot_convergence_u(
#     uarr_l2_stab, uarr_h1_stab, h;
#     uarr_l2_nostab = uarr_l2_nostab, uarr_h1_nostab = uarr_h1_nostab,
#     title_str="Convergence of CutFEM: Stokes Navier BC"
# )

# plot_convergence_p(
#     parr_l2_stab, parr_h1_stab, h;
#     parr_l2_nostab=parr_l2_nostab,
#     parr_h1_nostab=parr_h1_nostab,
#     title_str="Convergence of CutFEM: Stokes Navier BC"
# )

# # print eoc-verdier:
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

# ########### velocity convergence plot:
# plot(
#     0,
#     titlefont = 16,
#     guidefont = 14,
#     tickfont = 12
# )   

# start = 1
# #scatter!(arr_δ[idx], arr_l2u[idx], label=:"", marker=:circle, ms=4)
# plot!(arr_δ[start:end], arr_l2u[start:end],  yaxis=:log, lw=2, label="L2 stabilized")
# #scatter!(arr_δ[id2], arr_l2u_nostab[id2], label=:"", marker=:s, ms=4)
# plot!(arr_δ_nostab[start:end], arr_l2u_nostab[start:end], yaxis=:log, lw=2, label="L2 non-stabilized")
# #scatter!(arr_δ[idx], arr_h1u[idx], label=:"", marker=:circle, ms=4)
# plot!(arr_δ[start:end], arr_h1u[start:end], yaxis=:log, lw=2, label="H1 stabilized")
# #scatter!(arr_δ[id2], arr_h1u_nostab[id2], label=:"", marker=:s, ms=4)
# plot!(arr_δ_nostab[start:end], arr_h1u_nostab[start:end], yaxis=:log, lw=2, label="H1 non-stabilized")
# xlabel!("Perturbation δ")
# ylabel!("Velocity error")
# title!("Sensitivity CutFEM: Stokes with Navier BC")

# ########### pressure convergence plot:
# plot(
#     0,
#     titlefont = 16,
#     guidefont = 14,
#     tickfont = 12
# )   

# start = 1
# #scatter!(arr_δ[idx], arr_l2p[idx], label=:"", marker=:circle, ms=4)
# plot!(arr_δ[start:end], arr_l2p[start:end],  yaxis=:log, lw=2, label="L2 stabilized")
# #scatter!(arr_δ[id2], arr_l2p_nostab[id2], label=:"", marker=:s, ms=4)
# plot!(arr_δ_nostab[start:end], arr_l2p_nostab[start:end], yaxis=:log, lw=2, label="L2 non-stabilized")
# #scatter!(arr_δ[idx], arr_h1p[idx], label=:"", marker=:circle, ms=4)
# plot!(arr_δ[start:end], arr_h1p[start:end], yaxis=:log, lw=2, label="H1 stabilized")
# #scatter!(arr_δ[id2], arr_h1p_nostab[id2], label=:"", marker=:s, ms=4)
# plot!(arr_δ_nostab[start:end], arr_h1p_nostab[start:end], yaxis=:log, lw=2, label="H1 non-stabilized")
# xlabel!("Perturbation δ")
# ylabel!("Pressure error")
# title!("Sensitivity CutFEM: Stokes with Navier BC")