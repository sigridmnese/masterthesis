
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

"""Fitted FEM with weak boundary imposition using Nitsches method. Symmetric p-stokes equations with Navier boundary conditions. """

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
n = 32

# solver:
#fitted FEM stokes solver:
function stokes_FEM_navierBC_newton(;n, u_exact, p_exact, f, g, ud, order, geometry, βu0, γu1, γu2, γp, βp0, nu, stabilize, δ, save = false, calc_condition = false)
    """
    I solve the regular stokes equations with Nitsche's method for weak boundary imposition, using Newtons method to test.
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
    domain = (0,1,0,1)
    partition = (n,n)
    model = CartesianDiscreteModel(domain, partition)    
    h = 1/n
    γn = 0.05  # hentet fra Josefin sin artikkel.
    γt = 0.05  # hentet fra Josefin sin artikkel.
    e = 0        # e er det som er epsilon i artikkelen til Josefin. Denne er 0 for unit square. 

    # Embedded boundary
    # Dirichlet conditions on u
    Γd = BoundaryTriangulation(model)
    n_Γd = get_normal_vector(Γd)

    # lag tangent‐felt på Γd:
    t_Γd = CellField(x -> VectorValue(-n_Γd(x)[2],
                                   n_Γd(x)[1]),Γd)
    
    order = 2
    reffeᵤ = ReferenceFE(lagrangian,VectorValue{2,Float64},order)
    reffeₚ = ReferenceFE(lagrangian,Float64,order-1;space=:P)
  
    # definert de samme spacene som for stong BC implementering, bare fjernet at V på være lik ud på grensen
    V = TestFESpace(model,reffeᵤ,conformity=:H1) 
    Q = TestFESpace(model,reffeₚ,conformity=:L2,constraint=:zeromean)
    Y = MultiFieldFESpace([V,Q])
  
    U = TrialFESpace(V)         # har fjernet begrensningen til grensebetingelsen ud her..
    P = TrialFESpace(Q)
    X = MultiFieldFESpace([U,P])

    degree = order
    Ω = Triangulation(model)
    dΩ = Measure(Ω,degree)
    dΓd = Measure(Γd, degree)

    I₂ = one(TensorValue{2,2,Float64})

    # Projeksjons‐operatorene:
    # defining the Navier operators
    Pn(n) = n ⊗ n        # normal projection operator
    Pt(n) = I₂ - Pn(n)     # tangential projection operator

    ######################## kvadratish domain - men bruker nitsche implementasjon. Hvilke spaces skal jeg bruke da? #########################
    
    a0(u, v) = ∫( ε(v)⊙(flux(ε(u))))dΩ              # sjekket med Hanna
    a1(u, v) = ∫(((n_Γd ⋅ (flux(ε(u)))) ⋅ v))dΓd     # ok
    a2(u, v) = ∫((Pn(n_Γd)⋅ u) ⋅ (n_Γd ⋅ (flux(ε(v)))))dΓd      # OK
    a3(u, v) = ∫((2* nu0/(γn*h)*(n_Γd ⋅ u)) ⋅ (n_Γd ⋅ v)  )dΓd # linje 2, ledd 1                            OK
    a4(u, v) = ∫((e/(e + γt *h) * ((Pt(n_Γd) ⋅  (n_Γd ⋅ (flux(ε(u)))))) ⋅ v))dΓd # linje 2, ledd 2            # OK
    a5(u, v) = ∫((nu0/(e + γt*h)*(t_Γd ⋅ u)) ⋅ ( t_Γd ⋅ v))dΓd #linje 2, ledd 3                               # OK     
    a6(u, v) = ∫((e*γt * h/(e + γt*h) * (Pt(n_Γd) ⋅  (n_Γd ⋅ (flux∘ε(u)))) ) ⋅ (2*n_Γd ⋅ ε(v)))dΓd            # OK
    a7(u, v) = ∫((nu0 * γt * h / (e + γt *h) * (Pt(n_Γd)⋅ u)) ⋅ (2* n_Γd ⋅ ε(v)))dΓd                          # OK   
    
    b(p, v) = (∫(-1*(∇ ⋅ v*p))dΩ + ∫((n_Γd ⋅ v) * p)dΓd)

    l(v, q) = (∫(f ⋅ v)dΩ     # OK
    - ∫((n_Γd ⋅ ud) ⋅ ( n_Γd ⋅ (n_Γd ⋅ flux(ε(v)))))dΓd #OK
    + ∫(((2*nu0)/(γn *h) * (n_Γd ⋅ ud)) ⋅ (n_Γd ⋅ v) )dΓd #OK
    - ∫((n_Γd ⋅ ud) ⋅ q)dΓd  #OK
    + ∫((nu0/(e + γt * h) ⋅ (t_Γd ⋅ ud) )⋅ (t_Γd ⋅ v) )dΓd #OK
    - ∫((nu0 * γt * h/(e + γt * h) * Pt(n_Γd) ⋅ ud) ⋅ (2*n_Γd ⋅ ε(v)))dΓd # OK
    )

    # Har nå sjekket alle ledd i a og l sammen med Hanna, og de er alle OK.
    # til newton-løser:
    # res((u,p),(v,q)) = a0(u, v) - a1(u, v) - a2(u, v) + a3(u, v) + a4(u, v) + a5(u, v) - a6(u, v) - a7(u, v) + b(p, v) - b(q, u) -l(v, q)
    # jac((u, p), (du, dp), (v, q)) = a0(du, v) - a1(du, v) - a2(du, v) + a3(du, v) + a4(du, v) + a5(du, v) - a6(du, v) - a7(du, v) + b(dp, v) - b(q, du)
    
    # op = FEOperator(res, jac, X, Y)

    # # non-linear phase
    # nls = NLSolver(
    # show_trace=true, method=:newton, linesearch=BackTracking(), iterations=50)      #prøver å legge inn et max antall iterasjoner og en lav toleranse      
    # solver = FESolver(nls)

    # (uh, ph) = solve(solver, op)

    # hvis jeg vil teste uten newton-løser:
    res((u,p),(v,q)) = a0(u, v) - a1(u, v) - a2(u, v) + a3(u, v) + a4(u, v) + a5(u, v) - a6(u, v) - a7(u, v) + b(p, v) - b(q, u)
    jac((v, q)) = l(v, q)

    op = AffineFEOperator(res, jac, X, Y)
    (uh, ph) = solve(op)

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
        writevtk(Ω, "C:\\Users\\Sigri\\Documents\\Master\\report\\results\\stokes\\$n $geometry $order.vtu", cellfields=["u_ex" => u_exact, "uh"=>uh, "erru"=> erru, "p_ex" => p_exact, "ph"=>ph, "errp"=> errp, "nablau" => ∇(u_exact), "visk" => viskositet∘ε(u_exact)]) #, "erru" => erru]) 
    end
    return uh, u_exact, erru, l2_norm(uh - u_exact), h1_semi(uh - u_exact), ph, p_exact, errp, l2_norm(ph - p_exact), h1_semi(ph - p_exact), condition_numb, Ω
end

# solver parameters
stabilize = true
solver = stokes_FEM_navierBC_newton
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
numb_it = 6
uarr_l2, uarr_h1, parr_l2, parr_h1, h = convergence_stokes(;numb_it, u_exact, p_exact, f, g, ud, order, geometry, solver, δ, βu0, γu1, γu2, γp, βp0, nu, stabilize, save)

########## velocity convergence plot:

plot_convergence_u(uarr_l2,uarr_h1, h)

plot(
    0,
    titlefont = 16,
    guidefont = 14,
    tickfont = 12
)

plot!(h, uarr_l2, xaxis=:log, yaxis=:log, marker=:o, lw=2, label="L2")
plot!(h, uarr_h1, marker=:o, lw=2, label="H1")
xlabel!("Mesh size h")
ylabel!("Velocity error")
title!("Convergence of Stokes fitted FEM with Navier BC")

# # ########### pressure convergence plot:
plot(
    0,
    titlefont = 16,
    guidefont = 14,
    tickfont = 12
)

plot!(h, parr_l2, xaxis=:log, yaxis=:log, marker=:o, lw=2, label="L2")
plot!(h, parr_h1, marker=:o, lw=2, label="H1")
xlabel!("Mesh size h")
ylabel!("Pressure error")
title!("Convergence of Stokes fitted FEM with Navier BC")
