## In this programme, we use the STL mesh to define a 3-dimensional mesh illustraing a glacier domain. We solve the linear stokes equations using the cutFEM method. 

# Kan det være at jeg får store feil i hjørnene nå fordi det ikke eksisterer noen normalvektor til disse punktene???
# spørre André i morgen.

using Gridap
using GridapEmbedded
using STLCutters
using LinearAlgebra
using Plots
using LineSearches: BackTracking
import Random
using Logging
using LoggingExtras
using Test

include("C:\\Users\\Sigri\\Documents\\Master\\report\\code\\Julia_programmering\\Starting_up\\utils\\utils.jl")
include("C:\\Users\\Sigri\\Documents\\Master\\report\\code\\Julia_programmering\\Starting_up\\workfiles_stokes\\testing_non-linear_stokes.jl")

"""Fitted FEM with weak boundary imposition using Nitsches method. Symmetric stokes equations with Navier boundary conditions. """

# Defining constants
nu0 =  1   # klarer ikke mindre enn 0.01 og klarer heller ikke større enn 100 (størrelsesordener)
r = 4/3
A = nu0
ϵ_0 = 1e-6       

# Defining manufactured solutions
u_exact(x) =  VectorValue(2*x[1] + cos(2*π*x[2]), -2*x[2] + sin(2*π*x[1]), 0)#VectorValue(2*x[1] + exp(x[1]/2) * cos(2*π*x[2]), -2*x[2] + exp(x[2]/2) * sin(2*π*x[1])) #bytte til sin/cos-uttrykk  VectorValue(-x[2], x[1])
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

f(x) =  -divergence(flux∘(u_exact))(x) + ∇(p_exact)(x)      # prøver å endre f her...
u_D(x) = u_exact(x)
g(x) = tr(ε(u_exact)(x))

n = 32

# solver:
#fitted FEM stokes solver:
function stokes_cutFEM_3d(;n, u_exact, p_exact, f, g, ud, order, geometry, βu0, γu1, γu2, γp, βp0, nu, stabilize, δ, save = false, calc_condition = false)
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
    # Define background mesh and geometry
    p0 = Point(0.0, 0.0, 0.0)
    # Define background mesh
    n = 15
    partition = (n, n, n)
    dim = length(partition)
    pmin = p0+Point(-1.0, -1.0, 0.0)
    # pmin = p0+Point(-1.0, 0.0, 0.0)
    pmax = p0+Point(1.0, 1.0, 1.0)
    # mesh size
    h = (pmax-pmin)[1]/partition[1]

    # nitsche parameters 
    γn = 0.05  # hentet fra Josefin sin artikkel. Endre disse når du endrer geometrien!
    γt = 0.05  # hentet fra Josefin sin artikkel. Endre disse når du endrer geometrien!
    e = 0 

    bgmodel = CartesianDiscreteModel(pmin, pmax, partition)
    Ω_bg = Triangulation(bgmodel)
    writevtk(Ω_bg, "mesh_bg")

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

    # lag tangent‐felt på Γb:
    t_Γb = CellField(x -> VectorValue(-n_Γb(x)[2],
                                   n_Γb(x)[1]),Γb)

    # Outer mesh boundary 
    # Has Dirichlet for pF and no-stress for u
    Γs = BoundaryTriangulation(cutgeo_facets, PHYSICAL_IN)
    n_Γs = get_normal_vector(Γs)

    # lag tangent‐felt på Γd:
    t_Γs = CellField(x -> VectorValue(-n_Γs(x)[2],
                                   n_Γs(x)[1]),Γs)

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

    # Write out domains 
    # writevtk(Ω_act, "mesh_act")
    # writevtk(Ω, "mesh")
    # writevtk(Γb, "surface_gamma_d")
    # writevtk(Γs, "surface_gamma_s")
    # writevtk(Fg, "ghost_facets")

    # %% Define weak formulation
    # Define function spaces, Taylor Hood elements
    reffe_u  = ReferenceFE(lagrangian,VectorValue{dim, Float64},order_u)
    reffe_p = ReferenceFE(lagrangian,Float64, order_p)

    V = TestFESpace(Ω_act, reffe_u,  conformity=:H1)
    Q = TestFESpace(Ω_act, reffe_p, conformity=:H1, constraint=:zeromean)

    U = TrialFESpace(V)
    P = TrialFESpace(Q)

    X = MultiFieldFESpace([U, P])
    Y = MultiFieldFESpace([V, Q])

    # Nitsche and Ghost penalty stabilization parameter
    βu = 10.0*order_u^2
    γu1 = 1.0
    γu2 = 1.0
    γp = 0.1
    βp = 0.1

    I3 = one(TensorValue{3,3,Float64})

    # Projeksjons‐operatorene:
    # defining the Navier operators
    Pn(n) = n ⊗ n        # normal projection operator
    Pt(n) = I3 - Pn(n)     # tangential projection operator

    # TODO: For now we impose Dirichlet boundary conditions on the whole boundary, so need to change b.c. on Γs
    # Bilinear forms
    a(u,v) =  ( ∫( ε(u)⊙ε(v))dΩ )


    a_Γb(u, v) = (∫((-(n_Γb ⋅ (flux(ε(u)))) ⋅ v))dΓb     # ok
    - ∫((Pn(n_Γb)⋅ u) ⋅ (n_Γb ⋅ (flux(ε(v)))))dΓb      # OK
    + ∫((2* nu0/(γn*h)*(n_Γb ⋅ u)) ⋅ (n_Γb ⋅ v)  )dΓb # linje 2, ledd 1                            OK
    + ∫((e/(e + γt *h) * ((Pt(n_Γb) ⋅  (n_Γb ⋅ (flux(ε(u)))))) ⋅ v))dΓb # linje 2, ledd 2            # OK
    + ∫((nu0/(e + γt*h)*(t_Γb ⋅ u)) ⋅ ( t_Γb ⋅ v))dΓb #linje 2, ledd 3                               # OK     
    - ∫((e*γt * h/(e + γt*h) * (Pt(n_Γb) ⋅  (n_Γb ⋅ (flux∘ε(u)))) ) ⋅ (2*n_Γb ⋅ ε(v)))dΓb            # OK
    - ∫((nu0 * γt * h / (e + γt *h) * (Pt(n_Γb)⋅ u)) ⋅ (2* n_Γb ⋅ ε(v)))dΓb                          # OK   
    )

    a_Γs(u, v) = (∫(((-n_Γs ⋅ (flux(ε(u)))) ⋅ v))dΓs     # ok
    - ∫((Pn(n_Γs)⋅ u) ⋅ (n_Γs ⋅ (flux(ε(v)))))dΓs      # OK
    + ∫((2* nu0/(γn*h)*(n_Γs ⋅ u)) ⋅ (n_Γs ⋅ v)  )dΓs # linje 2, ledd 1                            OK
    + ∫((e/(e + γt *h) * ((Pt(n_Γs) ⋅  (n_Γs ⋅ (flux(ε(u)))))) ⋅ v))dΓs # linje 2, ledd 2            # OK
    + ∫((nu0/(e + γt*h)*(t_Γs ⋅ u)) ⋅ ( t_Γs ⋅ v))dΓs #linje 2, ledd 3                               # OK     
    - ∫((e*γt * h/(e + γt*h) * (Pt(n_Γs) ⋅  (n_Γs ⋅ (flux∘ε(u)))) ) ⋅ (2*n_Γs ⋅ ε(v)))dΓs            # OK
    - ∫((nu0 * γt * h / (e + γt *h) * (Pt(n_Γs)⋅ u)) ⋅ (2* n_Γs ⋅ ε(v)))dΓs                          # OK   
    )
    

    b(v, p) = ( ∫(-1*(∇⋅v*p))dΩ 
                + ∫((n_Γb⋅v)*p)dΓb
                + ∫(v⋅n_Γs*p)dΓs
    )

    g_u(u,v) = ( ∫( (γu1*h)*jump(n_Fg⋅∇(u))⋅jump(n_Fg⋅∇(v)) )dFg 
                +  
                ∫( (γu2*h^3)*jump_nn(u,n_Fg)⋅jump_nn(v,n_Fg) )dFg
    )

    g_p(p,q) = ∫( (γp*h^3)*jump(n_Fg⋅∇(p))*jump(n_Fg⋅∇(q)) )dFg
    s_p(p,q) = ∫( (βp*h^3)*jump(n_Fi⋅∇(p))*jump(n_Fi⋅∇(q)) )dFi
    


    # sjekker for stabilisering:
    if stabilize
      A((u, p), (v, q)) = a(u,v) + a_Γb(u,v) +  a_Γs(u,v) + b(v,p) + b(u,q) + g_u(u,v) - g_p(p,q) #- s_p(p,q) trenger ikke den siste her fordi jeg skal bruke taylor hood elementer
      # Linear forms
      Lb(v, q) = (∫(f ⋅ v)dΩ     # OK
    - ∫((n_Γb ⋅ ud) ⋅ ( n_Γb ⋅ (n_Γb ⋅ flux(ε(v)))))dΓb #OK
    + ∫(((2*nu0)/(γn *h) * (n_Γb ⋅ ud)) ⋅ (n_Γb ⋅ v) )dΓb #OK
    - ∫((n_Γb ⋅ ud) ⋅ q)dΓb  #OK
    + ∫((nu0/(e + γt * h) ⋅ (t_Γb ⋅ ud) )⋅ (t_Γb ⋅ v) )dΓb #OK
    - ∫((nu0 * γt * h/(e + γt * h) * Pt(n_Γb) ⋅ ud) ⋅ (2*n_Γb ⋅ ε(v)))dΓb # OK
    )

    #legger inn samme bc på Γb som på Γs
    Ls(v, q) = (∫(f ⋅ v)dΩ     # OK
    - ∫((n_Γs ⋅ ud) ⋅ ( n_Γs ⋅ (n_Γs ⋅ flux(ε(v)))))dΓs #OK
    + ∫(((2*nu0)/(γn *h) * (n_Γs ⋅ ud)) ⋅ (n_Γs ⋅ v) )dΓs #OK
    - ∫((n_Γs ⋅ ud) ⋅ q)dΓs  #OK
    + ∫((nu0/(e + γt * h) ⋅ (t_Γs ⋅ ud) )⋅ (t_Γs ⋅ v) )dΓs #OK
    - ∫((nu0 * γt * h/(e + γt * h) * Pt(n_Γs) ⋅ ud) ⋅ (2*n_Γs ⋅ ε(v)))dΓs # OK
    )

      L((v,q)) = ( ∫(f⋅v - g⋅q)dΩ) + Lb(v, q) + Ls(v, q)

      # %% Solve and write out results
      # Solve
      op = AffineFEOperator(A,L,X,Y)
      uh, ph = solve(op)

    else
      B((u, p), (v, q)) = a(u,v) + a_Γb(u,v) + b(v,p) + b(u,q) 
      # Linear forms
      M((v,q)) = ( ∫(f⋅v - g⋅q)dΩ
      # lh((v,q)) = ( ∫(fdiv⋅v - g⋅q)dΩ
                #  + ∫( -1*(n_Γb⋅∇(v))⋅u_D + βu/h*(u_D⋅v) + u_D⋅n_Γb*q)dΓb 
                + ∫( -1*(n_Γb⋅ε(v))⋅u_D + βu/h*(u_D⋅v) + u_D⋅n_Γb*q)dΓb 
                + ∫( -1*(n_Γs⋅ε(v))⋅u_D + βu/h*(u_D⋅v) + u_D⋅n_Γs*q)dΓs
    )

      # %% Solve and write out results
      # Solve
      op = AffineFEOperator(B, M, X, Y)
      uh, ph = solve(op)
    end


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
        outputfile = "stokes"
        if outputfile !== nothing
        writevtk(Ω,"C:\\Users\\Sigri\\Documents\\Master\\report\\results\\3d_stokes\\$n $geometry $order.vtu",cellfields=["uh"=>uh,
                                         "ph"=>ph, 
                                         "u_ex"=> u_exact, 
                                         "p_ex"=> p_exact,
                                         "erru" => erru,
                                         "errp" => errp,
                                         "f"=> f,
                                         "g"=> g,
                                            "nablau" => ∇(u_exact),
                                            "visk" => viskositet∘ε(u_exact)
                                         ])
        end
        #writevtk(Ω, "C:\\Users\\Sigri\\Documents\\Master\\report\\results\\stokes\\$n $geometry $order.vtu", cellfields=["u_ex" => u_exact, "uh"=>uh, "erru"=> erru, "p_ex" => p_exact, "ph"=>ph, "errp"=> errp, "nablau" => ∇(u_exact), "visk" => viskositet∘ε(u_exact)]) #, "erru" => erru]) 
    end
    return uh, u_exact, erru, l2_norm(uh - u_exact), h1_semi(uh - u_exact), ph, p_exact, errp, l2_norm(ph - p_exact), h1_semi(ph - p_exact), condition_numb, Ω
end

# solver parameters
stabilize = true
solver = stokes_cutFEM_3d
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



# plot_convergence_p(parr_l2,parr_h1, h)




# plot(
#     0,
#     titlefont = 16,
#     guidefont = 14,
#     tickfont = 12
# )

# plot!(h, uarr_l2, xaxis=:log, yaxis=:log, marker=:o, lw=2, label="L2 stabilized")
# plot!(h, uarr_h1, marker=:o, lw=2, label="H1 stabilized")
# xlabel!("Mesh size h")
# ylabel!("Velocity error")
# title!("Convergence of p-stokes FEM")

# # # ########### pressure convergence plot:
# plot(
#     0,
#     titlefont = 16,
#     guidefont = 14,
#     tickfont = 12
# )

# plot!(h, parr_l2, xaxis=:log, yaxis=:log, marker=:o, lw=2, label="L2 stabilized")
# plot!(h, parr_h1, marker=:o, lw=2, label="H1 stabilized")
# xlabel!("Mesh size h")
# ylabel!("Pressure error")
# title!("Convergence of p-stokes FEM")
