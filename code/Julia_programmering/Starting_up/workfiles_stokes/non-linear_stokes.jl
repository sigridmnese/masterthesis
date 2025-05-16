# prøver nå å bare bruke denne funksjonen som jeg har definert, kan være man må gjøre lit andre ing for at det skal fungere. 
#nu(du) = 1/2 * A^(1-r) * (1/2 * norm(du, 2))^((r-2)/2)      #for eksempel... Kan være man må ha noe annet her?
using Gridap
using GridapEmbedded
using STLCutters
using LinearAlgebra
using Plots
using LineSearches: BackTracking
import Random
using Symbolics


include("C:\\Users\\Sigri\\Documents\\Master\\report\\code\\Julia_programmering\\Starting_up\\utils\\utils.jl")

##### Divergensfri stokes løser ####
# hva hvis det ikke er divergensfritt?
nu = 1
r = 4/3
A = 1
# endrer nå alle symbolder før uder det er nabla til epsilon...

u_exact(x) = VectorValue(-x[2], x[1])   #VectorValue(2*x[1] + cos(2*π*x[2]), -2*x[2] + sin(2*π*x[1]))
p_exact(x) = sin(2*π*x[1])*cos(2*π*x[2])
flux(∇u) = norm(∇u)^(p-2) * ∇u
f(x) =  -divergence(flux∘ε(u_exact))(x) + ∇(p_exact)(x)      # prøver å endre f her...
ud(x) = u_exact(x)
#1/2 * A^(1-r) * (1/2 * norm(du, 2))^((r-2)/2) * ε(du)

domain = "circle"
n = 128
γ = 10* 2*2 
β_1 = 1
β_2 = 1
β_3 = 0.1
γ=10*2*2

g = VectorValue(0.0, 0.0)

function stokes_FEM(;n, u_exact, p_exact, f, g, ud, order, geometry, βu0, γu1, γu2, γp, βp0, nu, stabilize, δ, save = false, calc_condition = false)
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
    domain = (0,1,0,1)
    partition = (n,n)
    model = CartesianDiscreteModel(domain, partition)    

    labels = get_face_labeling(model)
    
    # alle grensetagsene får dirichlet = g
    add_tag_from_tags!(labels, "dirig", [1, 2, 3, 4, 5, 6, 7, 8])
    order = 2
    reffeᵤ = ReferenceFE(lagrangian,VectorValue{2,Float64},order)
    reffeₚ = ReferenceFE(lagrangian,Float64,order-1;space=:P)
  
    V = TestFESpace(model,reffeᵤ,labels=labels,dirichlet_tags="dirig",conformity=:H1)
    Q = TestFESpace(model,reffeₚ,conformity=:L2,constraint=:zeromean)
    Y = MultiFieldFESpace([V,Q])
  
    U = TrialFESpace(V,ud)
    P = TrialFESpace(Q)
    X = MultiFieldFESpace([U,P])

    degree = order
    Ω = Triangulation(model)
    dΩ = Measure(Ω,degree)
    Γ = Boundary(Ω)
    dΓ = Measure(Γ,degree)
     
    # Weak formulation components
    #f = VectorValue(0.0,0.0)
    a(u, v) = ∫( ∇(v)⊙∇(u))dΩ
    b(v, p) = ∫(-(∇⋅v)*p )dΩ
    #mean_p = integrate(p_exact,dΩ)     # skal bli ~ 0
    #@show mean_p
    println(sum( ∫( p_exact) * dΩ ))
    #println(integrate(p_exact, dΩ))
    #println(∫(∇(p_exact))*dΩ)

    l2_norm(u) = (sum( ∫( u ⋅ u )*dΩ ))
    h1_semi(u) = sum(∫(∇(u) ⊙ ∇(u))*dΩ)
      
    l1(v) = ∫(f ⋅ v)dΩ
    l2(v) = 0 #∫(-1* (n_Γd ⋅ ε(v)) ⋅ ud + γ/h ⋅ v ⋅ ud )dΓd
    l3(q) = 0 #∫(n_Γd ⋅ ud *q)dΓd

    if stabilize
        A((u,p),(v,q)) =(a(u, v) + b(v, p) - b(u, q))
        L((v, q)) = l1(v) #+ l2(v) + l3(q)
        op = AffineFEOperator(A,L,X,Y)
        uh, ph = solve(op)
    else
        B((u,p),(v,q)) = a(u,v) + b(v, p) - b(u, q) 
        M((v,q)) = l1(v) + l2(v) + l3(q)
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
    return uh, u_exact, erru, l2_norm(uh - u_exact), h1_semi(uh - u_exact), ph, p_exact, errp, l2_norm(ph - p_exact), h1_semi(ph - p_exact), condition_numb, Ω
end

function nonlinear_stokes_solver(;n, u_exact, p_exact, f, g, ud, order, geometry, βu0, γu1, γu2, γp, βp0, nu, stabilize, δ, save = false, calc_condition = false)
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
    domain = (0,1,0,1)
    partition = (n,n)
    model = CartesianDiscreteModel(domain, partition)    

    labels = get_face_labeling(model)
    
    # alle grensetagsene får dirichlet = g
    add_tag_from_tags!(labels, "dirig", [1, 2, 3, 4, 5, 6, 7, 8])
    order = 2
    reffeᵤ = ReferenceFE(lagrangian,VectorValue{2,Float64},order)
    reffeₚ = ReferenceFE(lagrangian,Float64,order-1;space=:P)
  
    V = TestFESpace(model,reffeᵤ,labels=labels,dirichlet_tags="dirig",conformity=:H1)
    Q = TestFESpace(model,reffeₚ,conformity=:L2,constraint=:zeromean)
    Y = MultiFieldFESpace([V,Q])
  
    U = TrialFESpace(V,ud)
    P = TrialFESpace(Q)
    X = MultiFieldFESpace([U,P])

    degree = order
    Ω = Triangulation(model)
    dΩ = Measure(Ω,degree)
    
    println(sum( ∫( p_exact) * dΩ ))

    l2_norm(u) = (sum( ∫( u ⋅ u )*dΩ ))
    h1_semi(u) = sum(∫(∇(u) ⊙ ∇(u))*dΩ)
    
    # fra klassisk Stokes FEM, så er det kun herfra og ned som er endret:)

    nu0 = 1
    ϵ_0 = 1e-6
    
    a(u, v) = ∫( ∇(v)⊙(flux∘∇(u)))dΩ
    b(v, p) = ∫(-(∇⋅v)*p )dΩ
    l(v) = ∫(f ⋅ v)dΩ

    # setting up non-linear system. Fluxen er bare en funksjon av du, men jakobianen må kanskje være en funksjon av 
    #flux(du) = nu0 *(ϵ_0 + norm(du)^2)^((r-2)/2) ⋅ du
    #res((u,v), (p, q)) = ∫( ∇(v)⊙(flux∘∇(u)) - v*f )*dΩ  + b(v, p) - b(u, q)       # bytte til epsilon her, og legge til uttrykkene for b(u, q), b(v, p)

    dflux(∇du,∇u) = nu0 *( (ϵ_0 + norm(∇u)^2)^((r-2)/2) * ∇du + (r-2) * (ϵ_0 + norm(∇u)^2)^((r-4)/2)*(∇u⊙∇du)⊙∇u)
    #jac(u, du, v) = ∫(∇(v)⊙(∇(du)⊙∇(u)) + ∇(v)⊙∇(flux∘∇(u))) * dΩ
    #jac((u, p),(du, dp),(v, q)) = ∫( ∇(v)⊙(dflux∘(∇(du),∇(u))) )*dΩ         # vet at jac egentlig ikke er riktig...
    
    res((u,v), (p, q)) = a(u, v)  + b(v, p) - b(u, q) - l(v)       # bytte til epsilon her, og legge til uttrykkene for b(u, q), b(v, p)
    
    # definer jacobian her:
    jac((u, p), (du, dp), (v, q)) = a(du, v) + b(v, dp) - b(du, q) + ∫(∇(v)⊙dflux∘∇(u))*dΩ

    op = FEOperator(res, jac, X, Y)

    # non-linear phase
    nls = NLSolver(
    show_trace=true, method=:newton, linesearch=BackTracking())
    solver = FESolver(nls)

    # if we need an iitial guess:
    #Random.seed!(1234)
    #x = rand(Float64,num_free_dofs(U))
    #y = rand(Float64,num_free_dofs(P))
    #uh0 = FEFunction(U,x)
    #ph0 = FEFunction(P,y)

    #(uh, ph), = solve!((uh0, ph0),solver,op)

    #op = AffineFEOperator(A,L,X,Y)
    #uh, ph = solve(op)
  
    uh, ph = solve(solver, op)

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
    return uh, u_exact, erru, l2_norm(uh - u_exact), h1_semi(uh - u_exact), ph, p_exact, errp, l2_norm(ph - p_exact), h1_semi(ph - p_exact), condition_numb, Ω
end

stabilize = true
δ = 0
save = true
calc_condition = false
order = 2
geometry = "circle"
βu0 = 1
γu1 = 10*2*2
γu2 = 10*2*2
γp = 10*2*2
βp0 = 1
n = 128
#uh, u_exact, erru, ul2_norm, uh1_semi, ph, p_exact, errp, pl2_norm, ph1_semi, condition_numb, Ω_act = stokes_FEM(;n, u_exact, p_exact, f, g, ud, order, geometry, βu0, γu1, γu2, γp, βp0, nu, stabilize, δ, save, calc_condition)
# numb_it = 6
# solver = stokes_FEM
#uarr_l2, uarr_h1, parr, parr_h1, h = convergence_stokes(;numb_it, u_exact, p_exact, f, g, ud, order, geometry, solver, δ, βu0, γu1, γu2, γp, βp0, nu, stabilize, save)

# plotter
# plot(
#     0,
#     title = "Convergence of Stokes FEM",
#     xlabel = "Mesh size h",
#     ylabel = "Velocity error",
#     titlefont = 16,
#     guidefont = 14,
#     tickfont = 12
# )
# plot!(h, uarr_l2, xaxis=:log, yaxis=:log, marker=:o, lw=2, label="L2 stabilized")
# plot!(h, uarr_h1, marker=:o, lw=2, label="H1 stabilized")
# xlabel!("Mesh size h")
# ylabel!("Error")
# title!("Convergence of Stokes FEM")

#plot!(h[start:end], uarr_l2_1_nostab[start:end], marker=:s, lw=2, label="L2 non-stabilized")
#plot!(h[start:end], uarr_h1_1_nostab[start:end], marker=:s, lw=2, label="H1 non-stabilized")

# # # Legger til aksetitler og tittel


##################### herfra prøver jeg å løse ikke-lineær stokes ######################
# med de samme parametrene som over
uh, u_exact, erru, ul2_norm, uh1_semi, ph, p_exact, errp, pl2_norm, ph1_semi, condition_numb, Ω_act = nonlinear_stokes_solver(;n, u_exact, p_exact, f, g, ud, order, geometry, βu0, γu1, γu2, γp, βp0, nu, stabilize, δ, save)

