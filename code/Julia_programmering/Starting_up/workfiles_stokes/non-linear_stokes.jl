# prøver nå å bare bruke denne funksjonen som jeg har definert, kan være man må gjøre lit andre ing for at det skal fungere. 
#nu(du) = 1/2 * A^(1-r) * (1/2 * norm(du, 2))^((r-2)/2)      #for eksempel... Kan være man må ha noe annet her?
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

##### Divergensfri stokes løser ####
# hva hvis det ikke er divergensfritt?
nu0 = 1
r = 4/3
A = 1
ϵ_0 = 1e-6          # hadde glemt å legge inn en epsilon0 i flux-termen. Newton-løseren divergerte
# endrer nå alle symbolder før uder det er nabla til epsilon...


# non-linear stokes problem
u_exact(x) = VectorValue(-x[2], x[1])   #VectorValue(2*x[1] + cos(2*π*x[2]), -2*x[2] + sin(2*π*x[1]))
p_exact(x) = sin(2*π*x[1])*cos(2*π*x[2])
flux(∇u) = nu0*(ϵ_0 + norm(∇u)^2)^((r-2)/2) * ∇u
f(x) =  -divergence(flux∘∇(u_exact))(x) + ∇(p_exact)(x)      # prøver å endre f her...
ud(x) = u_exact(x)
#1/2 * A^(1-r) * (1/2 * norm(du, 2))^((r-2)/2) * ε(du)

domain = "circle"
n = 16
γ = 10* 2*2 
β_1 = 1
β_2 = 1
β_3 = 0.1
γ=10*2*2

g = VectorValue(0.0, 0.0)

# function p_stokes_FEM(;n, u_exact, p_exact, f, g, ud, order, geometry, βu0, γu1, γu2, γp, βp0, nu, stabilize, δ, save = false, calc_condition = false)
#        """
#     Fitted FEM, non-linear stokes (p-stokes). Using P2-P1 Taylor-Hood elements.  
#     n: number of grid elements. Powers of 2 for simplicity and convergence estimates.
#     u_exact: exact solution for method of manufactured solutions
#     order: order of polynomial degree. 
#     f: lhs for first term, -Δ u_ex + ∇p = f
#     g: lhs for second term u = g
#     geometry: optional between "Circle", "Flower", "Heart", "Glacier".
#     stabilize: wheather to add the stabilization term or not
#     δ: perturbation of cut
#     """
#     # Define background mesh
#     domain = (0,1,0,1)
#     partition = (n,n)
#     model = CartesianDiscreteModel(domain, partition)    

#     labels = get_face_labeling(model)
    
#     # alle grensetagsene får dirichlet = g
#     add_tag_from_tags!(labels, "dirig", [1, 2, 3, 4, 5, 6, 7, 8])
#     order = 2
#     reffeᵤ = ReferenceFE(lagrangian,VectorValue{2,Float64},order)
#     reffeₚ = ReferenceFE(lagrangian,Float64,order-1;space=:P)

#     # taylor hood elements
#     V = TestFESpace(model,reffeᵤ,labels=labels,dirichlet_tags="dirig",conformity=:H1)
#     Q = TestFESpace(model,reffeₚ,conformity=:L2,constraint=:zeromean)
#     Y = MultiFieldFESpace([V,Q])
  
#     U = TrialFESpace(V,ud)
#     P = TrialFESpace(Q)
#     X = MultiFieldFESpace([U,P])

#     degree = order
#     Ω = Triangulation(model)
#     dΩ = Measure(Ω,degree)
    
#     println(sum( ∫( p_exact) * dΩ ))

#     l2_norm(u) = (sum( ∫( u ⋅ u )*dΩ ))
#     h1_semi(u) = sum(∫(∇(u) ⊙ ∇(u))*dΩ)
    
#     # fra klassisk Stokes FEM, så er det kun herfra og ned som er endret:)
#     nu0 = 1
#     ϵ_0 = 1e-6
    
#     a(u, v) = ∫( ∇(v)⊙(flux∘∇(u)))dΩ
#     b(v, p) = ∫(-(∇⋅v)*p )dΩ
#     l(v) = ∫(f ⋅ v)dΩ

#     # dflux calculated the same way as in the notebook p-Laplace...
#     dflux(∇du,∇u)=(r-2)*(ϵ_0 + norm(∇u)^2)^((r-4)/2)*(∇u⊙∇du) ⋅ ∇u + (ϵ_0 + norm(∇u)^2)^((r-2)/2)*∇du

#     # and introduced in the same way as in the notebook p-Laplace in the bilinear form a...
#     da(u, du, v) = ∫(∇(v)⊙(dflux∘(∇(du), ∇(u))))dΩ
    
#     # and then the Newton multifield system is assembled as in the Navier Stokes notebook...
#     res((u,p),(v,q)) = ∫( ∇(v)⊙(flux∘∇(u)))dΩ + ∫(-(∇⋅v)*p )dΩ - ∫(-(∇⋅u)*q )dΩ - ∫(f ⋅ v)dΩ    #a(u, v)  + b(v, p) - b(u, q) - l(v)       # bytte til epsilon her, og legge til uttrykkene for b(u, q), b(v, p)
#     jac((u, p), (du, dp), (v, q)) =  b(v, dp) - b(du, q) + da(u, du, v)

#     op = FEOperator(res, jac, X, Y)

#     # non-linear phase
#     nls = NLSolver(
#     show_trace=true, method=:newton, linesearch=BackTracking())      #prøver å legge inn et max antall iterasjoner og en lav toleranse      
#     solver = FESolver(nls)

#     (uh, ph) = solve(solver, op)

#     errp = p_exact - ph
#     erru = u_exact - uh
    
#     # condition number
#     if calc_condition
#       condition_numb= cond(Array(get_matrix(op)),2)   # kanskje bruke infinitynormen istedenfor
#     else
#       condition_numb = 1
#     end
  
#     if save
#         writevtk(Ω, "C:\\Users\\Sigri\\Documents\\Master\\report\\results\\stokes\\$n $geometry $order.vtu", cellfields=["u_ex" => u_exact, "uh"=>uh, "erru"=> erru, "p_ex" => p_exact, "ph"=>ph, "errp"=> errp, "nablau" => ∇(u_exact)]) #, "erru" => erru]) 
#     end
#     return uh, u_exact, erru, l2_norm(uh - u_exact), h1_semi(uh - u_exact), ph, p_exact, errp, l2_norm(ph - p_exact), h1_semi(ph - p_exact), condition_numb, Ω
# end

# function p_stokes_cutFEM(;n, u_exact, p_exact, f, g, ud, order, geometry, βu0, γu1, γu2, γp, βp0, nu, stabilize, δ, save = false, calc_condition = false)
#        """
#     Unfitted FEM with Nitsche boundary imposition, non-linear stokes (p-stokes).Using P2-P1 Taylor-Hood elements.  
#     n: number of grid elements. Powers of 2 for simplicity and convergence estimates.
#     u_exact: exact solution for method of manufactured solutions
#     order: order of polynomial degree. 
#     f: lhs for first term, -Δ u_ex + ∇p = f
#     g: lhs for second term u = g
#     geometry: optional between "Circle", "Flower", "Heart", "Glacier".
#     stabilize: wheather to add the stabilization term or not
#     δ: perturbation of cut
#     """
#     # Define background mesh
#       partition = (n, n)
#       dim = length(partition)
#       a = 1.2
#       pmin = Point(-a + δ, -a + δ)
#       pmax = Point(a + δ, a + δ)
#       bgmodel = CartesianDiscreteModel(pmin,pmax,partition)
#       # mesh size
#       h = (pmax-pmin)[1]/partition[1]
  
#       # defining ghost penalty constants
#       βu = βu0 *nu/(h^2)
#       βp = βp0/h
  
#       geo = create_geometry(geometry, n)
#       # Define active and physical mesh
#       cutgeo = cut(bgmodel,geo)
#       cutgeo_facets = cut_facets(bgmodel,geo)
#       Ω_bg = Triangulation(bgmodel)
#       Ω_act = Triangulation(cutgeo, ACTIVE)
#       Ω = Triangulation(cutgeo, PHYSICAL)
  
#       # Embedded boundary
#       # Dirichlet conditions on u
#       Γd = EmbeddedBoundary(cutgeo)
#       n_Γd = get_normal_vector(Γd)
  
#       # Get ghost penalty facets
#       Fg = GhostSkeleton(cutgeo)
#       n_Fg = get_normal_vector(Fg)
  
#       # Define measures
#       degree = 2*order
#       dΩ = Measure(Ω,degree)
#       dΓd = Measure(Γd, degree)
#       dFg = Measure(Fg, degree)
  
#       # Define function spaces 
#       reffe_u  = ReferenceFE(lagrangian,VectorValue{dim, Float64},order)
#       reffe_p = ReferenceFE(lagrangian,Float64, order - 1)
  
#       V = TestFESpace(Ω_act, reffe_u,  conformity=:H1)
#       Q = TestFESpace(Ω_act, reffe_p, conformity=:H1, constraint=:zeromean)
  
#       U = TrialFESpace(V)
#       P = TrialFESpace(Q)
  
#       X = MultiFieldFESpace([U, P])
#       Y = MultiFieldFESpace([V, Q])
    
#     println(sum( ∫( p_exact) * dΩ ))

#     l2_norm(u) = (sum( ∫( u ⋅ u )*dΩ ))
#     h1_semi(u) = sum(∫(∇(u) ⊙ ∇(u))*dΩ)
    
#     # fra klassisk Stokes FEM, så er det kun herfra og ned som er endret:)
#     nu0 = 1
#     ϵ_0 = 1e-6
#     γ = 10*2*2
#     # weak formulation components

    
#     a(u, v) = ∫( ∇(v)⊙(flux∘∇(u)))dΩ  + ∫(-((n_Γd ⋅ (flux∘∇(u))) ⋅ v) + (-(n_Γd ⋅ (flux∘∇(v))) ⋅ u)+(γ/h * (u ⋅ v)))dΓd      # denne må ha et ekstra boundary term. Finn ut hvordan det ser ut. 
#     b(v, p) = (∫(-1*(∇ ⋅ v*p))dΩ + ∫((n_Γd ⋅ v) * p)dΓd)   # b er den samme fom før. 
#     l1(v) = ∫(f ⋅ v)dΩ
#     l2(v) = ∫(-(n_Γd ⋅ (flux∘∇(v))) ⋅ ud)dΓd    # har brukt ud som dirichlet grense
#     l3(v) = ∫(γ/h * (ud ⋅ v))dΓd
#     l4(q) = ∫((n_Γd ⋅ ud) * q)dΓd
#     # dflux calculated the same way as in the notebook p-Laplace...
#     dflux(∇du,∇u)=(r-2)*(ϵ_0 + norm(∇u)^2)^((r-4)/2)*(∇u⊙∇du) ⋅ ∇u + (ϵ_0 + norm(∇u)^2)^((r-2)/2)*∇du
#     # and introduced in the same way as in the notebook p-Laplace in the bilinear form a...
#     #da(u, du, v) = ∫(∇(v)⊙(dflux∘(∇(du), ∇(u))))dΩ
   
#     da(u, du, v) = ∫( ∇(v)⊙(dflux∘(∇(du), ∇(u))))dΩ  + ∫(-((n_Γd ⋅ (dflux∘(∇(du), ∇(u)))) ⋅ v) + (-(n_Γd ⋅ (flux∘∇(v))) ⋅ du)+(γ/h * (du ⋅ v)))dΓd       
    
#     # and then the Newton multifield system is assembled as in the Navier Stokes notebook...
#     # når lineærformene differensieres med hensyn på u og p så forsvinner de
#     #res((u,p),(v,q)) = ∫( ∇(v)⊙(flux∘∇(u)))dΩ  + ∫(-((n_Γd ⋅ (flux∘∇(u))) ⋅ v) + (-(n_Γd ⋅ (flux∘∇(v))) ⋅ u)+(γ/h * (u ⋅ v)))dΓd + ∫(-1*(∇ ⋅ v*p))dΩ + ∫((n_Γd ⋅ v) * p)dΓd - ∫(-1*(∇ ⋅ u*q))dΩ - ∫((n_Γd ⋅ u) * q)dΓd - ∫(f ⋅ v)dΩ - ∫(-(n_Γd ⋅( flux∘∇(v))) ⋅ ud)dΓd - ∫(γ/h * (ud ⋅ v))dΓd - ∫((n_Γd ⋅ ud) * q)dΓd#a(u, v) + b(v, p) - b(u, q) -l1(v) -l2(v) -l3(v) -l4(q)
#     #jac((u, p), (du, dp), (v, q)) = ∫(-1*(∇ ⋅ v*dp))dΩ + ∫((n_Γd ⋅ v) * dp)dΓd - ∫(-1*(∇ ⋅ du*q))dΩ - ∫((n_Γd ⋅ du) * q)dΓd + ∫( ∇(v)⊙(dflux∘(∇(du), ∇(u))))dΩ  + ∫(-((n_Γd ⋅ (dflux∘(∇(du), ∇(u))) )⋅ v) + (-(n_Γd ⋅ (flux∘∇(v))) ⋅ du)+(γ/h *(du ⋅ v)))dΓd#b(v, dp) - b(du, q) + da(u, du, v)
    
#     gu(u,v) = ( ∫( (β_1*h)*jump(n_Fg ⋅ ∇(u))⋅jump(n_Fg⋅ ∇(v)) )dFg 
#               +  
#                  ∫( (β_2*h^3)*jump_nn(u,n_Fg)⋅jump_nn(v,n_Fg) )dFg)
  
#     gp(p, q) = (∫((β_3*h^3)*jump(n_Fg ⋅ ∇(p)) * jump(n_Fg ⋅ ∇(q)))dFg)

#     # res((u,p),(v,q)) = a(u, v) + b(v, p) + b(u, q) - l1(v) -l2(v) -l3(v) -l4(q)
#     # jac((u, p), (du, dp), (v, q)) = b(v, dp) + b(du, q) + da(u, du, v)

#     if stabilize
#       res((u,p),(v,q)) = a(u, v) + b(v, p) + b(u, q) + gu(u, v) - gp(p, q) -l1(v) -l2(v) -l3(v) -l4(q)
#       jac((u, p), (du, dp), (v, q)) = b(v, dp) + b(du, q) + da(u, du, v) + gu(du, v) - gp(dp, q) 
      
#       op = FEOperator(res, jac, X, Y)

#       # non-linear phase
#       nls = NLSolver(
#       show_trace=true, method=:newton, linesearch=BackTracking())      #prøver å legge inn et max antall iterasjoner og en lav toleranse      
#       solver = FESolver(nls)

#       (uh, ph) = solve(solver, op)
#     else
#       res_nostab((u,p),(v,q)) = a(u, v) + b(v, p) + b(u, q) - l1(v) -l2(v) -l3(v) -l4(q)
#       jac_nostab((u, p), (du, dp), (v, q)) = b(v, dp) + b(du, q) + da(u, du, v)
      
#       op = FEOperator(res_nostab, jac_nostab, X, Y)

#       # non-linear phase
#       nls = NLSolver(
#       show_trace=true, method=:newton, linesearch=BackTracking())      #prøver å legge inn et max antall iterasjoner og en lav toleranse      
#       solver = FESolver(nls)

#       (uh, ph) = solve(solver, op)
#     end

#     #op = FEOperator(res, jac, X, Y)

#     # non-linear phase
#     #nls = NLSolver(
#     #show_trace=true, method=:newton, linesearch=BackTracking())      #prøver å legge inn et max antall iterasjoner og en lav toleranse      
#     #solver = FESolver(nls)

#     #(uh, ph) = solve(solver, op)

#     errp = p_exact - ph
#     erru = u_exact - uh
    
#     # condition number
#     if calc_condition
#       condition_numb= cond(Array(get_matrix(op)),2)   # kanskje bruke infinitynormen istedenfor
#     else
#       condition_numb = 1
#     end
  
#     if save
#         writevtk(Ω, "C:\\Users\\Sigri\\Documents\\Master\\report\\results\\stokes\\$n $geometry $order.vtu", cellfields=["u_ex" => u_exact, "uh"=>uh, "erru"=> erru, "p_ex" => p_exact, "ph"=>ph, "errp"=> errp, "nablau" => ∇(u_exact)]) #, "erru" => erru]) 
#     end
#     return uh, u_exact, erru, l2_norm(uh - u_exact), h1_semi(uh - u_exact), ph, p_exact, errp, l2_norm(ph - p_exact), h1_semi(ph - p_exact), condition_numb, Ω
# end

stabilize = true
δ = 1/2
save = true
calc_condition = false
order = 2
geometry = "circle"
βu0 = 1
γu1 = 0.1
γu2 = 0.1
γp = 0.1
βp0 = 0.1
nu = 1
β_1 = 1
β_2 = 1
β_3 = 0.1
γ=10*2*2
#n = 128
#uh, u_exact, erru, ul2_norm, uh1_semi, ph, p_exact, errp, pl2_norm, ph1_semi, condition_numb, Ω_act = stokes_FEM(;n, u_exact, p_exact, f, g, ud, order, geometry, βu0, γu1, γu2, γp, βp0, nu, stabilize, δ, save, calc_condition)

################################# p stokes fitted FEM ##############################
# numb_it = 6
# solver = p_stokes_cutFEM
# uarr_l2, uarr_h1, parr_l2, parr_h1, h = convergence_stokes(;numb_it, u_exact, p_exact, f, g, ud, order, geometry, solver, δ, βu0, γu1, γu2, γp, βp0, nu, stabilize, save)


# plot(
#     0,
#     title = "Convergence of p-Stokes FEM",
#     xlabel = "Mesh size h",
#     ylabel = "Velocity error",
#     titlefont = 16,
#     guidefont = 14,
#     tickfont = 12
# )
# plot!(h, uarr_l2, xaxis=:log, yaxis=:log, marker=:o, lw=2, label="L2")
# plot!(h, uarr_h1, marker=:o, lw=2, label="H1")
# xlabel!("Mesh size h")
# ylabel!("Error")
# title!("Convergence of p-stokes cutFEM")

#plot!(h, uarr_l2_1_nostab, marker=:s, lw=2, label="L2 non-stabilized")
#plot!(h, uarr_h1_1_nostab, marker=:s, lw=2, label="H1 non-stabilized")

# # # Legger til aksetitler og tittel

##################### herfra prøver jeg å løse ikke-lineær stokes ######################
# med de samme parametrene som over
uh, u_exact, erru, ul2_norm, uh1_semi, ph, p_exact, errp, pl2_norm, ph1_semi, condition_numb, Ω_act = p_stokes_cutFEM(;n, u_exact, p_exact, f, g, ud, order, geometry, βu0, γu1, γu2, γp, βp0, nu, stabilize, δ, save)


################################# p stokes cut FEM ##############################

#uh, u_exact, erru, ul2_norm, uh1_semi, ph, p_exact, errp, pl2_norm, ph1_semi, condition_numb, Ω_act = p_stokes_cutFEM(;n, u_exact, p_exact, f, g, ud, order, geometry, βu0, γu1, γu2, γp, βp0, nu, stabilize, δ, save)

# now, doing a geometry robustness test

# ####### sensitivity_stokes test #######
# # kjører nå denne med n = 16, men bør nok kjøre for n = 32 eller n = 64, men det kan ta laaaag tid. 30 min per kjøring ved n = 64 -
n = 16           # øke denne
M = 100         #full kjøring med M = 2000 med 2000 så kjører det nok i 2 timer. 
order = 2
geometry = "circle"
solver = p_stokes_cutFEM
βu0 = 1
γu1 = 1
γu2 = 1
γp = 0.1
βp0 = 0.1
stabilize = true
save = false
calc_condition = false

arr_δ, arr_l2u, arr_h1u, arr_l2p, arr_h1p, arr_cond = sensitivity_stokes(;n, M, u_exact, p_exact, f, g, ud, order, geometry, solver, βu0, γu1, γu2, γp, βp0, nu, stabilize, save)
stabilize = false
start = 1
arr_δ_nostab, arr_l2u_nostab, arr_h1u_nostab, arr_l2p_nostab, arr_h1p_nostab, arr_cond_nostab = sensitivity_stokes(;n, M, u_exact, p_exact, f, g, ud, order, geometry, solver, βu0, γu1, γu2, γp, βp0, nu, stabilize, save)

plot(
    0,
    title = "Sensitivity of Stokes Solver",
    xlabel = "Perturbation δ",
    ylabel = "Velocity error",
    titlefont = 16,
    guidefont = 14,
    tickfont = 12
)
#Indeksene der vi vil ha markører (hver 100.)
idx = 1:100:1999
id2 = 51:100:1999
scatter!(arr_δ[idx], arr_l2p[idx], label=:"", marker=:circle, ms=4)
plot!(arr_δ[start:end], arr_l2p[start:end],  yaxis=:log, lw=2, label="L2 stabilized")
scatter!(arr_δ[id2], arr_l2p_nostab[id2], label=:"", marker=:s, ms=4)
plot!(arr_δ_nostab[start:end], arr_l2p_nostab[start:end], yaxis=:log, lw=2, label="L2 non-stabilized")
scatter!(arr_δ[idx], arr_h1p[idx], label=:"", marker=:circle, ms=4)
plot!(arr_δ[start:end], arr_h1p[start:end], yaxis=:log, lw=2, label="H1 stabilized")
scatter!(arr_δ[id2], arr_h1p_nostab[id2], label=:"", marker=:s, ms=4)
plot!(arr_δ_nostab[start:end], arr_h1p_nostab[start:end], yaxis=:log, lw=2, label="H1 non-stabilized")
xlabel!("Perturbation δ")
ylabel!("Pressure error")
title!("Sensitivity analysis of p-Stokes cutFEM")

#condition number plot:
plot(
    0,
    title = "Sensitivity of Stokes Solver",
    xlabel = "Perturbation δ",
    ylabel = "Condition number",
    titlefont = 16,
    guidefont = 14,
    tickfont = 12
)
#scatter!(arr_δ[idx], arr_cond[idx], label=:"", marker=:circle, ms=4)
#plot!(arr_δ[start:end], arr_cond[start:end], yaxis=:log, label = "Stabilized")
#scatter!(arr_δ[id2], arr_cond_nostab[id2], label=:"", marker=:s, ms=4)
#plot!(arr_δ[start:end], arr_cond_nostab[start:end],yaxis=:log, label = "Not stabilized")

