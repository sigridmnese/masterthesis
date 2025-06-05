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
include("C:\\Users\\Sigri\\Documents\\Master\\report\\code\\Julia_programmering\\Starting_up\\workfiles_stokes\\testing_non-linear_stokes.jl")

##### Divergensfri stokes løser ####
# hva hvis det ikke er divergensfritt?
nu0 =  1   # klarer ikke mindre enn 0.01 og klarer heller ikke større enn 100 (størrelsesordener)
r = 4/3
A = 1
ϵ_0 = 1e-6          # hadde glemt å legge inn en epsilon0 i flux-termen. Newton-løseren divergerte
# endrer nå alle symbolder før uder det er nabla til epsilon...

# non-linear stokes problem denne løsninen fungerer klink: VectorValue(x[2]^2, -x[1])# bør da for søren være mulig å kjøre denne likningen på en fitted mesh også.???
u_exact(x) =  VectorValue(2*x[1] + cos(2*π*x[2]), -2*x[2] + sin(2*π*x[1]))#VectorValue(2*x[1] + exp(x[1]/2) * cos(2*π*x[2]), -2*x[2] + exp(x[2]/2) * sin(2*π*x[1])) #bytte til sin/cos-uttrykk  VectorValue(-x[2], x[1])
p_exact(x) = sin(2*π*x[1])*cos(2*π*x[2])
flux(∇u) = nu0*(ϵ_0 + norm(∇u)^2)^((r-2)/2) * ∇u 
viskositet(∇u) = nu0^(1-r)*(ϵ_0^2 + norm(∇u)^2)^((r-2)/2)
dviskositet(∇u, ∇du) = nu0^(1-r) *(ϵ_0^2 + ∇u⋅ ∇u)^((r-4)/2)*(∇u⊙∇du) *(r-2)
f(x) =  -divergence(flux∘ε(u_exact))(x) + ∇(p_exact)(x)      # prøver å endre f her...
ud(x) = u_exact(x)
#1/2 * A^(1-r) * (1/2 * norm(du, 2))^((r-2)/2) * ε(du)

g = VectorValue(0.0, 0.0)
n = 32
stabilize = false
δ = 0 #(2000-1)/2000 *1.2/n            # perturbation of the cut. One element is 2 * 1.2/n
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
nu = 1      # denne brukes ikke i p_stokes_cutfem... Alt er definert i de funksjonene med nu0? sjekk opp dette...
uh, u_exact, erru, ul2_norm, uh1_semi, ph, p_exact, errp, pl2_norm, ph1_semi, condition_numb, Ω_act = stokes_cutFEM(;n, u_exact, p_exact, f, g, ud, order, geometry, βu0, γu1, γu2, γp, βp0, nu, stabilize, δ, save, calc_condition)

################################# p stokes fitted FEM ##############################
# numb_it = 6
# solver = p_stokes_FEM
# uarr_l2, uarr_h1, parr_l2, parr_h1, h = convergence_stokes(;numb_it, u_exact, p_exact, f, g, ud, order, geometry, solver, δ, βu0, γu1, γu2, γp, βp0, nu, stabilize, save)

# stabilize = false
# # uarr_l2_1_nostab, uarr_h1_1_nostab, parr_l2_1_nostab, parr_h1_1_nostab, h = convergence_stokes(;numb_it, u_exact, p_exact, f, g, ud, order, geometry, solver, δ, βu0, γu1, γu2, γp, βp0, nu, stabilize, save)
# plot(
#     0,
#     title = "Convergence of p-Stokes FEM",
#     xlabel = "Mesh size h",
#     ylabel = "Velocity error",
#     titlefont = 16,
#     guidefont = 14,
#     tickfont = 12
# )
# plot!(h, uarr_l2, xaxis=:log, yaxis=:log, marker=:o, lw=2, label="L2 stabilized")
# plot!(h,uarr_h1, marker=:o, lw=2, label="H1 stabilized")
# xlabel!("Mesh size h")
# ylabel!("Velocity error")
# title!("Convergence of p-stokes FEM")

#plot!(h, uarr_l2_1_nostab, marker=:s, lw=2, label="L2 non-stabilized")
#plot!(h, uarr_h1_1_nostab, marker=:s, lw=2, label="H1 non-stabilized")

# # # Legger til aksetitler og tittel

##################### herfra prøver jeg å løse ikke-lineær stokes ######################
# med de samme parametrene som over
#uh, u_exact, erru, ul2_norm, uh1_semi, ph, p_exact, errp, pl2_norm, ph1_semi, condition_numb, Ω_act = p_stokes_cutFEM_symmetric(;n, u_exact, p_exact, f, g, ud, order, geometry, βu0, γu1, γu2, γp, βp0, nu, stabilize, δ, save)

################################# p stokes cut FEM ##############################
#uh, u_exact, erru, ul2_norm, uh1_semi, ph, p_exact, errp, pl2_norm, ph1_semi, condition_numb, Ω_act = p_stokes_cutFEM(;n, u_exact, p_exact, f, g, ud, order, geometry, βu0, γu1, γu2, γp, βp0, nu, stabilize, δ, save)

# now, doing a geometry robustness test

# ####### sensitivity_stokes test #######
# # kjører nå denne med n = 16, men bør nok kjøre for n = 32 eller n = 64, men det kan ta laaaag tid. 30 min per kjøring ved n = 64 -
# n = 16           # øke denne
# M = 100         #full kjøring med M = 2000 med 2000 så kjører det nok i 2 timer. 
# order = 2
# geometry = "circle"
# solver = p_stokes_cutFEM
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
# arr_δ_nostab, arr_l2u_nostab, arr_h1u_nostab, arr_l2p_nostab, arr_h1p_nostab, arr_cond_nostab = sensitivity_stokes(;n, M, u_exact, p_exact, f, g, ud, order, geometry, solver, βu0, γu1, γu2, γp, βp0, nu, stabilize, save)

# plot(
#     0,
#     title = "Linear Basis Functions ϕ₃, ϕ₄, ϕ₅",
#     xlabel = "x",
#     ylabel = "ϕᵢ(x)",
#     titlefont = 16,
#     guidefont = 14,
#     tickfont = 12
# )
# #Indeksene der vi vil ha markører (hver 100.)
# idx = 1:100:1999
# id2 = 51:100:1999
# start = 1
# #scatter!(arr_δ[idx], arr_l2p[idx], label=:"", marker=:circle, ms=4)
# plot!(arr_δ[start:end], arr_l2u[start:end],  yaxis=:log, lw=2, label="L2 stabilized")
# #scatter!(arr_δ[id2], arr_l2p_nostab[id2], label=:"", marker=:s, ms=4)
# plot!(arr_δ_nostab[start:end], arr_l2u_nostab[start:end], yaxis=:log, lw=2, label="L2 non-stabilized")
# #scatter!(arr_δ[idx], arr_h1p[idx], label=:"", marker=:circle, ms=4)
# plot!(arr_δ[start:end], arr_h1u[start:end], yaxis=:log, lw=2, label="H1 stabilized")
# #scatter!(arr_δ[id2], arr_h1p_nostab[id2], label=:"", marker=:s, ms=4)
# plot!(arr_δ_nostab[start:end], arr_h1u_nostab[start:end], yaxis=:log, lw=2, label="H1 non-stabilized")
# xlabel!("Perturbation δ")
# ylabel!("Velocity error")
# title!("Sensitivity analysis of p-Stokes cutFEM")

# #condition number plot:
# plot(
#     0,
#     title = "Sensitivity of Stokes Solver",
#     xlabel = "Perturbation δ",
#     ylabel = "Condition number",
#     titlefont = 16,
#     guidefont = 14,
#     tickfont = 12
# )
# #scatter!(arr_δ[idx], arr_cond[idx], label=:"", marker=:circle, ms=4)
# plot!(arr_δ[start:end], arr_cond[start:end], yaxis=:log, label = "Stabilized")
# scatter!(arr_δ[id2], arr_cond_nostab[id2], label=:"", marker=:s, ms=4)
# plot!(arr_δ[start:end], arr_cond_nostab[start:end],yaxis=:log, label = "Not stabilized")