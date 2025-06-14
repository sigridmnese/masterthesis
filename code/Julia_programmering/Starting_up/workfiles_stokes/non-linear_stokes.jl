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

# Defining constants
nu0 =  0.01   # klarer ikke mindre enn 0.01 og klarer heller ikke større enn 100 (størrelsesordener)
r = 4/3
A = nu0
ϵ_0 = 1e-6       

# Defining manufactured solutions
u_exact(x) =  VectorValue(2*x[1] + cos(2*π*x[2]), -2*x[2] + sin(2*π*x[1]))#VectorValue(2*x[1] + exp(x[1]/2) * cos(2*π*x[2]), -2*x[2] + exp(x[2]/2) * sin(2*π*x[1])) #bytte til sin/cos-uttrykk  VectorValue(-x[2], x[1])
p_exact(x) = sin(2*π*x[1])*cos(2*π*x[2])

# Defining the problem 
flux(εu) = nu0*(ϵ_0 + norm(εu)^2)^((r-2)/2)⋅εu 
dflux(εdu,εu)=(r-2)*nu0*(ϵ_0 + norm(εu)^2)^((r-4)/2)*(εu⊙εdu) ⋅ εu + nu0*(ϵ_0 + norm(εu)^2)^((r-2)/2)*εdu
# Det var en feil i dflux!!!
viskositet(εu) = nu0^(1-r)*(ϵ_0^2 + norm(εu)^2)^((r-2)/2)
dviskositet(εu, εdu) = nu0^(1-r) *(ϵ_0^2 + εu ⋅ εu)^((r-4)/2)*(εu⊙εdu) *(r-2)

f(x) =  -divergence(flux∘ε(u_exact))(x) + ∇(p_exact)(x)      # prøver å endre f her...
ud(x) = u_exact(x)
g = VectorValue(0.0, 0.0)

# solver parameters
n = 32
stabilize = true
solver = p_stokes_cutFEM_symmetric
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

# stabilize = false
# uarr_l2_1_nostab, uarr_h1_1_nostab, parr_l2_1_nostab, parr_h1_1_nostab, h = convergence_stokes(;numb_it, u_exact, p_exact, f, g, ud, order, geometry, solver, δ, βu0, γu1, γu2, γp, βp0, nu, stabilize, save)

# ########## velocity convergence plot:
# plot(
#     0,
#     titlefont = 16,
#     guidefont = 14,
#     tickfont = 12
# )

# plot!(h, uarr_l2, xaxis=:log, yaxis=:log, marker=:o, lw=2, label="L2 stabilized")
# plot!(h, uarr_h1, marker=:o, lw=2, label="H1 stabilized")
# plot!(h, uarr_l2_1_nostab, marker=:s, lw=2, label="L2 non-stabilized")
# plot!(h, uarr_h1_1_nostab, marker=:s, lw=2, label="H1 non-stabilized")
# xlabel!("Mesh size h")
# ylabel!("Velocity error")
# title!("Convergence of p-stokes FEM")

# # ########### pressure convergence plot:
# plot(
#     0,
#     titlefont = 16,
#     guidefont = 14,
#     tickfont = 12
# )

# plot!(h, parr_l2, xaxis=:log, yaxis=:log, marker=:o, lw=2, label="L2 stabilized")
# plot!(h, parr_h1, marker=:o, lw=2, label="H1 stabilized")
# plot!(h, parr_l2_1_nostab, marker=:s, lw=2, label="L2 non-stabilized")
# plot!(h, parr_h1_1_nostab, marker=:s, lw=2, label="H1 non-stabilized")
# xlabel!("Mesh size h")
# ylabel!("Pressure error")
# title!("Convergence of p-stokes FEM")

# ###################################  sensitivity_stokes test ########################################## 
# n = 16                # øke denne?
# M = 4               #full kjøring med M = 2000 med 2000 så kjører det nok i 2 timer. 
# order = 2
# geometry = "heart"
# βu0 = 1
# γu1 = 1
# γu2 = 1
# γp = 0.1
# βp0 = 0.1
# stabilize = true
# save = true
# calc_condition = false

# arr_δ, arr_l2u, arr_h1u, arr_l2p, arr_h1p, arr_cond = sensitivity_stokes(;n, M, u_exact, p_exact, f, g, ud, order, geometry, solver, βu0, γu1, γu2, γp, βp0, nu, stabilize, save)
#stabilize = false
#start = 1
#save = false
#arr_δ_nostab, arr_l2u_nostab, arr_h1u_nostab, arr_l2p_nostab, arr_h1p_nostab, arr_cond_nostab = sensitivity_stokes(;n, M, u_exact, p_exact, f, g, ud, order, geometry, solver, βu0, γu1, γu2, γp, βp0, nu, stabilize, save)

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
# title!("Sensitivity analysis of p-Stokes cutFEM")

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
# title!("Sensitivity analysis of p-Stokes cutFEM")