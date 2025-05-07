using Gridap
using GridapEmbedded
using STLCutters
using LinearAlgebra
using Plots
include("C:\\Users\\Sigri\\Documents\\Master\\report\\code\\Julia_programmering\\Starting_up\\utils\\utils.jl")

##### Divergensfri stokes løser ####
# hva hvis det ikke er divergensfritt?
nu = 1
u_exact(x) = VectorValue(2*x[1] + cos(2*π*x[2]), -2*x[2] + sin(2*π*x[1]))
p_exact(x) = sin(2*π*x[1])
# Forcing term
#f(x)= -Δ(u_ex)(x)+ ∇(p_ex)(x)
f(x) = -divergence(∇(u_exact))(x) + ∇(p_exact)(x)
ud(x) = u_exact(x)
domain = "circle"
n = 128
γ = 10* 2*2 
β_1 = 1
β_2 = 1
β_3 = 0.1
γ=10*2*2

outputfolder = "C:\\Users\\Sigri\\Documents\\Master\\report\\results\\"
################ Løsning ################
# order = 2

n = 64
order = 2
geometry = "circle"
βu0 = 1
γu1 = 1
γu2 = 1.0
γp = 0.1
βp0 = 0.1
stabilize = true
δ = 0
save = true
g=1
calc_condition = false
#uh, u_exact, erru, l2_u, h1_semi_u, ph, p_exact, errp, l2_p, h1_semi_p, condition_numb, Ω_act  = stokes_solver(;n, u_exact, p_exact, f, g, ud, order, geometry, βu0, γu1, γu2, γp, βp0, nu, stabilize, δ, save, calc_condition)


# ####### convergence_stokes test #######
# #  stabilisering, order 1           # trenger ikke regne ut kondisjonstall når man gjør konvergestest...
numb_it = 6                         # Sånn koden er implementert nå så er det fra 2^1 til 2^{numb_it}. Tar kort tid å kjøre for numb_it = 6
order = 2                           # når jeg øker orden så øker kjøretid veeeeldig !! Bør vurdere å skru ned numb_it samtidig. 244 sekunder når jeg har på order = 2 for kjøringen 2^6
δ = 0                               # kan også se ut til at feilen havner på maskinnivå? vet ikke helt, men mulig å eksperimentere med dette. 
# med stabilisering, order 1
stabilize = true
solver = stokes_solver
geometry = "circle"
γu1 = 0.1
γu2 = 0.1
γp = 0.1
save = false
# fjernet disse: , eoc_l2, eoc_h1
uarr_l2_1_stab, uarr_h1_1_stab, parr_l2_1_stab, parr_h1_1_stab, h = convergence_stokes(;numb_it, u_exact, p_exact, f, g, ud, order, geometry, solver, δ, βu0, γu1, γu2, γp, βp0, nu, stabilize, save)
start = 1
plot(h[start:end], uarr_l2_1_stab[start:end], xaxis=:log, yaxis=:log, marker=:o, lw=2, label="L2 error")
plot!(h[start:end], uarr_h1_1_stab[start:end], marker=:s, lw=2, label="H1 error")

# # # # Uten stabilisering, order 1
# # stabilize = false
# # uarr_l2_1_nostab, uarr_h1_1_nostab, parr_l2_1_nostab, parr_h1_1_nostab, h = convergence_stokes(;numb_it, u_exact, p_exact, f, g, ud, order, geometry, solver, δ, βu0, γu1, γu2, γp, βp0, nu, stabilize, save)

# # plot!(h[start:end], uarr_l2_1_nostab[start:end], marker=:s, lw=2, label="L2 error, order $order (no Stab)")
# # plot!(h[start:end], uarr_h1_1_nostab[start:end], marker=:s, lw=2, label="H1 error, order $order (no Stab)")

# # # # Legger til aksetitler og tittel
# # xlabel!("Mesh size h")
# # ylabel!("Error")
# # title!("convergence_stokes of Stokes Solver")

# ####### sensitivity_stokes test #######
# # kjører nå denne med n = 16, men bør nok kjøre for n = 32 eller n = 64, men det kan ta laaaag tid. 30 min per kjøring ved n = 64 -
n = 16           # øke denne
M = 1000        #full kjøring med M = 2000
order = 2
geometry = "circle"
solver = stokes_solver
βu0 = 1
γu1 = 1
γu2 = 1
γp = 0.1
βp0 = 0.1
stabilize = true
save = false
calc_condition = true

arr_δ, arr_l2u, arr_h1u, arr_l2p, arr_h1p, arr_cond = sensitivity_stokes(;n, M, u_exact, p_exact, f, g, ud, order, geometry, solver, βu0, γu1, γu2, γp, βp0, nu, stabilize, save)
stabilize = false
start = 1
arr_δ_nostab, arr_l2u_nostab, arr_h1u_nostab, arr_l2p_nostab, arr_h1p_nostab, arr_cond_nostab = sensitivity_stokes(;n, M, u_exact, p_exact, f, g, ud, order, geometry, solver, βu0, γu1, γu2, γp, βp0, nu, stabilize, save)
plot(arr_δ[start:end], arr_l2u[start:end],  yaxis=:log, lw=2, label="L2 error, stabilized")
plot!(arr_δ_nostab[start:end], arr_l2u_nostab[start:end], yaxis=:log, lw=2, label="L2 error, non-stabilized")
plot!(arr_δ[start:end], arr_h1u[start:end], yaxis=:log, lw=2, label="h1 stabilized")
plot!(arr_δ_nostab[start:end], arr_h1u_nostab[start:end], yaxis=:log, lw=2, label="h1 non-stabilized")
xlabel!("Perturbation δ")
ylabel!("Condition number")
title!("sensitivity_stokes analysis of stokes solver, solution u")

plot!(arr_δ[start:end], arr_cond[start:end], yaxis=:log, label = "Stabilized")
plot!(arr_δ[start:end], arr_cond_nostab[start:end],yaxis=:log, label = "Not stabilized")

# # if I want the same plot for p
# plot(arr_δ[start:end], arr_l2p[start:end], xaxis=:log, yaxis=:log, lw=2, label="L2 error, stabilized")
# plot!(arr_δ_nostab[start:end], arr_l2p_nostab[start:end], xaxis=:log, yaxis=:log, lw=2, label="L2 error, non-stabilized")
# plot!(arr_δ[start:end], arr_h1p[start:end], xaxis=:log, yaxis=:log, lw=2, label="h1 stabilized")
# plot!(arr_δ_nostab[start:end], arr_h1p_nostab[start:end], xaxis=:log, yaxis=:log, lw=2, label="h1 non-stabilized")
# xlabel!("Perturbation δ")
# ylabel!("Condition number")
# title!("sensitivity_stokes analysis of stokes solver, solution p")

# # #### Varying geometry
# geometry_arr = ["circle", "flower", "heart"]
# s1 = plot()  # Lager et tomt plott
# s2 = plot()
# s3 = plot()
# start = 1
# for i = 1:3
#     geometry = geometry_arr[i]
#     stabilize = true
#     arr_δ, arr_l2u, arr_h1u, arr_l2p, arr_h1p, arr_cond = sensitivity_stokes(;n, M, u_exact, p_exact, f, g, ud, order, geometry, solver, βu0, γu1, γu2, γp, βp0, nu, stabilize, save)
#     stabilize = false
#     arr_δ, arr_l2u_nostab, arr_h1u_nostab, arr_l2p_nostab, arr_h1p_nostab, arr_cond_nostab = sensitivity_stokes(;n, M, u_exact, p_exact, f, g, ud, order, geometry, solver, βu0, γu1, γu2, γp, βp0, nu, stabilize, save)
#     plot!(s1, arr_δ[start:end], arr_cond[start:end], xaxis=:log, yaxis=:log, lw=2, label="Stabilized ($geometry)")
#     plot!(s1, arr_δ[start:end], arr_cond_nostab[start:end], xaxis=:log, yaxis=:log, lw=2, label="Non-stabilized ($geometry)")
#     plot!(s2, arr_δ[start:end], arr_l2u[start:end], xaxis=:log, yaxis=:log, lw=2, label="Stabilized ($geometry)")
#     plot!(s2, arr_δ[start:end], arr_l2u_nostab[start:end], xaxis=:log, yaxis=:log, lw=2, label="Non-stabilized ($geometry)")
#     plot!(s3, arr_δ[start:end], arr_h1u[start:end], xaxis=:log, yaxis=:log, lw=2, label="Stabilized ($geometry)")
#     plot!(s3, arr_δ[start:end], arr_h1u_nostab[start:end], xaxis=:log, yaxis=:log, lw=2, label="Non-stabilized ($geometry)")
# end
# title!(s1, "sensitivity_stokes analysis varing geometries")
# xlabel!(s1, "Perturbation δ")
# ylabel!(s1, "Condition number")
# title!(s2, "sensitivity_stokes analysis varing geometries")
# xlabel!(s2, "Perturbation δ")
# ylabel!(s2, "L2 norm")
# title!(s3, "sensitivity_stokes analysis varing geometries")
# xlabel!(s3, "Perturbation δ")
# ylabel!(s3, "H1 semi norm")

# display(s1)
# display(s2)
# display(s3)

# start = 1
# ### Sensitivitetsanalyse: Varying γu1 ###
# γu1_arr = [0.1, 1, 10]
# p1 = plot()
# p2 = plot()
# p3 = plot()
# for i = 1:3
#     γu1 = γu1_arr[i]
#     stabilize = true
#     arr_δ, arr_l2u, arr_h1u, arr_l2p, arr_h1p, arr_cond = sensitivity_stokes(;n, M, u_exact, p_exact, f, g, ud, order, geometry, solver, βu0, γu1, γu2, γp, βp0, nu, stabilize, save)
#     stabilize = false
#     arr_δ_nostab, arr_l2u_nostab, arr_h1u_nostab, arr_l2p_nostab, arr_h1p_nostab, arr_cond_nostab = sensitivity_stokes(;n, M, u_exact, p_exact, f, g, ud, order, geometry, solver, βu0, γu1, γu2, γp, βp0, nu, stabilize, save)
#     plot!(p1, arr_δ[start:end], arr_cond[start:end], xaxis=:log, yaxis=:log, lw=2, label="Stabilized γd = $γd")
#     plot!(p1, arr_δ[start:end], arr_cond_nostab[start:end], xaxis=:log, yaxis=:log, lw=2, label="Non-stabilized γd = $γd")
#     plot!(p2, arr_δ[start:end], arr_l2u[start:end], xaxis=:log, yaxis=:log, lw=2, label="Stabilized γd = $γd")
#     plot!(p2, arr_δ[start:end], arr_l2u_nostab[start:end], xaxis=:log, yaxis=:log, lw=2, label="Non-stabilized γd = $γd")
#     plot!(p3, arr_δ[start:end], arr_h1u[start:end], xaxis=:log, yaxis=:log, lw=2, label="Stabilized γd = $γd")
#     plot!(p3, arr_δ[start:end], arr_h1u_nostab[start:end], xaxis=:log, yaxis=:log, lw=2, label="Non-stabilized γd = $γd")
# end
# title!(p1, "sensitivity_stokes analysis varying γu1")
# xlabel!(p1, "Perturbation δ")
# ylabel!(p1, "Condition number")
# title!(p2, "sensitivity_stokes analysis varying γd")
# xlabel!(p2, "Perturbation δ")
# ylabel!(p2, "L2 norm (velocity)")
# title!(p3, "sensitivity_stokes analysis varying γd")
# xlabel!(p3, "Perturbation δ")
# ylabel!(p3, "H1 norm (velocity)")
# display(p1)
# display(p2)
# display(p3)

# # ### Sensitivitetsanalyse: Varying γg1 ###
# γg1_arr = [0.001, 0.01, 0.1, 1]
# γd = 0.1
# q1 = plot()
# q2 = plot()
# q3 = plot()
# for i = 1:4
#     γg1 = γg1_arr[i]
#     stabilize = true                   
#     arr_δ, arr_l2u, arr_h1u, arr_l2p, arr_h1p, arr_cond = sensitivity_stokes(;n, M, u_exact, p_exact, f, g, ud, order, geometry, solver, βu0, γu1, γu2, γp, βp0, nu, stabilize, save)
#     stabilize = false
#     arr_δ_nostab, arr_l2u_nostab, arr_h1u_nostab, arr_l2p_nostab, arr_h1p_nostab, arr_cond_nostab = sensitivity_stokes(;n, M, u_exact, p_exact, f, g, ud, order, geometry, solver, βu0, γu1, γu2, γp, βp0, nu, stabilize, save)
#     plot!(q1, arr_δ[start:end], arr_cond[start:end], xaxis=:log, yaxis=:log, lw=2, label="Stabilized γg1 = $γg1")
#     plot!(q1, arr_δ[start:end], arr_cond_nostab[start:end], xaxis=:log, yaxis=:log, lw=2, label="Non-stabilized γg1 = $γg1")
#     plot!(q2, arr_δ[start:end], arr_l2u[start:end], xaxis=:log, yaxis=:log, lw=2, label="Stabilized γg1 = $γg1")
#     plot!(q2, arr_δ[start:end], arr_l2u_nostab[start:end], xaxis=:log, yaxis=:log, lw=2, label="Non-stabilized γg1 = $γg1")
#     plot!(q3, arr_δ[start:end], arr_h1u[start:end], xaxis=:log, yaxis=:log, lw=2, label="Stabilized γg1 = $γg1")
#     plot!(q3, arr_δ[start:end], arr_h1u_nostab[start:end], xaxis=:log, yaxis=:log, lw=2, label="Non-stabilized γg1 = $γg1")
# end
# title!(q1, "sensitivity_stokes analysis varying γg1")
# xlabel!(q1, "Perturbation δ")
# ylabel!(q1, "Condition number")
# title!(q2, "sensitivity_stokes analysis varying γg1")
# xlabel!(q2, "Perturbation δ")
# ylabel!(q2, "L2 norm (velocity)")
# title!(q3, "sensitivity_stokes analysis varying γg1")
# xlabel!(q3, "Perturbation δ")
# ylabel!(q3, "H1 norm (velocity)")
# display(q1)
# display(q2)
# display(q3)

# ### Sensitivitetsanalyse: Varying γg3 ###
# γg3_arr = [0.1, 1, 10]
# γd = 0.1
# r1 = plot()
# r2 = plot()
# r3 = plot()
# for i = 1:3
#     γg3 = γg3_arr[i]
#     stabilize = true
#     arr_δ, arr_l2u, arr_h1u, arr_l2p, arr_h1p, arr_cond = sensitivity_stokes(;n, M, u_exact, p_exact, f, g, ud, order, geometry, solver, βu0, γu1, γu2, γp, βp0, nu, stabilize, save)
#     stabilize = false
#     arr_δ_nostab, arr_l2u_nostab, arr_h1u_nostab, arr_l2p_nostab, arr_h1p_nostab, arr_cond_nostab = sensitivity_stokes(;n, M, u_exact, p_exact, f, g, ud, order, geometry, solver, βu0, γu1, γu2, γp, βp0, nu, stabilize, save)
#     plot!(r1, arr_δ[start:end], arr_cond[start:end], xaxis=:log, yaxis=:log, lw=2, label="Stabilized γg3 = $γg3")
#     plot!(r1, arr_δ_nostab[start:end], arr_cond_nostab[start:end], xaxis=:log, yaxis=:log, lw=2, label="Non-stabilized γg3 = $γg3")
#     plot!(r2, arr_δ[start:end], arr_l2u[start:end], xaxis=:log, yaxis=:log, lw=2, label="Stabilized γg3 = $γg3")
#     plot!(r2, arr_δ_nostab[start:end], arr_l2u_nostab[start:end], xaxis=:log, yaxis=:log, lw=2, label="Non-stabilized γg3 = $γg3")
#     plot!(r3, arr_δ[start:end], arr_h1u[start:end], xaxis=:log, yaxis=:log, lw=2, label="Stabilized γg3 = $γg3")
#     plot!(r3, arr_δ_nostab[start:end], arr_h1u_nostab[start:end], xaxis=:log, yaxis=:log, lw=2, label="Non-stabilized γg3 = $γg3")
# end
# title!(r1, "sensitivity_stokes analysis varying γg3")
# xlabel!(r1, "Perturbation δ")
# ylabel!(r1, "Condition number")
# title!(r2, "sensitivity_stokes analysis varying γg3")
# xlabel!(r2, "Perturbation δ")
# ylabel!(r2, "L2 norm (velocity)")
# title!(r3, "sensitivity_stokes analysis varying γg3")
# xlabel!(r3, "Perturbation δ")
# ylabel!(r3, "H1 norm (velocity)")
# display(r1)
# display(r2)
# display(r3)


# start = 1

# ###########################
# ### sensitivity_stokes: βu
# ###########################
# βu_arr = [0.01, 0.1, 1.0]
# p1 = plot(); p2 = plot(); p3 = plot()
# for βu0 in βu_arr
#     stabilize = true
#     arr_δ, arr_l2u, arr_h1u, _, _, arr_cond = sensitivity_stokes(;n, M, u_exact, p_exact, f, g, ud, order, geometry, solver, βu0, γu1, γu2, γp, βp0, nu, stabilize, save)
#     stabilize = false
#     _, arr_l2u_ns, arr_h1u_ns, _, _, arr_cond_ns = sensitivity_stokes(;n, M, u_exact, p_exact, f, g, ud, order, geometry, solver, βu0, γu1, γu2, γp, βp0, nu, stabilize, save)

#     plot!(p1, arr_δ[start:end], arr_cond[start:end], xaxis=:log, yaxis=:log, label="Stabilized βu=$βu")
#     plot!(p1, arr_δ[start:end], arr_cond_ns[start:end], xaxis=:log, yaxis=:log, label="Non-stabilized βu=$βu")

#     plot!(p2, arr_δ[start:end], arr_l2u[start:end], xaxis=:log, yaxis=:log, label="Stabilized βu=$βu")
#     plot!(p2, arr_δ[start:end], arr_l2u_ns[start:end], xaxis=:log, yaxis=:log, label="Non-stabilized βu=$βu")

#     plot!(p3, arr_δ[start:end], arr_h1u[start:end], xaxis=:log, yaxis=:log, label="Stabilized βu=$βu")
#     plot!(p3, arr_δ[start:end], arr_h1u_ns[start:end], xaxis=:log, yaxis=:log, label="Non-stabilized βu=$βu")
# end
# title!(p1, "sensitivity_stokes: varying βu"); xlabel!(p1, "δ"); ylabel!(p1, "Condition number")
# title!(p2, "L2 error (velocity), varying βu"); xlabel!(p2, "δ"); ylabel!(p2, "L2 norm")
# title!(p3, "H1 error (velocity), varying βu"); xlabel!(p3, "δ"); ylabel!(p3, "H1 norm")
# display(p1); 
# display(p2); 
# display(p3)


# ###########################
# ### sensitivity_stokes: γu1
# ###########################
# γu1_arr = [0.1, 1.0, 10.0]
# q1 = plot(); q2 = plot(); q3 = plot()
# for γu1 in γu1_arr
#     stabilize = true
#     arr_δ, arr_l2u, arr_h1u, _, _, arr_cond = sensitivity_stokes(;n, M, u_exact, p_exact, f, g, ud, order, geometry, solver, βu0, γu1, γu2, γp, βp0, nu, stabilize, save)
#     stabilize = false
#     _, arr_l2u_ns, arr_h1u_ns, _, _, arr_cond_ns = sensitivity_stokes(;n, M, u_exact, p_exact, f, g, ud, order, geometry, solver, βu0, γu1, γu2, γp, βp0, nu, stabilize, save)

#     plot!(q1, arr_δ[start:end], arr_cond[start:end], xaxis=:log, yaxis=:log, label="Stabilized γu1=$γu1")
#     plot!(q1, arr_δ[start:end], arr_cond_ns[start:end], xaxis=:log, yaxis=:log, label="Non-stabilized γu1=$γu1")

#     plot!(q2, arr_δ[start:end], arr_l2u[start:end], xaxis=:log, yaxis=:log, label="Stabilized γu1=$γu1")
#     plot!(q2, arr_δ[start:end], arr_l2u_ns[start:end], xaxis=:log, yaxis=:log, label="Non-stabilized γu1=$γu1")

#     plot!(q3, arr_δ[start:end], arr_h1u[start:end], xaxis=:log, yaxis=:log, label="Stabilized γu1=$γu1")
#     plot!(q3, arr_δ[start:end], arr_h1u_ns[start:end], xaxis=:log, yaxis=:log, label="Non-stabilized γu1=$γu1")
# end
# title!(q1, "sensitivity_stokes: varying γu1"); xlabel!(q1, "δ"); ylabel!(q1, "Condition number")
# title!(q2, "L2 error (velocity), varying γu1"); xlabel!(q2, "δ"); ylabel!(q2, "L2 norm")
# title!(q3, "H1 error (velocity), varying γu1"); xlabel!(q3, "δ"); ylabel!(q3, "H1 norm")
# display(q1); 
# display(q2); 
# display(q3)


# ###########################
# ### sensitivity_stokes: γu2
# ###########################
# γu2_arr = [0.01, 0.1, 1.0]
# r1 = plot(); r2 = plot(); r3 = plot()
# for γu2 in γu2_arr
#     stabilize = true
#     arr_δ, arr_l2u, arr_h1u, _, _, arr_cond = sensitivity_stokes(;n, M, u_exact, p_exact, f, g, ud, order, geometry, solver, βu0, γu1, γu2, γp, βp0, nu, stabilize, save)
#     _, arr_l2u_ns, arr_h1u_ns, _, _, arr_cond_ns = sensitivity_stokes(;n, M, u_exact, p_exact, f, g, ud, order, geometry, solver, βu0, γu1, γu2, γp, βp0, nu, stabilize, save)

#     plot!(r1, arr_δ[start:end], arr_cond[start:end], xaxis=:log, yaxis=:log, label="Stabilized γu2=$γu2")
#     plot!(r1, arr_δ[start:end], arr_cond_ns[start:end], xaxis=:log, yaxis=:log, label="Non-stabilized γu2=$γu2")

#     plot!(r2, arr_δ[start:end], arr_l2u[start:end], xaxis=:log, yaxis=:log, label="Stabilized γu2=$γu2")
#     plot!(r2, arr_δ[start:end], arr_l2u_ns[start:end], xaxis=:log, yaxis=:log, label="Non-stabilized γu2=$γu2")

#     plot!(r3, arr_δ[start:end], arr_h1u[start:end], xaxis=:log, yaxis=:log, label="Stabilized γu2=$γu2")
#     plot!(r3, arr_δ[start:end], arr_h1u_ns[start:end], xaxis=:log, yaxis=:log, label="Non-stabilized γu2=$γu2")
# end
# title!(r1, "sensitivity_stokes: varying γu2"); xlabel!(r1, "δ"); ylabel!(r1, "Condition number")
# title!(r2, "L2 error (velocity), varying γu2"); xlabel!(r2, "δ"); ylabel!(r2, "L2 norm")
# title!(r3, "H1 error (velocity), varying γu2"); xlabel!(r3, "δ"); ylabel!(r3, "H1 norm")
# display(r1); 
# display(r2); 
# display(r3)


# ###########################
# ### sensitivity_stokes: γp
# ###########################
# γp_arr = [0.01, 0.1, 1.0]
# s1 = plot(); s2 = plot(); s3 = plot()
# for γp in γp_arr
#     stabilize = true
#     arr_δ, arr_l2u, arr_h1u, _, _, arr_cond = sensitivity_stokes(;n, M, u_exact, p_exact, f, g, ud, order, geometry, solver, βu0, γu1, γu2, γp, βp0, nu, stabilize, save)
#     stabilize = false
#     _, arr_l2u_ns, arr_h1u_ns, _, _, arr_cond_ns = sensitivity_stokes(;n, M, u_exact, p_exact, f, g, ud, order, geometry, solver, βu0, γu1, γu2, γp, βp0, nu, stabilize, save)

#     plot!(s1, arr_δ[start:end], arr_cond[start:end], xaxis=:log, yaxis=:log, label="Stabilized γp=$γp")
#     plot!(s1, arr_δ[start:end], arr_cond_ns[start:end], xaxis=:log, yaxis=:log, label="Non-stabilized γp=$γp")

#     plot!(s2, arr_δ[start:end], arr_l2u[start:end], xaxis=:log, yaxis=:log, label="Stabilized γp=$γp")
#     plot!(s2, arr_δ[start:end], arr_l2u_ns[start:end], xaxis=:log, yaxis=:log, label="Non-stabilized γp=$γp")

#     plot!(s3, arr_δ[start:end], arr_h1u[start:end], xaxis=:log, yaxis=:log, label="Stabilized γp=$γp")
#     plot!(s3, arr_δ[start:end], arr_h1u_ns[start:end], xaxis=:log, yaxis=:log, label="Non-stabilized γp=$γp")
# end
# title!(s1, "sensitivity_stokes: varying γp"); xlabel!(s1, "δ"); ylabel!(s1, "Condition number")
# title!(s2, "L2 error (velocity), varying γp"); xlabel!(s2, "δ"); ylabel!(s2, "L2 norm")
# title!(s3, "H1 error (velocity), varying γp"); xlabel!(s3, "δ"); ylabel!(s3, "H1 norm")
# display(s1); 
# display(s2); 
# display(s3)


# ###########################
# ### sensitivity_stokes: βp
# ###########################
# βp_arr = [0.001, 0.01, 0.1]
# t1 = plot(); t2 = plot(); t3 = plot()
# for βp0 in βp_arr
#     stabilize = true
#     arr_δ, arr_l2u, arr_h1u, _, _, arr_cond = sensitivity_stokes(;n, M, u_exact, p_exact, f, g, ud, order, geometry, solver, βu0, γu1, γu2, γp, βp0, nu, stabilize, save)
#     stabilize = false
#     _, arr_l2u_ns, arr_h1u_ns, _, _, arr_cond_ns = sensitivity_stokes(;n, M, u_exact, p_exact, f, g, ud, order, geometry, solver, βu0, γu1, γu2, γp, βp0, nu, stabilize, save)

#     plot!(t1, arr_δ[start:end], arr_cond[start:end], xaxis=:log, yaxis=:log, label="Stabilized βp=$βp")
#     plot!(t1, arr_δ[start:end], arr_cond_ns[start:end], xaxis=:log, yaxis=:log, label="Non-stabilized βp=$βp")

#     plot!(t2, arr_δ[start:end], arr_l2u[start:end], xaxis=:log, yaxis=:log, label="Stabilized βp=$βp")
#     plot!(t2, arr_δ[start:end], arr_l2u_ns[start:end], xaxis=:log, yaxis=:log, label="Non-stabilized βp=$βp")

#     plot!(t3, arr_δ[start:end], arr_h1u[start:end], xaxis=:log, yaxis=:log, label="Stabilized βp=$βp")
#     plot!(t3, arr_δ[start:end], arr_h1u_ns[start:end], xaxis=:log, yaxis=:log, label="Non-stabilized βp=$βp")
# end
# title!(t1, "sensitivity_stokes: varying βp"); xlabel!(t1, "δ"); ylabel!(t1, "Condition number")
# title!(t2, "L2 error (velocity), varying βp"); xlabel!(t2, "δ"); ylabel!(t2, "L2 norm")
# title!(t3, "H1 error (velocity), varying βp"); xlabel!(t3, "δ"); ylabel!(t3, "H1 norm")
# display(t1); 
# display(t2); 
# display(t3)
