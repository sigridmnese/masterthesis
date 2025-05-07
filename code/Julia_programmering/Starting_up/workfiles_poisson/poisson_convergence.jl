using Gridap
using LinearAlgebra
using GridapEmbedded
using Plots
include("C:\\Users\\Sigri\\Documents\\Master\\report\\code\\Julia_programmering\\Starting_up\\utils\\utils.jl")

u_ex(x) = 0.01 - x[1]^2 - x[2]^2
f(x) = 4
∇u_ex(x) = VectorValue(-2*x[1], -2*x[2])

import Gridap: ∇
∇(::typeof(u_ex)) = ∇u_ex
∇(u_ex) === ∇u_ex

# defining parameters
γd = 10
γg1 = 10
γg3 = 0.1
stabilize = true
n = 16
uh, u_exact, erru, _, _, _, Ω_active, Ω = poisson_solver(n, u_ex, f, 1, "heart", γd, γg1, γg3, stabilize, 0, false)
#writevtk(Ω, "testing_stabilized_poisson_circle.vtu", cellfields=["uh"=>uh, "u_ex"=>u_exact, "erru" => erru])

####### convergence_poisson test #######
# med stabilisering, order 1
# numb_it = 7                         # Sånn koden er implementert nå så er det fra 2^1 til 2^{numb_it}. Tar kort tid å kjøre for numb_it = 6
# order = 1                           # når jeg øker orden så øker kjøretid veeeeldig !! Bør vurdere å skru ned numb_it samtidig. 244 sekunder når jeg har på order = 2 for kjøringen 2^6
# δ = 0                               # kan også se ut til at feilen havner på maskinnivå? vet ikke helt, men mulig å eksperimentere med dette. 
# # med stabilisering, order 1
# stabilization = true
# solver = poisson_solver
# arr_l2_1_stab, arr_h1_1_stab, h = convergence_poisson(numb_it, u_ex, f, order, "circle",solver, δ, γd, γg1, γg3, stabilization, false)
# start = 2
# plot(h[start:end], arr_l2_1_stab[start:end], xaxis=:log, yaxis=:log, marker=:o, lw=2, label="L2 error, order 1 (stab)")
# plot!(h[start:end], arr_h1_1_stab[start:end], marker=:s, lw=2, label="H1 error, order 1 (Stab)")

# # Uten stabilisering, order 1
# stabilization = false
# solver = poisson_solver
# arr_l2_1_nostab, arr_h1_1_nostab, h = convergence_poisson(numb_it, u_ex, f, order, "circle", solver, δ, γd, γg1, γg3, stabilization, true)

# plot!(h[start:end], arr_l2_1_nostab[start:end], marker=:s, lw=2, label="L2 error, order 1 (no Stab)")
# plot!(h[start:end], arr_h1_1_nostab[start:end], marker=:s, lw=2, label="H1 error, order 1 (no Stab)")

# # Legger til aksetitler og tittel
# xlabel!("Mesh size h")
# ylabel!("Error")
# title!("convergence_poisson of Poisson Solver")

####### sensitivity_poisson test #######
M = 2000
n = 16
order = 1
stabilize = true
arr_δ, arr_l2, arr_h1, arr_cond = sensitivity_poisson(n, M, u_ex, f, order, "circle", poisson_solver, 0, γd, γg1, γg3, stabilize, false)
stabilize = false
arr_δ_nostab, arr_l2_nostab, arr_h1_nostab, arr_cond_nostab = sensitivity_poisson(n, M, u_ex, f, order, "circle", poisson_solver, 0,γd, γg1, γg3, stabilize, false)
start = 1
using Plots
e = 1999
#condition numbers
plot(arr_δ[start:e], arr_cond[start:e], yaxis=:log, lw=2, label="Stabilized")
plot!(arr_δ_nostab[start:e], arr_cond_nostab[start:e], yaxis=:log, lw=2, label="Non-stabilized")
xlabel!("Perturbation δ")
ylabel!("Condition number")
title!("sensitivity analysis of cutFEM poisson")
savefig("C:\\Users\\Sigri\\Documents\\Master\\report\\results\\poisson\\sensitivity_n16_order1_M2000_condition_number")

#errors
plot(arr_δ[start:end], arr_l2[start:end], xaxis=:log, yaxis=:log, marker=:o, lw=2, label="L^2 norm")
plot!(arr_δ[start:end], arr_h1[start:end], xaxis=:log, yaxis=:log, marker=:o, lw=2, label="H^1 norm")
xlabel!("Perturbation δ")
ylabel!("Error")
title!("sensitivity analysis of cutFEM poisson")
savefig("C:\\Users\\Sigri\\Documents\\Master\\report\\results\\poisson\\sensitivity_n16_order1_M2000_errors")


# #### Varying geometry
# geometry_arr = ["circle", "flower", "heart"]
# s = plot()  # Lager et tomt plott
# for i = 1:3
#     geo = geometry_arr[i]
#     stabilize = true
#     arr_δ, arr_l2, arr_h1, arr_cond = sensitivity_poisson(n, M, u_ex, f, order, geometry_arr[i], poisson_solver, 0, γd, γg1, γg3, stabilize, false)
#     stabilize = false
#     arr_δ_nostab, arr_l2_nostab, arr_h1_nostab, arr_cond_nostab = sensitivity_poisson(n, M, u_ex, f, order, geometry_arr[i], poisson_solver, 0,γd_arr[i], γg1, γg3, stabilize, false)
#     plot!(s, arr_δ[start:end], arr_cond[start:end], xaxis=:log, yaxis=:log, lw=2, label="Stabilized ($geo)")
#     plot!(s, arr_δ_nostab[start:end], arr_cond_nostab[start:end], xaxis=:log, yaxis=:log, lw=2, label="Non-stabilized ($geo)")
# end
# title!("sensitivity_poisson analysis varing geometries")

# display(s)

# start = 1
# # ### Varying γd
# γd_arr = [0.1, 1, 10]
# p = plot()  # Lager et tomt plott
# for i = 1:3
#     γd = γd_arr[i]
#     stabilize = true
#     arr_δ, arr_l2, arr_h1, arr_cond = sensitivity_poisson(n, M, u_ex, f, order, "circle", poisson_solver, 0, γd_arr[i], γg1, γg3, stabilize, false)
#     stabilize = false
#     arr_δ_nostab, arr_l2_nostab, arr_h1_nostab, arr_cond_nostab = sensitivity_poisson(n, M, u_ex, f, order, "circle", poisson_solver, 0,γd_arr[i], γg1, γg3, stabilize, false)
#     plot!(p, arr_δ[start:end], arr_cond[start:end], xaxis=:log, yaxis=:log, lw=2, label="Stabilized with γd = $γd")
#     plot!(p, arr_δ_nostab[start:end], arr_cond_nostab[start:end], xaxis=:log, yaxis=:log, lw=2, label="Non-stabilized with γd = $γd")
# end
# title!("sensitivity_poisson analysis varing γd")

# display(p)

# # ### Varying γg1
# γg1_arr = [0.001, 0.01, 0.1, 1]
# γd = 0.1
# q = plot()  # Lager et tomt plott
# for i = 1:3
#     γg1 = γg1_arr[i]
#     stabilize = true
#     arr_δ, arr_l2, arr_h1, arr_cond = sensitivity_poisson(n, M, u_ex, f, order, "circle", poisson_solver, 0, γd, γg1_arr[i], γg3, stabilize, false)
#     stabilize = false
#     arr_δ_nostab, arr_l2_nostab, arr_h1_nostab, arr_cond_nostab = sensitivity_poisson(n, M, u_ex, f, order, "circle", poisson_solver, 0,γd, γg1_arr[i], γg3, stabilize, false)
#     plot!(q, arr_δ[start:end], arr_cond[start:end], xaxis=:log, yaxis=:log, lw=2, label="Stabilized with γg1 = $γg1")
#     plot!(q, arr_δ_nostab[start:end], arr_cond_nostab[start:end], xaxis=:log, yaxis=:log, lw=2, label="Non-stabilized with γg1 = $γg1")
# end
# title!("sensitivity_poisson analysis varying γg1")

# display(q)

### varying γg3
# γg3_arr = [0.1, 1, 10]
# γd = 0.1
# r = plot()  #lager et tomt plott
# for i = 1:3
#     γg3 = γg3_arr[i]
#     stabilize = true
#     arr_δ, arr_l2, arr_h1, arr_cond = sensitivity_poisson(n, M, u_ex, f, order, "circle", poisson_solver, 0, γd, γg1, γg3_arr[i], stabilize, false)
#     stabilize = false
#     arr_δ_nostab, arr_l2_nostab, arr_h1_nostab, arr_cond_nostab = sensitivity_poisson(n, M, u_ex, f, order, "circle", poisson_solver, 0,γd, γg1, γg3_arr[i], stabilize, false)
#     plot!(r, arr_δ[start:end], arr_cond[start:end], xaxis=:log, yaxis=:log, lw=2, label="Stabilized with γg3 = $γg3")
#     plot!(r, arr_δ_nostab[start:end], arr_cond_nostab[start:end], xaxis=:log, yaxis=:log, lw=2, label="Non-stabilized with γg3 = $γg3")
# end
# title!("sensitivity_poisson analysis varing γg3")

# display(r)