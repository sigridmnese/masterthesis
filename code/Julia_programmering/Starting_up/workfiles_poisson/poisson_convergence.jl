using Gridap
using LinearAlgebra
using GridapEmbedded
using Plots
u_ex(x) = 0.01 - x[1]^2 - x[2]^2
f(x) = 4
∇u_ex(x) = VectorValue(-2*x[1], -2*x[2])

import Gridap: ∇
∇(::typeof(u_ex)) = ∇u_ex
∇(u_ex) === ∇u_ex

function create_geometry_2(name, n, δ)
    try
        if lowercase(name) == "circle"
            R  = 1.0
            geo = AnalyticalGeometry(x->(x[1]-δ/n)^2+(x[2]^2-δ/n)-R^2)
            return geo
        elseif lowercase(name) =="flower"
            LL= 2.1
            r0, r1 = 0.35*LL, 0.09*LL
            
            # using ! operator to define the interior
            geo = AnalyticalGeometry(x-> -r0 - r1*cos(5.0*atan(x[1], x[2]) - δ/n) +(x[1]^2 + x[2]^2)^0.5)   #lagt inn delta/n, kan fjernes?
            return geo
        elseif lowercase(name) =="heart"
            R  = 0.5
            geo = AnalyticalGeometry(x -> (x[1]-δ/n) ^2 + (5 * (x[2]-δ/n)/4 - √(abs(x[1] - δ/n )))^2 - R)
            return geo
        else
            error("$name er ikke en definert geometri!")
        end    
    catch e
        if e isa MethodError
            println("Ugyldig input for geometri!")
        else
            println("Ukjent feil:", e)
        end
    
    end
end

function create_geometry(name, n)
    try
        if lowercase(name) == "circle"
            R  = 1.0
            geo = AnalyticalGeometry(x->(x[1])^2+(x[2]^2)-R^2)
            return geo
        elseif lowercase(name) =="flower"
            LL= 2.1
            r0, r1 = 0.35*LL, 0.09*LL
            
            # using ! operator to define the interior
            geo = AnalyticalGeometry(x-> -r0 - r1*cos(5.0*atan(x[1], x[2]) ) +(x[1]^2 + x[2]^2)^0.5)   #lagt inn delta/n, kan fjernes?
            return geo
        elseif lowercase(name) =="heart"
            R  = 0.5
            geo = AnalyticalGeometry(x -> (x[1]) ^2 + (5 * (x[2])/4 - √(abs(x[1] )))^2 - R)
            return geo
        else
            error("$name er ikke en definert geometri!")
        end    
    catch e
        if e isa MethodError
            println("Ugyldig input for geometri!")
        else
            println("Ukjent feil:", e)
        end
    
    end
end

function jump_nn(u,n)
    return ( (n.plus⋅∇∇(u).plus)⋅ n.plus - (n.minus ⋅ ∇∇(u).minus) ⋅ n.minus )
end

function poisson_solver(n, u_exact, lhs, order, geometry, γd, γg1, γg3, stabilize, δ, save = false)
    """
    n: number of grid elements. Powers of 2 for simplicity and convergence estimates.
    u_exact: exact solution for method of manufactured solutions
    order: order of polynomial degree. 
    lhs: for poisson, -Δ u_ex 
    geometry: optional between "Circle", "Flower", "Heart"
    δ: perturbation of cut
    """
    #Construct geometries and mesh and cut
    # flytte bgmodel istedenfor!
    # sende nye figurer til andre etterhvert

    # Background model
    a = 1.11
    domain = (-a + δ, a + δ, -a + δ, a + δ)
    pmin = Point(-a + δ, -a + δ)            # kanskje fjerne /n?
    pmax = Point(a + δ, a + δ)
    partition = (n,n)
    h = (pmax - pmin)[1]/partition[1]
    bgmodel = CartesianDiscreteModel(pmin, pmax, partition)

    # Creating desired geometry for active mesh. Possible to perturb the cut by perturbing δ. 
    geo = create_geometry(geometry, n)
    # Cut the background model using active geometry
    cutgeo = cut(bgmodel, geo)

    # Set up interpolation mesh and function spaces
    Ω_act = Triangulation(cutgeo, ACTIVE)

    # Construct function spaces
    V = TestFESpace(Ω_act, ReferenceFE(lagrangian, Float64, order),conformity=:H1)
    U = TrialFESpace(V)

    # Set up integration meshes, measures and normals
    Ω = Triangulation(cutgeo, PHYSICAL)
    Γ = EmbeddedBoundary(cutgeo)
    Fg = GhostSkeleton(cutgeo)

    # Set up integration measures
    degree = 2*order
    dΩ   = Measure(Ω, degree)
    dΓ   = Measure(Γ, degree)
    dFg  = Measure(Fg, degree)

    # Set up normal vectors
    n_Γ = get_normal_vector(Γ)
    n_Fg = get_normal_vector(Fg)

    # Define weak form
    # define bilinear form
    a(u, v) = ∫( ∇(u)⋅∇(v) ) * dΩ  +
    ∫( (γd/h)*u*v  - u*(n_Γ⋅∇(v)) - (n_Γ⋅∇(u))*v ) * dΓ

    g(u, v) =  ∫( (γg1*h)*jump(n_Fg⋅∇(u))*jump(n_Fg⋅∇(v))  #adding stabilization term
               + (γg3*h^3)*jump_nn(u,n_Fg)⋅jump_nn(v,n_Fg))*dFg   #h^3 stabilization for second order... should maybe not ave anything to say?
    # linear form
    l(v) =
        ∫( lhs*v ) * dΩ +
        ∫( u_exact*( (γd/h)*v - (n_Γ⋅∇(v)) )  ) * dΓ

    if stabilize
        A(u, v) = a(u, v) + g(u, v)
       
        # FE problem
        op = AffineFEOperator(A,l,U,V)
        uh = solve(op)
    else
        # FE problem
        op = AffineFEOperator(a,l,U,V)
        uh = solve(op)
    end
    
    # error of u over entire domain
    erru = uh - u_exact

    # l2 and h1 normal
    l2(u) = √(∑( ∫( u*u )dΩ ))
    h1(u) = √(∑( ∫( u*u + ∇(u)⋅∇(u) )dΩ ))

    # condition number
    condition_numb= cond(Array(get_matrix(op)),2)

    if save
        writevtk(Ω_act, "mesh_act_$geometry $δ.vtu")
        writevtk(Ω, "mesh_$geometry $δ.vtu")
        writevtk(Γ, "surface_gamma_d_$geometry $δ.vtu")
        writevtk(Fg, "ghost_facets_$geometry $δ.vtu")
        writevtk(Ω, "$n $geometry $order $δ.vtu", cellfields=["u_ex" => u_exact, "uh"=>uh, "erru"=> erru]) #, "erru" => erru])
    end

return uh, u_exact, erru, l2(uh - u_exact), h1(uh - u_exact), condition_numb, Ω_act, Ω
end

function convergence(numb_it, u_exact, lhs, order, geometry, solver, δ, γd, γg1, γg3, stabilize, save = false)
    "function to calculate convergence of the poisson solver, or the stokes solver, with or without stabilization"
    "function to calculate convergence of the poisson solver, or the stokes solver, with or without stabilization"
    arr_l2 = zeros(Float64, numb_it)
    arr_h1 = zeros(Float64, numb_it)

    n = 2 .^ (1:(numb_it))
    h = 1.0 ./ n

    for i = 1:numb_it
        elapsed_time, solver_result = let
            t, val = @timed solver(n[i], u_exact, lhs, order, geometry, γd, γg1, γg3, stabilize, δ, save)
            (val, t) 
        end
        uh, u_exact, erru, l2u, h1u, condition_numb, Ω_act = solver_result

        arr_l2[i] = l2u  #l2 error
        arr_h1[i] = h1u  #h1 error

        println("$i: Solved system in $elapsed_time seconds.")
    end
    return arr_l2, arr_h1, h
end

function sensitivity(n, M, u_exact, lhs, order, geometry, solver, δ, γd, γ1, γ3, stabilize, save = false)
    arr_δ = zeros(Float64, M-1)
    arr_l2 = zeros(Float64, M-1)
    arr_h1 = zeros(Float64, M-1)
    arr_cond = zeros(Float64, M-1)
    #loop to perturb the active geometry
    for i = 1:(M-1)
        #delta_n = n*h/N_max*(1,1)/\sqrt{2}
        δ = i/n/M *1.1/ sqrt(2)
        #δ = i/M

        elapsed_time, solver_result = let
            t, val = @timed solver(n, u_exact, lhs, order, geometry, γd, γg1, γg3, stabilize, δ, save)
            (val, t)
        end
        uh, u_exact, erru, l2u, h1u, condition_numb = solver_result

        arr_δ[i] = δ
        arr_l2[i] = l2u
        arr_h1[i] = h1u
        arr_cond[i] = condition_numb

        if (i-1) % 100 == 0
            println("$i: Solved system in $elapsed_time seconds.")
            save = true
        else
            save = false
        end
    end
    return arr_δ, arr_l2, arr_h1, arr_cond
end

function logplot(x, yarr, start, stop, title, xlabel, ylabel, labels)
    plot(x, yarr[1], xaxis=:log, yaxis=:log, marker=:o, lw=2, label=labels[1])
    for i = 2:lastindex(x)
        plot(x, yarr[i], xaxis=:log, yaxis=:log, marker=:o, lw=2, label=labels[i])
    end
    # Legger til aksetitler og tittel
    xlabel!("Mesh size h")
    ylabel!("Error")
    title!("Convergence of Poisson Solver")
end
# defining parameters
γd = 10
γg1 = 10
γg3 = 0.1
stabilize = true
n = 16
uh, u_exact, erru, l2, h1, condition_numb, Ω_active, Ω = poisson_solver(n, u_ex, f, 1, "heart", γd, γg1, γg3, stabilize, 0, false)
writevtk(Ω, "testing_stabilized_poisson_circle.vtu", cellfields=["uh"=>uh, "u_ex"=>u_exact, "erru" => erru])

####### Convergence test #######
# med stabilisering, order 1
# numb_it = 6                         # Sånn koden er implementert nå så er det fra 2^1 til 2^{numb_it}. Tar kort tid å kjøre for numb_it = 6
# order = 2                           # når jeg øker orden så øker kjøretid veeeeldig !! Bør vurdere å skru ned numb_it samtidig. 244 sekunder når jeg har på order = 2 for kjøringen 2^6
# δ = 0                               # kan også se ut til at feilen havner på maskinnivå? vet ikke helt, men mulig å eksperimentere med dette. 
# # med stabilisering, order 1
# stabilization = true
# solver = poisson_solver
# arr_l2_1_stab, arr_h1_1_stab, h = convergence(numb_it,  order, "circle",solver, δ, γd, γg1, γg3, stabilization, false)
# start = 2
# plot(h[start:end], arr_l2_1_stab[start:end], xaxis=:log, yaxis=:log, marker=:o, lw=2, label="L2 error, order 1 (stab)")
# plot!(h[start:end], arr_h1_1_stab[start:end], marker=:s, lw=2, label="H1 error, order 1 (Stab)")

# # Uten stabilisering, order 1
# stabilization = false
# solver = poisson_solver
# arr_l2_1_nostab, arr_h1_1_nostab, h = convergence(numb_it, order, "circle", solver, δ, γd, γg1, γg3, stabilization, true)

# plot!(h[start:end], arr_l2_1_nostab[start:end], marker=:s, lw=2, label="L2 error, order 1 (no Stab)")
# plot!(h[start:end], arr_h1_1_nostab[start:end], marker=:s, lw=2, label="H1 error, order 1 (no Stab)")

# # Legger til aksetitler og tittel
# xlabel!("Mesh size h")
# ylabel!("Error")
# title!("Convergence of Poisson Solver")

####### Sensitivity test #######
# M = 2000
# n = 32
# order = 1
# stabilize = true
# arr_δ, arr_l2, arr_h1, arr_cond = sensitivity(n, M, u_ex, f, order, "circle", poisson_solver, 0, γd, γg1, γg3, stabilize, false)
# stabilize = false
# arr_δ_nostab, arr_l2_nostab, arr_h1_nostab, arr_cond_nostab = sensitivity(n, M, u_ex, f, order, "circle", poisson_solver, 0,γd, γg1, γg3, stabilize, false)
# start = 1
# using Plots
# e = 1999
# plot(arr_δ[start:e], arr_cond[start:e], yaxis=:log, lw=2, label="Stabilized")
# plot!(arr_δ_nostab[start:e], arr_cond_nostab[start:e], yaxis=:log, lw=2, label="Non-stabilized")
# #plot(arr_δ[start:end], arr_l2[start:end], xaxis=:log, yaxis=:log, marker=:o, lw=2, label="L^2 norm")
# #plot!(arr_δ[start:end], arr_h1[start:end], xaxis=:log, yaxis=:log, marker=:o, lw=2, label="H^1 norm")
# xlabel!("Perturbation δ")
# ylabel!("Condition number")
# title!("Sensitivity analysis of non-stabilized poisson")

# #### Varying geometry
# geometry_arr = ["circle", "flower", "heart"]
# s = plot()  # Lager et tomt plott
# for i = 1:3
#     geo = geometry_arr[i]
#     stabilize = true
#     arr_δ, arr_l2, arr_h1, arr_cond = sensitivity(n, M, u_ex, f, order, geometry_arr[i], poisson_solver, 0, γd, γg1, γg3, stabilize, false)
#     stabilize = false
#     arr_δ_nostab, arr_l2_nostab, arr_h1_nostab, arr_cond_nostab = sensitivity(n, M, u_ex, f, order, geometry_arr[i], poisson_solver, 0,γd_arr[i], γg1, γg3, stabilize, false)
#     plot!(s, arr_δ[start:end], arr_cond[start:end], xaxis=:log, yaxis=:log, lw=2, label="Stabilized ($geo)")
#     plot!(s, arr_δ_nostab[start:end], arr_cond_nostab[start:end], xaxis=:log, yaxis=:log, lw=2, label="Non-stabilized ($geo)")
# end
# title!("Sensitivity analysis varing geometries")

# display(s)

# start = 1
# # ### Varying γd
# γd_arr = [0.1, 1, 10]
# p = plot()  # Lager et tomt plott
# for i = 1:3
#     γd = γd_arr[i]
#     stabilize = true
#     arr_δ, arr_l2, arr_h1, arr_cond = sensitivity(n, M, u_ex, f, order, "circle", poisson_solver, 0, γd_arr[i], γg1, γg3, stabilize, false)
#     stabilize = false
#     arr_δ_nostab, arr_l2_nostab, arr_h1_nostab, arr_cond_nostab = sensitivity(n, M, u_ex, f, order, "circle", poisson_solver, 0,γd_arr[i], γg1, γg3, stabilize, false)
#     plot!(p, arr_δ[start:end], arr_cond[start:end], xaxis=:log, yaxis=:log, lw=2, label="Stabilized with γd = $γd")
#     plot!(p, arr_δ_nostab[start:end], arr_cond_nostab[start:end], xaxis=:log, yaxis=:log, lw=2, label="Non-stabilized with γd = $γd")
# end
# title!("Sensitivity analysis varing γd")

# display(p)

# # ### Varying γg1
# γg1_arr = [0.001, 0.01, 0.1, 1]
# γd = 0.1
# q = plot()  # Lager et tomt plott
# for i = 1:3
#     γg1 = γg1_arr[i]
#     stabilize = true
#     arr_δ, arr_l2, arr_h1, arr_cond = sensitivity(n, M, u_ex, f, order, "circle", poisson_solver, 0, γd, γg1_arr[i], γg3, stabilize, false)
#     stabilize = false
#     arr_δ_nostab, arr_l2_nostab, arr_h1_nostab, arr_cond_nostab = sensitivity(n, M, u_ex, f, order, "circle", poisson_solver, 0,γd, γg1_arr[i], γg3, stabilize, false)
#     plot!(q, arr_δ[start:end], arr_cond[start:end], xaxis=:log, yaxis=:log, lw=2, label="Stabilized with γg1 = $γg1")
#     plot!(q, arr_δ_nostab[start:end], arr_cond_nostab[start:end], xaxis=:log, yaxis=:log, lw=2, label="Non-stabilized with γg1 = $γg1")
# end
# title!("Sensitivity analysis varying γg1")

# display(q)

### varying γg3
# γg3_arr = [0.1, 1, 10]
# γd = 0.1
# r = plot()  #lager et tomt plott
# for i = 1:3
#     γg3 = γg3_arr[i]
#     stabilize = true
#     arr_δ, arr_l2, arr_h1, arr_cond = sensitivity(n, M, u_ex, f, order, "circle", poisson_solver, 0, γd, γg1, γg3_arr[i], stabilize, false)
#     stabilize = false
#     arr_δ_nostab, arr_l2_nostab, arr_h1_nostab, arr_cond_nostab = sensitivity(n, M, u_ex, f, order, "circle", poisson_solver, 0,γd, γg1, γg3_arr[i], stabilize, false)
#     plot!(r, arr_δ[start:end], arr_cond[start:end], xaxis=:log, yaxis=:log, lw=2, label="Stabilized with γg3 = $γg3")
#     plot!(r, arr_δ_nostab[start:end], arr_cond_nostab[start:end], xaxis=:log, yaxis=:log, lw=2, label="Non-stabilized with γg3 = $γg3")
# end
# title!("Sensitivity analysis varing γg3")

# display(r)