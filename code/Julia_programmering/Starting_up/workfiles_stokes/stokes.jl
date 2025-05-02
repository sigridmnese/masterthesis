using Gridap
using LinearAlgebra
using GridapEmbedded
using Plots
u_ex(x) = VectorValue(2*x[1] + cos(2*π*x[2]), -2*x[2] + sin(2*π*x[1]))
p_ex(x) = sin(2*π*x[1])
f(x) = -divergence(ε(u_ex))(x) + ∇(p_ex)(x)                       # skal fungere å lagre den sånn her...
g(x) = tr(∇(u_ex)(x))
ud(x) = u_ex(x)

function create_geometry(name, n, δ)
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

function jump_nn(u,n)
    return ( (n.plus⋅∇∇(u).plus)⋅ n.plus - (n.minus ⋅ ∇∇(u).minus) ⋅ n.minus )
end

function stokes_solver(n, u_exact, p_exact, f, g, ud, order, geometry, βu, γu1, γu2, γp, βp, stabilize, δ, save = false)
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
    n = 60
    partition = (n, n)
    dim = length(partition)
    pmin = Point(-1.11, -1.11)
    pmax = Point(1.11, 1.11)
    bgmodel = CartesianDiscreteModel(pmin,pmax,partition)
    
    geo = create_geometry(geometry, n, δ)

    # Define active and physical mesh
    cutgeo = cut(bgmodel,geo)
    # Needed to extract boundary facets from active mesh
    cutgeo_facets = cut_facets(bgmodel,geo)
    Ω_bg = Triangulation(bgmodel)
    Ω_act = Triangulation(cutgeo, ACTIVE)
    Ω = Triangulation(cutgeo, PHYSICAL)

    # Embedded boundary
    # Has Dirichlet for u and Neumann for pF
    Γd = EmbeddedBoundary(cutgeo)
    n_Γd = get_normal_vector(Γd)

    # Outer mesh boundary 
    # Has Dirichlet for pF and no-stress for u
    Γs = BoundaryTriangulation(cutgeo_facets, ACTIVE_IN)
    n_Γs = get_normal_vector(Γs)

    # Get ghost penalty facets
    Fg = GhostSkeleton(cutgeo)
    n_Fg = get_normal_vector(Fg)

    # All interior facets of active mesh
    Fi = SkeletonTriangulation(Ω_act)
    n_Fi = get_normal_vector(Fi)

    # Function spaces orders
    order_p = order - 1

    # Define measures
    degree = 2*order
    dΩ = Measure(Ω,degree)
    dΓd = Measure(Γd, degree)
    dΓs = Measure(Γs, degree)
    dFg = Measure(Fg, degree)
    dFi = Measure(Fi, degree)

    # %% Define weak formulation
    # Define function spaces
    reffe_u  = ReferenceFE(lagrangian,VectorValue{dim, Float64},order)
    reffe_p = ReferenceFE(lagrangian,Float64, order_p)

    V = TestFESpace(Ω_act, reffe_u,  conformity=:H1)
    Q = TestFESpace(Ω_act, reffe_p, conformity=:H1, constraint=:zeromean)

    U = TrialFESpace(V)
    P = TrialFESpace(Q)

    X = MultiFieldFESpace([U, P])
    Y = MultiFieldFESpace([V, Q])

    # Physical parameters (potentially rescaled after time-step discretization)

    # Weak formulation

    # mesh size
    h = (pmax-pmin)[1]/partition[1]

    # TODO: For now we impose Dirichlet boundary conditions on the whole boundary, so need to change b.c. on Γs
    # Bilinear forms
    # a(u,v) =  ( ∫( ∇(u)⊙∇(v))dΩ 
            #  + ∫( -(n_Γd⋅∇(u))⋅v - u⋅(n_Γd⋅∇(v)) + (βu/h)*u⋅v)dΓd
    a(u,v) =  ( ∫( ε(u)⊙ε(v))dΩ 
            + ∫( -(n_Γd⋅ε(u))⋅v -(n_Γd⋅ε(v))⋅u + βu/h*u⋅v)dΓd
            + ∫( -(n_Γs⋅ε(u))⋅v -(n_Γs⋅ε(v))⋅u + βu/h*u⋅v)dΓs
    )

    b(v, p) = ( ∫(-1*(∇⋅v*p))dΩ 
                + ∫((n_Γd⋅v)*p)dΓd
                + ∫(v⋅n_Γs*p)dΓs
    )

    # Define ghost penalty
    function jump_nn(u,n)
        return ( (n.plus⋅∇∇(u).plus)⋅ n.plus - (n.minus ⋅ ∇∇(u).minus) ⋅ n.minus )
    end

    g_u(u,v) = ( ∫( (γu1*h)*jump(n_Fg⋅∇(u))⋅jump(n_Fg⋅∇(v)) )dFg 
                +  
                ∫( (γu2*h^3)*jump_nn(u,n_Fg)⋅jump_nn(v,n_Fg) )dFg
    )

    g_p(p,q) = ∫( (γp*h^3)*jump(n_Fg⋅∇(p))*jump(n_Fg⋅∇(q)) )dFg
    s_p(p,q) = ∫( (βp*h^3)*jump(n_Fi⋅∇(p))*jump(n_Fi⋅∇(q)) )dFi

    # Total bilinear form
    A((u,p),(v,q)) = a(u,v) + g_u(u,v) + b(v,p) + b(u,q) - g_p(p,q)
    # Ah((u,p),(v,q)) = a(u,v) + g_u(u,v) + b(v,p) + b(u,q) - g_p(p,q) - s_p(p,q)

    # Linear forms
    L((v,q)) = ( ∫(f⋅v - g⋅q)dΩ
    # lh((v,q)) = ( ∫(fdiv⋅v - g⋅q)dΩ
                #  + ∫( -1*(n_Γd⋅∇(v))⋅u_D + βu/h*(u_D⋅v) + u_D⋅n_Γd*q)dΓd 
                + ∫( -1*(n_Γd⋅ε(v))⋅ud + βu/h*(ud⋅v) + ud⋅n_Γd*q)dΓd 
                + ∫( -1*(n_Γs⋅ε(v))⋅ud + βu/h*(ud⋅v) + ud⋅n_Γs*q)dΓs
    )

    # %% Solve and write out results
    # Solve
    op = AffineFEOperator(A,L,X,Y)
    uh, ph = solve(op)

    # error of u over entire domain
    erru = uh - u_exact
    errp = ph - p_exact

    # l2 and h1 normal
    #l2(u) = √(∑( ∫( u ⋅ u )dΩ ))
    #h1(u) = √(∑( ∫( u ⋅ u + ∇(u) ⋅ ∇(u) )dΩ ))

    l2(u) = (sum( ∫( u ⋅ u )*dΩ ))
    h1_semi(u) = sum(∫(∇(u) ⊙ ∇(u))*dΩ)             #kan være jeg må ha en annen norm for p enn for u ???

    # condition number
    condition_numb= cond(Array(get_matrix(op)),2)

    if save
        writevtk(Ω, "stokes $n $geometry $order", cellfields=["u_ex" => u_exact, "uh"=>uh, "erru"=> erru, "p_ex" => p_exact, "ph"=>ph, "errp"=> errp]) #, "erru" => erru])
    end
return uh, u_exact, erru, l2(uh - u_exact), h1_semi(uh - u_exact), ph, p_exact, errp, l2(ph - p_exact), h1_semi(ph - p_exact), condition_numb, Ω_act
end

function convergence(numb_it, u_exact, p_exact, f, g, ud, order, geometry, solver, δ, γd, γg1, γg3, stabilize, save = false)
    "function to calculate convergence of the poisson solver, or the stokes solver, with or without stabilization"
    "function to calculate convergence of the poisson solver, or the stokes solver, with or without stabilization"
    uarr_l2 = zeros(Float64, numb_it)
    uarr_h1 = zeros(Float64, numb_it)
    parr_l2 = zeros(Float64, numb_it)
    parr_h1 = zeros(Float64, numb_it)

    n = 2 .^ (1:(numb_it))
    h = 1.0 ./ n

    for i = 1:numb_it
        elapsed_time, solver_result = let
            t, val = @timed solver(n[i], u_exact, p_exact, f, g, ud, order, geometry, γd, γg1, γg3, stabilize, δ, save)
            (val, t) 
        end

        uh, u_exact, erru, l2u, h1u, ph, p_exact, errp, l2p, h1p, condition_numb, Ω_act = solver_result

        uarr_l2[i] = l2u  #l2 error
        uarr_h1[i] = h1u  #h1 error
        parr_l2[i] = l2p  #l2 error
        parr_h1[i] = h1p  #h1 error

        println("$i: Solved system in $elapsed_time seconds.")
    end
    return uarr_l2, uarr_h1, parr_l2, parr_l2, h
end

function sensitivity(n, M, u_exact, lhs, order, geometry, solver, δ, γd, γ1, γ3, stabilize, save = false)
    arr_δ = zeros(Float64, M-1)
    uarr_l2 = zeros(Float64, M-1)
    uarr_h1 = zeros(Float64, M-1)
    parr_l2 = zeros(Float64, M-1)
    parr_h1 = zeros(Float64, M-1)
    arr_cond = zeros(Float64, M-1)
    #loop to perturb the active geometry
    for i = 1:(M-1)
        δ = i/M

        elapsed_time, solver_result = let
            t, val = @timed solver(n, u_exact, lhs, order, geometry, γd, γg1, γg3, stabilize, δ, save)
            (val, t)
        end
        uh, u_exact, erru, l2u, h1u, ph, p_exact, errp, l2p, h1p, condition_numb, Ω_act = solver_result

        arr_δ[i] = δ
        uarr_l2[i] = l2u
        uarr_h1[i] = h1u
        parr_l2[i] = l2p
        parr_h1[i] = h1p
        arr_cond[i] = condition_numb

        println("$i: Solved system in $elapsed_time seconds.")
    end
    return arr_δ, uarr_l2, uarr_h1, parr_l2, parr_h1, arr_cond
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
# Nitsche and Ghost penalty stabilization parameter
βu = 10.0*2^2
γu1 = 1.0
γu2 = 1.0
γp = 0.1
βp = 0.1

# defining parameters
γd = 1
γg1 = 1
γg3 = 1
β0 = 0.1              # deffo change these!
β1 = 0.1
β2 = 0.1
δ = 0
order = 2           # taylor hood elements. P2 and P1 elements
stabilize = true
n = 16
uh, u_exact, erru, l2, h1, condition_numb, Ω_active = stokes_solver(n, u_ex, p_ex, f, g, ud, order, "circle", βu, γu1, γu2, γp, βp, stabilize, δ, true)
#writevtk(Ω_active, "testing_stabilized_poisson_circle_order3.vtu", cellfields=["uh"=>uh, "u_ex"=>u_exact, "erru" => erru])

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
# M = 50
# n = 25
# order = 1
# stabilize = true
# arr_δ, arr_l2, arr_h1, arr_cond = sensitivity(n, M, u_ex, f, order, "circle", poisson_solver, 0, γd, γg1, γg3, stabilize, false)
# stabilize = false
# arr_δ_nostab, arr_l2_nostab, arr_h1_nostab, arr_cond_nostab = sensitivity(n, M, u_ex, f, order, "circle", poisson_solver, 0,γd, γg1, γg3, stabilize, false)
# start = 1
# using Plots
# plot(arr_δ[start:end], arr_cond[start:end], xaxis=:log, yaxis=:log, lw=2, label="Stabilized")
# plot!(arr_δ_nostab[start:end], arr_cond_nostab[start:end], xaxis=:log, yaxis=:log, lw=2, label="Non-stabilized")
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
