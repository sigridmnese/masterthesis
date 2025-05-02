using Gridap
using GridapEmbedded
using LinearAlgebra
using Plots
using STLCutters

import Gridap: ∇
∇(::typeof(u)) = ∇u

∇(u) === ∇u

# Manufactured solution
u_ex(x) = VectorValue(2*x[1]+cos(2*π*x[2]),-2*x[2]+sin(2*π*x[1]))
p_ex(x) = sin(2*π*x[1])

# %% Compute PDE data
# Viscosity
μ = 1.0 
f(x) = -divergence(ε(u_ex))(x) + ∇(p_ex)(x)
g(x) = tr(∇(u_ex)(x))

# Dirichlet data
u_D(x) = u_ex(x)

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
            geo = AnalyticalGeometry(x -> x[1]^2 + (5 * x[2]/4 - √(abs(x[1])))^2 - R)
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

function convergence(numb_it, geometry, order, solver, δ, γd, γg1, γg3, stabilization)
    "function to calculate convergence of the poisson solver, or the stokes solver, with or without stabilization"
    arr_l2 = zeros(Float64, numb_it)
    arr_h1 = zeros(Float64, numb_it)

    n = 2 .^ (2:(numb_it+1))
    h = 1.0 ./ n

    for i = 1:numb_it
        elapsed_time, solver_result = let
            t, val = @timed solver(n[i], geometry, order, δ, γd, γg1, γg3, stabilization)
            (val, t)
        end
        uh, u_exact, erru, l2u, h1u, condition_numb, Ω  = solver_result

        arr_l2[i] = l2u  #l2 error
        arr_h1[i] = h1u  #h1 error

        println("$i: Solved system in $elapsed_time seconds.")
    end
    return arr_l2, arr_h1, h
end

function sensitivity(n, M, u_exact, geometry, solver, stabilization)
    arr_δ = zeros(Float64, M-1)
    arr_l2 = zeros(Float64, M-1)
    arr_h1 = zeros(Float64, M-1)
    arr_cond = zeros(Float64, M-1)
    #loop to perturb the active geometry
    for i = 1:(M-1)
        δ = i/M

        elapsed_time, solver_result = let
            t, val = @timed solver(n, geometry, order, δ, γd, γg1, γg3, stabilization)
            (val, t)
        end
        uh, u_exact, erru, l2u, h1u, condition_numb, Ω  = solver_result

        arr_δ[i] = δ
        arr_l2[i] = l2u
        arr_h1[i] = h1u
        arr_cond[i] = condition_numb

        println("$i: Solved system in $elapsed_time seconds.")
    end
    return arr_δ, arr_l2, arr_h1, arr_cond
end

function stokes(n, geometry, δ, order, γd, γg1, γg3, γp, βu, βp, h, stabilization)
    p0 = Point(0.0, 0.0)
    # Define background mesh
    partition = (n, n)
    dim = length(partition)
    pmin = p0-1.0
    pmax = p0+1.0
    bgmodel = CartesianDiscreteModel(pmin,pmax,partition)

    # Define active and physical mesh
    geo = create_geometry(geometry, n, δ)
    cutgeo = cut(bgmodel,geo)

    # Defining active and physical meshes and inspecting
    Ω_act = Triangulation(cutgeo, ACTIVE)
    Ω = Triangulation(cutgeo, PHYSICAL)
    dΩ   = Measure(Ω, 2*order)

    # solving using nitsches metod. Establish function spaces, a, b and solve according to nitsches method with stabilization
    uh, ph, op = nitsche_stokes(Ω_act, Ω, cutgeo, order, γd, γg1, γg3, γp, βu, βp, h, stabilization)
    
    perr = p_ex - ph
    uerr = u_ex - uh
    # l2 and h1 normal
    l2(u) = √(∑( ∫( u*u )dΩ ))
    h1(u) = √(∑( ∫( u*u + ∇(u)⋅∇(u) )dΩ ))

    # condition number
    condition_numb= cond(Array(get_matrix(op)), 2)

    return uh, u, ph, p, uerr, perr, l2(uh - u), h1(uh - u), condition_numb, Ω_act
end

function nitsche_stokes(Ω_act, Ω, cutgeo, order, γd, γg1, γg3, γp, βu, βp, h, stabilization)
    ## Construct function spaces that are not restricted to zero boundary conditions. Will apply boundary conditions later...
    # First reference spaces
    reffe_u  = ReferenceFE(lagrangian,VectorValue{2, Float64},order)
    reffe_p = ReferenceFE(lagrangian,Float64, order - 1)

    V = TestFESpace(Ω_act, reffe_u,  conformity=:H1)
    Q = TestFESpace(Ω_act, reffe_p, conformity=:H1, constraint=:zeromean)

    U = TrialFESpace(V)
    P = TrialFESpace(Q)

    X = MultiFieldFESpace([U, P])
    Y = MultiFieldFESpace([V, Q])

    # Set up integration meshes, measures and normals
    Ω = Triangulation(cutgeo, PHYSICAL)
    Γ = EmbeddedBoundary(cutgeo)
    Fg = GhostSkeleton(cutgeo)

    # maybe this needs to be different...?
    degree = 2*order
    dΩ = Measure(Ω,degree)
    n_Γ = get_normal_vector(Γ)

    a((u,p),(v,q)) =
    ∫( ∇(v)⋅∇(u) - p * ∇(v) + ∇(u) *q )dΩ -
    ∫(v*(n_Γ⋅∇(u)) + (p ⋅ n_Γ) * v)dΓ   + (γd/h)*v*u  - (n_Γ⋅∇(v))*u 

    l(v) = ∫( (γd/h)*v*u - (n_Γ⋅∇(v))*u )dΓ

    # Solving
    op = AffineFEOperator(a, l, X, Y)
    uh, ph = solve(op)

    # %% Solve and write out results
   return uh, ph, op
end

# parameters
 # Nitsche and Ghost penalty stabilization parameter
 order = 2
 βu = 10.0*order^2
 γg1 = 1.0
 γg3 = 1.0
 γp = 0.1
 βp = 0.1
 n = 2^4 

 uh, u, ph, p, uerr, perr, l2, h1, condition_numb, Ω_act = stokes(40, "circle", 0, order, γd, γg1, γg3, γp, βu, βp, h, stabilization)
