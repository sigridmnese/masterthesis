using Gridap
using GridapEmbedded
using STLCutters
using LinearAlgebra
using Plots

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

function create_geometry_old(name, n, δ = 0)
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
          geo = AnalyticalGeometry(x-> -r0 - r1*cos(5.0*atan(x[1], x[2])) +(x[1]^2 + x[2]^2)^0.5)   #lagt inn delta/n, kan fjernes?
          return geo
      elseif lowercase(name) =="heart"
          R  = 0.7
          geo = AnalyticalGeometry(x -> (-R     +   (x[1])^2   +     ( 5/4 * x[2]   -     sqrt(abs(x[1]))  )^2   ))
          #geo = AnalyticalGeometry(x -> (x[1])^2 + (5/4 * (x[2]) - √(abs(x[1])))^2 - R)
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

# Define ghost penalty
function jump_nn(u,n)
  return ( n.plus ⋅ (n.plus⋅∇∇(u).plus) - n.minus ⋅ (n.minus ⋅ ∇∇(u).minus) )       # andre ordens hopp... Forklare dette skikkelig. 
end

function stokes_solver(;n, u_exact, p_exact, f, g, ud, order, geometry, βu0, γu1, γu2, γp, βp0, nu, stabilize, δ, save = false, calc_condition = false)
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
    partition = (n, n)
    dim = length(partition)
    a = 1.2
    pmin = Point(-a + δ, -a + δ)
    pmax = Point(a + δ, a + δ)
    bgmodel = CartesianDiscreteModel(pmin,pmax,partition)
    # mesh size
    h = (pmax-pmin)[1]/partition[1]

    # defining ghost penalty constants
    βu = βu0 *nu/(h^2)
    βp = βp0/h

    geo = create_geometry(geometry, n)
    # Define active and physical mesh
    cutgeo = cut(bgmodel,geo)
    cutgeo_facets = cut_facets(bgmodel,geo)
    Ω_bg = Triangulation(bgmodel)
    Ω_act = Triangulation(cutgeo, ACTIVE)
    Ω = Triangulation(cutgeo, PHYSICAL)

    # Embedded boundary
    # Dirichlet conditions on u
    Γd = EmbeddedBoundary(cutgeo)
    n_Γd = get_normal_vector(Γd)

    # Get ghost penalty facets
    Fg = GhostSkeleton(cutgeo)
    n_Fg = get_normal_vector(Fg)

    # Define measures
    degree = 2*order
    dΩ = Measure(Ω,degree)
    dΓd = Measure(Γd, degree)
    dFg = Measure(Fg, degree)

    # Define function spaces 
    reffe_u  = ReferenceFE(lagrangian,VectorValue{dim, Float64},order)
    reffe_p = ReferenceFE(lagrangian,Float64, order - 1)

    V = TestFESpace(Ω_act, reffe_u,  conformity=:H1)
    Q = TestFESpace(Ω_act, reffe_p, conformity=:H1, constraint=:zeromean)

    U = TrialFESpace(V)
    P = TrialFESpace(Q)

    X = MultiFieldFESpace([U, P])
    Y = MultiFieldFESpace([V, Q])

    γ = 10*2*2
    u_D = ud
  | # Weak formulation components
    a(u, v) = (∫(∇(u) ⊙ ∇(v))dΩ + ∫(-((n_Γd ⋅ ∇(u)) ⋅ v) + (-(n_Γd ⋅ ∇(v)) ⋅ u)+(γ/h * u ⋅ v))dΓd )
    b(v, p) = (∫(-1*(∇ ⋅ v*p))dΩ
                + ∫((n_Γd ⋅ v) * p)dΓd)#
   

    gu(u,v) = ( ∫( (β_1*h)*jump(n_Fg⋅∇(u))⋅jump(n_Fg⋅∇(v)) )dFg 
            +  
               ∫( (β_2*h^3)*jump_nn(u,n_Fg)⋅jump_nn(v,n_Fg) )dFg)


    gp(p, q) = (∫((β_3*h^3)*jump(n_Fg ⋅ ∇(p)) * jump(n_Fg ⋅ ∇(q)))dFg)

    l2_norm(u) = (sum( ∫( u ⋅ u )*dΩ ))
    h1_semi(u) = sum(∫(∇(u) ⊙ ∇(u))*dΩ)
    
    l1(v) = ∫(f ⋅ v)dΩ
    l2(v) = ∫(-1* (n_Γd ⋅ ∇(v)) ⋅ ud + γ/h ⋅ v ⋅ ud )dΓd   #∫)dΓd # + γ/h ⋅ ud ⋅ v)dΓd
    l3(q) = ∫(n_Γd ⋅ ud *q)dΓd

    if stabilize
        A((u,p),(v,q)) =(a(u, v) + b(v, p) + b(u, q) 
        + gu(u,v)
        - gp(p, q)
        )

        L((v, q)) = l1(v) + l2(v) + l3(q)
    
        op = AffineFEOperator(A,L,X,Y)
        uh, ph = solve(op)

    else
        B((u,p),(v,q)) = a(u,v) + b(v, p) + b(u, q) 
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
    # writevtk(Ω_act, "mesh_act_$geometry")
    # writevtk(Ω_bg, "mesh_bg_$geometry")
    # writevtk(Ω, "mesh_$geometry")
    # writevtk(Γd, "surface_gamma_d_$geometry")
    # writevtk(Γs , "outer_$geometry")
    # writevtk(Fg, "ghost_facets_$geometry")
    println("lagret")
    writevtk(Ω,joinpath(outputfolder," $n $geometry $order"), cellfields=["u_ex" => u_exact, "uh"=>uh, "erru"=> erru, "p_ex" => p_exact, "ph"=>ph, "errp"=> errp, "nablau" => ∇(u_exact)]) #, "erru" => erru])
  end
return uh, u_exact, erru, l2_norm(uh - u_exact), h1_semi(uh - u_exact), ph, p_exact, errp, l2_norm(ph - p_exact), h1_semi(ph - p_exact), condition_numb, Ω_act
end

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
#uh, u_exact, erru, l2_u, h1_semi_u, ph, p_exact, errp, l2_p, h1_semi_p, condition_numb, Ω_act  = stokes_solver(;n, u_exact, p_exact, f, g, ud, order, geometry, βu0, γu1, γu2, γp, βp0, nu, stabilize, δ, save, calc_condition)

function convergence(;numb_it, u_exact, p_exact, f, g, ud, order, geometry, solver, δ, βu0, γu1, γu2, γp, βp0, nu, stabilize, save = false)
  #"""function to calculate convergence of the poisson solver, or the stokes solver, with or without stabilization"""
  calc_condition = false 
  uarr_l2 = zeros(Float64, numb_it)
  uarr_h1 = zeros(Float64, numb_it)
  parr_l2 = zeros(Float64, numb_it)
  parr_h1 = zeros(Float64, numb_it)

  n_arr = 2 .^ (2:(numb_it + 1))
  h = 1.0 ./ n_arr

  for i = 1:numb_it
    n = n_arr[i]
    elapsed_time, solver_result = let
        t, val = @timed solver(;n, u_exact, p_exact, f, g, ud, order, geometry, βu0, γu1, γu2, γp, βp0, nu, stabilize, δ, save, calc_condition)
        (val, t) 
    end
    uh, u_exact, erru, l2_u, h1_semi_u, ph, p_exact, errp, l2_p, h1_semi_p, condition_numb, Ω_act = solver_result

    uarr_l2[i] = l2_u  #l2 error
    uarr_h1[i] = h1_semi_u  #h1 error
    parr_l2[i] = l2_p  #l2 error
    parr_h1[i] = h1_semi_p  #h1 error

      println("$i: Solved system in $elapsed_time seconds.")
  end
  #println(length(uarr_h1))
  #println(length(h))
  #EOC_l2 = log.(uarr_l2[1:end-1] ./ uarr_l2[2:end]) ./ log.(h[1:end-1] ./ h[2:end])
  #EOC_h1 = log.(uarr_h1[1:end-1] ./ uarr_h1[2:end]) ./ log.(h[1:end-1] ./ h[2:end])

  return uarr_l2, uarr_h1, parr_l2, parr_l2, h#, EOC_l2, EOC_h1
end

function sensitivity(;n, M, u_exact, p_exact, f, g, ud, order, geometry, solver, βu0, γu1, γu2, γp, βp0, nu, stabilize, save)
    calc_condition = true
    arr_δ = zeros(Float64, M-1)
    arr_l2u = zeros(Float64, M-1)
    arr_h1u = zeros(Float64, M-1)
    arr_l2p = zeros(Float64, M-1)
    arr_h1p = zeros(Float64, M-1)
    arr_cond = zeros(Float64, M-1)
    #loop to perturb the active geometry
    for i = 1:(M-1)
        δ = i/n/M *1.1/ sqrt(2)
        elapsed_time, solver_result = let
            t, val = @timed solver(;n, u_exact, p_exact, f, g, ud, order, geometry, βu0, γu1, γu2, γp, βp0, nu, stabilize, δ, save, calc_condition)
            (val, t)
        end
        uh, u_exact, erru, l2_u, h1_semi_u, ph, p_exact, errp, l2_p, h1_semi_p, condition_numb, Ω_act = solver_result
        save = false
        arr_δ[i] = δ
        arr_l2u[i] = l2_u
        arr_h1u[i] = h1_semi_u
        arr_l2p[i] = l2_p
        arr_h1p[i] = h1_semi_p
        arr_cond[i] = condition_numb

        if i % 100 == 0
            println("$i: Solved system in $elapsed_time seconds.")
            save = true     #lagrer løsningen hver 100nde gang
        end
    end
    return arr_δ, arr_l2u, arr_h1u, arr_l2p, arr_h1p, arr_cond
end

####### Convergence test #######
#  stabilisering, order 1           # trenger ikke regne ut kondisjonstall når man gjør konvergestest...
# numb_it = 6                         # Sånn koden er implementert nå så er det fra 2^1 til 2^{numb_it}. Tar kort tid å kjøre for numb_it = 6
# order = 2                           # når jeg øker orden så øker kjøretid veeeeldig !! Bør vurdere å skru ned numb_it samtidig. 244 sekunder når jeg har på order = 2 for kjøringen 2^6
# δ = 0                               # kan også se ut til at feilen havner på maskinnivå? vet ikke helt, men mulig å eksperimentere med dette. 
# # med stabilisering, order 1
# stabilize = true
# solver = stokes_solver
# geometry = "circle"
# γu1 = 0.1
# γu2 = 0.1
# γp = 0.1
# save = false
# # fjernet disse: , eoc_l2, eoc_h1
# uarr_l2_1_stab, uarr_h1_1_stab, parr_l2_1_stab, parr_h1_1_stab, h = convergence(;numb_it, u_exact, p_exact, f, g, ud, order, geometry, solver, δ, βu0, γu1, γu2, γp, βp0, nu, stabilize, save)
# start = 1
# plot(h[start:end], uarr_l2_1_stab[start:end], xaxis=:log, yaxis=:log, marker=:o, lw=2, label="L2 error")
# plot!(h[start:end], uarr_h1_1_stab[start:end], marker=:s, lw=2, label="H1 error")

# # # Uten stabilisering, order 1
# stabilize = false
# uarr_l2_1_nostab, uarr_h1_1_nostab, parr_l2_1_nostab, parr_h1_1_nostab, h = convergence(;numb_it, u_exact, p_exact, f, g, ud, order, geometry, solver, δ, βu0, γu1, γu2, γp, βp0, nu, stabilize, save)

# plot!(h[start:end], uarr_l2_1_nostab[start:end], marker=:s, lw=2, label="L2 error, order $order (no Stab)")
# plot!(h[start:end], uarr_h1_1_nostab[start:end], marker=:s, lw=2, label="H1 error, order $order (no Stab)")

# # # Legger til aksetitler og tittel
# xlabel!("Mesh size h")
# ylabel!("Error")
# title!("Convergence of Stokes Solver")

####### Sensitivity test #######
# kjører nå denne med n = 16, men bør nok kjøre for n = 32 eller n = 64, men det kan ta laaaag tid. 30 min per kjøring ved n = 64 -
n = 16           # øke denne
M = 2000        #full kjøring med M = 2000
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

arr_δ, arr_l2u, arr_h1u, arr_l2p, arr_h1p, arr_cond = sensitivity(;n, M, u_exact, p_exact, f, g, ud, order, geometry, solver, βu0, γu1, γu2, γp, βp0, nu, stabilize, save)
stabilize = false
start = 1
arr_δ_nostab, arr_l2u_nostab, arr_h1u_nostab, arr_l2p_nostab, arr_h1p_nostab, arr_cond_nostab = sensitivity(;n, M, u_exact, p_exact, f, g, ud, order, geometry, solver, βu0, γu1, γu2, γp, βp0, nu, stabilize, save)
plot(arr_δ[start:end], arr_l2u[start:end],  yaxis=:log, lw=2, label="L2 error, stabilized")
plot!(arr_δ_nostab[start:end], arr_l2u_nostab[start:end], yaxis=:log, lw=2, label="L2 error, non-stabilized")
plot!(arr_δ[start:end], arr_h1u[start:end], yaxis=:log, lw=2, label="h1 stabilized")
plot!(arr_δ_nostab[start:end], arr_h1u_nostab[start:end], yaxis=:log, lw=2, label="h1 non-stabilized")
xlabel!("Perturbation δ")
ylabel!("Condition number")
title!("Sensitivity analysis of stokes solver, solution u")

# plot(arr_δ[start:end], arr_cond[start:end], yaxis=:log, label = "Stabilized")
# plot!(arr_δ[start:end], arr_cond_nostab[start:end],yaxis=:log, label = "Not stabilized")

# # # if I want the same plot for p
# plot(arr_δ[start:end], arr_l2p[start:end], xaxis=:log, yaxis=:log, lw=2, label="L2 error, stabilized")
# plot!(arr_δ_nostab[start:end], arr_l2p_nostab[start:end], xaxis=:log, yaxis=:log, lw=2, label="L2 error, non-stabilized")
# plot!(arr_δ[start:end], arr_h1p[start:end], xaxis=:log, yaxis=:log, lw=2, label="h1 stabilized")
# plot!(arr_δ_nostab[start:end], arr_h1p_nostab[start:end], xaxis=:log, yaxis=:log, lw=2, label="h1 non-stabilized")
# xlabel!("Perturbation δ")
# ylabel!("Condition number")
# title!("Sensitivity analysis of stokes solver, solution p")

#### Varying geometry
# geometry_arr = ["circle", "flower", "heart"]
# s1 = plot()  # Lager et tomt plott
# s2 = plot()
# s3 = plot()
# start = 1
# for i = 1:3
#     geometry = geometry_arr[i]
#     stabilize = true
#     arr_δ, arr_l2u, arr_h1u, arr_l2p, arr_h1p, arr_cond = sensitivity(;n, M, u_exact, p_exact, f, g, ud, order, geometry, solver, βu0, γu1, γu2, γp, βp0, nu, stabilize, save)
#     stabilize = false
#     arr_δ, arr_l2u_nostab, arr_h1u_nostab, arr_l2p_nostab, arr_h1p_nostab, arr_cond_nostab = sensitivity(;n, M, u_exact, p_exact, f, g, ud, order, geometry, solver, βu0, γu1, γu2, γp, βp0, nu, stabilize, save)
#     plot!(s1, arr_δ[start:end], arr_cond[start:end], xaxis=:log, yaxis=:log, lw=2, label="Stabilized ($geometry)")
#     plot!(s1, arr_δ[start:end], arr_cond_nostab[start:end], xaxis=:log, yaxis=:log, lw=2, label="Non-stabilized ($geometry)")
#     plot!(s2, arr_δ[start:end], arr_l2u[start:end], xaxis=:log, yaxis=:log, lw=2, label="Stabilized ($geometry)")
#     plot!(s2, arr_δ[start:end], arr_l2u_nostab[start:end], xaxis=:log, yaxis=:log, lw=2, label="Non-stabilized ($geometry)")
#     plot!(s3, arr_δ[start:end], arr_h1u[start:end], xaxis=:log, yaxis=:log, lw=2, label="Stabilized ($geometry)")
#     plot!(s3, arr_δ[start:end], arr_h1u_nostab[start:end], xaxis=:log, yaxis=:log, lw=2, label="Non-stabilized ($geometry)")
# end
# title!(s1, "Sensitivity analysis varing geometries")
# xlabel!(s1, "Perturbation δ")
# ylabel!(s1, "Condition number")
# title!(s2, "Sensitivity analysis varing geometries")
# xlabel!(s2, "Perturbation δ")
# ylabel!(s2, "L2 norm")
# title!(s3, "Sensitivity analysis varing geometries")
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
#     arr_δ, arr_l2u, arr_h1u, arr_l2p, arr_h1p, arr_cond = sensitivity(;n, M, u_exact, p_exact, f, g, ud, order, geometry, solver, βu0, γu1, γu2, γp, βp0, nu, stabilize, save)
#     stabilize = false
#     arr_δ_nostab, arr_l2u_nostab, arr_h1u_nostab, arr_l2p_nostab, arr_h1p_nostab, arr_cond_nostab = sensitivity(;n, M, u_exact, p_exact, f, g, ud, order, geometry, solver, βu0, γu1, γu2, γp, βp0, nu, stabilize, save)
#     plot!(p1, arr_δ[start:end], arr_cond[start:end], xaxis=:log, yaxis=:log, lw=2, label="Stabilized γd = $γd")
#     plot!(p1, arr_δ[start:end], arr_cond_nostab[start:end], xaxis=:log, yaxis=:log, lw=2, label="Non-stabilized γd = $γd")
#     plot!(p2, arr_δ[start:end], arr_l2u[start:end], xaxis=:log, yaxis=:log, lw=2, label="Stabilized γd = $γd")
#     plot!(p2, arr_δ[start:end], arr_l2u_nostab[start:end], xaxis=:log, yaxis=:log, lw=2, label="Non-stabilized γd = $γd")
#     plot!(p3, arr_δ[start:end], arr_h1u[start:end], xaxis=:log, yaxis=:log, lw=2, label="Stabilized γd = $γd")
#     plot!(p3, arr_δ[start:end], arr_h1u_nostab[start:end], xaxis=:log, yaxis=:log, lw=2, label="Non-stabilized γd = $γd")
# end
# title!(p1, "Sensitivity analysis varying γu1")
# xlabel!(p1, "Perturbation δ")
# ylabel!(p1, "Condition number")
# title!(p2, "Sensitivity analysis varying γd")
# xlabel!(p2, "Perturbation δ")
# ylabel!(p2, "L2 norm (velocity)")
# title!(p3, "Sensitivity analysis varying γd")
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
#     arr_δ, arr_l2u, arr_h1u, arr_l2p, arr_h1p, arr_cond = sensitivity(;n, M, u_exact, p_exact, f, g, ud, order, geometry, solver, βu0, γu1, γu2, γp, βp0, nu, stabilize, save)
#     stabilize = false
#     arr_δ_nostab, arr_l2u_nostab, arr_h1u_nostab, arr_l2p_nostab, arr_h1p_nostab, arr_cond_nostab = sensitivity(;n, M, u_exact, p_exact, f, g, ud, order, geometry, solver, βu0, γu1, γu2, γp, βp0, nu, stabilize, save)
#     plot!(q1, arr_δ[start:end], arr_cond[start:end], xaxis=:log, yaxis=:log, lw=2, label="Stabilized γg1 = $γg1")
#     plot!(q1, arr_δ[start:end], arr_cond_nostab[start:end], xaxis=:log, yaxis=:log, lw=2, label="Non-stabilized γg1 = $γg1")
#     plot!(q2, arr_δ[start:end], arr_l2u[start:end], xaxis=:log, yaxis=:log, lw=2, label="Stabilized γg1 = $γg1")
#     plot!(q2, arr_δ[start:end], arr_l2u_nostab[start:end], xaxis=:log, yaxis=:log, lw=2, label="Non-stabilized γg1 = $γg1")
#     plot!(q3, arr_δ[start:end], arr_h1u[start:end], xaxis=:log, yaxis=:log, lw=2, label="Stabilized γg1 = $γg1")
#     plot!(q3, arr_δ[start:end], arr_h1u_nostab[start:end], xaxis=:log, yaxis=:log, lw=2, label="Non-stabilized γg1 = $γg1")
# end
# title!(q1, "Sensitivity analysis varying γg1")
# xlabel!(q1, "Perturbation δ")
# ylabel!(q1, "Condition number")
# title!(q2, "Sensitivity analysis varying γg1")
# xlabel!(q2, "Perturbation δ")
# ylabel!(q2, "L2 norm (velocity)")
# title!(q3, "Sensitivity analysis varying γg1")
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
#     arr_δ, arr_l2u, arr_h1u, arr_l2p, arr_h1p, arr_cond = sensitivity(;n, M, u_exact, p_exact, f, g, ud, order, geometry, solver, βu0, γu1, γu2, γp, βp0, nu, stabilize, save)
#     stabilize = false
#     arr_δ_nostab, arr_l2u_nostab, arr_h1u_nostab, arr_l2p_nostab, arr_h1p_nostab, arr_cond_nostab = sensitivity(;n, M, u_exact, p_exact, f, g, ud, order, geometry, solver, βu0, γu1, γu2, γp, βp0, nu, stabilize, save)
#     plot!(r1, arr_δ[start:end], arr_cond[start:end], xaxis=:log, yaxis=:log, lw=2, label="Stabilized γg3 = $γg3")
#     plot!(r1, arr_δ_nostab[start:end], arr_cond_nostab[start:end], xaxis=:log, yaxis=:log, lw=2, label="Non-stabilized γg3 = $γg3")
#     plot!(r2, arr_δ[start:end], arr_l2u[start:end], xaxis=:log, yaxis=:log, lw=2, label="Stabilized γg3 = $γg3")
#     plot!(r2, arr_δ_nostab[start:end], arr_l2u_nostab[start:end], xaxis=:log, yaxis=:log, lw=2, label="Non-stabilized γg3 = $γg3")
#     plot!(r3, arr_δ[start:end], arr_h1u[start:end], xaxis=:log, yaxis=:log, lw=2, label="Stabilized γg3 = $γg3")
#     plot!(r3, arr_δ_nostab[start:end], arr_h1u_nostab[start:end], xaxis=:log, yaxis=:log, lw=2, label="Non-stabilized γg3 = $γg3")
# end
# title!(r1, "Sensitivity analysis varying γg3")
# xlabel!(r1, "Perturbation δ")
# ylabel!(r1, "Condition number")
# title!(r2, "Sensitivity analysis varying γg3")
# xlabel!(r2, "Perturbation δ")
# ylabel!(r2, "L2 norm (velocity)")
# title!(r3, "Sensitivity analysis varying γg3")
# xlabel!(r3, "Perturbation δ")
# ylabel!(r3, "H1 norm (velocity)")
# display(r1)
# display(r2)
# display(r3)


# start = 1

# ###########################
# ### Sensitivity: βu
# ###########################
# βu_arr = [0.01, 0.1, 1.0]
# p1 = plot(); p2 = plot(); p3 = plot()
# for βu0 in βu_arr
#     stabilize = true
#     arr_δ, arr_l2u, arr_h1u, _, _, arr_cond = sensitivity(;n, M, u_exact, p_exact, f, g, ud, order, geometry, solver, βu0, γu1, γu2, γp, βp0, nu, stabilize, save)
#     stabilize = false
#     _, arr_l2u_ns, arr_h1u_ns, _, _, arr_cond_ns = sensitivity(;n, M, u_exact, p_exact, f, g, ud, order, geometry, solver, βu0, γu1, γu2, γp, βp0, nu, stabilize, save)

#     plot!(p1, arr_δ[start:end], arr_cond[start:end], xaxis=:log, yaxis=:log, label="Stabilized βu=$βu")
#     plot!(p1, arr_δ[start:end], arr_cond_ns[start:end], xaxis=:log, yaxis=:log, label="Non-stabilized βu=$βu")

#     plot!(p2, arr_δ[start:end], arr_l2u[start:end], xaxis=:log, yaxis=:log, label="Stabilized βu=$βu")
#     plot!(p2, arr_δ[start:end], arr_l2u_ns[start:end], xaxis=:log, yaxis=:log, label="Non-stabilized βu=$βu")

#     plot!(p3, arr_δ[start:end], arr_h1u[start:end], xaxis=:log, yaxis=:log, label="Stabilized βu=$βu")
#     plot!(p3, arr_δ[start:end], arr_h1u_ns[start:end], xaxis=:log, yaxis=:log, label="Non-stabilized βu=$βu")
# end
# title!(p1, "Sensitivity: varying βu"); xlabel!(p1, "δ"); ylabel!(p1, "Condition number")
# title!(p2, "L2 error (velocity), varying βu"); xlabel!(p2, "δ"); ylabel!(p2, "L2 norm")
# title!(p3, "H1 error (velocity), varying βu"); xlabel!(p3, "δ"); ylabel!(p3, "H1 norm")
# display(p1); 
# display(p2); 
# display(p3)


# ###########################
# ### Sensitivity: γu1
# ###########################
# γu1_arr = [0.1, 1.0, 10.0]
# q1 = plot(); q2 = plot(); q3 = plot()
# for γu1 in γu1_arr
#     stabilize = true
#     arr_δ, arr_l2u, arr_h1u, _, _, arr_cond = sensitivity(;n, M, u_exact, p_exact, f, g, ud, order, geometry, solver, βu0, γu1, γu2, γp, βp0, nu, stabilize, save)
#     stabilize = false
#     _, arr_l2u_ns, arr_h1u_ns, _, _, arr_cond_ns = sensitivity(;n, M, u_exact, p_exact, f, g, ud, order, geometry, solver, βu0, γu1, γu2, γp, βp0, nu, stabilize, save)

#     plot!(q1, arr_δ[start:end], arr_cond[start:end], xaxis=:log, yaxis=:log, label="Stabilized γu1=$γu1")
#     plot!(q1, arr_δ[start:end], arr_cond_ns[start:end], xaxis=:log, yaxis=:log, label="Non-stabilized γu1=$γu1")

#     plot!(q2, arr_δ[start:end], arr_l2u[start:end], xaxis=:log, yaxis=:log, label="Stabilized γu1=$γu1")
#     plot!(q2, arr_δ[start:end], arr_l2u_ns[start:end], xaxis=:log, yaxis=:log, label="Non-stabilized γu1=$γu1")

#     plot!(q3, arr_δ[start:end], arr_h1u[start:end], xaxis=:log, yaxis=:log, label="Stabilized γu1=$γu1")
#     plot!(q3, arr_δ[start:end], arr_h1u_ns[start:end], xaxis=:log, yaxis=:log, label="Non-stabilized γu1=$γu1")
# end
# title!(q1, "Sensitivity: varying γu1"); xlabel!(q1, "δ"); ylabel!(q1, "Condition number")
# title!(q2, "L2 error (velocity), varying γu1"); xlabel!(q2, "δ"); ylabel!(q2, "L2 norm")
# title!(q3, "H1 error (velocity), varying γu1"); xlabel!(q3, "δ"); ylabel!(q3, "H1 norm")
# display(q1); 
# display(q2); 
# display(q3)


# ###########################
# ### Sensitivity: γu2
# ###########################
# γu2_arr = [0.01, 0.1, 1.0]
# r1 = plot(); r2 = plot(); r3 = plot()
# for γu2 in γu2_arr
#     stabilize = true
#     arr_δ, arr_l2u, arr_h1u, _, _, arr_cond = sensitivity(;n, M, u_exact, p_exact, f, g, ud, order, geometry, solver, βu0, γu1, γu2, γp, βp0, nu, stabilize, save)
#     _, arr_l2u_ns, arr_h1u_ns, _, _, arr_cond_ns = sensitivity(;n, M, u_exact, p_exact, f, g, ud, order, geometry, solver, βu0, γu1, γu2, γp, βp0, nu, stabilize, save)

#     plot!(r1, arr_δ[start:end], arr_cond[start:end], xaxis=:log, yaxis=:log, label="Stabilized γu2=$γu2")
#     plot!(r1, arr_δ[start:end], arr_cond_ns[start:end], xaxis=:log, yaxis=:log, label="Non-stabilized γu2=$γu2")

#     plot!(r2, arr_δ[start:end], arr_l2u[start:end], xaxis=:log, yaxis=:log, label="Stabilized γu2=$γu2")
#     plot!(r2, arr_δ[start:end], arr_l2u_ns[start:end], xaxis=:log, yaxis=:log, label="Non-stabilized γu2=$γu2")

#     plot!(r3, arr_δ[start:end], arr_h1u[start:end], xaxis=:log, yaxis=:log, label="Stabilized γu2=$γu2")
#     plot!(r3, arr_δ[start:end], arr_h1u_ns[start:end], xaxis=:log, yaxis=:log, label="Non-stabilized γu2=$γu2")
# end
# title!(r1, "Sensitivity: varying γu2"); xlabel!(r1, "δ"); ylabel!(r1, "Condition number")
# title!(r2, "L2 error (velocity), varying γu2"); xlabel!(r2, "δ"); ylabel!(r2, "L2 norm")
# title!(r3, "H1 error (velocity), varying γu2"); xlabel!(r3, "δ"); ylabel!(r3, "H1 norm")
# display(r1); 
# display(r2); 
# display(r3)


# ###########################
# ### Sensitivity: γp
# ###########################
# γp_arr = [0.01, 0.1, 1.0]
# s1 = plot(); s2 = plot(); s3 = plot()
# for γp in γp_arr
#     stabilize = true
#     arr_δ, arr_l2u, arr_h1u, _, _, arr_cond = sensitivity(;n, M, u_exact, p_exact, f, g, ud, order, geometry, solver, βu0, γu1, γu2, γp, βp0, nu, stabilize, save)
#     stabilize = false
#     _, arr_l2u_ns, arr_h1u_ns, _, _, arr_cond_ns = sensitivity(;n, M, u_exact, p_exact, f, g, ud, order, geometry, solver, βu0, γu1, γu2, γp, βp0, nu, stabilize, save)

#     plot!(s1, arr_δ[start:end], arr_cond[start:end], xaxis=:log, yaxis=:log, label="Stabilized γp=$γp")
#     plot!(s1, arr_δ[start:end], arr_cond_ns[start:end], xaxis=:log, yaxis=:log, label="Non-stabilized γp=$γp")

#     plot!(s2, arr_δ[start:end], arr_l2u[start:end], xaxis=:log, yaxis=:log, label="Stabilized γp=$γp")
#     plot!(s2, arr_δ[start:end], arr_l2u_ns[start:end], xaxis=:log, yaxis=:log, label="Non-stabilized γp=$γp")

#     plot!(s3, arr_δ[start:end], arr_h1u[start:end], xaxis=:log, yaxis=:log, label="Stabilized γp=$γp")
#     plot!(s3, arr_δ[start:end], arr_h1u_ns[start:end], xaxis=:log, yaxis=:log, label="Non-stabilized γp=$γp")
# end
# title!(s1, "Sensitivity: varying γp"); xlabel!(s1, "δ"); ylabel!(s1, "Condition number")
# title!(s2, "L2 error (velocity), varying γp"); xlabel!(s2, "δ"); ylabel!(s2, "L2 norm")
# title!(s3, "H1 error (velocity), varying γp"); xlabel!(s3, "δ"); ylabel!(s3, "H1 norm")
# display(s1); 
# display(s2); 
# display(s3)


# ###########################
# ### Sensitivity: βp
# ###########################
# βp_arr = [0.001, 0.01, 0.1]
# t1 = plot(); t2 = plot(); t3 = plot()
# for βp0 in βp_arr
#     stabilize = true
#     arr_δ, arr_l2u, arr_h1u, _, _, arr_cond = sensitivity(;n, M, u_exact, p_exact, f, g, ud, order, geometry, solver, βu0, γu1, γu2, γp, βp0, nu, stabilize, save)
#     stabilize = false
#     _, arr_l2u_ns, arr_h1u_ns, _, _, arr_cond_ns = sensitivity(;n, M, u_exact, p_exact, f, g, ud, order, geometry, solver, βu0, γu1, γu2, γp, βp0, nu, stabilize, save)

#     plot!(t1, arr_δ[start:end], arr_cond[start:end], xaxis=:log, yaxis=:log, label="Stabilized βp=$βp")
#     plot!(t1, arr_δ[start:end], arr_cond_ns[start:end], xaxis=:log, yaxis=:log, label="Non-stabilized βp=$βp")

#     plot!(t2, arr_δ[start:end], arr_l2u[start:end], xaxis=:log, yaxis=:log, label="Stabilized βp=$βp")
#     plot!(t2, arr_δ[start:end], arr_l2u_ns[start:end], xaxis=:log, yaxis=:log, label="Non-stabilized βp=$βp")

#     plot!(t3, arr_δ[start:end], arr_h1u[start:end], xaxis=:log, yaxis=:log, label="Stabilized βp=$βp")
#     plot!(t3, arr_δ[start:end], arr_h1u_ns[start:end], xaxis=:log, yaxis=:log, label="Non-stabilized βp=$βp")
# end
# title!(t1, "Sensitivity: varying βp"); xlabel!(t1, "δ"); ylabel!(t1, "Condition number")
# title!(t2, "L2 error (velocity), varying βp"); xlabel!(t2, "δ"); ylabel!(t2, "L2 norm")
# title!(t3, "H1 error (velocity), varying βp"); xlabel!(t3, "δ"); ylabel!(t3, "H1 norm")
# display(t1); 
# display(t2); 
# display(t3)
