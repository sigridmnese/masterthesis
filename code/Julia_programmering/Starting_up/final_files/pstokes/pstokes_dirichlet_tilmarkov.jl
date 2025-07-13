using Gridap
using GridapEmbedded
using STLCutters
using LinearAlgebra
using Plots
using LineSearches: BackTracking
import Random
using Logging
using LoggingExtras
using PGFPlotsX
using DataFrames
using CSV

"""Fitted FEM with weak boundary imposition using Nitsches method. Symmetric p-stokes equations, using Newtons method to linearize
denne funker for nu0 = 100, 1 og 0.01!!"""
#""" Lagt inn konvergensplott for denne i teams"""

function jump_nn(u,n)
    return ( n.plus ⋅ (n.plus⋅∇∇(u).plus) - n.minus ⋅ (n.minus ⋅ ∇∇(u).minus) )       # andre ordens hopp... Forklare dette skikkelig. 
end


function create_geometry(name, n)
    try
        if lowercase(name) == "circle"
            R  = 1.0
            geo = AnalyticalGeometry(x->(x[1])^2+(x[2]^2)-R^2)
            return geo

        elseif lowercase(name) == "circle2"
            R = 2
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
            geo = AnalyticalGeometry(x -> (-R     +   (x[1])^2   +     ( 5/4 * (x[2]+0.2)   -     sqrt(abs(x[1]))  )^2   ))
            #geo = AnalyticalGeometry(x -> (x[1])^2 + (5/4 * (x[2]) - √(abs(x[1])))^2 - R)
            return geo
        elseif lowercase(name) == "twoheartsline"
            # Vi lager to hjerter med sentre horisontalt separert med 1.0, 
            # begge med samme "R" og "y‐shift" som over, 
            # og under dem en horisontal linje 3.0 enheter under hjertenes bunnspiss.
            R      = 0.7
            yshift = 0.2    # samme y‐forskyvning som i enkelthjerter‐funksjonen
            d      = 0.5    # halve avstanden mellom hjertesentrene i x‐retning
            
            # Implicit‐funksjon for venstre hjerte (senter i x = -0.5):
            heart1(x) = -R
                        + (x[1] + d)^2
                        + ( (5/4)*(x[2] + yshift) - sqrt(abs(x[1] + d)) )^2

            # Implicit‐funksjon for høyre hjerte (senter i x = +0.5):
            heart2(x) = -R
                        + (x[1] - d)^2
                        + ( (5/4)*(x[2] + yshift) - sqrt(abs(x[1] - d)) )^2

            # Beregn hvor hjerte‐bunnspissen ligger i y‐koordinat (med R=0.7 og yshift=0.2
            # blir bunnspissen y_bunn = -0.2).
            y_bunn = -yshift

            # Horisontal linje 3.0 enheter under hjertenes bunnspiss:
            yline = y_bunn - 3.0   # = -0.2 - 3.0 = -3.2

            # Nå ønsker vi at et punkt x=(x[1],x[2]) skal være "inni" geometrien 
            # dersom det ligger *i minst ett* av hjertene OG *over* linjen y = yline.
            # Union av hjertene:   f_hearts(x) = min( heart1(x), heart2(x) )
            # Deretter klippe av mot linjen:   f_total(x) = max( f_hearts(x), yline - x[2] )
            # Et punkt er inne akkurat når f_total(x) <= 0.

            return AnalyticalGeometry(x -> heart1(x) + heart2(x) + yline - x[2])
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


function convergence_stokes(;numb_it, u_exact, p_exact, f, g, ud, order, geometry, solver, δ, βu0, γu1, γu2, γp, βp0, nu, stabilize, save = false)
  #"""function to calculate convergence of the poisson solver, or the stokes solver, with or without stabilization"""
  calc_condition = false 
  uarr_l2 = zeros(Float64, numb_it)
  uarr_h1 = zeros(Float64, numb_it)
  parr_l2 = zeros(Float64, numb_it)
  parr_h1 = zeros(Float64, numb_it)

  n_arr = 2 .^ (3:(numb_it + 2))
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

  return uarr_l2, uarr_h1, parr_l2, parr_h1, h#, EOC_l2, EOC_h1
end



function sensitivity_stokes(;n, M, u_exact, p_exact, f, g, ud, order, geometry, solver, βu0, γu1, γu2, γp, βp0, nu, stabilize, save)
    calc_condition = false
    arr_δ = zeros(Float64, M-1)
    arr_l2u = zeros(Float64, M-1)
    arr_h1u = zeros(Float64, M-1)
    arr_l2p = zeros(Float64, M-1)
    arr_h1p = zeros(Float64, M-1)
    arr_cond = zeros(Float64, M-1)
    #loop to perturb the active geometry
    for i = 1:(M-1)
        δ = i/n/M *1.1/ sqrt(2)
        
        #elapsed_time, solver_result = let
        #    t, val = @timed solver(;n, u_exact, p_exact, f, g, ud, order, geometry, βu0, γu1, γu2, γp, βp0, nu, stabilize, δ, save, calc_condition)
        #    (val, t)
        #end
        #uh, u_exact, erru, l2_u, h1_semi_u, ph, p_exact, errp, l2_p, h1_semi_p, condition_numb, Ω_act = solver_result
        uh, u_exact, erru, l2_u, h1_semi_u, ph, p_exact, errp, l2_p, h1_semi_p, condition_numb, Ω_act = solver(;n, u_exact, p_exact, f, g, ud, order, geometry, βu0, γu1, γu2, γp, βp0, nu, stabilize, δ, save, calc_condition)
        save = false
        arr_δ[i] = δ
        arr_l2u[i] = l2_u
        arr_h1u[i] = h1_semi_u
        arr_l2p[i] = l2_p
        arr_h1p[i] = h1_semi_p
        arr_cond[i] = condition_numb

        if i % 100 == 0
            println("$i") #: Solved system in $elapsed_time seconds.")
            save = true     #lagrer løsningen hver 100nde gang
        end
    end
    return arr_δ, arr_l2u, arr_h1u, arr_l2p, arr_h1p, arr_cond
end




  # Defining constants
   # test ulike verdier for nu0!!! Da synes jeg du skal kjøre konvergenstester men også bare generelt se på hva som skjer med løsningen når du endrer nu0.



r = 4/3
A = 2^(1/(1-r))
ϵ_0 = 1e-6
nu0 =  A

# Defining manufactured solutions
#test case 1
u_exact(x) =  VectorValue(2*x[1] + cos(2*π*x[2])*10, -2*x[2] + sin(2*π*x[1]))#VectorValue(2*x[1] + exp(x[1]/2) * cos(2*π*x[2]), -2*x[2] + exp(x[2]/2) * sin(2*π*x[1])) #bytte til sin/cos-uttrykk  VectorValue(-x[2], x[1])
p_exact(x) = sin(2*π*x[1])*cos(2*π*x[2])

#test case 2
# u_exact(x) = VectorValue( sin(π*x[1]) * cos(π*x[2]), -cos(π*x[1]) * sin(π*x[2]) )
# p_exact(x) = cos(2π*x[1]) * sin(2π*x[2])

# test case 3
# u_exact(x) = VectorValue( exp(-100*((x[1]-0.5)^2 + (x[2]-0.5)^2)), 0.0 )
# p_exact(x) = x[1] * (1 - x[1]) * x[2] * (1 - x[2])


# Defining the problem 
flux(εu) = nu0^(1-r)*(ϵ_0^2 + tr(εu' ⋅ εu))^((r-2)/2)⋅εu 
dflux(εu, εdu)=(r-2)*nu0^(1-r)*(ϵ_0^2 + tr(εu' ⋅ εu))^((r-4)/2)*(εu⊙εdu) ⋅ εu + nu0^(1-r)*(ϵ_0^2 + tr(εu' ⋅ εu))^((r-2)/2)*εdu
viskositet(εu) = nu0^(1-r)*(ϵ_0^2 + tr(εu' ⋅ εu))^((r-2)/2)
dviskositet(εu, εdu) = nu0^(1-r) *(ϵ_0^2 + tr(εu' ⋅ εu))^((r-4)/2)*(εu⊙εdu) *(r-2)


# viskositet(εu) = 1/2*A^(1-r) * (ϵ_0^2 + 1/2 * norm(εu)^2)^((r-2)/2)
# dviskositet(εu, εdu) =  (r-2)/4*A^(1-r) *(ϵ_0^2 + 1/2* norm(εu)^2)^((r-4)/2) * (εu⊙εdu)  # eller delt på 4?
# flux(εu) = A^(1-r) * (ϵ_0^2 + 1/2 *norm(εu)^2)^((r-2)/2) ⋅ εu
# dflux(εu, εdu)= (r-2)*A^(1-r)*(ϵ_0^2 + 1/2 * norm(εu)^2)^((r-4)/2)*(εu⊙εdu) ⋅ εu + A^(1-r)*(ϵ_0^2 + 1/2 * norm(εu)^2)^((r-2)/2)*εdu


f(x) =  -divergence(flux∘ε(u_exact))(x) + ∇(p_exact)(x)      # prøver å endre f her...
ud(x) = u_exact(x)
g(x) = -tr(ε(u_exact)(x))

# solver:
#fitted FEM stokes solver:
function pstokes_CutFEM_newton(;n, u_exact, p_exact, f, g, ud, order, geometry, βu0, γu1, γu2, γp, βp0, nu, stabilize, δ, save = false, calc_condition = false)
    """    
    Case A er mye bedre enn Case B
    Fitted FEM with weak boundary imposition using Nitsches method. Symmetric p-stokes equations, using Newtons method to linearize.
    Should test for different values of nu0.
    Symmetric case.
    n: number of grid elements. Powers of 2 for simplicity and convergence estimates.
    u_exact: exact solution for method of manufactured solutions
    order: order of polynomial degree. 
    f: lhs for first term, -Δ u_ex + ∇p = f
    g: lhs for second term u = g
    geometry: optional between "Circle", "Flower", "Heart", "Glacier". Will not affect anything - just added so that the function takes similar function arguments as the stokes_solver
    βu0, γu1, γu2, γp, βp0: parameters for the stokes solver. Will not affect anything - just added so that the function takes similar function arguments as the stokes_solver
    stabilize: will not affect anything - just added so that the function takes similar function arguments as the stokes_solver
    δ: perturbation of cut, will not affect anything - just added so that the function takes similar function arguments as the stokes_solver
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
    
    order = 2
    # Define function spaces 
    reffe_u  = ReferenceFE(lagrangian,VectorValue{2, Float64}, order)
    reffe_p = ReferenceFE(lagrangian,Float64, order - 1)
  
    V = TestFESpace(Ω_act, reffe_u,  conformity=:H1)
    Q = TestFESpace(Ω_act, reffe_p, conformity=:H1, constraint=:zeromean)
    Y = MultiFieldFESpace([V,Q])
  
    U = TrialFESpace(V)         # har fjernet begrensningen til grensebetingelsen ud her..
    P = TrialFESpace(Q)
    X = MultiFieldFESpace([U,P])

    ######################## kvadratish domain - men bruker nitsche implementasjon. Hvilke spaces skal jeg bruke da? #########################
    # weak formulation components   

    a(u, v) = ∫(ε(v)⊙(flux∘ε(u)))dΩ  + ∫(-((n_Γd ⋅ (flux∘ε(u))) ⋅ v) + (-(n_Γd ⋅ (flux∘ε(v))) ⋅ u) + 2*γ/h * (viskositet∘ε(u) ⋅ (u ⋅ v)))dΓd      # denne må ha et ekstra boundary term. Finn ut hvordan det ser ut. 
    b(v, p) = (∫(-1*(∇ ⋅ v*p))dΩ + ∫((n_Γd ⋅ v) * p)dΓd)   # b er den samme som før. 
    l1((u, p), (v, q)) = ∫(f ⋅ v - g ⋅ q)dΩ + ∫(-(n_Γd ⋅ (flux∘ε(v))) ⋅ ud)dΓd + ∫( 2*γ/h* (viskositet∘ε(u) ⋅ (ud ⋅ v)))dΓd
    l2(q) = ∫((n_Γd ⋅ ud) * q)dΓd

    dl1((u, du), (v, q)) = ∫( 2*γ/h*(dviskositet∘(ε(u), ε(du))⋅ (ud ⋅ v)))dΓd 
  
    # a(u, v) = ∫( ε(v)⊙(flux∘ε(u)))dΩ  + ∫(-((n_Γd ⋅ (flux∘ε(u))) ⋅ v) + (-(n_Γd ⋅ (flux∘ε(v))) ⋅ u) + 2*nu0*γ/h * (u ⋅ v))dΓd      # denne må ha et ekstra boundary term. Finn ut hvordan det ser ut. 
    # b(v, p) = (∫(-1*(∇ ⋅ v*p))dΩ + ∫((n_Γd ⋅ v) * p)dΓd)   # b er den samme som før. 
    # l1((v, q)) = ∫(f ⋅ v + g ⋅ q)dΩ
    # l2(v) = ∫(-(n_Γd ⋅ (flux∘ε(v))) ⋅ ud)dΓd    
    # l3((v, q)) = ∫( 2*nu0*γ/h* (ud ⋅ v))dΓd
    # l4(q) = ∫(-1*(n_Γd ⋅ ud) * q)dΓd
    
    da(u, du, v) = ∫( ε(v)⊙(dflux∘(ε(u), ε(du))))dΩ  + ∫(-((n_Γd ⋅ (dflux∘(ε(u), ε(du)))) ⋅ v) + (-(n_Γd ⋅ (flux∘ε(v))) ⋅ du) + (2*γ/h*(dviskositet∘(ε(u), ε(du))⋅ (u ⋅ v)))   +   (2*γ/h*(viskositet∘ε(u)⋅ (du ⋅ v))))dΓd    

    g_u(u,v) = ( ∫( (γu1*h)*  (jump(n_Fg⋅∇(u))⋅jump(n_Fg⋅∇(v))))dFg 
            +  
      ∫( (γu2*h^3)*  (jump_nn(u,n_Fg)⋅jump_nn(v,n_Fg)))dFg
)

    g_p(p, q) = ∫( (γp*h^3)* (jump(n_Fg⋅∇(p))*jump(n_Fg⋅∇(q))))dFg

    
    if stabilize
      res((u,p),(v,q)) = a(u, v) + b(v, p) - b(u, q) - l1((u,p), (v, q)) + l2(q) + g_u(u, v) + g_p(p, q)
      jac((u, p), (du, dp), (v, q)) = b(v, dp) - b(du, q) + da(u, du, v) + g_u(du, v) + g_p(dp, q) - dl1((u, du), (v, q))
    
      op = FEOperator(res, jac, X, Y)

      # non-linear phase
      nls = NLSolver(
      show_trace=true, method=:newton, linesearch=BackTracking(), iterations=100)      #prøver å legge inn et max antall iterasjoner og en lav toleranse      
      solver = FESolver(nls)

      (uh, ph) = solve(solver, op)

    else
      res2((u,p),(v,q)) = a(u, v) + b(v, p) - b(u, q) - l1((u,p), (v, q)) + l2(q) 
      jac2((u, p), (du, dp), (v, q)) = b(v, dp) - b(du, q) + da(u, du, v)   - dl1((u, du), (v, q))
    
      op = FEOperator(res2, jac2, X, Y)

      # non-linear phase
      nls = NLSolver(
      show_trace=true, method=:newton, linesearch=BackTracking(), iterations=100)      #prøver å legge inn et max antall iterasjoner og en lav toleranse      
      solver = FESolver(nls)

      (uh, ph) = solve(solver, op)
    end

    errp = p_exact - ph
    erru = u_exact - uh

    l2_norm(u) = (sum( ∫( u ⋅ u )*dΩ ))^(1/2)
    h1_semi(u) = sum(∫(∇(u) ⊙ ∇(u))*dΩ)^(1/2)
    
    if save
      writevtk(Ω, "test_case_1 $n.vtu", cellfields=["u_ex" => u_exact, "uh"=>uh, "erru"=> erru, "p_ex" => p_exact, "ph"=>ph, "errp"=> errp, "nablau" => ∇(u_exact), "viskositet" => viskositet∘ε(u_exact)]) #, "erru" => erru]) 
    end
    # condition number
    if calc_condition
      condition_numb= cond(Array(get_matrix(op)),2)   # kanskje bruke infinitynormen istedenfor
    else
      condition_numb = 1
    end
    return uh, u_exact, erru, l2_norm(uh - u_exact), h1_semi(uh - u_exact), ph, p_exact, errp, l2_norm(ph - p_exact), h1_semi(ph - p_exact), condition_numb, Ω
end

# solver parameters
n = 128
stabilize = false
solver = pstokes_CutFEM_newton
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

# # ################################################## convergence ##########################################################
numb_it = 8
stabilize = true
uarr_l2, uarr_h1, parr_l2, parr_h1, harr = convergence_stokes(;numb_it, u_exact, p_exact, f, g, ud, order, geometry, solver, δ, βu0, γu1, γu2, γp, βp0, nu, stabilize, save)

stabilize = false
uarr_l2_nostab, uarr_h1_nostab, parr_l2_nostab, parr_h1_nostab, harr = convergence_stokes(;numb_it, u_exact, p_exact, f, g, ud, order, geometry, solver, δ, βu0, γu1, γu2, γp, βp0, nu, stabilize, save)


########################################## lagrer konvergensresultater ##########################################################
# Konvergensresultater lagres til CSV
df_convergence = DataFrame(
    h = harr,
    u_L2_stab = uarr_l2,
    u_H1_stab = uarr_h1,
    u_L2_nostab = uarr_l2_nostab,
    u_H1_nostab = uarr_h1_nostab,
    p_L2_stab = parr_l2,
    p_H1_stab = parr_h1,
    p_L2_nostab = parr_l2_nostab,
    p_H1_nostab = parr_h1_nostab
)

CSV.write("convergence_results_1.csv", df_convergence)

###################################### sensitivitetstester ##########################################################

n = 64                # øke denne?
M = 1000               #full kjøring med M = 2000 med 2000 så kjører det nok i 2 timer. 
# βu0 = 1
# γu1 = 1
# γu2 = 1
# γp = 0.1
# βp0 = 0.1
stabilize = true
save = false
calc_condition = false

arr_δ, arr_l2u, arr_h1u, arr_l2p, arr_h1p, arr_cond = sensitivity_stokes(;n, M, u_exact, p_exact, f, g, ud, order, geometry, solver, βu0, γu1, γu2, γp, βp0, nu, stabilize, save)
stabilize = false
start = 1
save = false
arr_δ_nostab, arr_l2u_nostab, arr_h1u_nostab, arr_l2p_nostab, arr_h1p_nostab, arr_cond_nostab = sensitivity_stokes(;n, M, u_exact, p_exact, f, g, ud, order, geometry, solver, βu0, γu1, γu2, γp, βp0, nu, stabilize, save)


############################################ lagrer sensitivitetstester til CSV ##########################################################

# Sensitivitetsresultater lagres til CSV
df_sensitivity = DataFrame(
    δ = arr_δ,
    u_L2_stab = arr_l2u,
    u_H1_stab = arr_h1u,
    p_L2_stab = arr_l2p,
    p_H1_stab = arr_h1p,
    u_L2_nostab = arr_l2u_nostab,
    u_H1_nostab = arr_h1u_nostab,
    p_L2_nostab = arr_l2p_nostab,
    p_H1_nostab = arr_h1p_nostab
)

CSV.write("sensitivity_results_1.csv", df_sensitivity)
