outputfolderpoisson = "C:\\Users\\Sigri\\Documents\\Master\\report\\results\\poisson"
outputfolderstokes = "C:\\Users\\Sigri\\Documents\\Master\\report\\results\\stokes"
  
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
            geo = AnalyticalGeometry(x -> (-R     +   (x[1])^2   +     ( 5/4 * (x[2]+0.2)   -     sqrt(abs(x[1]))  )^2   ))
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

function jump_nn(u,n)
    return ( (n.plus⋅∇∇(u).plus)⋅ n.plus - (n.minus ⋅ ∇∇(u).minus) ⋅ n.minus )
end

# cut FEM poisson solver:
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
    pmin = Point(-a + δ, -a + δ)            
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
        writevtk(bgmodel, "C:\\Users\\Sigri\\Documents\\Master\\report\\figures\\Domenefigurer\\mesh_bg$geometry $δ.vtu")
        writevtk(Ω_act, "C:\\Users\\Sigri\\Documents\\Master\\report\\figures\\Domenefigurer\\mesh_act_$geometry $δ.vtu")
        writevtk(Ω, "C:\\Users\\Sigri\\Documents\\Master\\report\\figures\\Domenefigurer\\mesh_$geometry $δ.vtu")
        writevtk(Γ, "C:\\Users\\Sigri\\Documents\\Master\\report\\figures\\Domenefigurer\\surface_gamma_d_$geometry $δ.vtu")
        writevtk(Fg, "C:\\Users\\Sigri\\Documents\\Master\\report\\figures\\Domenefigurer\\ghost_facets_$geometry $δ.vtu")
        
        writevtk(Ω, "C:\\Users\\Sigri\\Documents\\Master\\report\\results\\poisson\\$n $geometry $order $δ.vtu", cellfields=["u_ex" => u_exact, "uh"=>uh, "erru"=> erru]) #, "erru" => erru])
    end

return uh, u_exact, erru, l2(uh - u_exact), h1(uh - u_exact), condition_numb, Ω_act, Ω
end

# cut FEM poisson solver:
function convergence_poisson(numb_it, u_exact, lhs, order, geometry, solver, δ, γd, γg1, γg3, stabilize, save = false)
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

# function to calculate the sensitivity of the poisson solver
function sensitivity_poisson(n, M, u_exact, lhs, order, geometry, solver, δ, γd, γ1, γ3, stabilize, save = false)
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
  
# cut FEM stokes solver:
function stokes_cutFEM(;n, u_exact, p_exact, f, g, ud, order, geometry, βu0, γu1, γu2, γp, βp0, nu, stabilize, δ, save = false, calc_condition = false)
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
      a(u, v) = (∫(ε(u) ⊙ ε(v))dΩ + ∫(-((n_Γd ⋅ ε(u)) ⋅ v) + (-(n_Γd ⋅ ε(v)) ⋅ u)+(γ/h * u ⋅ v))dΓd )
      b(v, p) = (∫(-1*(∇ ⋅ v*p))dΩ
                  + ∫((n_Γd ⋅ v) * p)dΓd)#
     
      gu(u,v) = ( ∫( (β_1*h)*jump(n_Fg ⋅ ∇(u))⋅jump(n_Fg⋅ ∇(v)) )dFg 
              +  
                 ∫( (β_2*h^3)*jump_nn(u,n_Fg)⋅jump_nn(v,n_Fg) )dFg)
  
  
      gp(p, q) = (∫((β_3*h^3)*jump(n_Fg ⋅ ∇(p)) * jump(n_Fg ⋅ ∇(q)))dFg)
  
      l2_norm(u) = (sum( ∫( u ⋅ u )*dΩ ))
      h1_semi(u) = sum(∫(∇(u) ⊙ ∇(u))*dΩ)
      
      l1(v) = ∫(f ⋅ v)dΩ
      l2(v) = ∫(-1* (n_Γd ⋅ ε(v)) ⋅ ud)dΓd 
      l3(v) = ∫(γ/h ⋅ v ⋅ ud )dΓd
      l4(q) = ∫(n_Γd ⋅ ud *q)dΓd
  
      if stabilize
          A((u,p),(v,q)) =(a(u, v) + b(v, p) + b(u, q) 
          + gu(u,v)
          - gp(p, q)
          )
          L((v, q)) = l1(v) + l2(v) + l3(v) + l4(q)
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
      #println("lagret")
      writevtk(Ω, "C:\\Users\\Sigri\\Documents\\Master\\report\\results\\stokes\\$n $geometry $order.vtu", cellfields=["u_ex" => u_exact, "uh"=>uh, "erru"=> erru, "p_ex" => p_exact, "ph"=>ph, "errp"=> errp, "nablau" => ∇(u_exact)]) #, "erru" => erru]) 
    end
    return uh, u_exact, erru, l2_norm(uh - u_exact), h1_semi(uh - u_exact), ph, p_exact, errp, l2_norm(ph - p_exact), h1_semi(ph - p_exact), condition_numb, Ω_act
end

#fitted FEM stokes solver:
function stokes_FEM(;n, u_exact, p_exact, f, g, ud, order, geometry, βu0, γu1, γu2, γp, βp0, nu, stabilize, δ, save = false, calc_condition = false)
    """
    Regular fitted FEM solver    
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
    domain = (0,1,0,1)
    partition = (n,n)
    model = CartesianDiscreteModel(domain, partition)    

    labels = get_face_labeling(model)
    
    # alle grensetagsene får dirichlet = g
    add_tag_from_tags!(labels, "dirig", [1, 2, 3, 4, 5, 6, 7, 8])
    order = 2
    reffeᵤ = ReferenceFE(lagrangian,VectorValue{2,Float64},order)
    reffeₚ = ReferenceFE(lagrangian,Float64,order-1;space=:P)
  
    V = TestFESpace(model,reffeᵤ,labels=labels,dirichlet_tags="dirig",conformity=:H1)
    Q = TestFESpace(model,reffeₚ,conformity=:L2,constraint=:zeromean)
    Y = MultiFieldFESpace([V,Q])
  
    U = TrialFESpace(V,ud)
    P = TrialFESpace(Q)
    X = MultiFieldFESpace([U,P])

    degree = order
    Ω = Triangulation(model)
    dΩ = Measure(Ω,degree)
     
    # Weak formulation components
    a(u, v) = ∫( ∇(v)⊙∇(u))dΩ
    b(v, p) = ∫(-(∇⋅v)*p )dΩ
    

    l2_norm(u) = (sum( ∫( u ⋅ u )*dΩ ))
    h1_semi(u) = sum(∫(∇(u) ⊙ ∇(u))*dΩ)
      
    l1(v) = ∫(f ⋅ v)dΩ
    
    A((u,p),(v,q)) =(a(u, v) + b(v, p) - b(u, q))
    L((v, q)) = l1(v)
    
    op = AffineFEOperator(A,L,X,Y)
    uh, ph = solve(op)
  
    errp = p_exact - ph
    erru = u_exact - uh
    
    # condition number
    if calc_condition
      condition_numb= cond(Array(get_matrix(op)),2)   # kanskje bruke infinitynormen istedenfor
    else
      condition_numb = 1
    end
  
    if save
        writevtk(Ω, "C:\\Users\\Sigri\\Documents\\Master\\report\\results\\stokes\\$n $geometry $order.vtu", cellfields=["u_ex" => u_exact, "uh"=>uh, "erru"=> erru, "p_ex" => p_exact, "ph"=>ph, "errp"=> errp, "nablau" => ∇(u_exact)]) #, "erru" => erru]) 
    end
    return uh, u_exact, erru, l2_norm(uh - u_exact), h1_semi(uh - u_exact), ph, p_exact, errp, l2_norm(ph - p_exact), h1_semi(ph - p_exact), condition_numb, Ω
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
      a(u, v) = (∫(ε(u) ⊙ ε(v))dΩ + ∫(-((n_Γd ⋅ ε(u)) ⋅ v) + (-(n_Γd ⋅ ε(v)) ⋅ u)+(γ/h * u ⋅ v))dΓd )
      b(v, p) = (∫(-1*(∇ ⋅ v*p))dΩ
                  + ∫((n_Γd ⋅ v) * p)dΓd)#
     
      gu(u,v) = ( ∫( (β_1*h)*jump(n_Fg ⋅ ∇(u))⋅jump(n_Fg⋅ ∇(v)) )dFg 
              +  
                 ∫( (β_2*h^3)*jump_nn(u,n_Fg)⋅jump_nn(v,n_Fg) )dFg)
  
      gp(p, q) = (∫((β_3*h^3)*jump(n_Fg ⋅ ∇(p)) * jump(n_Fg ⋅ ∇(q)))dFg)
  
      l2_norm(u) = (sum( ∫( u ⋅ u )*dΩ ))
      h1_semi(u) = sum(∫(∇(u) ⊙ ∇(u))*dΩ)
      
      l1(v) = ∫(f ⋅ v)dΩ
      l2(v) = ∫(-1* (n_Γd ⋅ ε(v)) ⋅ ud + γ/h ⋅ v ⋅ ud )dΓd
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
      writevtk(Ω, "C:\\Users\\Sigri\\Documents\\Master\\report\\results\\stokes\\$n $geometry $order.vtu", cellfields=["u_ex" => u_exact, "uh"=>uh, "erru"=> erru, "p_ex" => p_exact, "ph"=>ph, "errp"=> errp, "nablau" => ∇(u_exact)]) #, "erru" => erru]) 
    end
    return uh, u_exact, erru, l2_norm(uh - u_exact), h1_semi(uh - u_exact), ph, p_exact, errp, l2_norm(ph - p_exact), h1_semi(ph - p_exact), condition_numb, Ω_act
end

#fitted FEM for p-stokes:
function nonlinear_stokes_FEM(;n, u_exact, p_exact, f, g, ud, order, geometry, βu0, γu1, γu2, γp, βp0, nu, stabilize, δ, save = false, calc_condition = false)
       """
    Fitted FEM, non-linear stokes (p-stokes).Using P2-P1 Taylor-Hood elements.  
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
    domain = (0,1,0,1)
    partition = (n,n)
    model = CartesianDiscreteModel(domain, partition)    

    labels = get_face_labeling(model)
    
    # alle grensetagsene får dirichlet = g
    add_tag_from_tags!(labels, "dirig", [1, 2, 3, 4, 5, 6, 7, 8])
    order = 2
    reffeᵤ = ReferenceFE(lagrangian,VectorValue{2,Float64},order)
    reffeₚ = ReferenceFE(lagrangian,Float64,order-1;space=:P)
  
    V = TestFESpace(model,reffeᵤ,labels=labels,dirichlet_tags="dirig",conformity=:H1)
    Q = TestFESpace(model,reffeₚ,conformity=:L2,constraint=:zeromean)
    Y = MultiFieldFESpace([V,Q])
  
    U = TrialFESpace(V,ud)
    P = TrialFESpace(Q)
    X = MultiFieldFESpace([U,P])

    degree = order
    Ω = Triangulation(model)
    dΩ = Measure(Ω,degree)
    
    println(sum( ∫( p_exact) * dΩ ))

    l2_norm(u) = (sum( ∫( u ⋅ u )*dΩ ))
    h1_semi(u) = sum(∫(∇(u) ⊙ ∇(u))*dΩ)
    
    # fra klassisk Stokes FEM, så er det kun herfra og ned som er endret:)
    nu0 = 1
    ϵ_0 = 1e-6
    
    a(u, v) = ∫( ∇(v)⊙(flux∘∇(u)))dΩ
    b(v, p) = ∫(-(∇⋅v)*p )dΩ
    l(v) = ∫(f ⋅ v)dΩ

    # dflux calculated the same way as in the notebook p-Laplace...
    dflux(∇du,∇u)=(r-2)*(ϵ_0 + norm(∇u)^2)^((r-4)/2)*(∇u⊙∇du) ⋅ ∇u + (ϵ_0 + norm(∇u)^2)^((r-2)/2)*∇du

    # and introduced in the same way as in the notebook p-Laplace in the bilinear form a...
    da(u, du, v) = ∫(∇(v)⊙(dflux∘(∇(du), ∇(u))))dΩ
    
    # and then the Newton multifield system is assembled as in the Navier Stokes notebook...
    res((u,p),(v,q)) = ∫( ∇(v)⊙(flux∘∇(u)))dΩ + ∫(-(∇⋅v)*p )dΩ - ∫(-(∇⋅u)*q )dΩ - ∫(f ⋅ v)dΩ    #a(u, v)  + b(v, p) - b(u, q) - l(v)       # bytte til epsilon her, og legge til uttrykkene for b(u, q), b(v, p)
    jac((u, p), (du, dp), (v, q)) =  b(v, dp) - b(du, q) + da(u, du, v)
    # her har det ikke noe å si om jeg plusser eller trekker fra b(u, q) siden det ikke er lagt til noe annet ledd (pga strong dirichlet)
    op = FEOperator(res, jac, X, Y)

    # non-linear phase
    nls = NLSolver(
    show_trace=true, method=:newton, linesearch=BackTracking())      #prøver å legge inn et max antall iterasjoner og en lav toleranse      
    solver = FESolver(nls)

    (uh, ph) = solve(solver, op)

    errp = p_exact - ph
    erru = u_exact - uh
    
    # condition number
    if calc_condition
      condition_numb= cond(Array(get_matrix(op)),2)   # kanskje bruke infinitynormen istedenfor
    else
      condition_numb = 1
    end
  
    if save
        writevtk(Ω, "C:\\Users\\Sigri\\Documents\\Master\\report\\results\\stokes\\$n $geometry $order.vtu", cellfields=["u_ex" => u_exact, "uh"=>uh, "erru"=> erru, "p_ex" => p_exact, "ph"=>ph, "errp"=> errp, "nablau" => ∇(u_exact)]) #, "erru" => erru]) 
    end
    return uh, u_exact, erru, l2_norm(uh - u_exact), h1_semi(uh - u_exact), ph, p_exact, errp, l2_norm(ph - p_exact), h1_semi(ph - p_exact), condition_numb, Ω
end

function convergence_stokes(;numb_it, u_exact, p_exact, f, g, ud, order, geometry, solver, δ, βu0, γu1, γu2, γp, βp0, nu, stabilize, save = false)
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

  return uarr_l2, uarr_h1, parr_l2, parr_h1, h#, EOC_l2, EOC_h1
end

function sensitivity_stokes(;n, M, u_exact, p_exact, f, g, ud, order, geometry, solver, βu0, γu1, γu2, γp, βp0, nu, stabilize, save)
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



function stokes_plot(
    δ_stab, δ_nostab,                          # δ-verdier (x-verdier)
    datasets::Vector{Tuple},                   # Liste med: (y_stab, y_nostab, label, ylabel)
    idx = 1:100:1999,                          # Markørindekser for stabiliserte
    id2 = 51:100:1999,                         # Markørindekser for ustabiliserte
    start::Int = 1, end_::Int = 1999,          # Plot-utvalg
    plot_title = "Stokes Solver Analysis"
)

    plt = plot(0, title=plot_title,
               xlabel="Perturbation δ",
               titlefont=16,
               guidefont=14,
               tickfont=12)

    for (y_stab, y_nostab, label, ylabel) in datasets
        # Markører for stabiliserte
        scatter!(δ_stab[idx], y_stab[idx], label="", marker=:circle, ms=4)
        # Linje for stabiliserte
        plot!(δ_stab[start:end_], y_stab[start:end_], yaxis=:log, lw=2, label="$label stabilized")
        # Markører for ikke-stabiliserte
        scatter!(δ_nostab[id2], y_nostab[id2], label="", marker=:s, ms=4)
        # Linje for ikke-stabiliserte
        plot!(δ_nostab[start:end_], y_nostab[start:end_], yaxis=:log, lw=2, label="$label non-stabilized")
        ylabel!(ylabel)
    end

    return plt
end