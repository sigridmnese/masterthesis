# Denne fungerer bra for nu = 1, MEN skal nå få lagt inn viskositetsavhengighet ijacobi-determinanten. Prøver først å kjøre noen tester der jeg fjerner gu på residualen.
function slettes(;n, u_exact, p_exact, f, g, ud, order, geometry, βu0, γu1, γu2, γp, βp0, nu, stabilize, δ, save = false, calc_condition = false)
    """
    Unfitted FEM with Nitsche boundary imposition, non-linear stokes (p-stokes).Using P2-P1 Taylor-Hood elements.  
    n: number of grid elements. Powers of 2 for simplicity and convergence estimates.
    u_exact: exact solution for method of manufactured solutions
    order: order of polynomial degree. 
    f: lhs for first term, -∇ (ν ∇(u_ex) + ∇p = f
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
  
    geo = create_geometry(geometry, n)
    # Define active and physical mesh
    cutgeo = cut(bgmodel,geo)
    cutgeo_facets = cut_facets(bgmodel,geo)
    Ω_bg = Triangulation(bgmodel)
    Ω_act = Triangulation(cutgeo, ACTIVE)
    Ω = Triangulation(cutgeo, PHYSICAL)
    γ = 2*2*order
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
    
    #println(sum( ∫( p_exact) * dΩ ))
    l2_norm(u) = (sum( ∫( u ⋅ u )*dΩ ))
    h1_semi(u) = sum(∫(∇(u) ⊙ ∇(u))*dΩ)
    
    # fra klassisk Stokes FEM, så er det kun og ned som er endret:)
    ϵ_0 = 1e-6
    γ = 10*2*2
    #γ(∇u) = viskositet(∇u)/h
    #γ(u) = viskositet∘∇(u)
    # weak formulation components    nu0^(1-r)*(ϵ_0^2 + (norm∘ε(u))^2)^((r-2)/2)
    a(u, v) = ∫( ε(v)⊙(flux∘ε(u)))dΩ  + ∫(-((n_Γd ⋅ (flux∘ε(u))) ⋅ v) + (-(n_Γd ⋅ (flux∘ε(v))) ⋅ u)+ 1/h * (u ⋅ v))dΓd      # denne må ha et ekstra boundary term. Finn ut hvordan det ser ut. 
    b(v, p) = (∫(-1*(∇ ⋅ v*p))dΩ + ∫((n_Γd ⋅ v) * p)dΓd)   # b er den samme som før. 
    l1(v) = ∫(f ⋅ v)dΩ
    l2(v) = ∫(-(n_Γd ⋅ (flux∘ε(v))) ⋅ ud)dΓd    # har brukt ud som dirichlet grense
    l3(v) = ∫(γ/h * (ud ⋅ v))dΓd
    l4(q) = ∫((n_Γd ⋅ ud) * q)dΓd               # kan det være at denne skal være annerledes? mtp flux * ... og ikke n_Γd ⋅ ud?
    # dflux calculated the same way as in the notebook p-Laplace...
    dflux(∇du,∇u)=(r-2)*(ϵ_0 + norm(∇u)^2)^((r-4)/2)*(∇u⊙∇du) ⋅ ∇u + (ϵ_0 + norm(∇u)^2)^((r-2)/2)*∇du
    
    da(u, du, v) = ∫( ε(v)⊙(dflux∘(ε(du), ε(u))))dΩ  + ∫(-((n_Γd ⋅ (dflux∘(ε(du), ε(u)))) ⋅ v) + (-(n_Γd ⋅ (flux∘ε(v))) ⋅ du)+(γ/h*(du ⋅ v)))dΓd       

    gu(u, v) = (∫((β_1 * h)*jump(n_Fg ⋅ ε(u))⋅jump(n_Fg⋅ ε(v)) )dFg  # prøver å multiplisere med viskositet her
              +  
                 ∫( (β_2*h^3 )*jump_nn(u,n_Fg)⋅jump_nn(v,n_Fg) )dFg)

    gp(p, q) = (∫((β_3*h)*jump(n_Fg ⋅ ∇(p)) * jump(n_Fg ⋅ ∇(q)))dFg)

    if stabilize
      res((u,p),(v,q)) = a(u, v) + b(v, p) + b(u, q) - l1(v) - l2(v) - l3(v) - l4(q) #- gp(p, q)+ gu(u, v)
      jac((u, p), (du, dp), (v, q)) = b(v, dp) + b(du, q) + da(u, du, v) + gu(u, v) - gp(p, q)  #nå er ghos penaltiene lineære, men prøver å heller få inn ikke-linearite i konstane foran...
      
      op = FEOperator(res, jac, X, Y)

      # non-linear phase
      nls = NLSolver(
      show_trace=true, method=:newton, linesearch=BackTracking(), iterations=40)      #prøver å legge inn et max antall iterasjoner og en lav toleranse      
      solver = FESolver(nls)

      (uh, ph) = solve(solver, op)

    else
      res_nostab((u,p),(v,q)) = a(u, v) + b(v, p) + b(u, q) - l1(v) - l2(v) - l3(v) - l4(q)
      jac_nostab((u, p), (du, dp), (v, q)) = b(v, dp) + b(du, q) + da(u, du, v)
      
      op = FEOperator(res_nostab, jac_nostab, X, Y)

      # non-linear phase
      nls = NLSolver(
      show_trace=true, method=:newton, linesearch=BackTracking(), iterations=40)      #prøver å legge inn et max antall iterasjoner og en lav toleranse      
      solver = FESolver(nls)

      (uh, ph) = solve(solver, op)
    end

    errp = p_exact - ph
    erru = u_exact - uh
    
    condition_numb = 1
    if save
        #writevtk(bgmodel, "C:\\Users\\Sigri\\Documents\\Master\\report\\figures\\Domenefigurer\\mesh_bg$geometry $δ.vtu")
        #writevtk(Γd, "C:\\Users\\Sigri\\Documents\\Master\\report\\figures\\Domenefigurer\\surface_gamma_d_$geometry $δ.vtu")
        #writevtk(n_Γd, "C:\\Users\\Sigri\\Documents\\Master\\report\\figures\\Domenefigurer\\surface_gamma_d_$geometry $δ.vtu")
        writevtk(Ω, "C:\\Users\\Sigri\\Documents\\Master\\report\\results\\stokes\\$n $geometry $order $δ.vtu", cellfields=["u_ex" => u_exact, "uh"=>uh, "erru"=> erru, "p_ex" => p_exact, "ph"=>ph, "errp"=> errp, "nablau" => ∇(u_exact), "flux" => flux∘ε(u_exact), "viskositet" => viskositet∘ε(u_exact)]) #, "erru" => erru]) 
    end 
    return uh, u_exact, erru, l2_norm(uh - u_exact), h1_semi(uh - u_exact), ph, p_exact, errp, l2_norm(ph - p_exact), h1_semi(ph - p_exact), condition_numb, Ω
end

# p stokes cut FEM fungerer veldig bra, men må få lagt inn rikig viskositetsavhengighet.
function blablabla(;n, u_exact, p_exact, f, g, ud, order, geometry, βu0, γu1, γu2, γp, βp0, nu, stabilize, δ, save = false, calc_condition = false)
    """
    cutFEM solver for the non-linear stokes (p-stokes) equations. Using P2-P1 Taylor-Hood elements.  
    n: number of grid elements. Powers of 2 for simplicity and convergence estimates.
    u_exact: exact solution for method of manufactured solutions
    order: order of polynomial degree. 
    f: lhs for first term, -∇ (ν ∇(u_ex) + ∇p = f
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
    
    #println(sum( ∫( p_exact) * dΩ ))
    l2_norm(u) = (sum( ∫( u ⋅ u )*dΩ ))
    h1_semi(u) = sum(∫(∇(u) ⊙ ∇(u))*dΩ)
    
    ϵ_0 = 1e-6

    #defining viskosity parameter:
    # You have to chage the call for γ in three functions: a - last term, l3 - first term, da - last term
    γ = 10*order*order
    #γ(∇u) = viskositet(∇u)
    #γ(u) = viskositet∘ε(u)

    # weak formulation components    
    a(u, v) = ∫( ε(v)⊙(flux∘ε(u)))dΩ  + ∫(-((n_Γd ⋅ (flux∘ε(u))) ⋅ v) + (-(n_Γd ⋅ (flux∘ε(v))) ⋅ u)+ γ/h* (u ⋅ v))dΓd      # denne må ha et ekstra boundary term. Finn ut hvordan det ser ut. 
    b(v, p) = (∫(-1*(∇ ⋅ v*p))dΩ + ∫((n_Γd ⋅ v) * p)dΓd)   # b er den samme som før. 
    l1(v) = ∫(f ⋅ v)dΩ
    l2(v) = ∫(-(n_Γd ⋅ (flux∘ε(v))) ⋅ ud)dΓd    # har brukt ud som dirichlet grense
    l3(v) = ∫( γ/h* (ud ⋅ v))dΓd
    l4(q) = ∫((n_Γd ⋅ ud) * q)dΓd               # kan det være at denne skal være annerledes? mtp flux * ... og ikke n_Γd ⋅ ud?
    # dflux calculated the same way as in the notebook p-Laplace...
    dflux(∇du,∇u)=(r-2)*(ϵ_0 + norm(∇u)^2)^((r-4)/2)*(∇u⊙∇du) ⋅ ∇u + (ϵ_0 + norm(∇u)^2)^((r-2)/2)*∇du
    
    da(u, du, v) = ∫( ε(v)⊙(dflux∘(ε(du), ε(u))))dΩ  + ∫(-((n_Γd ⋅ (dflux∘(ε(du), ε(u)))) ⋅ v) + (-(n_Γd ⋅ (flux∘ε(v))) ⋅ du)+(γ/h*(du ⋅ v)))dΓd       

    gu(u, v) = (∫((β_1 * h)*jump(n_Fg ⋅ ε(u))⋅jump(n_Fg⋅ ε(v)) )dFg  # prøver å multiplisere med viskositet her
              +  
                 ∫((β_2*h^3 )*jump_nn(u,n_Fg)⋅jump_nn(v,n_Fg) )dFg)

    gp(p, q) = (∫((β_3*h)*jump(n_Fg ⋅ ∇(p)) * jump(n_Fg ⋅ ∇(q)))dFg)  #eventuelt h^3 her, men synes h^1 fungerer tilsnelatende bedre...


    if stabilize # have tested with different combinations of calling the gp in the residual and jacobian, but this one seems to work best.
      res((u,p),(v,q)) = a(u, v) + b(v, p) + b(u, q)+ gu(u,v) - gp(p, q) -l1(v) -l2(v) -l3(v) -l4(q) 
      jac((u, p), (du, dp), (v, q)) = b(v, dp) + b(du, q) + da(u, du, v) + gu(du,v) - gp(dp, q)  
      op = FEOperator(res, jac, X, Y)

      # non-linear phase
      nls = NLSolver(
      show_trace=true, method=:newton, linesearch=BackTracking(), iterations=40)      #prøver å legge inn et max antall iterasjoner og en lav toleranse      
      solver = FESolver(nls)

      (uh, ph) = solve(solver, op)

    else
      res_nostab((u,p),(v,q)) = a(u, v) + b(v, p) + b(u, q) - l1(v) -l2(v) -l3(v) -l4(q)
      jac_nostab((u, p), (du, dp), (v, q)) = b(v, dp) + b(du, q) + da(u, du, v)
      
      op = FEOperator(res_nostab, jac_nostab, X, Y)

      # non-linear phase
      nls = NLSolver(
      show_trace=true, method=:newton, linesearch=BackTracking(), iterations=40)      #prøver å legge inn et max antall iterasjoner og en lav toleranse      
      solver = FESolver(nls)

      (uh, ph) = solve(solver, op)
    end

    errp = p_exact - ph
    erru = u_exact - uh
    
    condition_numb = 1
    if save
        #writevtk(bgmodel, "C:\\Users\\Sigri\\Documents\\Master\\report\\figures\\Domenefigurer\\mesh_bg$geometry $δ.vtu")
        #writevtk(Γd, "C:\\Users\\Sigri\\Documents\\Master\\report\\figures\\Domenefigurer\\surface_gamma_d_$geometry $δ.vtu")
        #writevtk(n_Γd, "C:\\Users\\Sigri\\Documents\\Master\\report\\figures\\Domenefigurer\\surface_gamma_d_$geometry $δ.vtu")
        writevtk(Ω, "C:\\Users\\Sigri\\Documents\\Master\\report\\results\\stokes\\$n $geometry $order $δ.vtu", cellfields=["u_ex" => u_exact, "uh"=>uh, "erru"=> erru, "p_ex" => p_exact, "ph"=>ph, "errp"=> errp, "nablau" => ∇(u_exact), "flux" => flux∘ε(u_exact), "viskositet" => viskositet∘ε(u_exact)]) #, "erru" => erru]) 
    end 
    return uh, u_exact, erru, l2_norm(uh - u_exact), h1_semi(uh - u_exact), ph, p_exact, errp, l2_norm(ph - p_exact), h1_semi(ph - p_exact), condition_numb, Ω
end

############################################## ALLE UNDER HER ER PICARD LØSERE###############################################
function picard_iteration(;make_a, make_l, nu_func, tol=1e-5, maxiter = 100)
    while true
        nu = nu_func(∇(uh))
        A = make_a(nu)
        L = make_l(nu)
        
        # Solve the linear system A * u = L
        uh = solve(A, L)

        # Calculate the residual
        res = L - A * uh

        # Check convergence
        if norm(res) < tol
            break
        end
        
        # Check for maximum iterations
        if maxiter <= 0
            error("Picard iteration did not converge within the maximum number of iterations.")
        end
        maxiter -= 1
    end
end

# I denne prøvde jeg å bruke cutfem med newton-løser, men det går selvfølgelig ikke, for det er ikke en lineær form.
# FUNGERER IKKE:
function teste(;n, u_exact, p_exact, f, g, ud, order, geometry, βu0, γu1, γu2, γp, βp0, nu, stabilize, δ, save = false, calc_condition = false)
    """
    Unfitted FEM with Nitsche boundary imposition, non-linear stokes (p-stokes).Using P2-P1 Taylor-Hood elements.  
    n: number of grid elements. Powers of 2 for simplicity and convergence estimates.
    u_exact: exact solution for method of manufactured solutions
    order: order of polynomial degree. 
    f: lhs for first term, -∇ (ν ∇(u_ex) + ∇p = f
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
  
    geo = create_geometry(geometry, n)
    # Define active and physical mesh
    cutgeo = cut(bgmodel,geo)
    cutgeo_facets = cut_facets(bgmodel,geo)
    Ω_bg = Triangulation(bgmodel)
    Ω_act = Triangulation(cutgeo, ACTIVE)
    Ω = Triangulation(cutgeo, PHYSICAL)
    γ = 2*2*order
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
    
    #println(sum( ∫( p_exact) * dΩ ))
    l2_norm(u) = (sum( ∫( u ⋅ u )*dΩ ))
    h1_semi(u) = sum(∫(∇(u) ⊙ ∇(u))*dΩ)
    
    # fra klassisk Stokes FEM, så er det kun herfra og ned som er endret:)
    ϵ_0 = 1e-6
    γ = 10*2*2
    # weak formulation components    
    a(u, v) = ∫( ε(v)⊙(flux∘ε(u)))dΩ  + ∫(-((n_Γd ⋅ (flux∘ε(u))) ⋅ v) + (-(n_Γd ⋅ (flux∘ε(v))) ⋅ u)+(γ/h * (u ⋅ v)))dΓd      # denne må ha et ekstra boundary term. Finn ut hvordan det ser ut. 
    b(v, p) = (∫(-1*(∇ ⋅ v*p))dΩ + ∫((n_Γd ⋅ v) * p)dΓd)   # b er den samme fom før. 
    l1(v) = ∫(f ⋅ v)dΩ
    l2(v) = ∫(-(n_Γd ⋅ (flux∘ε(v))) ⋅ ud)dΓd    # har brukt ud som dirichlet grense
    l3(v) = ∫(γ/h * (ud ⋅ v))dΓd
    l4(q) = ∫((n_Γd ⋅ ud) * q)dΓd
    # dflux calculated the same way as in the notebook p-Laplace...
    dflux(∇du,∇u)=(r-2)*(ϵ_0 + norm(∇u)^2)^((r-4)/2)*(∇u⊙∇du) ⋅ ∇u + (ϵ_0 + norm(∇u)^2)^((r-2)/2)*∇du
    # and introduced in the same way as in the notebook p-Laplace in the bilinear form a...
    #da(u, du, v) = ∫(∇(v)⊙(dflux∘(∇(du), ∇(u))))dΩ
   
    da(u, du, v) = ∫( ε(v)⊙(dflux∘(ε(du), ε(u))))dΩ  + ∫(-((n_Γd ⋅ (dflux∘(ε(du), ε(u)))) ⋅ v) + (-(n_Γd ⋅ (flux∘ε(v))) ⋅ du)+(γ/h * (du ⋅ v)))dΓd       
    
    gu(u,v) = ((β_1 * h)*∫(jump(n_Fg ⋅ ε(u))⋅jump(n_Fg ⋅ ε(v)) )dFg)  # prøver å bare straffe med første orden - bør hjelpe litt tror jeg?
              #+  
              #   ∫( (β_2*h^3 )*jump_nn(u,n_Fg)⋅jump_nn(v,n_Fg) )dFg)
    
    gp(p, q) = (∫((β_3*h^3)*jump(n_Fg ⋅ ∇(p)) * jump(n_Fg ⋅ ∇(q)))dFg)

    nit = 20
    if stabilize
      res((u,p),(v,q)) = a(u, v) + b(v, p) + b(u, q) + (β_1 * h) *gu(u, v) - gp(p, q) -l1(v) -l2(v) -l3(v) -l4(q)
      jac((u, p), (du, dp), (v, q)) = b(v, dp) + b(du, q) + da(u, du, v)    # kan ikke ha med ghost penalties i jacobien, siden disse ikke er differensierbare
      
      op = FEOperator(res, jac, X, Y)

      # non-linear phase
      nls = NLSolver(
      show_trace=true, method=:newton, linesearch=BackTracking(), iterations=nit)      #prøver å legge inn et max antall iterasjoner og en lav toleranse      
      solver = FESolver(nls)

      (uh, ph) = solve(solver, op)

    else
      res_nostab((u,p),(v,q)) = a(u, v) + b(v, p) + b(u, q) - l1(v) -l2(v) -l3(v) -l4(q)
      jac_nostab((u, p), (du, dp), (v, q)) = b(v, dp) + b(du, q) + da(u, du, v)
      
      op = FEOperator(res_nostab, jac_nostab, X, Y)

      # non-linear phase
      nls = NLSolver(
      show_trace=false, method=:newton, linesearch=BackTracking(), iterations=nit)      #prøver å legge inn et max antall iterasjoner og en lav toleranse      
      solver = FESolver(nls)

      (uh, ph) = solve(solver, op)
    end

    errp = p_exact - ph
    erru = u_exact - uh

    condition_numb = 1
    if save
        writevtk(Ω, "C:\\Users\\Sigri\\Documents\\Master\\report\\results\\stokes\\$n $geometry $order $δ.vtu", cellfields=["u_ex" => u_exact, "uh"=>uh, "erru"=> erru, "p_ex" => p_exact, "ph"=>ph, "errp"=> errp, "nablau" => ∇(u_exact), "flux" => flux∘ε(u_exact), "viskositet" => viskositet∘ε(u_exact)]) #, "erru" => erru]) 
    end 
    return uh, u_exact, erru, l2_norm(uh - u_exact), h1_semi(uh - u_exact), ph, p_exact, errp, l2_norm(ph - p_exact), h1_semi(ph - p_exact), condition_numb, Ω
end

# teste 2 er med picard-iterasjon på fitted FEM. 
function teste2(;n, u_exact, p_exact, f, g, ud, order, geometry, βu0, γu1, γu2, γp, βp0, nu, stabilize, δ, save = false, calc_condition = false)
    """
    Forbedringspunkter:
    Lage egen Picard-iterasjon-funksjon hvis det er mulig
    Får konstant konvergensrate med Picard-iterasjon, hva kan være årsaken til det?
    Får fortsatt konstant konvergensrate med Picard-iterasjonen, og klarer ikke å redusere feilen...
    Når jeg kjører for nu = 100, så fluktuerer 2-steg feilen rundt 10^2
    Og 2-steg feilen er omtrent 10^7 fra starten av. 
    Spørre Andre om dette. Kan være det er noe feil med løseren? Den var kanskje ikke helt riktig i utgangspuntket?

    Fitted FEM for the p-stokes equation.
    n: number of grid elements. Powers of 2 for simplicity and convergence estimates.
    u_exact: exact solution for method of manufactured solutions
    order: order of polynomial degree. 
    f: lhs for first term, -Δ u_ex + ∇p = f
    g: lhs for second term u = g
    geometry: optional between "Circle", "Flower", "Heart", "Glacier".
    stabilize: wheather to add the stabilization term or not
    δ: perturbation of cut

    Defining Picard iteration constants and variables further down in the code.
    """
    # Define background mesh
    domain = (0,1,0,1)
    partition = (n,n)
    model = CartesianDiscreteModel(domain, partition)    
    
    # dirichlet bc på alle rander
    labels = get_face_labeling(model)
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
    
    # classical stokes fem
    nu0 = 1
    ϵ_0 = 1e-6

    # initializing the picard iteration - consider putting this into a separate function
    maxiter = 100
    tol = 1e-8

    Random.seed!(1234)
    x = rand(Float64,num_free_dofs(U))
    y = rand(Float64,num_free_dofs(P))
    uh = FEFunction(U,x)
    ph = FEFunction(P,y)

    while true
      order_ν = 1
      reffe_ν  = ReferenceFE(lagrangian, Float64, order_ν)
      V0       = TestFESpace(model, reffe_ν, conformity=:H1)

      ν_expr(x) = begin
        εu = ε(uh)(x)                   # TensorValue{2,2}
        mag2 = sum(εu .* εu)            # Float64
        return nu0 * (ϵ_0^2 + mag2)^((r-2)/2)
      end
      nuh = interpolate(ν_expr, V0)     # nuh er nå et skalar‐felt over hele meshet

      println("Picard iteration: ", maxiter)
      a(u, v) = ∫( ε(v)⊙(nuh*ε(u)))dΩ # + ∫(-((∇(u) ⋅ v) + (∇(v) ⋅ u)) + (γ/h * (u ⋅ v)))dΩ      # denne må ha et ekstra boundary term. Finn ut hvordan det ser ut.
      b(v, p) = (∫(-1*(∇ ⋅ v*p))dΩ)   # b er den samme fom før.
      l(v) = ∫(f ⋅ v)dΩ # + ∫(-(∇(uh) ⋅ v) * ud)dΩ + ∫(γ/h * (ud ⋅ v))dΩ  # har brukt ud som dirichlet grense
      A((u,p),(v,q)) = a(u, v) + b(v, p) - b(u, q)
      L((v,q)) = l(v)

      # assembling system matrices and solving
      op = AffineFEOperator(A,L,X,Y)
      u_new, p_new = solve(op)

      # change in the Picard iteration
      Δu = u_new - uh
      Δp = p_new - ph
      println("Iterasjon $maxiter: ‖Δu‖ = ", l2_norm(Δu), ", ‖Δp‖ = ", l2_norm(Δp))
      # Check convergence
      if l2_norm(Δu) < tol
          break
      end
        
      # Check for maximum iterations
      if maxiter <= 0
        println("Picard iteration did not converge within the maximum number of iterations.")
        break
      end

      # updating with the new solutions
      maxiter -= 1
      uh = u_new
      ph = p_new
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
        writevtk(Ω, "C:\\Users\\Sigri\\Documents\\Master\\report\\results\\stokes\\$n $geometry $order.vtu", cellfields=["u_ex" => u_exact, "uh"=>uh, "erru"=> erru, "p_ex" => p_exact, "ph"=>ph, "errp"=> errp, "nablau" => ∇(u_exact), "viskositet" => viskositet∘ε(u_exact)]) #, "erru" => erru]) 
    end
    return uh, u_exact, erru, l2_norm(uh - u_exact), h1_semi(uh - u_exact), ph, p_exact, errp, l2_norm(ph - p_exact), h1_semi(ph - p_exact), condition_numb, Ω
end


function teste3(;n, u_exact, p_exact, f, g, ud, order, geometry, βu0, γu1, γu2, γp, βp0, nu, stabilize, δ, save = false, calc_condition = false)
       """
    Fitted FEM, non-linear stokes (p-stokes). Using P2-P1 Taylor-Hood elements.  
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

    # taylor hood elements
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
    

    ############################################ from here I comment out################################
    # a(u, v) = ∫( ε(v)⊙(flux∘ε(u)))dΩ
    # b(v, p) = ∫(-(∇⋅v)*p )dΩ
    # l(v) = ∫(f ⋅ v)dΩ

    # # dflux calculated the same way as in the notebook p-Laplace...
    # dflux(∇du,∇u)=(r-2)*(ϵ_0 + norm(∇u)^2)^((r-4)/2)*(∇u⊙∇du) ⋅ ∇u + (ϵ_0 + norm(∇u)^2)^((r-2)/2)*∇du

    # # and introduced in the same way as in the notebook p-Laplace in the bilinear form a...
    # da(u, du, v) = ∫(ε(v)⊙(dflux∘(ε(du), ε(u))))dΩ
    
    # # and then the Newton multifield system is assembled as in the Navier Stokes notebook...
    # res((u,p),(v,q)) = ∫( ε(v)⊙(flux∘ε(u)))dΩ + ∫(-(∇⋅v)*p )dΩ - ∫(-(∇⋅u)*q )dΩ - ∫(f ⋅ v)dΩ    #a(u, v)  + b(v, p) - b(u, q) - l(v)       # bytte til epsilon her, og legge til uttrykkene for b(u, q), b(v, p)
    # jac((u, p), (du, dp), (v, q)) =  b(v, dp) - b(du, q) + da(u, du, v)

    # op = FEOperator(res, jac, X, Y)

    # # non-linear phase
    # nls = NLSolver(
    # show_trace=true, method=:newton, linesearch=BackTracking(), iterations=20)      #prøver å legge inn et max antall iterasjoner og en lav toleranse      
    # solver = FESolver(nls)

    # (uh, ph) = solve(solver, op)

    ############################################To here!! ################################
    # initializing the picard iteration - consider putting this into a separate function
    maxiter = 10
    tol = 1e-8

    Random.seed!(1234)
    x = rand(Float64,num_free_dofs(U))
    y = rand(Float64,num_free_dofs(P))
    uh = FEFunction(U,x)
    ph = FEFunction(P,y)

    while true
      order_ν = 2
      reffe_ν  = ReferenceFE(lagrangian, Float64, order_ν)
      V0       = TestFESpace(model, reffe_ν, conformity=:H1)

      ν_expr(x) = begin
        εu = ε(uh)(x)                   # TensorValue{2,2}
        mag2 = sum(εu .* εu)            # Float64
        return nu0 * (ϵ_0^2 + mag2)^((r-2)/2)
      end
      nuh = interpolate(ν_expr, V0)     # nuh er nå et skalar‐felt over hele meshet
      
      a(u, v) = ∫( nuh * (ε(v)⊙(ε(u))))dΩ # + ∫(-((∇(u) ⋅ v) + (∇(v) ⋅ u)) + (γ/h * (u ⋅ v)))dΩ      # denne må ha et ekstra boundary term. Finn ut hvordan det ser ut.
      b(v, p) = ∫(-(∇⋅v)*p )dΩ   # b er den samme fom før.
      l(v) = ∫(f ⋅ v)dΩ # + ∫(-(∇(uh) ⋅ v) * ud)dΩ + ∫(γ/h * (ud ⋅ v))dΩ  # har brukt ud som dirichlet grense
      A((u,p),(v,q)) = a(u, v) + b(v, p) + b(u, q)
      L((v,q)) = l(v)

      # assembling system matrices and solving
      op = AffineFEOperator(A,L,X,Y)
      u_new, p_new = solve(op)

      # change in the Picard iteration
      Δu = u_new - uh
      Δp = p_new - ph
      println("Picard iteration: $maxiter: ‖Δu‖ = ", l2_norm(Δu), ", ‖Δp‖ = ", l2_norm(Δp))
      # Check convergence
      if l2_norm(Δu + Δp) < tol
          break
      end
        
      # Check for maximum iterations
      if maxiter <= 0
        println("Picard iteration did not converge within the maximum number of iterations.")
        break
      end

      # updating with the new solutions
      maxiter -= 1
      uh = u_new
      ph = p_new
    end

    errp = p_exact - ph
    erru = u_exact - uh
    
    # condition number
    #if calc_condition
    #  condition_numb= cond(Array(get_matrix(op)),2)   # kanskje bruke infinitynormen istedenfor
    #else
    condition_numb = 1
    
  
    if save
        writevtk(Ω, "C:\\Users\\Sigri\\Documents\\Master\\report\\results\\stokes\\$n $geometry $order.vtu", cellfields=["u_ex" => u_exact, "uh"=>uh, "erru"=> erru, "p_ex" => p_exact, "ph"=>ph, "errp"=> errp, "nablau" => ∇(u_exact), "viskositet" => viskositet∘ε(u_exact)]) #, "erru" => erru]) 
    end
    return uh, u_exact, erru, l2_norm(uh - u_exact), h1_semi(uh - u_exact), ph, p_exact, errp, l2_norm(ph - p_exact), h1_semi(ph - p_exact), condition_numb, Ω
end