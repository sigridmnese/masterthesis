#using STLCutters
using Gridap
using GridapEmbedded
#using Test
#using Plots
using LaTeXStrings
using Printf

#Domain 
#Bilin forms 
#R.H.s

function Stokes(;f,u_D, domain = "circle",n=20, γ = 10, β_1=10, β_2=10, β_3=10)
    if domain == "flower"
        p0 =Point(0,0)
        LL= 2.1
        r0, r1 = 0.35*LL, 0.09*LL
        
        #defining flower function
        flower(x) =  -r0 - r1*cos(5.0*atan(x[1], x[2])) +(x[1]^2 + x[2]^2)^0.5
        
        # using ! operator to define the interior
        geo = AnalyticalGeometry(x-> flower(x))
        
    else
        domain = "circle"
        R = 0.712
        # Disc center
        p0 = Point(0.0, 0.0)
        geo = disk(R, x0=p0)
    end

    partition = (n, n)
    pmin = p0 - 1.0
    pmax = p0 + 1.0
    bgmodel = CartesianDiscreteModel(pmin, pmax, partition)

    # Cut geometry
    cutgeo = cut(bgmodel, geo)
    cutgeo_facets = cut_facets(bgmodel, geo)
    Ω_bg = Triangulation(bgmodel)
    Ω_act = Triangulation(cutgeo, ACTIVE)
    Ω = Triangulation(cutgeo, PHYSICAL)

    # Boundaries and normals
    Γd = EmbeddedBoundary(cutgeo)
    n_Γd = get_normal_vector(Γd)

    Γs = BoundaryTriangulation(cutgeo_facets, ACTIVE_IN)
    n_Γs = get_normal_vector(Γs)


    Fg = GhostSkeleton(cutgeo)
    n_Fg = get_normal_vector(Fg)

    Fi = SkeletonTriangulation(Ω_act)
    n_Fi = get_normal_vector(Fi)

    #writevtk(Ω_act, "CutFEMimplementation/data/mesh/mesh_act_$domain")
    #writevtk(Ω_bg, "CutFEMimplementation/data/mesh/mesh_bg_$domain")
    #writevtk(Ω, "CutFEMimplementation/data/mesh/mesh_$domain")
    #writevtk(Γd, "CutFEMimplementation/data/mesh/surface_gamma_d_$domain")
    #writevtk(Fg, "CutFEMimplementation/data/mesh/ghost_facets_$domain")

    dim = 2
    order_u = 2
    order_p = 1

      # Define measures
      degree = 2 * order_u
      dΩ = Measure(Ω, degree)
      dΓd = Measure(Γd, degree)
      dFg = Measure(Fg, degree)
      dΓs = Measure(Γs, degree)  # Define measures
      dFi = Measure(Fi, degree)
   
    # Define function spaces
    reffe_u = ReferenceFE(lagrangian, VectorValue{dim, Float64}, order_u)
    reffe_p = ReferenceFE(lagrangian, Float64, order_p)

    V = TestFESpace(Ω_act, reffe_u, conformity=:H1)
    Q = TestFESpace(Ω_act, reffe_p, conformity=:H1, constraint=:zeromean)

    
    U = TrialFESpace(V)
    P = TrialFESpace(Q)
    X = MultiFieldFESpace([U, P])
    Y = MultiFieldFESpace([V, Q])

    # Mesh size
    h = 2/n
    #h = (pmax - pmin)[1] / n


    function jump_nn(u,n)
        return ( n.plus⋅ (n.plus⋅∇∇(u).plus) - n.minus ⋅ (n.minus ⋅ ∇∇(u).minus))
    end
    
        # Weak formulation components
    a(u, v) = (∫(∇(u) ⊙ ∇(v))dΩ +
                ∫(-((n_Γd ⋅ ∇(u)) ⋅ v) - ((n_Γd ⋅ ∇(v)) ⋅ u)+(γ/h * u ⋅ v))dΓd )
    b(v, p) = (∫(-1*(∇ ⋅ v*p))dΩ
                + ∫((n_Γd ⋅ v) * p)dΓd)

    # i(u, v) = β_1*∫(h*jump(n_Fg ⋅ ∇(u)) ⋅ jump(n_Fg ⋅ ∇(v)))dFg
    #         + β_2*∫(h^3*jump_nn(u,n_Fg)⊙jump_nn(v,n_Fg))dFg # missing parameter?¨
    # i(u, v) = (β_1*∫(h*jump(n_Fg ⋅ ∇(u)) ⋅ jump(n_Fg ⋅ ∇(v)))dFg
    # + β_2*∫(h^3*jump_nn(u,n_Fg)⋅jump_nn(v,n_Fg))dFg) # missing parameter?

    i(u,v) = ( ∫( (β_1*h)*jump(n_Fg⋅∇(u))⋅jump(n_Fg⋅∇(v)) )dFg 
            +  
               ∫( (β_2*h^3)*jump_nn(u,n_Fg)⋅jump_nn(v,n_Fg) )dFg)


    j(p, q) = (∫((β_3*h^3)*jump(n_Fg ⋅ ∇(p)) * jump(n_Fg ⋅ ∇(q)))dFg)
    #s_p(p,q) = ∫( (β_3*h^3)*jump(n_Fi⋅∇(p))*jump(n_Fi⋅∇(q)) )dFi

    Ah((u, p), (v, q)) = (a(u, v) + b(v, p) + b(u, q) 
                        + i(u,v)
                        - j(p, q)
                        )

    lh((v, q)) = (∫(f ⋅ v)dΩ 
                    + ∫(γ /h *(v ⋅ u_D))dΓd # Nitsche term
                    - ∫((n_Γd ⋅ ∇(v)) ⋅ u_D)dΓd # symetri
                    + ∫((n_Γd ⋅ u_D)*q)dΓd) # b(u,q) term

    op = AffineFEOperator(Ah, lh, X, Y)
    uh, ph = solve(op)
    
    return(uh,ph,dΩ,dΓd,n_Γd,Ω)
end 

function slope(h::Vector{Float64}, E::Vector{Float64})
    # Ensure the arrays have the same length and at least two elements
    if length(E) != length(h)
        throw(ArgumentError("Arrays E and h must be of the same length."))
    elseif length(E) < 2
        throw(ArgumentError("Arrays E and h must have at least two elements."))
    end
    
    # Compute EOC for each k and store in an array
    eoc_values = [log(E[i-1] / E[i]) / log(h[i-1] / h[i]) for i in 2:length(E)]
    # Return the average EOC
    return eoc_values
end

# function save_to_latex(filename, ns, errors_u_l2, eoc_u_l2, errors_u_h1, eoc_u_h1 ,errors_p_l2, eoc_p_l2)
#     open(filename, "w") do file
#         for i in 1:length(ns)
#             # Format the errors in scientific notation with 3 significant digits
#             error_u_l2_str = @sprintf("%.3e", errors_u_l2[i])  # x.xxE±x
#             error_u_h1_str = @sprintf("%.3e", errors_u_h1[i])  # x.xxE±x
#             error_p_l2_str = @sprintf("%.3e", errors_p_l2[i])  # x.xxE±x
            
#             # Convert the scientific notation from `e` to `\times 10^` for LaTeX
#             error_u_l2_latex = replace(error_u_l2_str, "e" => "\\cdot 10^{") * "}"
#             error_u_h1_latex = replace(error_u_h1_str, "e" => "\\cdot 10^{") * "}"
#             error_p_l2_latex = replace(error_p_l2_str, "e" => "\\cdot 10^{") * "}"

#             # Get EOC values from the array (eoc arrays have one less element)
#             eoc_u_l2_value = i == 1 ? "-" : string(round(eoc_u_l2[i-1], sigdigits=3))
#             eoc_u_h1_value = i == 1 ? "-" : string(round(eoc_u_h1[i-1], sigdigits=3))
#             eoc_p_l2_value = i == 1 ? "-" : string(round(eoc_p_l2[i-1], sigdigits=3))
#             # Create the LaTeX row
#             row = @sprintf("%d & %s & %s & %s & %s & %s & %s \\\\", 
#                 ns[i], 
#                 error_u_l2_latex, eoc_u_l2_value, 
#                 error_u_h1_latex, eoc_u_h1_value,
#                 error_p_l2_latex, eoc_p_l2_value
#             )
            
#             # Write the row to the file
#             println(file, row)
#         end
#     end
# end

function errors(;ns, domain = "circle",γ =10, β_1=0.1, β_2=0.01, β_3=0.1)

    
u_ex(x) = VectorValue(2*x[1] + cos(2*π*x[2]), -2*x[2] + sin(2*π*x[1]))
p_ex(x) = sin(2*π*x[1])
# Forcing term
#f(x)= -Δ(u_ex)(x)+ ∇(p_ex)(x)
f(x) = -divergence(∇(u_ex))(x) + ∇(p_ex)(x)
u_D(x) = u_ex(x)

    e_u_l2_norm_vec = Float64[]

    e_u_h1_semi_norm_vec = Float64[]

    e_p_l2_norm_vec = Float64[]
 
    hs =Float64[]

    
    for n in ns
        if domain == "flower"
            uh,ph,dΩ,dΓd,n_Γd = Stokes(f =f ,u_D=u_D, domain = "flower",n=n,γ = γ, β_1=β_1, β_2=β_2, β_3=β_3) 
        else
            uh,ph,dΩ,dΓd,n_Γd = Stokes(f =f ,u_D=u_D, domain = "circle",n=n,γ = γ, β_1=β_1, β_2=β_2, β_3=β_3)
        end

        h = 2/n 

        println(n)
        l2_norm(e) = (sum( ∫( e ⋅ e )*dΩ ))
        h1_semi_norm(e) = sum(∫(∇(e) ⊙ ∇(e))*dΩ)
  
        e_p = ph-p_ex
        e_u = uh-u_ex

        push!(e_u_l2_norm_vec, sqrt(l2_norm(e_u)))

        push!(e_u_h1_semi_norm_vec, sqrt(h1_semi_norm(e_u)))

        push!(e_p_l2_norm_vec,  sqrt(l2_norm(e_p)))
      
        push!(hs, h)
    end 
    return (e_u_l2_norm_vec, e_u_h1_semi_norm_vec, e_p_l2_norm_vec, hs)
end

#parameters 
ns = [16,32,64,128,256]
β_1 = 1
β_2 = 1
β_3 = 0.1

u_ex(x) = VectorValue(2 * x[1] + cos(2 * π * x[2]), -2 * x[2] + sin(2 * π *x[1]))
p_ex(x) = sin(2 * π *x[1])
f(x) = -divergence(∇(u_ex))(x) + ∇(p_ex)(x)
u_D(x) = u_ex(x)
n = 32
γ = 100
uh, ph, dΩ, dΓd, n_Γd, Ω = Stokes(f = f, u_D = u_D, domain = "circle", n=n, γ = γ, β_1 = β_1, β_2 = β_2, β_3 = β_3)
writevtk(Ω, "Stokes", cellfields=["uh"=>uh, "u_ex"=>u_ex, "ph"=> ph, "p_ex" => p_ex]) #, "erru" => erru])

"""
errors_u_l2, errors_u_h1, errors_p_l2, hs= errors(ns= ns,domain = "flower" ,γ = γ , β_1 = β_1, β_2 = β_2, β_3 = β_3)

eoc_u_l2 = slope(hs, errors_u_l2)
eoc_u_h1 = slope(hs, errors_u_h1)
eoc_p_l2 = slope(hs, errors_p_l2)

print(errors_u_l2)
print(errors_u_h1)
print(errors_p_l2)

println("β_1, β_2, β_3: ", β_1," ", β_2, " ", β_3)
println(eoc_u_l2)
println(eoc_u_h1)
println(eoc_p_l2)

save_to_latex("convergence_rates_Stokes_flower.txt", ns, errors_u_l2, eoc_u_l2, errors_u_h1, eoc_u_h1 ,errors_p_l2, eoc_p_l2)
 """

function Plain_stokes() 
    """#    p0 = Point(0.0, 0.0)
    #     model = disk(R, x0=p0)
    
    #     reffe_u = ReferenceFE(lagrangian, VectorValue{dim, Float64}, order_u)
    #     reffe_p = ReferenceFE(lagrangian, Float64, order_p)
    
    #     V = TestFESpace(Ω_act, reffe_u, conformity=:H1)
    #     Q = TestFESpace(Ω_act, reffe_p, conformity=:L2, constraint=:zeromean)
    
        
    #     U = TrialFESpace(V)
    #     P = TrialFESpace(Q)
    #     X = MultiFieldFESpace([U, P])
    #     Y = MultiFieldFESpace([V, Q])
    
    #     # Define measures
    #     degree = 2 * order_u
    #     Ω = Triangulation(model)
    #     dΩ = Measure(Ω, degree)
    
    #     # Mesh size
    #     h = (pmax - pmin)[1] / n
    #     a(u, v) = (-∫(∇(u) ⊙ ∇(v))dΩ)
        
    #         b(v, p) = (∫(-(∇ ⋅ v) * p)dΩ)
                        
        
    #         Ah((u, p), (v, q)) = (a(u, v) + b(v, p) + b(u, q) 
    #                             + i(u,v)
    #                             - j(p, q))
        
    #         lh((v, q)) = (∫(f ⋅ v)dΩ 
    #                         + ∫(γ * h^(-1) * v ⋅ u_D)dΓd # Nitsche term
    #                         - ∫((n_Γd ⋅ ∇(v)) ⋅ u_D)dΓd # symetri
    #                         + ∫((n_Γd ⋅ u_D)*q)dΓd # b(u,q) term 
    #                         )
        
    #         op = AffineFEOperator(Ah, lh, X, Y)
    #         uh, ph = solve(op)
            
             return(uh,ph,dΩ,dΓd,n_Γd)"""
    end
    