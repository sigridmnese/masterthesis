using STLCutters
using Gridap
using GridapEmbedded
using Test
using Plots
using Printf
# Manufactured solution
u(x) = exp(-10*(x[1]^2 + 0.5*x[2]^2)) 
#(x[1]^2 + x[2]^2)*sin(x[2]*x[1])^3  #bruk sin + cos /exp(helst ikke)(for ikke få for god konv)
f(x) = - Δ(u)(x)
ud(x) = u(x)
 
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


# %% Define geometry
#R = 0.712
 
# Disc center
p0 = Point(0.0, 0.0)
#geo = disk(R, x0=p0)

# Define background mesh
for n in [64] # [16,32,64,128,256,512] #
    println(n)
    geo = create_geometry("heart", n)
 
    partition = (n, n)
    dim = length(partition)
    pmin = p0-1.11
    pmax = p0+1.11
    bgmodel = CartesianDiscreteModel(pmin,pmax,partition)
    # Cut the background model
 
    cutgeo = cut(bgmodel,geo)
 
 
    # Setup integration meshes
    Ω_bg = Triangulation(bgmodel)
    Ω = Triangulation(cutgeo,PHYSICAL)
 
    # Setup boundary and normals
    Γd = EmbeddedBoundary(cutgeo)
    n_Γd = get_normal_vector(Γd)
 
    # Get ghost penalty facets and normals
    Fg = GhostSkeleton(cutgeo)
    n_Fg = get_normal_vector(Fg)
 
    # Setup Lebesgue measures
    order = 1
    degree = 2*order
    dΩ = Measure(Ω,degree)
    dΓd = Measure(Γd,degree)
    dFg = Measure(Fg, degree)
 
    # Setup FESpace
    Ω_act = Triangulation(cutgeo,ACTIVE)
 
    V = FESpace(Ω_act,ReferenceFE(lagrangian,Float64,order),conformity=:H1)
    U = TrialFESpace(V)
 
    # Weak form
    γd = 10.0 * order*order
    γu1 = 0.1 #test, kjenn på effekten, se hva som skjer
    γu2 = 0.1
    h = (pmax - pmin)[1] / partition[1]
 
    a(u,v) =
    ∫( ∇(v)⋅∇(u) ) * dΩ  +
    ∫( (γd/h)*v*u  - v*(n_Γd⋅∇(u)) - (n_Γd⋅∇(v))*u ) * dΓd
   
    #legge til for p2
    function jump_nn(u,n)
        return ( (n.plus⋅∇∇(u).plus)⋅ n.plus - (n.minus ⋅ ∇∇(u).minus) ⋅ n.minus )
    end
   
    g_u(u,v) =  (∫( (γu1*h)*jump(n_Fg⋅∇(u))⋅jump(n_Fg⋅∇(v)) )dFg
               +  
                  ∫( (γu2*h^3)*jump_nn(u,n_Fg)⋅jump_nn(v,n_Fg) )dFg)
 
    l(v) = ∫( v*f ) * dΩ + ∫( (γd/h)*v*ud - (n_Γd⋅∇(v))*ud ) * dΓd
    A(u,v)= a(u,v) + g_u(u,v)
 
    # FE problem
    op = AffineFEOperator(A,l,U,V)
    uh = solve(op)
    
    e = u - uh
   
    outputfile = true

    if outputfile !== nothing
      writevtk(Ω,"hjerte_med_ghost",cellfields=[
                                        "uh"=>uh,
                                        "u_ex"=> u,
                                        "erru" => e
                                        ])
    end

    # Postprocess
    #l2_norm(e) = sqrt(sum( ∫( e*e )*dΩ ))
    #semi_norm(e) = sqrt(sum(∫( ∇(e)⋅∇(e) )*dΩ))
 
    #e_ghost = ghost_semi_norm(e)
    #e_semi= semi_norm(e)
    #e_l2 = l2_norm(e)
 
    #push!(e_semi_vec, e_semi)
    #push!(e_l2_vec, e_l2)
    #push!(hs, h)
end
 
# println("e_semi_vec: ", e_semi_vec)
# println("e_l2_vec:", e_l2_vec)"