outputfolderpoisson = "C:\\Users\\Sigri\\Documents\\Master\\report\\results\\poisson"
outputfolderstokes = "C:\\Users\\Sigri\\Documents\\Master\\report\\results\\stokes"
  
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

function plot_convergence_u(uarr_l2, uarr_h1, h;
                            uarr_l2_nostab=nothing, uarr_h1_nostab=nothing,
                            title_str="Convergence of p-stokes FEM")
    plot(
        0,
        titlefont = 16,
        guidefont = 14,
        tickfont = 12
    )

    has_nostab = !(uarr_l2_nostab === nothing || uarr_h1_nostab === nothing)

    label_l2 = has_nostab ? "L² stabilized" : "L²"
    label_h1 = has_nostab ? "H¹ stabilized" : "H¹"

    plot!(h, uarr_l2, xaxis=:log, yaxis=:log, marker=:o, lw=2, label=label_l2)
    plot!(h, uarr_h1, marker=:o, lw=2, label=label_h1)

    if has_nostab
        plot!(h, uarr_l2_nostab, marker=:x, lw=2, label="L² unstabilized")
        plot!(h, uarr_h1_nostab, marker=:x, lw=2, label="H¹ unstabilized")
    end

    xlabel!("Mesh size h")
    ylabel!("Velocity error")
    title!(title_str)
end

function plot_convergence_p(parr_l2, parr_h1, h;
                            parr_l2_nostab=nothing, parr_h1_nostab=nothing,
                            title_str="Convergence of p-stokes FEM")
    plot(
        0,
        titlefont = 16,
        guidefont = 14,
        tickfont = 12
    )

    has_nostab = !(parr_l2_nostab === nothing || parr_h1_nostab === nothing)

    label_l2 = has_nostab ? "L² stabilized" : "L²"
    label_h1 = has_nostab ? "H¹ stabilized" : "H¹"

    plot!(h, parr_l2, xaxis=:log, yaxis=:log, marker=:o, lw=2, label=label_l2)
    plot!(h, parr_h1, marker=:o, lw=2, label=label_h1)

    if has_nostab
        plot!(h, parr_l2_nostab, marker=:x, lw=2, label="L² unstabilized")
        plot!(h, parr_h1_nostab, marker=:x, lw=2, label="H¹ unstabilized")
    end

    xlabel!("Mesh size h")
    ylabel!("Pressure error")
    title!(title_str)
end

using Plots

function plot_sensitivity_velocity(
    δ_stab::Vector{<:Real}, l2_stab::Vector{<:Real}, h1_stab::Vector{<:Real},
    δ_nostab::Vector{<:Real}, l2_nostab::Vector{<:Real}, h1_nostab::Vector{<:Real};
    start::Int = 1,
    title_str::String = "Sensitivity CutFEM: Stokes with Navier BC"
)
    plot(
        xlabel = "Perturbation δ",
        ylabel = "Velocity error",
        title = title_str,
        yaxis = :log,
        lw = 2,
        titlefont = 16,
        guidefont = 14,
        tickfont = 12,
        legend = :topright
    )

    plot!(δ_stab[start:end], l2_stab[start:end], label = "L² stabilized")
    plot!(δ_nostab[start:end], l2_nostab[start:end], label = "L² non-stabilized")
    plot!(δ_stab[start:end], h1_stab[start:end], label = "H¹ stabilized")#, marker = :diamond)
    plot!(δ_nostab[start:end], h1_nostab[start:end], label = "H¹ non-stabilized")#, marker = :diamond)
end

function plot_sensitivity_pressure(
    δ_stab::Vector{<:Real}, l2_stab::Vector{<:Real}, h1_stab::Vector{<:Real},
    δ_nostab::Vector{<:Real}, l2_nostab::Vector{<:Real}, h1_nostab::Vector{<:Real};
    start::Int = 1,
    title_str::String = "Sensitivity CutFEM: Stokes with Navier BC"
)
    plot(
        xlabel = "Perturbation δ",
        ylabel = "Pressure error",
        title = title_str,
        yaxis = :log,
        lw = 2,
        titlefont = 16,
        guidefont = 14,
        tickfont = 12,
        legend = :topright
    )

    plot!(δ_stab[start:end], l2_stab[start:end], label = "L² stabilized")#, marker = :circle)
    plot!(δ_nostab[start:end], l2_nostab[start:end], label = "L² non-stabilized")#, marker = :circle)
    plot!(δ_stab[start:end], h1_stab[start:end], label = "H¹ stabilized")#, marker = :diamond)
    plot!(δ_nostab[start:end], h1_nostab[start:end], label = "H¹ non-stabilized")#, marker = :diamond)
end
function plot_sensitivity_poisson(
    δ_stab::Vector{<:Real}, l2_stab::Vector{<:Real}, h1_stab::Vector{<:Real},
    δ_nostab::Vector{<:Real}, l2_nostab::Vector{<:Real}, h1_nostab::Vector{<:Real};
    start::Int = 1,
    title_str::String = "Sensitivity CutFEM: Poisson",
    marker_l2_stab = :circle,
    marker_l2_nostab = :diamond,
    marker_h1_stab = :utriangle,
    marker_h1_nostab = :square,
    markstep::Int = 50
)
    function sparse_markers(x, y, step)
        x_marker = [i % step == 0 ? xi : missing for (i, xi) in enumerate(x)]
        y_marker = [i % step == 0 ? yi : missing for (i, yi) in enumerate(y)]
        return x_marker, y_marker
    end

    plot(
        xlabel = "Perturbation δ",
        ylabel = "Error",
        title = title_str,
        yaxis = :log,
        lw = 2,
        titlefont = 16,
        guidefont = 14,
        tickfont = 12,
        legend = :topright,
        markersize = 6
    )

    # L² stabilized
    plot!(δ_stab[start:end], l2_stab[start:end], label = "L² stabilized", marker = :none)
    if marker_l2_stab !== nothing
        x, y = sparse_markers(δ_stab[start:end], l2_stab[start:end], markstep)
        scatter!(x, y, label = "", marker = marker_l2_stab)
    end

    # L² non-stabilized
    plot!(δ_nostab[start:end], l2_nostab[start:end], label = "L² non-stabilized", marker = :none)
    if marker_l2_nostab !== nothing
        x, y = sparse_markers(δ_nostab[start:end], l2_nostab[start:end], markstep)
        scatter!(x, y, label = "", marker = marker_l2_nostab)
    end

    # H¹ stabilized
    plot!(δ_stab[start:end], h1_stab[start:end], label = "H¹ stabilized", marker = :none)
    if marker_h1_stab !== nothing
        x, y = sparse_markers(δ_stab[start:end], h1_stab[start:end], markstep)
        scatter!(x, y, label = "", marker = marker_h1_stab)
    end

    # H¹ non-stabilized
    plot!(δ_nostab[start:end], h1_nostab[start:end], label = "H¹ non-stabilized", marker = :none)
    if marker_h1_nostab !== nothing
        x, y = sparse_markers(δ_nostab[start:end], h1_nostab[start:end], markstep)
        scatter!(x, y, label = "", marker = marker_h1_nostab)
    end
end

function plot_condition_poisson(
    δ_stab::Vector{<:Real}, cond_stab::Vector{<:Real},
    δ_nostab::Vector{<:Real}, cond_nostab::Vector{<:Real};
    start::Int = 1,
    title_str::String = "Condition number: CutFEM Poisson",
    marker_stab = :circle,
    marker_nostab = :diamond,
    markstep::Int = 50
)
    function sparse_markers(x, y, step)
        x_marker = [i % step == 0 ? xi : missing for (i, xi) in enumerate(x)]
        y_marker = [i % step == 0 ? yi : missing for (i, yi) in enumerate(y)]
        return x_marker, y_marker
    end

    plot(
        xlabel = "Perturbation δ",
        ylabel = "Condition number",
        title = title_str,
        yaxis = :log,
        lw = 2,
        titlefont = 16,
        guidefont = 14,
        tickfont = 12,
        legend = :topleft,
        markersize = 6
    )

    plot!(δ_stab[start:end], cond_stab[start:end], label = "Stabilized", marker = :none)
    if marker_stab !== nothing
        x, y = sparse_markers(δ_stab[start:end], cond_stab[start:end], markstep)
        scatter!(x, y, label = "", marker = marker_stab)
    end

    plot!(δ_nostab[start:end], cond_nostab[start:end], label = "Non-stabilized", marker = :none)
    if marker_nostab !== nothing
        x, y = sparse_markers(δ_nostab[start:end], cond_nostab[start:end], markstep)
        scatter!(x, y, label = "", marker = marker_nostab)
    end
end



function compute_eoc(h, error)
    return [log(error[i+1]/error[i]) / log(h[i+1]/h[i]) for i in 1:length(h)-1]
end

function plot_eoc(h, error; label="EOC", color=:blue)
    eocs = compute_eoc(h, error)
    h_mid = [sqrt(h[i]*h[i+1]) for i in 1:length(eocs)]  # log-midtpunkt
    plot(h_mid, eocs, xaxis=:log, marker=:circle, lw=2, label=label, xlabel="Mesh size h", ylabel="EOC", color=color)
end


function plot_convergence_generic(
    h::Vector{<:Real},
    u_l2_stab::Vector{<:Real}, u_h1_stab::Vector{<:Real};
    u_l2_nostab::Union{Nothing, Vector{<:Real}} = nothing,
    u_h1_nostab::Union{Nothing, Vector{<:Real}} = nothing,
    start::Int = 1,
    title_str::String = "Convergence of CutFEM Stokes Solver",
    annotate_mode::Symbol = :last # :last eller :all
)
    plot(
        xlabel = "Mesh size h",
        ylabel = "Velocity error",
        title = title_str,
        xaxis = :log,
        yaxis = :log,
        markerstrokewidth = 0,
        titlefont = 16,
        guidefont = 14,
        tickfont = 12,
        legend = :bottomleft
    )

    hplot = h[start:end]
    curves = [
        (u_l2_stab[start:end], "L2 stabilized", :o),
        (u_h1_stab[start:end], "H1 stabilized", :o)
    ]
    if u_l2_nostab !== nothing && u_h1_nostab !== nothing
        push!(curves, (u_l2_nostab[start:end], "L2 non-stabilized", :s))
        push!(curves, (u_h1_nostab[start:end], "H1 non-stabilized", :s))
    end

    for (yvals, labelname, markertype) in curves
        plot!(hplot, yvals, marker = markertype, lw = 2, label = labelname)

        eocs = compute_eoc(hplot, yvals)

        if annotate_mode == :last
            xpos = sqrt(hplot[end-1]*hplot[end])
            ypos = sqrt(yvals[end-1]*yvals[end])
            annotate!(xpos, ypos*0.1, text("slope ≈ $(round(eocs[end], digits=2))", 9))
        elseif annotate_mode == :all
            for i = 1:length(eocs)
                xpos = sqrt(hplot[i]*hplot[i+1])
                ypos = sqrt(yvals[i]*yvals[i+1]) * 1.1^(length(curves) - 1) * 0.95^i
                annotate!(xpos, ypos, text("$(round(eocs[i], digits=2))", 8))
            end
        end
    end
    return current()
end
  # Define ghost penalty
function jump_nn(u,n)
    return ( n.plus ⋅ (n.plus⋅∇∇(u).plus) - n.minus ⋅ (n.minus ⋅ ∇∇(u).minus) )       # andre ordens hopp... Forklare dette skikkelig. 
end

function jump_nn(u,n)
    return ( (n.plus⋅∇∇(u).plus)⋅ n.plus - (n.minus ⋅ ∇∇(u).minus) ⋅ n.minus )
end

function jump_nn_symm(u,n)
    εu_minus = ε(u).minus
    εεu_minus = ∇(εu_minus)
    εu_plus = ε(u).plus
    εεu_plus = ∇(εu_plus)

    return ( n.plus ⋅ (n.plus⋅εεu_plus) - n.minus ⋅ (n.minus ⋅ εεu_minus))       # andre ordens hopp... Forklare dette skikkelig. 
end

function print_eoc_latex_combined_pressure(h::Vector{<:Real};
    parr_l2_stab::Union{Nothing, Vector{<:Real}} = nothing,
    parr_h1_stab::Union{Nothing, Vector{<:Real}} = nothing,
    parr_l2_nostab::Union{Nothing, Vector{<:Real}} = nothing,
    parr_h1_nostab::Union{Nothing, Vector{<:Real}} = nothing,
    start::Int = 1
)
    hsub = h[start:end]
    n_rows = length(hsub) - 1

    datasets = []
    if parr_l2_stab !== nothing
        push!(datasets, ("\\(L^2\\) stab", parr_l2_stab[start:end]))
    end
    if parr_h1_stab !== nothing
        push!(datasets, ("\\(H^1\\) stab", parr_h1_stab[start:end]))
    end
    if parr_l2_nostab !== nothing
        push!(datasets, ("\\(L^2\\) no stab", parr_l2_nostab[start:end]))
    end
    if parr_h1_nostab !== nothing
        push!(datasets, ("\\(H^1\\) no stab", parr_h1_nostab[start:end]))
    end

    n_cols = 1 + 2 * length(datasets)
    println("\\begin{tabular}{", "c"^n_cols, "}")
    println("\\toprule")

    print("\$h_1 \\rightarrow h_2\$")
    for (label, _) in datasets
        print(" & \\multicolumn{2}{c}{", label, "} ")
    end
    println(" \\\\")

    print(" ")
    for _ in 1:length(datasets)
        print(" & Error & EOC")
    end
    println(" \\\\")
    println("\\midrule")

    for i in 1:n_rows
        k1 = round(Int, log2(1/hsub[i]))
        k2 = round(Int, log2(1/hsub[i+1]))
        print("\$ 2^{-$k1} \\rightarrow 2^{-$k2} \$")
        for (_, errors) in datasets
            err = round(errors[i+1], sigdigits=3)
            eoc = round(log(errors[i+1]/errors[i]) / log(hsub[i+1]/hsub[i]), digits=2)
            print(" & \$ $err \$ & \$ $eoc \$")
        end
        println(" \\\\")
    end

    println("\\bottomrule")
    println("\\end{tabular}")
end


function print_eoc_latex_combined(h::Vector{<:Real};
    uarr_l2_stab::Union{Nothing, Vector{<:Real}} = nothing,
    uarr_h1_stab::Union{Nothing, Vector{<:Real}} = nothing,
    uarr_l2_nostab::Union{Nothing, Vector{<:Real}} = nothing,
    uarr_h1_nostab::Union{Nothing, Vector{<:Real}} = nothing,
    start::Int = 1
)
    hsub = h[start:end]
    n_rows = length(hsub) - 1

    # Samle alle datasett
    datasets = []
    if uarr_l2_stab !== nothing
        push!(datasets, ("\\(L^2\\) stab", uarr_l2_stab[start:end]))
    end
    if uarr_h1_stab !== nothing
        push!(datasets, ("\\(H^1\\) stab", uarr_h1_stab[start:end]))
    end
    if uarr_l2_nostab !== nothing
        push!(datasets, ("\\(L^2\\) no stab", uarr_l2_nostab[start:end]))
    end
    if uarr_h1_nostab !== nothing
        push!(datasets, ("\\(H^1\\) no stab", uarr_h1_nostab[start:end]))
    end

    n_cols = 1 + 2 * length(datasets)
    println("\\begin{tabular}{", "c"^n_cols, "}")
    println("\\toprule")

    print("\$h_1 \\rightarrow h_2\$")
    for (label, _) in datasets
        print(" & \\multicolumn{2}{c}{", label, "} ")
    end
    println(" \\\\")

    print(" ")
    for _ in 1:length(datasets)
        print(" & Error & EOC")
    end
    println(" \\\\")
    println("\\midrule")

    for i in 1:n_rows
        k1 = round(Int, log2(1/hsub[i]))
        k2 = round(Int, log2(1/hsub[i+1]))
        print("\$ 2^{-$k1} \\rightarrow 2^{-$k2} \$")
        for (_, errors) in datasets
            err = round(errors[i+1], sigdigits=3)
            eoc = round(log(errors[i+1]/errors[i]) / log(hsub[i+1]/hsub[i]), digits=2)
            print(" & \$ $err \$ & \$ $eoc \$")
        end
        println(" \\\\")
    end

    println("\\bottomrule")
    println("\\end{tabular}")
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
    l2(u) = √(∑( ∫( u*u )dΩ ))^(1/2)
    h1(u) = √(∑( ∫( u*u + ∇(u)⋅∇(u) )dΩ ))^(1/2)

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
      Symmetric gradient. 
      Constant viscosity
      Working?
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
     
    #   gu(u,v) = ( ∫( (β_1*h)*jump(n_Fg ⋅ ∇(u))⋅jump(n_Fg⋅ ∇(v)) )dFg 
    #           +  
    #              ∫( (β_2*h^3)*jump_nn(u,n_Fg)⋅jump_nn(v,n_Fg) )dFg)
  
  
    #   gp(p, q) = (∫((β_3*h^3)*jump(n_Fg ⋅ ∇(p)) * jump(n_Fg ⋅ ∇(q)))dFg)

    gu(u,v) = ( ∫( (γu1*h)*jump(n_Fg⋅∇(u))⋅jump(n_Fg⋅∇(v)) )dFg 
            +  
               ∫( (γu2*h^3)*jump_nn(u,n_Fg)⋅jump_nn(v,n_Fg) )dFg
)

    gp(p,q) = ∫( (γp*h^3)*jump(n_Fg⋅∇(p))*jump(n_Fg⋅∇(q)) )dFg
  
      l2_norm(u) = (sum( ∫( u ⋅ u )*dΩ ))^(1/2)
      h1_semi(u) = sum(∫(∇(u) ⊙ ∇(u))*dΩ)^(1/2)
      
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
          B((u,p),(v,q)) = (a(u,v) + b(v, p) + b(u, q) )
          M((v,q)) = l1(v) + l2(v) + l3(v) + l4(q)
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

function stokes_CutFEM(;n, u_exact, p_exact, f, g, ud, order, geometry, βu0, γu1, γu2, γp, βp0, nu, stabilize, δ, save = false, calc_condition = false)
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
  
      l2_norm(u) = (sum( ∫( u ⋅ u )*dΩ ))^(1/2)
      h1_semi(u) = sum(∫(∇(u) ⊙ ∇(u))*dΩ)^(1/2)
      
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

########################### p stokes solver ###########################
# fitted FEM p-stokes solver:
function p_stokes_FEM(;n, u_exact, p_exact, f, g, ud, order, geometry, βu0, γu1, γu2, γp, βp0, nu, stabilize, δ, save = false, calc_condition = false)
    """
    Denne fungerer veldig bra. Tror ikke den fungerer for variasjoner i nu0, som gir store variasjoner i viskositeten generelt
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
    
    a(u, v) = ∫( ε(v)⊙(flux∘ε(u)))dΩ
    b(v, p) = ∫(-(∇⋅v)*p )dΩ
    l(v) = ∫(f ⋅ v)dΩ

    # dflux calculated the same way as in the notebook p-Laplace...
    # and introduced in the same way as in the notebook p-Laplace in the bilinear form a...
    da(u, du, v) = ∫(ε(v)⊙(dflux∘(ε(du), ε(u))))dΩ
    
    # and then the Newton multifield system is assembled as in the Navier Stokes notebook...
    res((u,p),(v,q)) = a(u, v) + b(v, p) - b(u, q) - l(v)
    jac((u, p), (du, dp), (v, q)) =  b(v, dp) - b(du, q) + da(u, du, v)

    op = FEOperator(res, jac, X, Y)

    # non-linear phase
    nls = NLSolver(
    show_trace=true, method=:newton, linesearch=BackTracking(), iterations=100)      #prøver å legge inn et max antall iterasjoner og en lav toleranse      
    solver = FESolver(nls)

    (uh, ph) = solve(solver, op)

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

# p stokes cut FEM fungerer veldig bra, men må få lagt inn rikig viskositetsavhengighet.
function p_stokes_cutFEM_symmetric(;n, u_exact, p_exact, f, g, ud, order, geometry, βu0, γu1, γu2, γp, βp0, nu, stabilize, δ, save = false, calc_condition = false)
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
    γ = nu0*order*order
    #γ(∇u) = viskositet(∇u)
    #γ(u) = viskositet∘ε(u)

    # weak formulation components    
    # a(u, v) = ∫( ε(v)⊙(flux∘ε(u)))dΩ  + ∫(-((n_Γd ⋅ (flux∘ε(u))) ⋅ v) + (-(n_Γd ⋅ (flux∘ε(v))) ⋅ u) + viskositet∘ε(u)/h * (u ⋅ v))dΓd      # denne må ha et ekstra boundary term. Finn ut hvordan det ser ut. 
    # b(v, p) = (∫(-1*(∇ ⋅ v*p))dΩ + ∫((n_Γd ⋅ v) * p)dΓd)   # b er den samme som før. 
    # l1((v, q)) = ∫(f ⋅ v)dΩ
    # l2(v) = ∫(-(n_Γd ⋅ (flux∘ε(v))) ⋅ ud)dΓd    # har brukt ud som dirichlet grense
    # l3((u, p), (v, q)) = ∫( viskositet∘ε(u)/h* (ud ⋅ v))dΓd
    # l4(q) = ∫((n_Γd ⋅ ud) * q)dΓd               # kan det være at denne skal være annerledes? mtp flux * ... og ikke n_Γd ⋅ ud?
    
    # da(u, du, v) = ∫( ε(v)⊙(dflux∘(ε(du), ε(u))))dΩ  + ∫(-((n_Γd ⋅ (dflux∘(ε(du), ε(u)))) ⋅ v) + (-(n_Γd ⋅ (flux∘ε(v))) ⋅ du)+(viskositet∘ε(u)/h*(du ⋅ v)))dΓd       
    a(u, v) = ∫( ε(v)⊙(flux∘ε(u)))dΩ  + ∫(-((n_Γd ⋅ (flux∘ε(u))) ⋅ v) + (-(n_Γd ⋅ (flux∘ε(v))) ⋅ u) + γ/h * (u ⋅ v))dΓd      # denne må ha et ekstra boundary term. Finn ut hvordan det ser ut. 
    b(v, p) = (∫(-1*(∇ ⋅ v*p))dΩ + ∫((n_Γd ⋅ v) * p)dΓd)   # b er den samme som før. 
    l1((v, q)) = ∫(f ⋅ v - g ⋅ q)dΩ
    l2(v) = ∫(-(n_Γd ⋅ (flux∘ε(v))) ⋅ ud)dΓd    
    l3((u, p), (v, q)) = ∫( γ/h* (ud ⋅ v))dΓd
    l4(q) = ∫((n_Γd ⋅ ud) * q)dΓd
    
    da(u, du, v) = ∫( ε(v)⊙(dflux∘(ε(u), ε(du))))dΩ  + ∫(-((n_Γd ⋅ (dflux∘(ε(u), ε(du)))) ⋅ v) + (-(n_Γd ⋅ (flux∘ε(v))) ⋅ du)+(γ/h*(du ⋅ v)))dΓd    

    gu((u,p), (v, q)) = (∫(( β_1 * h*γ)*jump(n_Fg ⋅ ε(u))⋅jump(n_Fg⋅ ε(v)) )dFg  # prøver å multiplisere med viskositet her
              +  
                ∫((β_2*h^3*γ)*jump_nn(u,n_Fg)⋅jump_nn(v,n_Fg) )dFg)

    gp((u,p),(v, q)) = (∫((β_3 * h^3/γ)*jump(n_Fg ⋅ ∇(p)) * jump(n_Fg ⋅ ∇(q)))dFg)  #eventuelt h^3 her, men synes h^1 fungerer tilsnelatende bedre...

    if stabilize # have tested with different combinations of calling the gp in the residual and jacobian, but this one seems to work best.
      res((u,p),(v,q)) = a(u, v) + b(v, p) + b(u, q) + gu((u,p), (v, q)) - gp((u,p), (v, q)) - l1((v, q)) -l2(v) - l3((u, p), (v, q)) -l4(q) 
      jac((u, p), (du, dp), (v, q)) = b(v, dp) + b(du, q) + da(u, du, v) + gu((du, p), (v, q)) - gp((u,dp), (v, q))
      
      op = FEOperator(res, jac, X, Y)

      # non-linear phase
      nls = NLSolver(
      show_trace=true, method=:newton, linesearch=BackTracking(), iterations=50)      #prøver å legge inn et max antall iterasjoner og en lav toleranse      
      solver = FESolver(nls)

      (uh, ph) = solve(solver, op)

    else
      res_nostab((u,p),(v,q)) = a(u, v) + b(v, p) + b(u, q) - l1(v) -l2(v) -l3(v) -l4(q)
      jac_nostab((u, p), (du, dp), (v, q)) = b(v, dp) + b(du, q) + da(u, du, v)
      
      op = FEOperator(res_nostab, jac_nostab, X, Y)

      # non-linear phase
      nls = NLSolver(
      show_trace=true, method=:newton, linesearch=BackTracking(), iterations=50)      #prøver å legge inn et max antall iterasjoner og en lav toleranse      
      solver = FESolver(nls)

      (uh, ph) = solve(solver, op)
    end

    errp = p_exact - ph
    erru = u_exact - uh
    
    condition_numb = 1
    if save
        writevtk(bgmodel, "C:\\Users\\Sigri\\Documents\\Master\\report\\figures\\Domenefigurer\\mesh_bg$geometry $δ.vtu")
        # writevtk(Fg, "C:\\Users\\Sigri\\Documents\\Master\\report\\figures\\Domenefigurer\\Fg$geometry $δ.vtu")
        # writevtk(Γd, "C:\\Users\\Sigri\\Documents\\Master\\report\\figures\\Domenefigurer\\surface_gamma_d_$geometry $δ.vtu")
        # writevtk(Ω_act, "C:\\Users\\Sigri\\Documents\\Master\\report\\figures\\Domenefigurer\\Ω_act$geometry $δ.vtu")
        # writevtk(Ω_act, "C:\\Users\\Sigri\\Documents\\Master\\report\\figures\\Domenefigurer\\Ω_act$geometry $δ.vtu")
        writevtk(Ω, "C:\\Users\\Sigri\\Documents\\Master\\report\\results\\stokes\\$n $geometry $order $δ.vtu", cellfields=["u_ex" => u_exact, "uh"=>uh, "erru"=> erru, "p_ex" => p_exact, "ph"=>ph, "errp"=> errp, "nablau" => ∇(u_exact), "flux" => flux∘ε(u_exact), "viskositet" => viskositet∘ε(u_exact)]) #, "erru" => erru]) 
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


function convergence_stokes_weird_domain(;numb_it, u_exact, p_exact, f, g, ud, order, geometry, solver, δ, βu0, γu1, γu2, γp, βp0, nu, stabilize, save = false)
  #"""function to calculate convergence of the poisson solver, or the stokes solver, with or without stabilization"""
  calc_condition = false 
  uarr_l2 = zeros(Float64, numb_it)
  uarr_h1 = zeros(Float64, numb_it)
  parr_l2 = zeros(Float64, numb_it)
  parr_h1 = zeros(Float64, numb_it)

  n_arr = 2 .^ (4:(numb_it + 3))
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

