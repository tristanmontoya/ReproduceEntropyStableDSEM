struct AdvectionPRefinementDriver{d, 
    ApproximationType <: AbstractApproximationType,
    ElemShape <: AbstractElemShape, ResidualForm <: AbstractResidualForm,
    Strategy <: AbstractStrategy, 
    OperatorAlgorithm <: AbstractOperatorAlgorithm,
    ODEAlgorithm <: OrdinaryDiffEqAlgorithm}
    p_min::Int
    p_max::Int
    l::Int
    C_t::Float64
    n_s::Int
    scheme::String
    element_type::ElemShape
    form::ResidualForm
    strategy::Strategy
    operator_algorithm::OperatorAlgorithm
    ode_algorithm::ODEAlgorithm
    path::String
    M0::Int
    λ::Float64
    L::Float64
    a::NTuple{d,Float64}
    T::Float64
    mesh_perturb::Float64
    load_from_file::Bool
    overwrite::Bool
    run::Bool
    spectral_radius::Bool
    eigensolver::String
    flops::Bool
    r::Int
    tol::Float64
end

function rhs_flops(
    reference_approximation::ReferenceApproximation{<:RefElemData{d}};
    L::Float64=1.0, M::Int=1, mesh_perturb=1.0/16.0,
    strategy::AbstractStrategy=ReferenceOperator(),
    operator_algorithm::AbstractOperatorAlgorithm=DefaultOperatorAlgorithm(),
    assume_orthonormal::Bool=true) where {d}

    conservation_law = LinearAdvectionEquation(Tuple(1.0 for m in 1:d))

    original_mesh = uniform_periodic_mesh(reference_approximation, 
        Tuple((0.0,L) for m in 1:d), Tuple(M for m in 1:d))
        
    mesh = warp_mesh(original_mesh, reference_approximation, 
        ChanWarping(mesh_perturb, Tuple(L for m in 1:d)))

    form = StandardForm(mapping_form=SkewSymmetricMapping(), 
    inviscid_numerical_flux=LaxFriedrichsNumericalFlux(1.0))

    spatial_discretization = SpatialDiscretization(mesh, 
        reference_approximation, project_jacobian=isa(reference_approximation,
        Union{ModalTensor,ModalMulti}))

    solver = Solver(conservation_law, spatial_discretization, form,
        strategy, operator_algorithm, WeightAdjustedSolver(
        spatial_discretization, operator_algorithm,
        assume_orthonormal=assume_orthonormal), Serial())

    (N_p, N_c, N_e) = size(solver)
    u = rand(N_p, N_c, N_e)
    dudt = similar(u)
    cnt = @count_ops StableSpectralElements.Analysis.rhs_benchmark!(
        dudt, u, solver, 0.0)
    return cnt
end

function AdvectionPRefinementDriver(p_min, p_max; l=2, C_t=0.1, n_s=50,
    element_type="Tri", scheme="ModalMulti",
    mapping_form="SkewSymmetricMapping", 
    strategy="ReferenceOperator",
    operator_algorithm="DefaultOperatorAlgorithm",
    ode_algorithm="CarpenterKennedy2N54", 
    path="../results/20230501_flops/", M0 = 2,
    λ=1.0, L=1.0, a = nothing, T = 1.0, mesh_perturb = 1.0/16.0,load_from_file=true, overwrite=false, run=true, spectral_radius=false, 
    eigensolver="KrylovKit",
    flops=false, r=2, tol=1.0e-2)

    if Int(round(λ)) == 0 
        path = string(path, scheme, "_", element_type, "_",
            mapping_form, "/central/")
    elseif Int(round(λ)) == 1 
        path = string(path, scheme, "_", element_type,
            "_", mapping_form, "/upwind/")
    else 
        path = string(path, scheme, "_", element_type, "_",
            mapping_form, "/lambda", replace(string(lambda), 
            "." => "_"), "/")
    end

    element_type = eval(Symbol(element_type))()
    form = StandardForm(mapping_form=eval(Symbol(mapping_form))(),
        inviscid_numerical_flux=LaxFriedrichsNumericalFlux(λ))
    operator_algorithm = eval(Symbol(operator_algorithm))()
    ode_algorithm = eval(Symbol(ode_algorithm))()
    strategy = eval(Symbol(strategy))()

    if isnothing(a)
        a = Tuple(1.0 for m in 1:dim(element_type))
    end

    return AdvectionPRefinementDriver(p_min, p_max,
        l,C_t,n_s,scheme,element_type,
        form,strategy, operator_algorithm, ode_algorithm,path,M0,λ,L,
        a, T, mesh_perturb, load_from_file, overwrite, run,
        spectral_radius, eigensolver, flops, r, tol)
end

function run_driver(driver::AdvectionPRefinementDriver{d}) where {d}

    @unpack p_min, p_max,l,C_t,n_s,scheme,element_type,form,strategy, operator_algorithm, ode_algorithm,path,M0,λ,L,a,T,mesh_perturb, load_from_file, overwrite, run, spectral_radius, eigensolver, flops, r, tol = driver

    if (!load_from_file || !isdir(path)) path = new_path(
        path,overwrite,overwrite) end
    if !isdir(string(path, "/p", p_min, "/")) p_start = p_min else
        for i in p_min:p_max
            if !isdir(string(path, "p", i + 1, "/"))
                p_start = i
                break
            end
        end
    end
    open(string(path,"screen.txt"), "a") do io
        println(io, "Starting refinement from p = ", p_start)
    end

    initial_data = InitialDataSine(1.0,Tuple(2*π/L for m in 1:d))

    conservation_law = LinearAdvectionEquation(a)

    for p in p_start:p_max
        
        if !isfile(string(path, "poly_degrees.jld2"))
            save_object(string(path, "poly_degrees.jld2"), Int64[])
        end

        poly_degrees=load_object(string(path, "poly_degrees.jld2"))
        save_object(string(path, "poly_degrees.jld2"),
            push!(poly_degrees, p))

        M = M0

        reference_approximation =ReferenceApproximation(
            eval(Symbol(scheme))(p), element_type, mapping_degree=l)

        original_mesh = uniform_periodic_mesh(reference_approximation, 
            Tuple((0.0,L) for m in 1:d), Tuple(M for m in 1:d))
        
        mesh = warp_mesh(original_mesh, reference_approximation, 
            ChanWarping(mesh_perturb,Tuple(L for m in 1:d)))

        spatial_discretization = SpatialDiscretization(mesh, 
            reference_approximation, 
            project_jacobian=!isa(reference_approximation.V,UniformScalingMap))

        solver = Solver(conservation_law,spatial_discretization, form,
            strategy, operator_algorithm, default_mass_matrix_solver(
            spatial_discretization, operator_algorithm), Serial())

        results_path = string(path, "p", p, "/")
        if !isdir(results_path)
            save_project(conservation_law,
                spatial_discretization, initial_data, form, 
                (0.0, T), results_path, overwrite=true, clear=true)
            open(string(results_path,"screen.txt"), "a") do io
                println(io, "Number of Julia threads: ", Threads.nthreads())
                println(io, "Number of BLAS threads: ", 
                    BLAS.get_num_threads(),"\n")
                println(io, "Results Path: ", "\"", results_path, "\"\n")
            end
        end

        if run
            time_steps = load_time_steps(results_path)
            if !isempty(time_steps)
                restart_step = last(time_steps)
                u0, t0 = load_solution(results_path, restart_step)
                open(string(results_path,"screen.txt"), "a") do io
                    println(io, "\nRestarting from time step ", restart_step,
                        "  t = ", t0)
                end
            else
                restart_step = 0
                u0, t0 = initialize(initial_data, spatial_discretization), 0.0
            end
            ode_problem = ODEProblem(semi_discrete_residual!, u0, (t0,T),
                solver)

            h = L/M
            dt = C_t*h/(norm(a)*p^2)

            open(string(results_path,"screen.txt"), "a") do io
                println(io, "Using ", ode_algorithm, " with  dt = ", dt)
            end
            
            reset_timer!()
            sol = solve(ode_problem, ode_algorithm, adaptive=false,
                dt=dt, save_everystep=false, callback=save_callback(
                    results_path, (t0,T), ceil(Int, T/(dt*n_s)), restart_step))

            if sol.retcode != :Success 
                open(string(results_path,"screen.txt"), "a") do io
                    println(io, "Solver failed! Retcode: ", string(sol.retcode))
                end
                continue
            end

            error_analysis = ErrorAnalysis(results_path, conservation_law, 
                spatial_discretization)
            
            error = analyze(error_analysis, last(sol.u), initial_data)

            open(string(results_path,"screen.txt"), "a") do io
                println(io, "Solver successfully finished!\n")
                println(io, @capture_out print_timer(), "\n")
                println(io,"L2 error:\n", error)
            end
            open(string(path,"screen.txt"), "a") do io
                println(io, "p = ", p, " L2 error: ", error)
            end

            if !isfile(string(path, "errors.jld2"))
                save_object(string(path, "errors.jld2"), Float64[])
            end

            errors=load_object(string(path, "errors.jld2"))
            save_object(string(path, "errors.jld2"),
                push!(errors, error[1]))
        end

        if spectral_radius

            if !isfile(string(path, "spectral_radii.jld2"))
                save_object(string(path, "spectral_radii.jld2"), Float64[])
                save_object(string(path, "max_real.jld2"), Float64[])
                save_object(string(path, "max_abs_real.jld2"), Float64[])
            end

            if eigensolver == "ArPack"

                open(string(results_path,"screen.txt"), "a") do io
                    println(io, "starting arpack eigenvalue solve")
                end

                eigenvalues, _ = eigs(LinearResidual(solver), 
                    nev=r, which=:LM, maxiter=1000)
                specr=maximum(abs.(eigenvalues))

                open(string(path,"screen.txt"), "a") do io
                    println(io, "p = ", p, ", spectral radius = ",specr)
                end

            elseif eigensolver == "KrylovKit"

                open(string(results_path,"screen.txt"), "a") do io
                    println(io, "starting krylovkit eigenvalue solve")
                end

                solver_map = LinearResidual(solver)
                f(x) = solver_map*x
                eigenvalues, _, info = eigsolve(f, size(solver_map,1), 1, :LM,
                    tol=tol, verbosity=1, ishermitian=false,
                     issymmetric=false)
                specr=maximum(abs.(eigenvalues))

                open(string(results_path,"screen.txt"), "a") do io
                    println(io, "spectral radius = ", specr)
                    println(io, @capture_out display(info))
                end

                open(string(path,"screen.txt"), "a") do io
                    println(io, "p = ", p, ", spectral radius = ",specr)
                end

            else
                solver_map = Matrix(LinearResidual(solver))

                open(string(results_path,"screen.txt"), "a") do io
                    println(io, "starting direct eigenvalue solve")
                end
                F = eigen(solver_map)
                specr=maximum(abs.(F.values))
                max_real = maximum(real.(F.values))
                max_abs_real = maximum(abs.(real.(F.values)))

                save_object(string(results_path, "spectrum.jld2"), F.values)

                open(string(path,"screen.txt"), "a") do io
                    println(io, "p = ", p, ", spectral radius = ",specr,
                        ", max real = ", max_real, 
                        ", max abs real = ", max_abs_real)
                end
            
                max_real_array=load_object(string(path, "max_real.jld2"))
                save_object(string(path, "max_real.jld2"),
                    push!(max_real_array, max_real))

                max_abs_real_array=load_object(string(path, "max_abs_real.jld2"))
                save_object(string(path, "max_abs_real.jld2"),
                    push!(max_abs_real_array, max_abs_real))
            end

            spectral_radii=load_object(string(path, "spectral_radii.jld2"))
            save_object(string(path, "spectral_radii.jld2"),
                push!(spectral_radii, specr))
        end

        if flops
            if !isfile(string(path, "flops.jld2"))
                save_object(string(path, "flops.jld2"), Float64[])
            end

            count = rhs_flops(reference_approximation, strategy=strategy,
                operator_algorithm=operator_algorithm)
            total = 2*count.muladd64 + count.add64 + count.mul64
            flop_count_array =load_object(string(path, "flops.jld2"))
            save_object(string(path, "flops.jld2"),
                push!(flop_count_array, total))

            open(string(path,"screen.txt"), "a") do io
                println(io, "p = ", p, ", flops = ", total)
            end
        end
    end
end
