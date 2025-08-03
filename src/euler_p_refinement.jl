struct EulerPRefinementDriver{ElemShape,ResidualForm,Strategy,Algorithm}
    p_min::Int
    p_max::Int
    p_quad::Int
    C_t::Float64
    n_s::Int
    scheme::String
    element_type::ElemShape
    form::ResidualForm
    strategy::Strategy
    ode_algorithm::Algorithm
    path::String
    M0::Int
    L::Float64
    mach_number::Float64
    base_density::Float64
    γ::Float64
    T::Float64
    mesh_perturb::Float64
    load_from_file::Bool
    overwrite::Bool
    test_type::Int
end

function EulerPRefinementDriver(
    p_min::Int,
    p_max::Int;
    p_quad = 35,
    C_t = 0.1,
    n_s = 500,
    element_type = "Tri",
    scheme = "ModalMulti",
    form = "FluxDifferencingForm",
    numerical_flux = "LaxFriedrichsNumericalFlux",
    strategy = "ReferenceOperator",
    ode_algorithm = "DP8",
    path = "../results/euler_p_test/",
    M0 = 16,
    L = 1.0,
    mach_number = 0.4,
    base_density = 0.5,
    γ = 1.4,
    T = 15.0,
    mesh_perturb = 1 / 16,
    load_from_file = true,
    overwrite = false,
    test_type = 1,
)

    path = string(path, scheme, "_", element_type, "_", form, "/", numerical_flux, "/")

    element_type = eval(Symbol(element_type))()

    form = eval(Symbol(form))(inviscid_numerical_flux = eval(Symbol(numerical_flux))())
    ode_algorithm = eval(Symbol(ode_algorithm))()
    strategy = eval(Symbol(strategy))()

    return EulerPRefinementDriver(
        p_min,
        p_max,
        p_quad,
        C_t,
        n_s,
        scheme,
        element_type,
        form,
        strategy,
        ode_algorithm,
        path,
        M0,
        L,
        mach_number,
        base_density,
        γ,
        T,
        mesh_perturb,
        load_from_file,
        overwrite,
        test_type,
    )
end

function run_driver(driver::EulerPRefinementDriver)

    (;
        p_min,
        p_max,
        p_quad,
        C_t,
        n_s,
        scheme,
        element_type,
        form,
        strategy,
        ode_algorithm,
        path,
        M0,
        L,
        mach_number,
        base_density,
        γ,
        T,
        mesh_perturb,
        load_from_file,
        overwrite,
        test_type,
    ) = driver

    if (!load_from_file || !isdir(path))
        path = new_path(path, overwrite, overwrite)
    end
    if !isdir(string(path, "/p", p_min, "/"))
        p_start = p_min
    else
        for i = p_min:p_max
            if !isdir(string(path, "p", i + 1, "/"))
                p_start = i
                break
            end
        end
    end
    open(string(path, "screen.txt"), "a") do io
        println(io, "Starting refinement from p = ", p_start)
    end

    d = dim(element_type)
    conservation_law = EulerEquations{d}(γ)

    if test_type == 1 # accuracy
        L = 2.0
        initial_data = EulerPeriodicTest(conservation_law, 0.2, L)
        a = sqrt(d)
        T = L # one period
    elseif test_type == 2 # robustness
        if d == 1
            L = 2.0
            initial_data = EulerPeriodicTest(conservation_law, 0.2, L)
            a = 1.0
            T = L
        elseif d == 2
            L = 2.0
            initial_data = KelvinHelmholtzInstability(conservation_law, base_density)
            a = 0.5
        elseif d == 3
            L = 2π
            initial_data = TaylorGreenVortex(conservation_law, mach_number)
            a = 1.0
        end
        @error "Test type should be 1 for accuracy or 2 for robustness"
    end

    alg = DefaultOperatorAlgorithm()

    for p = p_start:p_max

        if !isfile(string(path, "poly_degrees.jld2"))
            save_object(string(path, "poly_degrees.jld2"), Int64[])
        end

        poly_degrees = load_object(string(path, "poly_degrees.jld2"))
        save_object(string(path, "poly_degrees.jld2"), push!(poly_degrees, p))

        M = M0

        approx_type = eval(Symbol(scheme))(p)

        reference_approximation =
            ReferenceApproximation(approx_type, element_type, mapping_degree = p)

        original_mesh = uniform_periodic_mesh(
            reference_approximation,
            Tuple((0.0, L) for m = 1:d),
            Tuple(M for m = 1:d),
        )

        if d == 1
            spatial_discretization =
                SpatialDiscretization(original_mesh, reference_approximation)
        else
            mesh = warp_mesh(
                original_mesh,
                reference_approximation,
                ChanWarping(mesh_perturb, Tuple(L for m = 1:d)),
            )
            spatial_discretization =
                SpatialDiscretization(mesh, reference_approximation, ChanWilcoxMetrics())
        end

        solver = Solver(
            conservation_law,
            spatial_discretization,
            form,
            strategy,
            alg,
            default_mass_matrix_solver(spatial_discretization, alg),
            Threaded(),
        )

        results_path = string(path, "p", p, "/")
        if !isdir(results_path)
            save_project(
                conservation_law,
                spatial_discretization,
                initial_data,
                form,
                (0.0, T),
                results_path,
                overwrite = true,
                clear = true,
            )
            open(string(results_path, "screen.txt"), "a") do io
                println(io, "Number of Julia threads: ", Threads.nthreads())
                println(io, "Number of BLAS threads: ", BLAS.get_num_threads(), "\n")
                println(io, "Results Path: ", "\"", results_path, "\"\n")
            end
        end

        time_steps = load_time_steps(results_path)
        if !isempty(time_steps)
            restart_step = last(time_steps)
            u0, t0 = load_solution(results_path, restart_step)
            open(string(results_path, "screen.txt"), "a") do io
                println(io, "\nRestarting from time step ", restart_step, " t = ", t0)
            end
        else
            restart_step = 0
            u0, t0 = initialize(initial_data, spatial_discretization), 0.0
        end

        ode_problem = ODEProblem(semi_discrete_residual!, u0, (t0, T), solver)
        h = L / M
        dt = C_t * h / (a * p^2)

        open(string(results_path, "screen.txt"), "a") do io
            println(io, "Using ", ode_algorithm, " with  dt = ", dt, ", T = ", T)
        end

        reset_timer!()
        integrator = init(
            ode_problem,
            ode_algorithm,
            adaptive = false,
            dt = dt,
            save_everystep = false,
            callback = save_callback(
                results_path,
                (t0, T),
                ceil(Int, T / (dt * n_s)),
                restart_step,
            ),
        )
        try
            solve!(integrator)
        catch ex
            open(string(results_path, "screen.txt"), "a") do io
                println(io, "Solver failed! ", ex)
            end
        end

        open(string(path, "screen.txt"), "a") do io
            println(io, "p = ", p)
        end

        if test_type == 1
            error_analysis = ErrorAnalysis(
                results_path,
                conservation_law,
                spatial_discretization,
                DefaultQuadrature(p_quad),
            )

            error = analyze(error_analysis, last(integrator.sol.u), initial_data)

            open(string(results_path, "screen.txt"), "a") do io
                println(io, @capture_out print_timer(), "\n")
                println(io, "L2 error:\n", error)
            end
            open(string(path, "screen.txt"), "a") do io
                println(io, ", L2 error: ", error)
            end

            if !isfile(string(path, "errors.jld2"))
                save_object(string(path, "errors.jld2"), Vector{Float64}[])
            end
            errors = load_object(string(path, "errors.jld2"))
            save_object(string(path, "errors.jld2"), push!(errors, error))
        end

        open(string(path, "screen.txt"), "a") do io
            println(io, "end time: ", integrator.t)
        end
        if !isfile(string(path, "end_times.jld2"))
            save_object(string(path, "end_times.jld2"), Float64[])
        end

        end_times = load_object(string(path, "end_times.jld2"))
        save_object(string(path, "end_times.jld2"), push!(end_times, integrator.t))

        time_steps = load_time_steps(results_path)
        conservation = analyze(
            PrimaryConservationAnalysis(
                results_path,
                conservation_law,
                spatial_discretization,
            ),
            time_steps,
        )
        entropy = analyze(
            EntropyConservationAnalysis(
                results_path,
                conservation_law,
                spatial_discretization,
            ),
            time_steps,
        )
        save_object(string(results_path, "conservation.jld2"), conservation.E)
        save_object(string(results_path, "conservation_residual.jld2"), conservation.dEdt)
        save_object(string(results_path, "entropy.jld2"), entropy.E)
        save_object(string(results_path, "entropy_residual.jld2"), entropy.dEdt)
    end
end
