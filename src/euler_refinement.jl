struct EulerDriver{ApproximationType,ElemShape,ResidualForm,Strategy,Algorithm}
    l::Int
    p_quad::Int
    C_t::Float64
    scheme::ApproximationType
    element_type::ElemShape
    form::ResidualForm
    strategy::Strategy
    ode_algorithm::Algorithm
    path::String
    M0::Int
    L::Float64
    mach_number::Float64
    γ::Float64
    n_periods::Int
    mesh_perturb::Float64
    n_grids::Int
    load_from_file::Bool
    overwrite::Bool
    test_type::Int
end

function EulerDriver(
    p;
    l = p,
    p_quad = 35,
    C_t = 0.01,
    element_type = "Tri",
    scheme = "ModalMulti",
    form = "FluxDifferencingForm",
    numerical_flux = "LaxFriedrichsNumericalFlux",
    strategy = "ReferenceOperator",
    ode_algorithm = "DP8",
    path = "../results/euler_20230918/",
    M0 = 2,
    L = 1.0,
    mach_number = 0.4,
    γ = 1.4,
    n_periods = 1,
    mesh_perturb = 1 / 16,
    n_grids = 6,
    load_from_file = true,
    overwrite = false,
    test_type = 1,
)

    path = string(
        path,
        scheme,
        "_",
        element_type,
        "_",
        form,
        "_p",
        p,
        "/",
        numerical_flux,
        "/",
    )

    element_type = eval(Symbol(element_type))()
    scheme = eval(Symbol(scheme))(p)

    form = eval(Symbol(form))(inviscid_numerical_flux = eval(Symbol(numerical_flux))())
    ode_algorithm = eval(Symbol(ode_algorithm))()
    strategy = eval(Symbol(strategy))()

    return EulerDriver(
        l,
        p_quad,
        C_t,
        scheme,
        element_type,
        form,
        strategy,
        ode_algorithm,
        path,
        M0,
        L,
        mach_number,
        γ,
        n_periods,
        mesh_perturb,
        n_grids,
        load_from_file,
        overwrite,
        test_type,
    )
end

function run_driver(driver::EulerDriver)

    (;
        l,
        p_quad,
        C_t,
        scheme,
        element_type,
        form,
        strategy,
        ode_algorithm,
        path,
        M0,
        L,
        mach_number,
        γ,
        n_periods,
        mesh_perturb,
        n_grids,
        load_from_file,
        overwrite,
        test_type,
    ) = driver

    if (!load_from_file || !isdir(path))
        path = new_path(path, overwrite, overwrite)
    end
    if !isdir(string(path, "grid_1/"))
        n_start = 1
    else
        for i = 1:n_grids
            if !isdir(string(path, "grid_", i + 1, "/"))
                n_start = i
                break
            end
        end
    end

    open(string(path, "screen.txt"), "a") do io
        println(io, "Starting refinement from grid level ", n_start)
    end

    d = dim(element_type)
    conservation_law = EulerEquations{d}(γ)

    if (d == 2) && (test_type == 1)
        initial_data = IsentropicVortex(
            conservation_law,
            θ = 0.0,
            Ma = mach_number,
            β = sqrt(2 / (γ - 1) * (1 - 0.75^(γ - 1))),
            R = 0.1,
            x_0 = (L / 2, L / 2),
        )
        T = n_periods * L / mach_number
        a = mach_number
    else
        L = 2.0
        initial_data = EulerPeriodicTest(conservation_law, 0.2, L)
        T = n_periods * L
        a = sqrt(d)
    end

    eoc = -1.0

    alg = DefaultOperatorAlgorithm()

    if scheme isa ModalTensor
        reference_approximation = ReferenceApproximation(
            scheme,
            element_type,
            mapping_degree = l,
            sum_factorize_vandermonde = false,
        )
    else
        reference_approximation =
            ReferenceApproximation(scheme, element_type, mapping_degree = l)
    end

    for n = n_start:n_grids

        M = M0 * 2^(n - 1)

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

        results_path = string(path, "grid_", n, "/")
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
                println(io, "\nRestarting from time step ", restart_step, "  t = ", t0)
            end
        else
            restart_step = 0
            u0, t0 = initialize(initial_data, spatial_discretization), 0.0
        end
        ode_problem = ODEProblem(semi_discrete_residual!, u0, (t0, T), solver)

        h = L / M
        dt = C_t * h / (a * scheme.p^2)
        open(string(results_path, "screen.txt"), "a") do io
            println(io, "Using ", ode_algorithm, " with  dt = ", dt)
        end

        reset_timer!()
        sol = solve(
            ode_problem,
            ode_algorithm,
            adaptive = false,
            dt = dt,
            save_everystep = false,
            callback = save_callback(
                results_path,
                (t0, T),
                ceil(Int, T / (dt * 100)),
                restart_step,
            ),
        )

        if sol.retcode != :Success
            open(string(results_path, "screen.txt"), "a") do io
                println(io, "Solver failed! Retcode: ", string(sol.retcode))
            end
            continue
        end

        error_analysis = ErrorAnalysis(
            results_path,
            conservation_law,
            spatial_discretization,
            DefaultQuadrature(p_quad),
        )

        open(string(results_path, "screen.txt"), "a") do io
            println(io, "Solver successfully finished!\n")
            println(io, @capture_out print_timer(), "\n")
            println(io, "L2 error:\n", analyze(error_analysis, last(sol.u), initial_data))
        end

        if n > 1
            refinement_results = analyze(
                RefinementAnalysis(initial_data, path, "./", "euler_test"),
                n,
                max_derivs = false,
            )
            open(string(path, "screen.txt"), "a") do io
                println(
                    io,
                    tabulate_analysis(refinement_results, e = 1, print_latex = false),
                )
            end
            eoc = refinement_results.eoc[end, 1]
        end
    end
    return eoc
end
