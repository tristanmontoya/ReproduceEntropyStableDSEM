struct AdvectionDriver{d, ApproximationType <: AbstractApproximationType,
    ElemShape <: AbstractElemShape, ResidualForm <: AbstractResidualForm,
    Strategy <: AbstractStrategy, Algorithm <: OrdinaryDiffEqAlgorithm}
    p::Int
    l::Int
    C_t::Float64
    n_s::Int
    scheme::ApproximationType
    element_type::ElemShape
    form::ResidualForm
    strategy::Strategy
    ode_algorithm::Algorithm
    path::String
    M0::Int
    λ::Float64
    L::Float64
    a::NTuple{d,Float64}
    T::Float64
    mesh_perturb::Float64
    n_grids::Int
    load_from_file::Bool
    overwrite::Bool
end

function AdvectionDriver(p; l=4, C_t=1.0e-4, n_s=50,
    element_type="Tri", scheme="ModalMulti",
    mapping_form="SkewSymmetricMapping", 
    strategy="ReferenceOperator", ode_algorithm="CarpenterKennedy2N54", 
    path="../results/20230421/", M0 = 2,
    λ=1.0, L=1.0, a = nothing, T = 1.0, mesh_perturb = 1/16,
    n_grids = 6, load_from_file=true, overwrite=false)

    if Int(round(λ)) == 0 
        path = string(path, scheme, "_", element_type, "_",
            mapping_form, "_p", p, "/central/")
    elseif Int(round(λ)) == 1 
        path = string(path, scheme, "_", element_type,
            "_", mapping_form, "_p", p, "/upwind/")
    else 
        path = string(path, scheme, "_", element_type, "_",
            mapping_form, "_p", p, "/lambda", replace(string(lambda), 
            "." => "_"), "/")
    end

    element_type = eval(Symbol(element_type))()
    scheme = eval(Symbol(scheme))(p)
    form = StandardForm(mapping_form=eval(Symbol(mapping_form))(),
        inviscid_numerical_flux=LaxFriedrichsNumericalFlux(λ))
    ode_algorithm = eval(Symbol(ode_algorithm))()
    strategy = eval(Symbol(strategy))()

    if isnothing(a)
        a = Tuple(1.0 for m in 1:dim(element_type))
    end

    return AdvectionDriver(p,l,C_t,n_s,scheme,element_type,
        form,strategy,ode_algorithm,path,M0,λ,L,
        a, T, mesh_perturb, n_grids, load_from_file, overwrite)
end

function run_driver(driver::AdvectionDriver{d}) where {d}

    @unpack p,l,C_t,n_s,scheme,element_type,form,strategy,ode_algorithm,path,M0,λ,L,a,T,mesh_perturb, n_grids, load_from_file, overwrite = driver

    if (!load_from_file || !isdir(path)) path = new_path(
        path,overwrite,overwrite) end
    if !isdir(string(path, "grid_1/")) n_start = 1 else
        for i in 1:n_grids
            if !isdir(string(path, "grid_", i + 1, "/"))
                n_start = i
                break
            end
        end
    end
    open(string(path,"screen.txt"), "a") do io
        println(io, "Starting refinement from grid level ", n_start)
    end

    initial_data = InitialDataSine(1.0,Tuple(2*π/L for m in 1:d))
    conservation_law = LinearAdvectionEquation(a)
    eoc = -1.0

    alg = DefaultOperatorAlgorithm()

    for n in n_start:n_grids

        M = M0*2^(n-1)

        reference_approximation =ReferenceApproximation(
            scheme, element_type, mapping_degree=l)

        original_mesh = uniform_periodic_mesh(reference_approximation, 
            Tuple((0.0,L) for m in 1:d), Tuple(M for m in 1:d))
        
        mesh = warp_mesh(original_mesh, reference_approximation, 
            ChanWarping(mesh_perturb,Tuple(L for m in 1:d)))

        spatial_discretization = SpatialDiscretization(mesh, 
            reference_approximation, 
            project_jacobian=!isa(reference_approximation.V,UniformScalingMap))

        solver = Solver(conservation_law,spatial_discretization, form,
            strategy, alg, default_mass_matrix_solver(spatial_discretization,
            alg), Threaded())
        
        results_path = string(path, "grid_", n, "/")
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
        ode_problem = ODEProblem(semi_discrete_residual!, u0, (t0,T), solver)

        h = L/M
        dt = C_t*h/norm(a)

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

        open(string(results_path,"screen.txt"), "a") do io
            println(io, "Solver successfully finished!\n")
            println(io, @capture_out print_timer(), "\n")
            println(io,"L2 error:\n", analyze(error_analysis,
                last(sol.u), initial_data))
        end

        if n > 1
            refinement_results = analyze(RefinementAnalysis(initial_data, path,
                "./", "advection test"), n, max_derivs=false)
            open(string(path,"screen.txt"), "a") do io
                println(io, tabulate_analysis(refinement_results, e=1,
                    print_latex=false))
            end
            eoc = refinement_results.eoc[end,1]
        end
    end
    return eoc
end