include("seirs_det.jl")
include("seirs_stoch.jl")
include("seirs_ibm.jl")

function run_all(p, do_plot = false)
	if !do_plot
		println("\nIndividual Based model:")
		@btime soln_ibm = seirs_ibm($p)
		nothing
	else
		#	println("Deterministic model:")
		#	@time soln_det = seirs_det(p)

		#	println("\nStochastic model:")
		#	@time soln_stoch = seirs_stoch(p)

		println("\nIndividual Based model:")
		@time soln_ibm = seirs_ibm(p)
		
		#	p1 = plot_det(soln_det)
		#	p2 = plot_stoch(soln_stoch)
		p3 = plot_ibm(soln_ibm)

		#	plot(p1, p2, p3, layout=(3, 1))
		plot(p3)
	end
end #fn
	
# parameters
# shared between det, stoch, and IBM models
# note: this is a Named Tuple, which is immutable, but very convenient for typing
N = 10_000
timestep = 1.0
track_days = 14.0
mat_type = UInt8
params = (
	# general parameters
	I0 = inv(100.0), # initial proportion of infected individuals
	G0 = inv(100.0), # initial proportion of infected groups
	T = 150.0,   # time to run simulation
	c = N,    # population size
	α = inv(5.0),    # incubation rate (1/period)
	β = 2.0/7.0,    # transmission coefficient
	γ = inv(7.0),    # recovery rate
	σ = 0.0,    # immunity waning rate
	μ = 0.0,    # mortality rate
	ν = 0.0,    # disease induced mortality rate
	κ = inv(365.0),  # movement rate

	# stochastic parameters
	method = :tau,  # :gillespie or :tau (stoch model)
	num_itns = 2,  # iteration
	num_recs = 301, # records per iteration
	num_groups = 1, # metapopulations
	contact_rate = 20.0, # contact rate (IBM)
	dt = timestep
)

run_all(params)
