using Statistics
using StatsBase
using Distributions
using Plots

function overshoot(X, dX)
	for i in eachindex(X)
		X[i] + dX[i] < 0 && return true
	end #for
	false
end #fn


function check_pop_error(X)
	for x in X
		x < 0 && return true
	end #for
	false
end #fn


function seirs_stoch(params)
	T = params.T
	num_recs = params.num_recs
	num_itns = params.num_itns
	num_groups = params.num_groups
	method = params.method

	# main parameters
	c = round(Int, params.c)
	α = params.α
	β = params.β
	γ = params.γ
	σ = params.σ
	μ = params.μ
	ν = params.ν
	κ = num_groups > 1 ? params.κ : 0.0

	# initial conditions
	X0 = [c, 0, 0, 0]
	E0 = ceil(Int, 0.01c)
	I0 = [c - E0, E0, 0, 0]
	G0 = params.G0
	initial_groups = 1:ceil(Int, G0 * num_groups) # number of groups seeded

	# event table
	# not sure if using collect() makes any difference, but it turns V' from an
	# LinearAlgebra.Adjoint into a regular array
	V = [
		#S  E  I  R
		+1  0  0  0 #  1 birth S
		-1  0  0  0 #  2 death S
		 0 -1  0  0 #  3 death E
		 0  0 -1  0 #  4 death I
		 0  0  0 -1 #  5 death R
		-1 +1  0  0 #  6 infection
		 0 -1 +1  0 #  7 incubation
		 0  0 -1 +1 #  8 recovery
		+1  0  0 -1 #  9 loss of immunity
		-1  0  0  0 # 10 movement S
		 0 -1  0  0 # 11 movement E
		 0  0 -1  0 # 12 movement I
		 0  0  0 -1 # 13 movement R
	]' |> collect

	num_vars, num_events = size(V)

	# data recording
	rec_width = T / (num_recs - 1)
	records = zeros(num_vars, num_recs, num_itns)

	# this is mostly pre-allocating arrays
	X         = zeros(Int, num_vars, num_groups)
	dX        = zeros(Int, size(X))
	rates     = zeros(num_events, num_groups)
	cum_rates = zeros(length(rates)) # this must be 1d
	update    = zeros(Bool, num_groups) # only calculates rates that have changed
	K         = zeros(Int, size(rates))
	moved     = zeros(Int, num_vars)
	move_events = 10:num_events


	for itn in 1:num_itns
		# set initial X and seed infections
		# (this could be randomised, so keep inside the for loop)
		X .= X0
		X[:, initial_groups] .= I0

		# at start of iteration, using Gillespie
		# all groups need to be updated
		update .= true

		t = 0.0
		dt = min(rec_width, 0.1)

		t_next_rec = 0.0
		rec = 0

		while t ≤ T
			while t ≥ t_next_rec && rec < num_recs
				records[:, rec += 1, itn] = mean(X, dims=2)
				t_next_rec += rec_width
			end #while

			for group in 1:num_groups
				update[group] || continue

				S, E, I, R = @view X[:, group]
				N = S + E + I + R

				rates[1, group] = μ * c # birth S
				rates[2, group] = μ * S # death S
				rates[3, group] = μ * E # death E
				rates[4, group] = μ * I # death I
				rates[5, group] = μ * R # death R
				rates[6, group] = β * S * I / N # infection
				rates[7, group] = α * E # incubation
				rates[8, group] = γ * I # recovery
				rates[9, group] = σ * R # loss of immunity
				rates[10, group] = κ * S # movement of S
				rates[11, group] = κ * E # movement of E
				rates[12, group] = κ * I # movement of I
				rates[13, group] = κ * R # movement of R
			end #for

			if method == :tau
				# this is my super unsafe "take longer steps, fix if things break" method
				# but it works
				dt = min(1.2dt, rec_width)

				dX .= 0
				moved .= 0
				while true
					K .= rand.(Poisson.(rates * dt))
					dX .= V * K

					# this checks if dX would send any X<0
					# if so, halve dt and try again
					if overshoot(X, dX)
						dt /= 2
						continue
					end #if

					break
				end #while


				X .+= dX

				# to handle movement, collect all individuals removed via movement
				# then randomly reassign them to other groups
				# this is safe because it's only adding
				moved .= vec(sum(@view(K[move_events, :]), dims=2))
				for (v, n) in enumerate(moved)
					X[v, :] .+= rand(Multinomial(n, num_groups))
				end #for
			elseif method == :gillespie
				cum_rates .= cumsum(@view rates[:])
				Rtot = sum(rates)

				if Rtot < 0
					@show X rates Rtot
					error("Rtot < 0 at t = $t")
				end #if

				dt = rand(Exponential(1 / Rtot))
				# this chooses event from the rates array
				event_idx = sample(1:length(rates), ProbabilityWeights(@view rates[:]))
				# but since we can only sample from a 1-dim array, we need to coerce
				# rates to be 1-dim then re-extract event and group from the idx.
				event = mod1(event_idx, num_events)
				group = ceil.(Int, event_idx / num_events)

				X[:, group] .+= @view(V[:, event])

				update[group] = true

				# handle movement
				if event ∈ move_events
					to_group = rand(1:num_groups)
					X[:, to_group] .-= @view(V[:, event])
					update[to_group] = true
				end #if

			else
				error("Method: choose either :tau or :gillespie")
			end #if

			# check population for errors
			if check_pop_error(X)
				@show X rates
				error("negative population at t = $t")
			end #if

			t += dt
		end #while

		# finish recording in case a leap is > rec_width
		# alternatively could cap dt at rec_width, I just prefer it this way
		while rec < num_recs
			records[:, rec += 1, itn] = mean(X, dims=2)
		end
	end #for

	(params = params, soln = records)
end #fn


function plot_stoch(results)
	# extract parameters to get time
	T = results.params.T
	num_recs = results.params.num_recs
	# t = (0:num_recs - 1) * T / num_recs
	t = range(0, T, length=num_recs)

	# simple mean of variables
	soln = results.soln
	X = dropdims(mean(soln, dims=3), dims=3)'

	plot(t, X, lw=3,
	xlim=(0, T),
	xlab="Time",
	ylab="Population",
	title="Stochastic SEIRS model",
	label=["S" "E" "I" "R"])
end #fn
