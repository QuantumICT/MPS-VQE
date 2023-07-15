using Distributed, SlurmClusterManager
addprocs(SlurmManager())

cprefix = "../H2/"*ARGS[1]

@everywhere push!(LOAD_PATH, "./QuantumSpin/src")

@everywhere using QuantumSpins
@everywhere using NLopt
@everywhere using Parsers
@everywhere using SharedArrays
@everywhere using DistributedArrays

@everywhere u3(theta, phi, lambda) = [cos(theta/2) -exp(im * lambda)*sin(theta/2); exp(im * phi) * sin(theta/2) exp(im*(phi+lambda)) * cos(theta/2)]

@everywhere function parse_gate(s::AbstractString)
	item = split(s, [' ', '(', ')', ',', ';'])
	item = [item_ for item_ in item if !isempty(item_)]
	gate = item[1]
	# println("gate is $gate")
	if gate == "x"
		pos = parse(Int, item[2][findfirst('[', item[2])+1:findfirst(']', item[2])-1])
		return [XGate(pos)]
	elseif gate == "cx"
		pos1 = parse(Int, item[2][findfirst('[', item[2])+1:findfirst(']', item[2])-1])
		pos2 = parse(Int, item[3][findfirst('[', item[3])+1:findfirst(']', item[3])-1])
		return [CNOTGate(pos1, pos2)]
	elseif gate == "cy"
		pos1 = parse(Int, item[2][findfirst('[', item[2])+1:findfirst(']', item[2])-1])
		pos2 = parse(Int, item[3][findfirst('[', item[3])+1:findfirst(']', item[3])-1])
		return [RzGate(pos2, -pi/2.0), CNOTGate(pos1, pos2), RzGate(pos2, pi/2.0)]
	elseif gate == "cz"
		pos1 = parse(Int, item[2][findfirst('[', item[2])+1:findfirst(']', item[2])-1])
		pos2 = parse(Int, item[3][findfirst('[', item[3])+1:findfirst(']', item[3])-1])
		return [CZGate(pos1, pos2)]
	elseif gate == "h"
		pos = parse(Int, item[2][findfirst('[', item[2])+1:findfirst(']', item[2])-1])
		return [HGate(pos)]
	elseif gate[1:2] == "rz"
		theta = parse(Float64, item[2])
		pos = parse(Int, item[3][findfirst('[', item[3])+1:findfirst(']', item[3])-1])
		if (length(item) > 3)
			amp = parse(Int, item[4][findfirst('[', item[4])+1:findfirst(']', item[4])-1]) + 1
			return [AmpRzGate(pos, theta, amp)]
		end
		return [RzGate(pos, theta)]
	elseif gate[1:2] == "u3"
		theta = parse(Float64, item[2])
		phi = parse(Float64, item[3])
		lbd = parse(Float64, item[4])
		pos = parse(Int, item[5][findfirst('[', item[5])+1:findfirst(']', item[5])-1])
		return [Gate(pos, u3(theta, phi, lbd))]
	elseif gate == "measure"
		return nothing
	else
		error("unknown operation $gate")
	end
end


@everywhere function p_range(n::Int, n_procs::Int, pid::Int)::Tuple{Int, Int}
	aprocs = n_procs - n % n_procs + 1
	q = n ÷ n_procs 
	if (pid < aprocs)
		pstart = (pid - 1) * q + 1
		pend = pstart + q - 1
	else
		pstart = (aprocs-1) * q + (pid - aprocs)*(q+1) + 1
		pend = pstart + q
	end
	return (pstart, pend)
end

@everywhere function parse_circuit(file_name::AbstractString)
    # println("read circuit from path $file_name")
    lines = readlines(file_name)[5:end]
    circuit = QCircuit()
    for line in lines
    	gate = parse_gate(line)
    	if !isnothing(gate)
		for gate_ in gate
    			push!(circuit, gate_)
		end
    	end
    end
    return circuit
end
@everywhere function read_circuit_local(cprefix::AbstractString, prange::Tuple{Int,Int}) :: Vector{QuantumCircuit}
	C = Vector{QuantumCircuit}()
	for i = prange[1]:prange[2]
		circuit=shift(parse_circuit(cprefix*"/circuit_qasm_$(i).txt"), 1)
		push!(C, circuit)
	end
	return C
end

function read_hamiltonian(file_name::AbstractString)
	hf = open(file_name,"r")
	lines = readlines(hf)
	h = Float64[0.0 for line in lines]
	i = 0
	for line in lines
		i = i + 1
		h[i] = Parsers.parse(Float64,line)
	end
	close(hf)
	return h
end

@everywhere function get_expec(n_qubits::Int, DC :: DArray{QuantumCircuit, 1, Vector{QuantumCircuit}}, r::Tuple{Int,Int}, amplitudes :: Vector{Float64}) :: Vector{Float64}
		len = r[2] - r[1] + 1
		expecs = fill(0.0, len)
		new_circuit = QCircuit()
		for j in 1:len
			#lcircuit = localpart(DC)[j]
			new_circuit = copy(localpart(DC)[j])
			for i = 1:length(new_circuit)
				g = getindex(new_circuit, i)
				if (typeof(g)==AmpRzGate{Float64})
					setindex!(new_circuit, RzGate(g.positions[1], g.parameter * amplitudes[g.amp]), i)	
				end	
			end


			# new_circuit = fuse_gates(new_circuit)
			# println("number of gates after gate fusion $(length(circuit)).")
			state = statevector_mps(n_qubits)
			mpstrunc = MPSTruncation(ϵ=1.0e-5, D=50)
			apply!(new_circuit, state, trunc=mpstrunc)
			measure_qubit = n_qubits
			outcome, prob_at_measure_qubit = measure!(state, measure_qubit, keep=true)
			# println("Probility of $(measure_qubit) qubit -> |$(outcome)> is $(prob_at_measure_qubit)")
			if outcome == 0
				expecs[j] = 2 * prob_at_measure_qubit - 1
			elseif outcome == 1
				expecs[j] = 2 * (1 - prob_at_measure_qubit) - 1
			else
				error("Outcome should be 0 or 1!")
			end
			new_circuit = QCircuit()
			GC.gc()
		end
		return expecs
end


# you may need the funciton expectation to compute expectation values
function main()

	n_qubits = 5
	n_amplitudes = 2
	n_circuits = 15
	
	n_procs = nprocs()
	println("working on",n_procs,"processings")		
	flush(stdout)	

	h = read_hamiltonian(cprefix*"/hamiltonian_coefficients.csv")


	DCircuit = DArray( [@spawnat p read_circuit_local(cprefix, p_range(n_circuits - 1, n_procs, p)) for p in 1:n_procs])


	function get_total_energy(amplitudes::Vector{Float64}) :: Float64
		expectation_vals = [1.0 for i in 1:n_circuits]
		futures = Array{Future}(undef, n_procs)
		for i = 1:n_procs
			range = p_range(n_circuits - 1, n_procs, i)
			if range[2] >= range[1]
				futures[i] = remotecall(get_expec, i, n_qubits, DCircuit, range, amplitudes)
			end
		end
		for i = 1:n_procs
			range = p_range(n_circuits - 1, n_procs, i)	
			if range[2] >= range[1]
				expectation_vals[range[1]+1:range[2]+1] = fetch(futures[i])
			end
		end
		energy=sum(h.*expectation_vals)
		println("-----end get_total_energy------- ",energy)
		flush(stdout)
		return energy	
	end
	
	
	for i in 1:2
		println("test ",i)
		@time println("energy at 0:",get_total_energy(Float64[0.0 for i in 1:n_amplitudes]))
		println("")
	end


	opt = Opt(:LN_COBYLA, n_amplitudes)
	opt.min_objective = (amplitudes::Vector{Float64}, grad::Vector{Float64}) -> @time get_total_energy(amplitudes)
	opt.maxeval = 42
	t1=time()	
	(minf, minx, ret) = optimize(opt, Float64[0.0 for i in 1:n_amplitudes])
	t2=time()
	println("time: ",t2-t1," seconds, got $minf at $minx (returned $ret)")

end

main()








