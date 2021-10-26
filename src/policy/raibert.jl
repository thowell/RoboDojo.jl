"""
    Raibert heuristic policy for Hopper
"""

mutable struct Raibert{T} <: Policy{T}
	kr_c::T
	kr_p::T
	kr_v_stance::T
	kr_v_flight::T
	kθ_c::T
	kθ_p::T
	kθ_v::T
	u::Vector{T}
	contact::Bool
	v0::T
	Tstance::T
	Tflight::T
	q0::Vector{T}
	q1::Vector{T}
    h::T
end

function raibert_policy(model::Hopper;
		v0=0.5, 
        Tstance=0.13, 
        Tflight=0.62, 
        h=0.1,
		kr_c=8.0e1,
		kr_p=-1.0e3,
		kr_v_stance=-1.0e-2,
		kr_v_flight=-1.0e1,
		kθ_c=0.0,
		kθ_p=-3.0e1,
		kθ_v=-1.0e1)

	u  = zeros(model.nu)
	q0 = zeros(model.nq)
	q1 = zeros(model.nq)
	contact = false

	Raibert(kr_c, kr_p, kr_v_stance, kr_v_flight, kθ_c, kθ_p, kθ_v,
		u, contact, v0, Tstance, Tflight, q0, q1, h)
end

function policy(p::Raibert, traj::Trajectory, t)
	# Initialization
	h = p.h
	p.q0 .= traj.q[t]
	p.q1 .= traj.q[t+1]

	# Detect contact
	if any(traj.γ[max(1,t-1)] .> 1.5e-2)
		p.contact = true
	else
		p.contact = false
	end

	# Velocities
	qv = (p.q1 - p.q0) / h
	θv = qv[3]
	rv = qv[4]

	# References
	θ1 = p.q1[3]
	r1 = p.q1[4]
	rref = 0.5
	# td = touchdown
	θtd  = 0.5 * asin(p.v0 * p.Tstance / (2.0 * rref))

	# Gains
	kr_c = p.kr_c
	kr_p = p.kr_p
	kr_v_stance = p.kr_v_stance
	kr_v_flight = p.kr_v_flight

	kθ_c = p.kθ_c
	kθ_p = p.kθ_p
	kθ_v = p.kθ_v

	if p.contact
		# regulate around θtd using same gain kθ_p but taking into account
		# that stance time is much lower than flight time
		p.u[1] = kθ_c + kθ_p*(θ1 + θtd) * p.Tflight / p.Tstance
		p.u[2] = kr_c + kr_p*(r1 - rref) + kr_v_stance * rv
	else
		# regulate around θtd
		p.u[1] = kθ_p * (θ1 - θtd)  + kθ_v * θv
		p.u[2] = kr_p * (r1 - rref) + kr_v_flight * rv
	end

	return p.u .* h
end