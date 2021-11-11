abstract type Model{T} end

floating_base_dim(::Model) = 3
nominal_configuration(::Model) = zeros(model.nq) 
