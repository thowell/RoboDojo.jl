abstract type Model{T} end

floating_base_dim(::Model) = 3
nominal_configuration(model::Model) = zeros(model.nq) 
