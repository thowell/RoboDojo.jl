function num_var(model) 
    nq = model.nq 
    nc = model.nc
    nq + nc + nc + nc + nc + nc + nc 
end 

function num_data(model; nf=1) 
    nq = model.nq 
    nu = model.nu 
    nw = model.nw
    2nq + nu + nw + nf + 1 
end