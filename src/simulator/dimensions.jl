function num_var(model) 
    nq = model.nq 
    nc = model.nc
    nq + nc + nc + nc + nc + nc + nc 
end 

function num_data(model; nf1=1, nf2=1) 
    nq = model.nq 
    nu = model.nu 
    nw = model.nw
    2nq + nu + nw + nf1 + nf2 + 1 
end