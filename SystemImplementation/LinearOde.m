function dx= LinearOde(time,x,u,theta)
    A = theta{1};
    B = theta{2};
    x_ss = theta{3};
    u_ss = theta{4};
    dx = A*(x-x_ss)+B*(u-u_ss);
    
end