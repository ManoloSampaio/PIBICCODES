function x=simulation(f,T,x_initial,u,theta)
x = zeros(length(x_initial),length(T));
for i=1:length(T)-1
    t_step = T(i):0.01:T(i+1);
    [t_step,x_aux] = ode45(@(t,x_aux)f(t_step,x_aux,u(:,i),theta),t_step,x(:,i));
    x(:,i+1)=x_aux(end,:);
end

end