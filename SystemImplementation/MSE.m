vector = [1:0.01:5];
error = ones(length(vector),1);
for i=1:length(vector)    
    y0= [vector(i) 0 0 0];
    T =0:0.01:40;
    [time,yode] = ode45(@(time,yode)invertedpendulumode(yode),T,y0);
    [A1,B1,A2,B2,C,D]=linearpend(yode);
    invertedpendulum = ss(A1,B1,C,D);
    T =0:0.01:40;
    U = 0*ones(1,length(T));
    UE = 0*ones(1,length(T));
    y = lsim(invertedpendulum,U-UE,T,y0-[pi 0 0 0]);
    error(i) = immse(yode(:,1),y(:,1)+pi);
end
plot(vector,error,'k')
xlabel('The initial angle (rad)','fontsize',15,'fontweight','bold')
ylabel('MSE','fontsize',15,'fontweight','bold')