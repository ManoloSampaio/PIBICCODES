close all;clear all;clc;
y0= [3 0 0 0];
T =0:0.01:40;
[time,yode] = ode45(@(time,yode)invertedpendulumode(yode),T,y0);
[A1,B1,A2,B2,C,D]=linearpend(yode);
invertedpendulum = ss(A1,B1,C,D);
T =0:0.01:40;
U = 0*ones(1,length(T));
UE = 0*ones(1,length(T));
y = lsim(invertedpendulum,U-UE,T,[3 0  0 0]-[pi 0 0 0]);
%% Inverted pendulum ODE.
figure()
plot(time,yode(:,1),'b')
xlabel('Time (s)','fontsize',15,'fontweight','bold')
ylabel('Pendulum angle (rad)','fontsize',15,'fontweight','bold')
%title("Pendulum angular position nonlinear model")
%print('NonLinearAngularPosition','-depsc')
hold on
plot(T,y(:,1)+pi,'r','LineStyle','--')
legend('Nonlinear Model','Linear Model')
figure()

%print('LinearCartPosition','-depsc')
plot(time,yode(:,3),'b')
xlabel('Time (s)','fontsize',15,'fontweight','bold')
ylabel('Cart position (m)','fontsize',15,'fontweight','bold')
%title("Cart postion nonlinear model")
hold on
plot(T, y(:,2),'r','LineStyle','--')
legend('Nonlinear Model','Linear Model')
%title("Cart position using the linear model approximation")

%print('NonLinearCartPosition','-depsc')
%print('LinearAngularPosition','-depsc')   