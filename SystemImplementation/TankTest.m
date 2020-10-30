%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%> Inteligent Control System 2020.1
%> Homework 1. 
%> Student: Emmanphil Sampaio
%> System choosed: Three Tank
%> This code contains:
%>      1. Linearization in the easy mode.
%>      2. Checking the transient response information.
%>      3. Compare the response of the systems models (linear and non linear).
%>      4. Test the system with other responses.
%>      4. Plotting the poles of the system.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setting the inputs and states.
    clc;
    clear;
    syms phi1 phi2 h1 h2 h3
    % Setting the parameters
    % mi2,mi2: Adjust the flow rates of the pipes.
    % a1,a2: Cross sectional areas: of the tanks, of the pipes.
    % g: Gravity
    % k: Constant for approximate the signal function.

    mi2 = 0.825;
    mi1 = 0.65;
    a1 = pi*((10/2)*10^-3)^2;
    a2 = pi*((7/2)*10^-3)^2;
    %a1 = 0.154;
    %a2 = 10^-4;
    g = 9.81;
    k = 91.01;

    % Setting the signal function approximation 
    % sign(x) = 2/(1+e^(-2(k)(h1-h2))
    signapprox1 = -1 + 2*((1+exp(-2*(k)*(h1-h2))).^-1);
    signapprox2 = -1 + 2*((1+exp(-2*(k)*(h2-h3))).^-1);

% Setting an approximation for the absolute function.
abs1 = ((h1-h2)^2+0.00001)^0.5;
abs2 = ((h2-h3)^2+0.00001)^0.5;

% The ODEs of the model.
f1 = (phi1-mi1*a2*signapprox1*(2*g*(abs1))^(1/2)-(0*a2*((2*g*h1)^(1/2))))/a1;
f2 = (mi1*a2*signapprox1*((2*g*abs1)^(1/2))-mi1*a2*signapprox2*((2*g*(abs2))^(1/2))-(0*a2*((2*g*h2)^(1/2))))/a1;
f3 = (phi2+mi1*a2*signapprox2*(2*g*(abs2))^(1/2)-(mi2*a2*(2*g*h3)^(1/2)))/a1;

% The matrices of the state space model.
A = jacobian([f1;f2;f3],[h1 h2 h3]);
B = jacobian([f1;f2;f3],[phi1 phi2]);
C = [1 0 0;0 0 1];
D = [0 0;0 0];
% Calculating the steady-state of the system.
theta = [mi1 mi2 a1 a2];
T = 0:0.1:20;
phi1 = 0.3*10^-4;
phi2 = 0.4*10^-4;
u = [phi1;phi2];
U = u.*ones(1,length(T));
x_initial = [0,0,0];
%[t,y]= ode45(@(t,y)threetankode(t,y,u,theta),T,[0,0,0]);
x=simulation(@threetankode,T,x_initial,U,theta);
h3 = ((phi1+phi2)^2)/(2*g*(mi2^2)*(a2^2));
h2 = h3 + (phi1^2/(2*g*(mi1^2)*(a2^2)));
h1 = h2 + (phi1^2/(2*g*(mi1^2)*(a2^2)));

linpoints = x(:,end);
%h1 = linpoints(1);
%h2 = linpoints(2);
%h3 = linpoints(3);

A = eval(A);
B = eval(B);


tank = ss(A,B,C,D);
theta_linear = {tank.A,tank.B,linpoints,u};
U = u.*ones(1,length(T));
x_linear=simulation(@LinearOde,T,x_initial,U,theta_linear);
erro_1 = immse(x,x_linear);

figure()
plot(T,x)
hold on
plot(T,x_linear,'--')
xlabel('Time (s)','fontsize',12,'fontweight','bold')
ylabel('Water Level (m)','fontsize',12,'fontweight','bold')
legend('Tank 1','Tank 2','Tank 3','Tank 1 linear','Tank 2 linear','Tank 3 linear')
figure()


transient_response_data = stepinfo(tank);

pzmap(tank)