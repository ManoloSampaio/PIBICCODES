%% This code is the ode for the three tank with the data of the report.
function dx = threetankode(t,y,u,theta)
    %% Gravity.
    g = 9.81;
   %% Set the level values.
    h1 = y(1);
    h2 = y(2);
    h3 = y(3);
    %% Set the leak coeficientes.
    mi1 = theta(1);
    mi2 = theta(2);
    %% Set input values.
    phi1 = u(1);
    phi2 = u(2);
    %% Set mi1 value.
    %% Set a1,a2 values.
    a1 = theta(3);
    a2 = theta(4);
    k= 91.1;
    %% Sign function approximation.
    sign1 = -1 + 2*((1+exp(-2*(k)*(h1-h2))).^-1);
    sign2 = -1 + 2*((1+exp(-2*(k)*(h2-h3))).^-1);
    %% Absolute value function approximation.
    abs1 = ((h1-h2)^2+0.0001)^0.5;
    abs2 = ((h2-h3)^2+0.0001)^0.5;
    %% ODEs.
    d1 = (phi1-mi1*a2*sign1*(2*g*(abs1))^(1/2))/a1;
    d2 = (mi1*a2*sign1*((2*g*abs1)^(1/2))-mi1*a2*sign2*((2*g*(abs2))^(1/2)))/a1;
    d3 = (phi2+mi1*a2*sign2*(2*g*(abs2))^(1/2)-(mi2*a2*(2*g*h3)^(1/2)))/a1;
    dx = [d1;d2;d3];
end