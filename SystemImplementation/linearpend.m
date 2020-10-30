function [A1,B1,A2,B2,C,D]=linearpend(yode)
  %% Variables of the system
  % m = Mass of the pendulum.
  % M = Mass of the cart.
  % u = Input.
  % l = Length of the pendulum.
  % g = Earth gravity.
  % x1 = Angular position.
  % x2 = Angular velocity.
  % x3 = Cart position.
  % x4 = Cart velocity.
  syms m M u L g b x1 x2 x3 x4 u
  %% ODEs of the system using the symbolic representation.
  term_1 = (M+m)/(M+m-m*(cos(x1))^2);
  term_2 = 1/(M+m-m*(cos(x1))^2);
  f1 = x2;
  f2 = term_1*((g/L)*sin(x1)+((-b*x4+u-m*L*x2^2*sin(x1)))/(L*(M+m))*cos(x1)-(b*x2/(m*L^2)));  
  %f2 = ((-b*x2/((l^2)*m))+(g*sin(x1)/l)+(-m*l*(x2^2)*sin(x1)*cos(x1)-b*x4*cos(x1))/((M+m)*l))/(M+m-m*(cos(x1)^2));
  f3 = x4;
  %f4 = (u+m*L*(cos(x1)*term_1*((g/L)*sin(x1)+((-b*x4+u-m*L*x2^2*sin(x1))/(L*(M+m)))*cos(x1)-b*x2/(m*L^2))-(x2^2)*sin(x1))-b*x4)/(M+m);
  f4 = term_2*(u+m*g*sin(x1)*cos(x1)-(b*x2/L)*cos(x1)-m*L*(x2^2)*sin(x1)-b*x4);
  %% Linearization of the system.
  A = jacobian([f1;f2;f3;f4],[x1 x2 x3 x4]);
  B = jacobian([f1;f2;f3;f4],[u]);
  C = [1 0 0 0;0 0 1 0];
  D = [0;0];
  %% Using the numerical values of the problem.
  m =2;
  M = 5;
  L=2;
  b=1;
  u=0;
  l= 2;
  g=9.81;
  %% Linearization point 1 [pi,0,0,0].
  x1 = pi;
  x2 = 0;
  x3 = 0;
  x4 = 0;
  linpoints1 = [x1 x2 x3 x4];
  A1 = eval(A);
  B1 = eval(B);
  %% Linearization point 2 [0,0,0,0].
  x1 = 0;
  x2 = 0;
  x3 = 0;
  x4 = 0;
  linpoints2 = [x1 x2 x3 x4];
  A2 = eval(A);
  B2 =eval(B);
  %invertedpendulum1 = ss(A1,B1,C,D);
  %invertedpendulum2 = ss(A2,B2,C,D);
  %figure()
  %linearpendsim(invertedpendulum1,linpoints1,0,0);
  %figure()
  %linearpendsim(invertedpendulum2,linpoints2,0,0);
end
function linearpendsim(invertedpend,linpoints,u1,ue1)
    T =0:0.01:10;
    U = u1*ones(1,length(T));
    UE = ue1*ones(1,length(T));
    y = lsim(invertedpend , U-UE , T , [5 0  0 0]-linpoints);
    figure()
    plot(T, y(:,2)+yode(end,3),'k')
    xlabel('Time','fontsize',16,'fontweight','bold')
    ylabel('Cart position ','fontsize',16,'fontweight','bold')
    title("Cart position using the linear model approximation")
    print('LinearCartPosition','-depsc')
    figure()
    plot(T,y(:,1)+pi,'k')
    xlabel('Time','fontsize',12,'fontweight','bold');
    ylabel('Pendulum angle','fontsize',12,'fontweight','bold')
    title("Pendulum angular position using the linear model approximation")
    print('LinearAngularPosition','-depsc')
    figure()   
end
