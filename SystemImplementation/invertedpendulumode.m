% this code implements the ode of the inverted pendulum.
function dx = invertedpendulumode(y)
    m =2;
    M = 5;
    %L=2;
    b=1;
    u=0;
    L= 2;
    g =9.81;
    x1=y(1);
    x2=y(2);
    x3=y(3);
    x4=y(4);
    term_1 = (M+m)/(M+m-m*(cos(x1))^2);
    term_2 = 1/(M+m-m*(cos(x1))^2);
    f1 = x2;
    f2 = term_1*((g/L)*sin(x1)+((-b*x4+u-m*L*x2^2*sin(x1)))/(L*(M+m))*cos(x1)-b*x2/(m*L^2));    
    %f2 = ((-b*x2/((l^2)*m))+(g*sin(x1)/l)+(-m*l*(x2^2)*sin(x1)*cos(x1)-b*x4*cos(x1))/((M+m)*l))/(M+m-m*(cos(x1)^2));
    f3 = x4;
    %f4 =  (u+m*L*(cos(x1)*term_1*((g/L)*sin(x1)+((-b*x4+u-m*L*x2^2*sin(x1))/(L*(M+m)))*cos(x1)-b*x2/(m*L^2))-(x2^2)*sin(x1))-b*x4)/(M+m);
    f4 = term_2*(u+m*g*sin(x1)*cos(x1)-(b*x2/L)*cos(x1)-m*L*(x2^2)*sin(x1)-b*x4);
    
    %f4 = (m*l*(-b*x2/((l^2)*m)+(((g*sin(x1)/l)+(-m*l*(x2^2)*sin(x1)*cos(x1)-b*x4*cos(x1))/((M+m)*l))/(M+m-m*(cos(x1)^2)))*cos(x1)-m*l*(x2^2)*sin(x1)-b*x4))/(m+M);
    dx =[f1;f2;f3;f4];
    end
    
