clc
clear
close all

theta1D=30*pi/180; %desired angle for theta1 
theta2D=30*pi/180; %desired angle for theta2

[t, y]=ode45(@EQTNS, [0 100],[0 0 0 0 0 0]); %calling the function for the equations

plot(t, y(:,1)) %theta1 plot
hold on
plot(t, y(:,2)) %theta2 plot
hold on
yline(theta1D,'--');
hold on
yline(theta2D,'--');
hold off 

function dydt = EQTNS(t, y)
theta1D=30*pi/180;
theta2D=30*pi/180;
m1=5; %mass 1
m2=4; %mass 2
g=9.8; %acceleration of gravity
L1=1; %length 1
L2=0.5; %length 2

%Trajectory equations for first angle: 

a_note1=0; %degrees
tf=50;
a21=(3/tf^2)*(theta1D-a_note1);
a31=(-2/tf^3)*(theta1D-a_note1);
position1=a_note1+a21*t^2+a31*t^3;
dposition1=2*a21*t+3*a31*t^2;
intposition1=a_note1*t+(a21*t^3)/3+(a31*t^4)/4;

%Trajectory equations for second angle :

a_note2=0; %degrees
a22=(3/tf^2)*(theta2D-a_note2);
a32=(-2/tf^3)*(theta2D-a_note2);
position2=a_note2+a22*t^2+a32*t^3;
dposition2=2*a22*t+3*a32*t^2; 
intposition2=a_note2*t+(a22*t^3)/3+(a32*t^4)/4;

if t > tf
    position1=a_note1+a21*tf^2+a31*tf^3;
    dposition1=2*a21*tf+3*a31*tf^2;
    intposition1=a_note1*tf+(a21*tf^3)/3+(a31*tf^4)/4 + theta1D*(t-tf);

    position2=a_note2+a22*tf^2+a32*tf^3;
    dposition2=2*a22*tf+3*a32*tf^2;
    intposition2=a_note2*tf+(a22*tf^3)/3+(a32*tf^4)/4 + theta2D*(t-tf);
end

% Without trajectory following

% position1=theta1D;
% dposition1=0; 
% intposition1=theta1D*t;
% 
% position2=theta2D;
% dposition2=0; 
% intposition2=theta2D*t;

%y(1) = theta1
%y(2) = theta2
%y(3) = theta1dot
%y(4) = theta2dot

Kp1 = 100; 
Kd1 = 500;
Kp2 = 100;
Kd2 = 500;
Ki1= 20;
Ki2= 20;
    
tau1=Kp1*(position1-y(1)) + Kd1*(dposition1-y(3))+ Ki1*(intposition1-y(5));
tau2=Kp2*(position2-y(2)) + Kd2*(dposition2-y(4))+ Ki2*(intposition2-y(6));

%In matrix form, the equations can be written as Mass*thetaddot = V

Mass11 = m1*L1^2+m2*(L1^2+2*L1*L2*cos(y(2))+L2^2);
Mass12 = m2*(L1*L2*cos(y(2))+L2^2);
Mass21 = Mass12;
Mass22 = m2*L2^2;

detMass = Mass11*Mass22 - Mass12*Mass21; % determinant of the mass matrix

V1 = tau1 + m2*L1*L2*sin(y(2))*(2*y(3)*y(4)+y(4)^2)... 
    - (m1+m2)*L1*g*cos(y(1))...
    - m2*g*L2*cos(y(1)+y(2));
V2 = tau2 - m2*L1*L2*y(3)^2*sin(y(2)) - m2*g*L2*cos(y(1)+y(2));

dydt(1)=y(3); %first derivative of theta1
dydt(2)=y(4); %first derivative of theta2
dydt(3)=(Mass22*V1 - Mass12*V2)/detMass; %second derivative of theta1
dydt(4)=(-Mass21*V1 + Mass11*V2)/detMass; %second derivative of theta2
dydt(5)=y(1);
dydt(6)=y(2);

dydt=dydt';
end

