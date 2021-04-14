clear all;

##function Xdot = fpap(t, X)
##  %A = 0.15e-3;
##  A = 0;
##  theta0 = [ -pi/3 -pi/12 0 pi/12 pi/2 ]';
##  a = [ 1.2 -5 30 -7.5 0.75 ]';
##  b = [ 0.25 0.1 0.1 0.1 0.4 ]';
##  f2 = 0.5;
##  z0 = A*sin(2*pi*f2*t); %baseline wander breath
##  alpha = 1 - sqrt(X(1)^2 + X(2)^2);
##  f0 = 1;
##  omega = 2*pi*f0;
##  theta = atan2(X(2), X(1));
##  dth = theta - theta0;
##  
##  Xdot = [ alpha*X(1) - omega*X(2);
##           alpha*X(2) + omega*X(1);
##          -sum(a.*dth.*exp(-dth.^2./(2*b.^2))) - (X(3) - z0) ];
##end

function Xdot = f2nd(t, X)
  %A = 0.15e-3;
  A = 0;
  theta0 = [ -pi/3 -pi/12 0 pi/12 pi/2 ]';
  a = [ 1.2 -5 30 -7.5 0.75 ]';
  b = [ 0.25 0.1 0.1 0.1 0.4 ]';
  f2 = 0.5;
  z0 = A*sin(2*pi*f2*t); %baseline wander breath
  f0 = 1;
  omega = 2*pi*f0;
  theta = atan2(X(2), X(1));
  dth = theta - theta0;
  
  Xdot = [ -omega*X(2);
           omega*X(1);
          -sum(a.*dth.*exp(-dth.^2./(2*b.^2))) - (X(3) - z0) 
          0];
end

function Xk1 = fd(X,T)
  %A = 0.15e-3;
  A = 0;
  theta0 = [ -pi/3 -pi/12 0 pi/12 pi/2 ]';
  a = [ 1.2 -5 30 -7.5 0.75 ]';
  b = [ 0.25 0.1 0.1 0.1 0.4 ]';
  f2 = 0.5;
  z0 = 0; %baseline wander breath

  theta = atan2(X(2), X(1));
  dth = theta - theta0;
  
  Xdot = [ -X(4)*X(2);
           X(4)*X(1);
          -sum(a.*dth.*exp(-dth.^2./(2*b.^2))) - (X(3) - z0) 
          0];
          
  Xk1 = T*Xdot + X;
end

Amtx;

function A=Adf(X,dt)
  Ab = 0;
  theta0 = [ -pi/3 -pi/12 0 pi/12 pi/2 ]';
  a = [ 1.2 -5 30 -7.5 0.75 ]';
  b = [ 0.25 0.1 0.1 0.1 0.4 ]';
  f2 = 0.5;
  z0 = 0; %baseline wander breath
  
  theta = atan2(X(2), X(1));
  dth = theta - theta0;
  csi = exp(-dth.^2./(2*b.^2));
  m2 = X(1)^2 + X(2)^2;
  sdth1 = -X(2)/m2;
  sdth2 = X(1)/m2;
  pf = (dth./b).^2 - 1;
  X3dX1 = sum(a.*csi.*sdth1.*pf);
  X3dX2 = sum(a.*csi.*sdth2.*pf);
  Ac = [  0    -X(4)   0    -X(2) ;
         X(4)    0     0     X(1) ;
        X3dX1  X3dX2  -1      0   ;
          0      0     0      0  ];
         
  A = eye(4) + Ac*dt;
end


theta0 = [ -pi/3 -pi/12 0 pi/12 pi/2 ]';
a = [ 1.2 -5 30 -7.5 0.75 ]';
b = [ 0.25 0.1 0.1 0.1 0.4 ]';
f0 = 1;
omega = 2*pi*f0;



X0 = [ 0 1 0 omega]';
X0 = rand(4,1); 
SR = 100;
T = 10;
L = SR*T;
tc = (0:L-1)/SR;


Ad = @(X,a,b,th) eye(4) + (A(X,a,b,th)/SR);

[t2nd, X2nd] = ode45(@f2nd, tc, X0);

Ad(X0, a, b, theta0)
Adf(X0, 1/SR)

Xd = zeros(L,4);
Xd(1,:) = X0';
for i = 1:L-1
  Xd(i+1,:) = fd(Xd(i,:)', 1/SR);
end

%figure(1); plot(t2nd, X2nd(:,3), tc, Xd(:,3)); legend;


Xf = fft(X2nd(:,3));
%figure(2); stem(abs(Xf));
