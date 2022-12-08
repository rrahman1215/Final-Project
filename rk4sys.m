function [t,Y] = rk4sys(Fun,a,b,Y0,N)

h = (b-a)/N;
m = length(Y0); % number of equations
t = zeros(1,N+1); % initialize t
Y = zeros(m,N+1); % initialize Y
t(1) = a;
Y(:,1) = Y0; % the first row of Y is Y0

for i = 1:N
    t(i+1) = t(i) + h;
    K1 = h*Fun(t(i),Y(:,i));
    K2 = h*Fun(t(i)+h/2,Y(:,i)+K1/2);
    K3 = h*Fun(t(i)+h/2,Y(:,i)+K2/2);
    K4 = h*Fun(t(i)+h,Y(:,i)+K3);
    Y(:,i+1) = Y(:,i) + 1/6*(K1 + 2*K2 + 2*K3 + K4);
end
