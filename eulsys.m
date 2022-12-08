function [t,Y] = eulsys(Fun,a,b,Y0,N)

h = (b-a)/N;
m = length(Y0); % number of equations
t = zeros(1,N+1); % initialize t
Y = zeros(m,N+1); % initialize Y
t(1) = a;
Y(:,1) = Y0; % the first column of Y is Y0

for i = 1:N
    t(i+1) = t(i) + h;
    Y(:,i+1) = Y(:,i) + h*Fun(t(i),Y(:,i));
end
