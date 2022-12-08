function dy = model_1(y,k,N)                        % DE
           beta = k(1);                                % Assignes the parameters in the DE the current
           gamma = k(2);                                % value of the parameters 
           dy = zeros(3,1);                         % assigns zeros to dy 
           dy(1) = -beta*y(1)*y(2)/N;                  % RHS of first equation
           dy(2) =  beta*y(1)*y(2)/N- gamma*y(2);          % RHS of second equation
           dy(3) = gamma*y(2);
end 
