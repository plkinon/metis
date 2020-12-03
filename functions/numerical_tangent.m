function K_num = numerical_tangent(this_integrator,this_problem,zn1, zn)
%% FUNCTION: numerical tangent
% Computes numerical tangent matrix for a given residual defined by
% integrator and problem zn1 and zn

% Pre-allocate tangent matrix, initialized with zeros:
N       = length(zn1);
K_num   = zeros(N,N);

% Define epsilon which manipulates solution vector
epsilon = 1e-10;

% For-loop which computes column by column of numerical tangent matrix
for j=1:N
  
    % Save current entry of solution vector
    zsave=zn1(j);

    % Increment the jth component of the solution vector 
    delp=epsilon*(1.0+abs(zn1(j)));
    zn1(j)=zsave+delp;
    [R1,~] = this_integrator.compute_resi_tang(zn1, zn, this_problem);

    % Decrement the jth component of the solution vector
    zn1(j)=zsave-delp;
    [R2,~] = this_integrator.compute_resi_tang(zn1, zn, this_problem);

    % Restore the original vector 
    zn1(j)=zsave;

    % Compute the approximate tangent matrix entry
    K_num(:,j)=(R1-R2)/(2*delp);
  
end

end