syms I_1 I_3 theta p_theta p_psi p_phi

M= [I_1 0 0;
    0 I_3 I_3*cos(theta);
    0 I_3*cos(theta) I_3*cos(theta)^2 + I_1*sin(theta)^2];

Minv = eye(3) / M;

p = [p_theta; p_psi; p_phi];

T = 1/2*p'*Minv*p;

assume(p_psi,'real')
assume(p_phi,'real')
assume(p_theta,'real')

D_theta_T = diff(T,theta);