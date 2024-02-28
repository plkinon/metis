t=simulation.t';
x_1 = simulation.z(:,1)+1/2*simulation.z(:,3);
v_1 = simulation.z(:,61)+1/2*simulation.z(:,63);
H   = simulation.H-simulation.H(1);
g_q = vecnorm(simulation.constraint_position')';
g_v = vecnorm(simulation.constraint_velocity')';
results = [t, x_1, v_1, H, g_q, g_v];

fileID = fopen('results.txt','w');
fprintf(fileID,'%5.3f %12.8f %12.8f %12.8e %12.8e %12.8e -not provided- \n', results');
fclose(fileID);
