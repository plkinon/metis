%% Export results for comparison with pydykit
time = simulation.t';
coordinates = simulation.z(:,1:12);
velocities = simulation.z(:,13:24);
multiplier = simulation.z(:,25:26);
save("results","time","coordinates","velocities","multiplier")