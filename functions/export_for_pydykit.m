%% Export results for comparison with pydykit
time = simulation.t';
coordinates = simulation.z(:,1:4);
momenta = simulation.z(:,5:8);
multiplier = simulation.z(:,9);
save("results","time","coordinates","momenta","multiplier")