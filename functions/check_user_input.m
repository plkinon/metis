function check_user_input(CONFIG)

%% Check if user input is valid
directory_integrator = 'classes/Integrator';
directory_system     = 'classes/System';
directory_solver     = 'classes/Solver';
valid_integrator_classes = dir(fullfile(directory_integrator,'*.m'));
valid_system_classes     = dir(fullfile(directory_system,'*.m'));
valid_solver_classes     = dir(fullfile(directory_solver,'*.m'));

is_correct_integrator = false;
is_correct_system     = false;
is_correct_solver     = false;

for i = 1:numel(valid_integrator_classes)
    if strcmp(valid_integrator_classes(i).name,[CONFIG.INTEGRATOR,'.m'])
        is_correct_integrator = true;
    end
end

for j = 1:numel(valid_system_classes)
    if strcmp(valid_system_classes(j).name,[CONFIG.SYSTEM,'.m'])
        is_correct_system = true;
    end
end

for k = 1:numel(valid_solver_classes)
    if strcmp(valid_solver_classes(k).name,[CONFIG.SOLVER,'.m'])
        is_correct_solver = true;
    end
end

if ~is_correct_integrator 
    error('User input for integrator not available.');
elseif ~is_correct_system 
    error('User input for system not available.');
elseif ~is_correct_solver
    error('User input for solver not available.');
else
    fprintf('               Valid user input.                   \n');
    fprintf('************************************************** \n');
end

end