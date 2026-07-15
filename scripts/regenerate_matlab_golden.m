repo_root = fileparts(fileparts(mfilename('fullpath')));
drttools_root = fullfile(repo_root, 'DRTtools');
addpath(genpath(drttools_root));

if exist('quadprog', 'file') ~= 2 || ~license('test', 'Optimization_Toolbox')
    error('Optimization Toolbox with quadprog is required.');
end

input_relative = fullfile('tests', 'fixtures', 'eis_cleaned.csv');
input_path = fullfile(repo_root, input_relative);
output_root = fullfile(repo_root, 'tests', 'golden', 'drttools', 'end_to_end');
cases = {
    struct('name', 'real_order1_no_inductance', 'order', 1, 'fit_inductance', false), ...
    struct('name', 'real_order1_with_inductance', 'order', 1, 'fit_inductance', true), ...
    struct('name', 'real_order2_no_inductance', 'order', 2, 'fit_inductance', false)
};
lambda = 1.0e-3;

raw = readmatrix(input_path, 'FileType', 'text');
freq = raw(:, 1); b_re = raw(:, 2); b_im = raw(:, 3);
valid = isfinite(freq) & isfinite(b_re) & isfinite(b_im) & freq > 0;
freq = freq(valid); b_re = b_re(valid); b_im = b_im(valid);
[freq, order_index] = sort(freq, 'descend');
b_re = b_re(order_index); b_im = b_im(order_index);

for case_index = 1:numel(cases)
    config = cases{case_index};
    case_dir = fullfile(output_root, config.name);
    if ~isfolder(case_dir); mkdir(case_dir); end

    A_re = assemble_A_re(freq, 0, 'Piecewise linear');
    A_im = assemble_A_im(freq, 0, 'Piecewise linear');
    A_re(:, 2) = 1;
    if config.fit_inductance; A_im(:, 1) = 2*pi*freq; end
    if config.order == 1
        M = assemble_M_1(freq, 0, 'Piecewise linear');
    else
        M = assemble_M_2(freq, 0, 'Piecewise linear');
    end
    [H, c] = quad_format_combined(A_re, A_im, b_re, b_im, M, lambda);
    lb = zeros(numel(freq) + 2, 1);
    ub = Inf(numel(freq) + 2, 1);
    x0 = ones(numel(freq) + 2, 1);
    if ~config.fit_inductance; ub(1) = 0; x0(1) = 0; end
    options = optimoptions('quadprog', 'Display', 'off', ...
        'Algorithm', 'interior-point-convex', 'OptimalityTolerance', 1e-10);
    [x, fval, exitflag, output] = quadprog(H, c(:), [], [], [], [], lb, ub, x0, options);
    if exitflag <= 0; error('quadprog failed for %s with exitflag %d', config.name, exitflag); end

    tau = 1 ./ freq;
    gamma = x(3:end);
    z_fit = [A_re*x, A_im*x];
    writematrix(freq, fullfile(case_dir, 'frequency.csv'));
    writematrix(tau, fullfile(case_dir, 'tau.csv'));
    writematrix(x, fullfile(case_dir, 'coefficients.csv'));
    writematrix(gamma, fullfile(case_dir, 'gamma.csv'));
    writematrix(z_fit, fullfile(case_dir, 'reconstructed_impedance.csv'));

    summary.lambda = lambda;
    summary.regularization_order = config.order;
    summary.fit_inductance = logical(config.fit_inductance);
    summary.constraints = struct('gamma_nonnegative', true, 'r_inf_nonnegative', true, ...
        'inductance_mode', ternary(config.fit_inductance, 'nonnegative', 'fixed_zero'));
    summary.objective_value = fval;
    summary.R_inf = x(2);
    summary.inductance = x(1);
    summary.polarization_resistance = trapz(log(tau), gamma);
    summary.quadprog_exit_flag = exitflag;
    summary.quadprog_iterations = output.iterations;
    write_json(fullfile(case_dir, 'summary.json'), summary);

    metadata.MATLAB_version = version;
    metadata.Optimization_Toolbox_available = logical(license('test', 'Optimization_Toolbox'));
    metadata.DRTtools_commit = '034d9c4c4a4916a38a0e2f10381d931ffe1981b3';
    metadata.generation_script = 'scripts/regenerate_matlab_golden.m';
    metadata.generation_timestamp = char(datetime('now', 'TimeZone', 'UTC', ...
        'Format', 'yyyy-MM-dd''T''HH:mm:ssXXX'));
    metadata.input_fixture = strrep(input_relative, '\', '/');
    metadata.case_configuration = config;
    write_json(fullfile(case_dir, 'metadata.json'), metadata);
end
fprintf('Generated %d MATLAB golden cases under %s\n', numel(cases), output_root);

% DRTtools Simple Run defaults requested by the Rust Gaussian golden test:
% Gaussian RBF, FWHM coefficient 0.5, combined Re/Im, lambda 1e-3,
% first derivative regularization, nonnegative coefficients, and L fixed to 0.
gaussian_name = 'gaussian_simple_run';
gaussian_dir = fullfile(output_root, gaussian_name);
if ~isfolder(gaussian_dir); mkdir(gaussian_dir); end
rbf_type = 'Gaussian';
shape_control = 'FWHM Coefficient';
shape_coefficient = 0.5;
epsilon = compute_epsilon(freq, shape_coefficient, rbf_type, shape_control);
A_re = assemble_A_re(freq, epsilon, rbf_type);
A_im = assemble_A_im(freq, epsilon, rbf_type);
A_re(:, 2) = 1;
M = assemble_M_1(freq, epsilon, rbf_type);
[H, c] = quad_format_combined(A_re, A_im, b_re, b_im, M, lambda);
lb = zeros(numel(freq) + 2, 1);
ub = Inf(numel(freq) + 2, 1);
ub(1) = 0;
x0 = ones(numel(freq) + 2, 1);
x0(1) = 0;
options = optimoptions('quadprog', 'Display', 'off', ...
    'Algorithm', 'interior-point-convex', 'OptimalityTolerance', 1e-10);
[x, fval, exitflag, output] = quadprog(H, c(:), [], [], [], [], lb, ub, x0, options);
if exitflag <= 0; error('quadprog failed for %s with exitflag %d', gaussian_name, exitflag); end

tau = 1 ./ freq;
[gamma, ~] = map_array_to_gamma(freq, freq, x(3:end), epsilon, rbf_type);
gamma = gamma(:);
z_fit = [A_re*x, A_im*x];
writematrix(freq, fullfile(gaussian_dir, 'frequency.csv'));
writematrix(tau, fullfile(gaussian_dir, 'tau.csv'));
writematrix(A_re, fullfile(gaussian_dir, 'A_re.csv'));
writematrix(A_im, fullfile(gaussian_dir, 'A_im.csv'));
writematrix(M, fullfile(gaussian_dir, 'M_1.csv'));
writematrix(x, fullfile(gaussian_dir, 'coefficients.csv'));
writematrix(gamma, fullfile(gaussian_dir, 'gamma.csv'));
writematrix(z_fit, fullfile(gaussian_dir, 'reconstructed_impedance.csv'));

clear summary metadata;
summary.lambda = lambda;
summary.regularization_order = 1;
summary.fit_inductance = false;
summary.basis = rbf_type;
summary.shape_control = shape_control;
summary.shape_coefficient = shape_coefficient;
summary.epsilon = epsilon;
summary.constraints = struct('gamma_nonnegative', true, 'r_inf_nonnegative', true, ...
    'inductance_mode', 'fixed_zero');
summary.objective_value = fval;
summary.R_inf = x(2);
summary.inductance = x(1);
summary.polarization_resistance = trapz(log(tau), gamma);
summary.quadprog_exit_flag = exitflag;
summary.quadprog_iterations = output.iterations;
write_json(fullfile(gaussian_dir, 'summary.json'), summary);

metadata.MATLAB_version = version;
metadata.Optimization_Toolbox_available = logical(license('test', 'Optimization_Toolbox'));
metadata.DRTtools_commit = '034d9c4c4a4916a38a0e2f10381d931ffe1981b3';
metadata.generation_script = 'scripts/regenerate_matlab_golden.m';
metadata.generation_timestamp = char(datetime('now', 'TimeZone', 'UTC', ...
    'Format', 'yyyy-MM-dd''T''HH:mm:ssXXX'));
metadata.input_fixture = strrep(input_relative, '\', '/');
metadata.case_configuration = struct('basis', rbf_type, 'data_used', 'Combined Re-Im', ...
    'shape_control', shape_control, 'shape_coefficient', shape_coefficient, ...
    'lambda', lambda, 'regularization_order', 1, 'nonnegative', true, ...
    'fit_inductance', false);
write_json(fullfile(gaussian_dir, 'metadata.json'), metadata);
fprintf('Generated MATLAB Gaussian Simple Run golden under %s\n', gaussian_dir);

function write_json(path, value)
fid = fopen(path, 'w');
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, '%s', jsonencode(value, 'PrettyPrint', true));
end

function result = ternary(condition, yes, no)
if condition; result = yes; else; result = no; end
end
