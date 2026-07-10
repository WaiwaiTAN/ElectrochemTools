repo_root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(repo_root, 'DRTtools')));

if ~exist('DRT_INPUT_CSV', 'var')
    DRT_INPUT_CSV = fullfile(repo_root, 'examples', 'data', 'eis_cleaned.csv');
end
if ~exist('DRT_OUTPUT_DIR', 'var')
    DRT_OUTPUT_DIR = fullfile(repo_root, 'target', 'matlab_reference');
end
if ~exist('DRT_LAMBDA', 'var')
    DRT_LAMBDA = 1.0e-3;
end
if ~exist('DRT_REG_ORDER', 'var')
    DRT_REG_ORDER = 1;
end
if ~exist('DRT_FIT_INDUCTANCE', 'var')
    DRT_FIT_INDUCTANCE = false;
end

if exist('quadprog', 'file') ~= 2 || ~license('test', 'Optimization_Toolbox')
    error('Optimization Toolbox with quadprog is required for this DRTtools reference run.');
end

if ~isfolder(DRT_OUTPUT_DIR)
    mkdir(DRT_OUTPUT_DIR);
end

raw = readmatrix(DRT_INPUT_CSV, 'FileType', 'text');
if size(raw, 2) < 3
    error('Input CSV must contain at least frequency, Z_real, and Z_imag columns.');
end

freq = raw(:, 1);
b_re = raw(:, 2);
b_im = raw(:, 3);
valid = isfinite(freq) & isfinite(b_re) & isfinite(b_im) & freq > 0;
freq = freq(valid);
b_re = b_re(valid);
b_im = b_im(valid);

[freq, order] = sort(freq, 'descend');
b_re = b_re(order);
b_im = b_im(order);

epsilon = 0;
rbf_type = 'Piecewise linear';
A_re = assemble_A_re(freq, epsilon, rbf_type);
A_im = assemble_A_im(freq, epsilon, rbf_type);
A_re(:, 2) = 1;
if DRT_FIT_INDUCTANCE
    A_im(:, 1) = 2*pi*freq;
end

if DRT_REG_ORDER == 1
    M = assemble_M_1(freq, epsilon, rbf_type);
elseif DRT_REG_ORDER == 2
    M = assemble_M_2(freq, epsilon, rbf_type);
else
    error('DRT_REG_ORDER must be 1 or 2 for this reference script.');
end

[H, c] = quad_format_combined(A_re, A_im, b_re, b_im, M, DRT_LAMBDA);
lb = zeros(numel(freq) + 2, 1);
ub = Inf(numel(freq) + 2, 1);
x0 = ones(numel(freq) + 2, 1);
if ~DRT_FIT_INDUCTANCE
    ub(1) = 0;
    x0(1) = 0;
end

options = optimoptions('quadprog', 'Display', 'off');
[x, fval, exitflag, output] = quadprog(H, c(:), [], [], [], [], lb, ub, x0, options);
if exitflag <= 0
    error('quadprog failed with exitflag %d.', exitflag);
end

tau = 1 ./ freq;
gamma = x(3:end);
z_fit_re = A_re * x;
z_fit_im = A_im * x;

drt_path = fullfile(DRT_OUTPUT_DIR, 'drttools_reference_drt.csv');
fid = fopen(drt_path, 'w');
fprintf(fid, 'L, %.12e\n', x(1));
fprintf(fid, 'R, %.12e\n', x(2));
fprintf(fid, 'tau, gamma(tau)\n');
for i = 1:numel(tau)
    fprintf(fid, '%.12e, %.12e\n', tau(i), gamma(i));
end
fclose(fid);

regression_path = fullfile(DRT_OUTPUT_DIR, 'drttools_reference_eis_regression.txt');
fid = fopen(regression_path, 'w');
fprintf(fid, 'freq,mu_Z_re,mu_Z_im,Z_re,Z_im\n');
for i = 1:numel(freq)
    fprintf(fid, '%.12e,%.12e,%.12e,%.12e,%.12e\n', freq(i), z_fit_re(i), z_fit_im(i), b_re(i), b_im(i));
end
fclose(fid);

summary.lambda = DRT_LAMBDA;
summary.regularization_order = DRT_REG_ORDER;
summary.fit_inductance = logical(DRT_FIT_INDUCTANCE);
summary.n_points = numel(freq);
summary.inductance = x(1);
summary.r_inf = x(2);
summary.polarization_resistance = trapz(log(tau), gamma);
summary.quadprog_fval = fval;
summary.quadprog_iterations = output.iterations;
summary.drt_path = drt_path;
summary.regression_path = regression_path;

summary_path = fullfile(DRT_OUTPUT_DIR, 'drttools_reference_summary.json');
fid = fopen(summary_path, 'w');
fprintf(fid, '%s', jsonencode(summary, 'PrettyPrint', true));
fclose(fid);

fprintf('DRTtools reference written to %s\n', DRT_OUTPUT_DIR);
