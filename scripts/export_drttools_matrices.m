repo_root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(repo_root, 'DRTtools', 'src')));

if ~exist('DRT_MATRIX_OUTPUT_DIR', 'var')
    DRT_MATRIX_OUTPUT_DIR = fullfile(repo_root, 'tests', 'golden', 'drttools', 'matrix_piecewise_linear');
end
if ~isfolder(DRT_MATRIX_OUTPUT_DIR)
    mkdir(DRT_MATRIX_OUTPUT_DIR);
end

freq = [1000; 100; 10; 1; 0.1];
tau = 1 ./ freq;
r_inf = 0.8;
inductance = 2.0e-5;
resistance = 12.0;
tau_peak = 1.0e-2;
omega = 2*pi*freq;
z = r_inf + 1i*omega*inductance + resistance ./ (1 + 1i*omega*tau_peak);
b_re = real(z);
b_im = imag(z);
lambda = 1.0e-3;
epsilon = 0;
rbf_type = 'Piecewise linear';

A_re = assemble_A_re(freq, epsilon, rbf_type);
A_im = assemble_A_im(freq, epsilon, rbf_type);
A_re(:, 2) = 1;
A_im(:, 1) = omega;
M_1 = assemble_M_1(freq, epsilon, rbf_type);
M_2 = assemble_M_2(freq, epsilon, rbf_type);
[H, c] = quad_format_combined(A_re, A_im, b_re, b_im, M_1, lambda);
x_ridge = (H / 2) \ (-c(:) / 2);
gamma = x_ridge(3:end);
z_reconstructed = A_re*x_ridge + 1i*(A_im*x_ridge);

writematrix(freq, fullfile(DRT_MATRIX_OUTPUT_DIR, 'frequency.csv'));
writematrix(tau, fullfile(DRT_MATRIX_OUTPUT_DIR, 'tau.csv'));
writematrix(A_re, fullfile(DRT_MATRIX_OUTPUT_DIR, 'A_re.csv'));
writematrix(A_im, fullfile(DRT_MATRIX_OUTPUT_DIR, 'A_im.csv'));
writematrix(M_1, fullfile(DRT_MATRIX_OUTPUT_DIR, 'M_1.csv'));
writematrix(M_2, fullfile(DRT_MATRIX_OUTPUT_DIR, 'M_2.csv'));
writematrix(H, fullfile(DRT_MATRIX_OUTPUT_DIR, 'H.csv'));
writematrix(c(:), fullfile(DRT_MATRIX_OUTPUT_DIR, 'c.csv'));
writematrix(x_ridge, fullfile(DRT_MATRIX_OUTPUT_DIR, 'x_ridge.csv'));
writematrix(gamma, fullfile(DRT_MATRIX_OUTPUT_DIR, 'gamma.csv'));
writematrix([real(z_reconstructed), imag(z_reconstructed)], fullfile(DRT_MATRIX_OUTPUT_DIR, 'z_reconstructed.csv'));
writematrix([b_re, b_im], fullfile(DRT_MATRIX_OUTPUT_DIR, 'observations.csv'));

metadata.matlab_version = version;
metadata.drttools_commit = '034d9c4c4a4916a38a0e2f10381d931ffe1981b3';
metadata.rbf_type = rbf_type;
metadata.lambda = lambda;
metadata.fit_inductance = true;
metadata.regularization_order = 1;
metadata.generated_at = char(datetime('now', 'TimeZone', 'UTC', 'Format', 'yyyy-MM-dd''T''HH:mm:ssXXX'));
fid = fopen(fullfile(DRT_MATRIX_OUTPUT_DIR, 'metadata.json'), 'w');
fprintf(fid, '%s', jsonencode(metadata, 'PrettyPrint', true));
fclose(fid);

fprintf('DRTtools matrix golden references written to %s\n', DRT_MATRIX_OUTPUT_DIR);
