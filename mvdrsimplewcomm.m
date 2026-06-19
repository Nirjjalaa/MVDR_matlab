clc;
clear;
close all;

%% =========================
% ARRAY PARAMETERS
% =========================
N = 64;              % number of elements
d = 0.5;             % spacing (lambda/2)

theta_target = 30;
theta_interf  = 60;

num_snap = 2000;

%% =========================
% INDEX VECTOR
% =========================
n = (0:N-1)';

%% =========================
% STEERING VECTORS
% =========================
a_target = steering_vec(theta_target, n, d);
a_interf = steering_vec(theta_interf, n, d);

%% =========================
% SIGNAL GENERATION
% =========================
target_power = 1;
interf_power = 10;
noise_power = 0.1;

s_target = sqrt(target_power/2) * (randn(1,num_snap) + 1j*randn(1,num_snap));
s_interf = sqrt(interf_power/2) * (randn(1,num_snap) + 1j*randn(1,num_snap));
noise    = sqrt(noise_power/2) * (randn(N,num_snap) + 1j*randn(N,num_snap));

%% =========================
% RECEIVED SIGNAL
% =========================
X = a_target * s_target + a_interf * s_interf + noise;

%% =========================
% COVARIANCE MATRIX
% =========================
R = (X * X') / num_snap;

%% =========================
% DIAGONAL LOADING
% =========================
delta = 0.01 * trace(R) / N;
R = R + delta * eye(N);

%% =========================
% MVDR WEIGHT COMPUTATION
% Solve: Rw = a_target
% =========================
w = R \ a_target;

% normalization (distortionless constraint)
w = w / (a_target' * w);

%% =========================
% DIAGNOSTICS
% =========================
disp('Gain toward target:');
disp(abs(w' * a_target));

disp('Gain toward interference:');
disp(abs(w' * a_interf));

%% =========================
% BEAM PATTERN
% =========================
theta_scan = -90:0.2:90;
beam = zeros(size(theta_scan));

for k = 1:length(theta_scan)
    a = steering_vec(theta_scan(k), n, d);
    beam(k) = abs(w' * a);
end

beam = beam / max(beam);
beam_dB = 20*log10(beam + 1e-12);

%% =========================
% PLOT
% =========================
figure;
plot(theta_scan, beam_dB, 'b', 'LineWidth', 2);
grid on;

xlabel('Angle (deg)');
ylabel('Normalized Gain (dB)');
title('Clean MVDR Beamformer (Reference Model)');
ylim([-60 0]);

hold on;
xline(theta_target, 'g--', 'Target');
xline(theta_interf, 'r--', 'Interference');
legend('Beam Pattern','Target','Interference');

function a = steering_vec(theta, n, d)

% theta in degrees
theta = deg2rad(theta);

% spatial frequency
phase_shift = 2*pi*d*sin(theta);

% steering vector
a = exp(1j * phase_shift * n);

% normalization (important for FPGA consistency)
a = a / sqrt(length(n));

end
