clc;
clear;
close all;

%% Parameters
N = 64;                  % Number of array elements
d = 0.5;                 % Element spacing (lambda/2)

theta_target = 30;       % Desired direction
theta_interf = 60;       % Interference direction

num_snap = 1000;         % Number of snapshots

%% Array Steering Vectors
n = (0:N-1)';

a_target = exp(1j*2*pi*d*sind(theta_target)*n);
a_interf = exp(1j*2*pi*d*sind(theta_interf)*n);

%% Generate Signals

target_power = 1;
interf_power = 10;
noise_power = 0.1;

s_target = sqrt(target_power/2) * ...
    (randn(1,num_snap)+1j*randn(1,num_snap));

s_interf = sqrt(interf_power/2) * ...
    (randn(1,num_snap)+1j*randn(1,num_snap));

noise = sqrt(noise_power/2) * ...
    (randn(N,num_snap)+1j*randn(N,num_snap));

%% Received Data

X = a_target*s_target + a_interf*s_interf + noise;

%% Covariance Matrix

R = (X*X')/num_snap;

%% Diagonal Loading

delta = 0.01*trace(R)/N;
R = R + delta*eye(N);

%% MVDR Weight Computation

w_mvdr = R\a_target;
w_mvdr = w_mvdr/(a_target'*w_mvdr);

%% Apply Hamming Window

ham = hamming(N);

w_mvdr = w_mvdr .* ham;

% Renormalize
w_mvdr = w_mvdr/(a_target'*w_mvdr);

%% Diagnostics

disp('Gain towards Target')
disp(w_mvdr'*a_target)

disp('Gain towards Interference')
disp(abs(w_mvdr'*a_interf))

%% Beam Pattern

theta_scan = -90:0.1:90;

beam = zeros(size(theta_scan));

for k = 1:length(theta_scan)

    a_scan = exp(1j*2*pi*d*sind(theta_scan(k))*n);

    beam(k) = abs(w_mvdr'*a_scan);

end

beam = beam/max(beam);

beam_dB = 20*log10(beam);

beam_dB(beam_dB<-60) = -60;

%% Find Beam Peak

[~,idx] = max(beam);

fprintf('Maximum Beam at %.2f degrees\n',theta_scan(idx));

%% Plot

figure;

plot(theta_scan,beam_dB,'b','LineWidth',2);

hold on;

xline(theta_target,'g--','LineWidth',2);

xline(theta_interf,'r--','LineWidth',2);

grid on;

xlabel('Angle (degrees)');
ylabel('Normalized Magnitude (dB)');
title('64-Element MVDR Beamformer with Hamming Window');

legend('Beam Pattern','Target','Interference');

axis([-90 90 -60 0]);

%% Display Information

disp('------------------------------------');
disp(['Data Matrix Size      : ',mat2str(size(X))]);
disp(['Covariance Size       : ',mat2str(size(R))]);
disp(['Weight Vector Size    : ',mat2str(size(w_mvdr))]);
disp(['Weight Vector Norm    : ',num2str(norm(w_mvdr))]);
