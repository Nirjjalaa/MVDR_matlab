clear; clc; close all;

%% =========================
% PARAMETERS
% =========================
c = physconst('LightSpeed');
fc = 1e9;
lambda = c/fc;
d = lambda/2;

M = 8;
N = 8;
MN = M*N;

%% =========================
% ANGLES
% =========================
az_s = 40;   el_s = 0;
az_i1 = 30;  el_i1 = 10;
az_i2 = -60; el_i2 = -25;

%% =========================
% SNAPSHOTS
% =========================
K = 50000;
t = (0:K-1)/K;

s  = randn(1,K) + 1j*randn(1,K);
i1 = randn(1,K) + 1j*randn(1,K);
i2 = randn(1,K) + 1j*randn(1,K);

%% =========================
% NOISE
% =========================
noise = sqrt(0.1/2) * (randn(MN,K) + 1j*randn(MN,K));

%% =========================
% STEERING VECTOR (URA)
% =========================
steer = @(az,el) steering_vector_ura(M,N,az,el,d);

a_s  = steer(az_s,el_s);
a_i1 = steer(az_i1,el_i1);
a_i2 = steer(az_i2,el_i2);

%% =========================
% NORMALIZE SOURCE POWER
% =========================
s  = s  / std(s);
i1 = i1 / std(i1);
i2 = i2 / std(i2);

%% =========================
% SIGNAL MODEL
% =========================
Xtrain = a_i1*i1 + a_i2*i2 + noise;
Xdata  = a_s*s + a_i1*i1 + a_i2*i2 + noise;

%% =========================
% COVARIANCE MATRIX
% =========================
R = (Xtrain * Xtrain') / size(Xtrain,2);

disp('Condition before loading:');
disp(cond(R));

%% diagonal loading
delta = 0.1 * trace(R)/size(R,1);
R = R + delta*eye(MN);

disp('Condition after loading:');
disp(cond(R));

%% =========================
% QR DECOMPOSITION (GIVENS)
% =========================
[Q,Rq] = givens_qr(R);

%% =========================
% INVERSE
% =========================
R_inv = invert_with_qr(Q,Rq);

%% =========================
% MVDR WEIGHTS
% =========================
num = R_inv * a_s;
den = a_s' * num;

w_manual = num / (den + 1e-12);

%% =========================
% BEAM PATTERN (AZIMUTH CUT)
% =========================
az_scan = -90:1:90;
el_fixed = 0;

P = zeros(1,length(az_scan));

for k = 1:length(az_scan)

    a = steer(az_scan(k), el_fixed);
    P(k) = abs(w_manual' * a);

end

P = P / max(P);
P_dB = 20*log10(P + 1e-12);

%% =========================
% PLOT
% =========================
figure;
plot(az_scan, P_dB, 'b', 'LineWidth', 2);
grid on;

xlabel('Azimuth (deg)');
ylabel('Gain (dB)');
title('MVDR Beam Pattern (Azimuth Cut, El = 0°)');
ylim([-40 0]);

hold on;
xline(az_s,'g','Target');
xline(az_i1,'r','Interf 1');
xline(az_i2,'r','Interf 2');

%% =========================
% ELEVATION CUT
% =========================
el_scan = -90:1:90;
az_fixed = az_s;

P_el = zeros(1,length(el_scan));

for k = 1:length(el_scan)

    a = steer(az_fixed, el_scan(k));
    P_el(k) = abs(w_manual' * a);

end

P_el = P_el / max(P_el);
P_el_dB = 20*log10(P_el + 1e-12);

figure;
plot(el_scan, P_el_dB, 'm', 'LineWidth', 2);
grid on;

xlabel('Elevation (deg)');
ylabel('Gain (dB)');
title('MVDR Beam Pattern (Elevation Cut, Az fixed)');
ylim([-40 0]);

function a = steering_vector_ura(M,N,az,el,d)

k = 2*pi;

[x,y] = meshgrid(0:N-1,0:M-1);

% center array (IMPORTANT FIX)
x = x - mean(x(:));
y = y - mean(y(:));

x = x(:);
y = y(:);

az = deg2rad(az);
el = deg2rad(el);

ux = cos(el)*cos(az);
uy = cos(el)*sin(az);

phase = k*d*(x*ux + y*uy);

a = exp(1j*phase);

a = a / sqrt(length(a));
end


function [Q,R] = givens_qr(A)

[m,n] = size(A);
R = A;
Q = eye(m);

for j = 1:n
    for i = m:-1:(j+1)

        if R(i,j) ~= 0

            x = [R(i-1,j); R(i,j)];
            G = planerot(x);

            R([i-1 i], j:end) = G * R([i-1 i], j:end);
            Q(:,[i-1 i]) = Q(:,[i-1 i]) * G';

        end
    end
end
end


function A_inv = invert_with_qr(Q,R)

n = size(R,1);
R_inv = zeros(n);

for k = n:-1:1
    R_inv(k,k) = 1/R(k,k);

    for i = k-1:-1:1
        R_inv(i,k) = -R(i,i+1:k)*R_inv(i+1:k,k)/R(i,i);
    end
end

A_inv = R_inv * Q';
end

