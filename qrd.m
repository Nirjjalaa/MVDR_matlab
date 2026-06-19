clc;
clear;
close all;

N = 64;
d = 0.5;

theta_target = 30;
theta_interf = 60;

num_snap = 2000;

n = (0:N-1)';

a_target = steering_vec(theta_target,n,d);
a_interf = steering_vec(theta_interf,n,d);

target_power = 1;
interf_power = 10;
noise_power = 0.1;

s_target = sqrt(target_power/2)*(randn(1,num_snap)+1j*randn(1,num_snap));
s_interf = sqrt(interf_power/2)*(randn(1,num_snap)+1j*randn(1,num_snap));
noise = sqrt(noise_power/2)*(randn(N,num_snap)+1j*randn(N,num_snap));

X = a_target*s_target + a_interf*s_interf + noise;

R = (X*X')/num_snap;

delta = 0.01*trace(R)/N;
R = R + delta*eye(N);

[Q,Ru] = givens_qr(R);

y = Q' * a_target;

w = back_substitution(Ru,y);

w = w/(a_target'*w);

disp("Gain target:");
disp(abs(w'*a_target));

disp("Gain interferer:");
disp(abs(w'*a_interf));

theta_scan = -90:0.5:90;
beam = zeros(size(theta_scan));

for k = 1:length(theta_scan)
    a = steering_vec(theta_scan(k),n,d);
    beam(k) = abs(w'*a);
end

beam = beam/max(beam);
beam_dB = 20*log10(beam + 1e-12);

figure;
plot(theta_scan,beam_dB,'b','LineWidth',2);
grid on;
xlabel('Angle (deg)');
ylabel('Gain (dB)');
title('MVDR using Givens QRD');
ylim([-60 0]);

hold on;
xline(theta_target,'g--','Target');
xline(theta_interf,'r--','Interference');
function a = steering_vec(theta,n,d)

theta = deg2rad(theta);
psi = 2*pi*d*sin(theta);

a = exp(1j*psi*n);
a = a/sqrt(length(n));

end
function [Q,R] = givens_qr(A)

[m,n] = size(A);
Q = eye(m);
R = A;

for j = 1:n
    for i = m:-1:(j+1)

        if R(i,j) ~= 0

            a = R(i-1,j);
            b = R(i,j);

            r = sqrt(a^2 + b^2);

            c = a/r;
            s = -b/r;

            G = [c -s; s c];

            R([i-1 i], j:n) = G * R([i-1 i], j:n);
            Q(:,[i-1 i]) = Q(:,[i-1 i]) * G';

        end
    end
end

end
function x = back_substitution(R,b)

n = length(b);
x = zeros(n,1);

for i = n:-1:1
    x(i) = b(i);
    for j = i+1:n
        x(i) = x(i) - R(i,j)*x(j);
    end
    x(i) = x(i)/R(i,i);
end

end
