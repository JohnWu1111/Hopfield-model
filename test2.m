clear;
tic;
dt = 0.5;
T = 500;
t = 0:dt:T;
nt = length(t);
w = 0:2*pi/((nt-1)*dt):2*pi/dt;

load('result_N14_p2_No1.mat');
fm1 = abs(fft(m1));
fm1a = abs(fft(m1a));
fm2 = abs(fft(m2));
fm2a = abs(fft(m2a));

figure;
plot(w,fm1,w,fm1a,w,fm2,w,fm2a);
legend('m1','m1a','m2','m2a');

toc;