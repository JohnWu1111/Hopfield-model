clear;
% close all;
clc;
format long
tic;

% Definition of parameters
N = 16; %size
dt = 5;
T = 50000;
t = 0:dt:T;
nt = length(t);
J = 1*2^2;
h = 0.3*2;
ol = 0;
N1 = N/2 + ol;
N2 = N/2 - ol;

NN = (N1+1)*(N2+1);

S1 = N1/2;
S2 = N2/2;

S1_z = zeros(N1+1);
S1_p = zeros(N1+1);
S1_m = zeros(N1+1);
S2_z = zeros(N2+1);
S2_p = zeros(N2+1);
S2_m = zeros(N2+1);

for i = 1:N1+1
    S1_z(i,i) = N1/2 - (i-1);
end

for i = 1:N2+1
    S2_z(i,i) = N2/2 - (i-1);
end

M1 = diag(S1_z);
M2 = diag(S2_z);

for i = 1:N1
    S1_p(i,i+1) = sqrt(S1*(S1+1)-M1(i+1)*(M1(i+1)+1));
    S1_m(i+1,i) = sqrt(S1*(S1+1)-M1(i)*(M1(i)-1));
end

for i = 1:N2
    S2_p(i,i+1) = sqrt(S2*(S2+1)-M2(i+1)*(M2(i+1)+1));
    S2_m(i+1,i) = sqrt(S2*(S2+1)-M2(i)*(M2(i)-1));
end

S1_x = (S1_p + S1_m)/2;
S2_x = (S2_p + S2_m)/2;
H = -J*(kron(S1_z.^2,eye(N2+1))+kron(eye(N1+1),S2_z.^2))/(N) + h*(kron(S1_x,eye(N2+1))+kron(eye(N1+1),S2_x));

[V,D] = eig(H);
e = diag(D);

spin1 = zeros(NN,1);
spin1a = zeros(NN,1);
spin1(1) = 1;
spin1a(end) = 1;


spin0 = V'*spin1;
trans = exp(-1i*e*t);
spin = spin0.*trans;
spint = V*spin;

m1 = abs(spin1'*spint);
m1a = abs(spin1a'*spint);

dw = 2 * pi / ((nt - 1) * dt);
w = 0:dw:2 * pi / dt;

fm1 = abs(fft(m1));
fm1a = abs(fft(m1a));

figure;
plot(t,m1,t,m1a);
xlabel('t');
ylabel('m');
axis([0 T 0 1])

figure;
plot(w, fm1, w, fm1a);
legend('fm1', 'fm1a');
xlabel('\omega');
ylabel('fft of m');

toc;