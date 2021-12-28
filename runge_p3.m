clear;
% close all;
clc;
format long
tic;

% Definition of parameters
N = 80; %size
p = 3;
dt = 0.005;
T = 500;
t = 0:dt:T;
nt = length(t);
J = 1;
h = 0.3;
No = 1;
cut = 1;

load(strcat('p3_pattern\p3_nondeg_N',num2str(N),'.mat'));
N1 = ol(No,1);
N2 = ol(No,2);
N3 = ol(No,3);
N4 = ol(No,4);

NN = (N1+1)*(N2+1)*(N3+1)*(N4+1);

S1 = N1/2;
S2 = N2/2;
S3 = N3/2;
S4 = N4/2;

S1_z = zeros(N1+1,1);
S1_p = zeros(N1+1);
S1_m = zeros(N1+1);
S2_z = zeros(N2+1,1);
S2_p = zeros(N2+1);
S2_m = zeros(N2+1);
S3_z = zeros(N3+1,1);
S3_p = zeros(N3+1);
S3_m = zeros(N3+1);
S4_z = zeros(N4+1,1);
S4_p = zeros(N4+1);
S4_m = zeros(N4+1);

% construction of matrice
for m = 1:N1+1
    S1_z(m) = N1/2 - (m-1);
end
for m = 1:N2+1
    S2_z(m) = N2/2 - (m-1);
end
for m = 1:N3+1
    S3_z(m) = N3/2 - (m-1);
end
for m = 1:N4+1
    S4_z(m) = N4/2 - (m-1);
end

for m = 1:N1
    S1_p(m,m+1) = sqrt(S1*(S1+1)-S1_z(m+1)*(S1_z(m+1)+1));
    S1_m(m+1,m) = sqrt(S1*(S1+1)-S1_z(m)*(S1_z(m)-1));
end
for m = 1:N2
    S2_p(m,m+1) = sqrt(S2*(S2+1)-S2_z(m+1)*(S2_z(m+1)+1));
    S2_m(m+1,m) = sqrt(S2*(S2+1)-S2_z(m)*(S2_z(m)-1));
end
for m = 1:N3
    S3_p(m,m+1) = sqrt(S3*(S3+1)-S3_z(m+1)*(S3_z(m+1)+1));
    S3_m(m+1,m) = sqrt(S3*(S3+1)-S3_z(m)*(S3_z(m)-1));
end
for m = 1:N4
    S4_p(m,m+1) = sqrt(S4*(S4+1)-S4_z(m+1)*(S4_z(m+1)+1));
    S4_m(m+1,m) = sqrt(S4*(S4+1)-S4_z(m)*(S4_z(m)-1));
end

S1_x = (S1_p + S1_m)/2;
S2_x = (S2_p + S2_m)/2;
S3_x = (S3_p + S3_m)/2;
S4_x = (S4_p + S4_m)/2;

S1_z = sparse(S1_z);
S2_z = sparse(S2_z);
S3_z = sparse(S3_z);
S4_z = sparse(S4_z);
S1_x = sparse(S1_x);
S2_x = sparse(S2_x);
S3_x = sparse(S3_x);
S4_x = sparse(S4_x);

% construction of Hamiltonian
H1 = -J*(3*kron(S1_z.^2,ones((N2+1)*(N3+1)*(N4+1),1))...
        +3*kron3(ones(N1+1,1),S2_z.^2,ones((N3+1)*(N4+1),1))...
        +3*kron3(ones((N1+1)*(N2+1),1),S3_z.^2,ones(N4+1,1))...
        +3*kron(ones((N1+1)*(N2+1)*(N3+1),1),S4_z.^2)...%inner term
        +2*kron3(S1_z,kron(S2_z,ones(N3+1,1))+kron(ones(N3+1,1),S2_z),ones(N4+1,1))...
        +2*kron3(ones(N1+1,1),kron(S2_z,ones(N3+1,1))+kron(ones(N3+1,1),S2_z),S4_z)...
        -2*kron4(ones(N1+1,1),S2_z,S3_z,ones(N4+1,1))...
        -2*kron3(S1_z,ones((N2+1)*(N3+1),1),S4_z))/(2*N);
H1 = diag(H1);
H2 = h*(kron(S1_x,eye((N2+1)*(N3+1)*(N4+1)))...
    +kron3(eye(N1+1),S2_x,eye((N3+1)*(N4+1)))...
    +kron3(eye((N1+1)*(N2+1)),S3_x,eye(N4+1))...
    +kron(eye((N1+1)*(N2+1)*(N3+1)),S4_x));
H = H1 + H2;

spin1 = zeros(NN,1);
spin1a = zeros(NN,1);
spin1(1) = 1;
spin1a(end) = 1;

temp1 = zeros((N1+1)*(N2+1),1);
temp2 = zeros((N3+1)*(N4+1),1);
temp1(1) = 1;
temp2(end) = 1;
spin2 = kron(temp1,temp2);
temp1 = zeros((N1+1)*(N2+1),1);
temp2 = zeros((N3+1)*(N4+1),1);
temp1(end) = 1;
temp2(1) = 1;
spin2a = kron(temp1,temp2);

temp1 = zeros(N1+1,1);
temp2 = zeros(N2+1,1);
temp3 = zeros(N3+1,1);
temp4 = zeros(N4+1,1);
temp1(1) = 1;
temp2(end) = 1;
temp3(1) = 1;
temp4(end) = 1;
spin3 = kron4(temp1,temp2,temp3,temp4);
temp1 = zeros(N1+1,1);
temp2 = zeros(N2+1,1);
temp3 = zeros(N3+1,1);
temp4 = zeros(N4+1,1);
temp1(end) = 1;
temp2(1) = 1;
temp3(end) = 1;
temp4(1) = 1;
spin3a = kron4(temp1,temp2,temp3,temp4);

% time revolution

pro1 = zeros(1,nt);
pro1(1) = 1;
pro2 = zeros(1,nt);
pro2(1) = spin2'*spin1;
pro3 = zeros(1,nt);
pro3(1) = spin3'*spin1;
pro1a = zeros(1,nt);
pro1a(1) = 0;
pro2a = zeros(1,nt);
pro2a(1) = spin2a'*spin1;
pro3a = zeros(1,nt);
pro3a(1) = spin3a'*spin1;

M1 = kron_p4(S1_z,S2_z,S3_z,S4_z)/N;
M2 = kron_p4(S1_z,S2_z,-S3_z,-S4_z)/N;
M3 = kron_p4(S1_z,-S2_z,S3_z,-S4_z)/N;
m1 = zeros(1,nt);
m2 = zeros(1,nt);
m3 = zeros(1,nt);
m1(1) = 1/2;
m2(1) = M2'*spin1;
m3(1) = M3'*spin1;

spin = spin1;
H = sparse(H);
% H2 = H^2;
% H3 = H2*H;
% H4 = H3*H;

for i = 2:nt
%     spin = runge(spin);
    c1 = -1i*H*spin;
    c2 = -1i*H*(spin+c1.*dt/2);
    c3 = -1i*H*(spin+c2.*dt/2);
    c4 = -1i*H*(spin+c3.*dt);
    spin = spin + dt*(c1+2*c2+2*c3+c4)/6;

%     trans = -1i*dt*H-dt^2*H2/2+1i*dt^3*H3/6+dt^4*H4/24;
%     spin = trans*spin + spin;
    
    pro1(i) = spin'*spin1;
    pro2(i) = spin'*spin2;
    pro3(i) = spin'*spin3;
    
    spin_sq = abs(spin.^2);
    m1(i) = M1'*spin_sq;
    m2(i) = M2'*spin_sq;
    m3(i) = M3'*spin_sq;
end

pro1 = abs(pro1);
pro1a = abs(pro1a);
pro2 = abs(pro2);
pro2a = abs(pro2a);
pro3 = abs(pro3);
pro3a = abs(pro3a);

% Fourier analysis
W = 2 * pi / dt;
dw = 2 * pi / ((nt - 1) * dt);
w = 0:dw:2 * pi / dt;

fpro1 = abs(fft(pro1));
fpro1a = abs(fft(pro1a));
fpro2 = abs(fft(pro2));
fpro2a = abs(fft(pro2a));
fpro3 = abs(fft(pro3));
fpro3a = abs(fft(pro3a));

% delete the peak at zero
fpro1(1) = 0;
fpro1a(1) = 0;
fpro2(1) = 0;
fpro2a(1) = 0;
fpro3(1) = 0;
fpro3a(1) = 0;

mean_m1 = mean(m1);
mean_m2 = mean(m2);
mean_m3 = mean(m3);
std_m1 = std(m1);
std_m2 = std(m2);
std_m3 = std(m3);

fm1 = abs(fft(m1));
fm2 = abs(fft(m2));
fm3 = abs(fft(m3));
% delete the peak at zero
fm1(1) = 0;
fm2(1) = 0;
fm3(1) = 0;

%{
figure;
set(gcf, 'position', [250 70 1400 900]);
subplot(3,2,1)
plot(t,pro1,t,pro2,t,pro3);
xlabel('t');
ylabel('pro');
axis([0 T 0 1])
legend('pro1', 'pro2', 'pro3');

subplot(3,2,2)
plot(t,pro1a,t,pro2a,t,pro3a);
xlabel('t');
ylabel('pro');
axis([0 T 0 1])
legend('pro1a', 'pro2a', 'pro3a');

subplot(3,2,3)
plot(w(1:cut),fpro1(1:cut),w(1:cut),fpro2(1:cut),w(1:cut),fpro3(1:cut));
xlabel('w');
ylabel('fpro');
legend('fpro1', 'fpro2', 'fpro3');
% axis([0 1 0 inf])

subplot(3,2,4)
plot(w(1:cut),fpro1a(1:cut),w(1:cut),fpro2a(1:cut),w(1:cut),fpro3a(1:cut));
xlabel('w');
ylabel('fpro');
legend('fpro1a', 'fpro2a', 'fpro3a');
% axis([0 1 0 inf])

subplot(3,2,5)
plot(w,fpro1,w,fpro2,w,fpro3);
xlabel('w');
ylabel('fpro');
legend('fpro1', 'fpro2', 'fpro3');

subplot(3,2,6)
plot(t,m1,t,m2,t,m3)
xlabel('t');
ylabel('pro');
axis([0 T -1/2 1/2])
%}


figure;
set(gcf, 'position', [250 70 1400 900]);
subplot(2,2,1)
plot(t,pro1,t,pro2,t,pro3);
xlabel('t');
ylabel('pro');
axis([0 T 0 1])
legend('pro1', 'pro2', 'pro3');

subplot(2,2,2)
plot(t,m1,t,m2,t,m3);
xlabel('t');
ylabel('m');
axis([0 T -1/2 1/2])
legend('m1', 'm2', 'm3');

subplot(2,2,3)
plot(w(1:cut),fpro1(1:cut),w(1:cut),fpro2(1:cut),w(1:cut),fpro3(1:cut));
xlabel('w');
ylabel('fpro');
legend('fpro1', 'fpro2', 'fpro3');
% axis([0 1 0 inf])

subplot(2,2,4)
plot(w(1:cut),fm1(1:cut),w(1:cut),fm2(1:cut),w(1:cut),fm3(1:cut));
xlabel('w');
ylabel('fpro');
legend('fm1a', 'fm2a', 'fm3a');

toc;

function y = kron4(a,b,c,d)
y = kron(kron(kron(a,b),c),d);
end

function y = kron3(a,b,c)
y = (kron(kron(a,b),c));
end

function y = kron_p(a,b)
la = length(a);
lb = length(b);
y = zeros(la*lb,1);
for i = 1:la
    for j = 1:lb
        y((i-1)*lb+j) = a(i) + b(j);
    end
end
end

function y = kron_p4(a,b,c,d)
y = kron_p(kron_p(kron_p(a,b),c),d);
end