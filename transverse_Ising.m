clear;
% close all;
clc;
format long
tic;

% Definition of parameters
N = 4; %size
dt = 1;
T = 100;
t = 0:dt:T;
nt = length(t);
J = 1;
len = 2^N;

h = 0.1;

% sigma_x = sparse([0,1;1,0])./2;
% sigma_z = sparse([1,0;0,-1])./2;
sigma_x = sparse([0,1;1,0]);
sigma_z = [1;-1];
I2 = sparse(eye(2));

Hz = zeros(len,1);

for i = 1:N-1
    for j = i+1:N
        H1 = 1;
        H1 = kron(H1,ones(2^(i-1),1));
        H1 = kron(H1,sigma_z/2);
        H1 = kron(H1,ones(2^(j-i-1),1));
        H1 = kron(H1,sigma_z/2);
        H1 = kron(H1,ones(2^(N-j),1));
        Hz = Hz - H1*J/N;
    end
end

Hx = sparse(len,len);
for i = 1:N
    H2 = sparse(1);
    H2 = kron(H2,sparse(eye(2^(i-1))));
    H2 = kron(H2,sigma_x/2);
    H2 = kron(H2,sparse(eye(2^(N-i))));
    Hx = Hx + h*H2;
end

H = diag(Hz) + full(Hx);
[V,D] = eig(H);
e = diag(D);
NN = length(e);

phi0 = zeros(NN,1);
phi0(1) = 1;

M_Sz = zeros(NN,N);
for i = 1:N
    temp = ones(2^(i-1),1);
    temp = kron(temp,[1;-1]/2);
    temp = kron(temp,ones(2^(N-i),1));
    M_Sz(:,i) = temp;
end

temp = V'*phi0;
trans = exp(-1i*e*t);
temp = trans.*temp;
phit = V*temp;

Sz = zeros(N,nt);
for i = 1:N
    Sz(i,:) = sum(conj(phit).*(M_Sz(:,i).*phit));
end

Sz_mean = sum(Sz)/N;

figure;
plot(t,Sz_mean)

toc;

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