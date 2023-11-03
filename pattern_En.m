clear;
% close all;
clc;
format long
tic;

% Definition of parameters
N = 24; %size
p = 4;
dt = 0.005;
T = 500;
t = 0:dt:T;
nt = length(t);
J = 1;
h = 0.1;
No = 1;
cut = 100;

load(strcat('mem_N',num2str(N),'_p',num2str(p),'_No',num2str(No),'.mat'));

NN = 1;
for i = 1:len_p
    NN = NN*(N_s(i)+1);
end

S = N_s/2;

% construction of matrice
S_z = cell(len_p,1);
for i = 1:len_p
    S_z{i} = zeros(N_s(i)+1,1);
    for m = 1:N_s(i)+1
        S_z{i}(m) = N_s(i)/2 - (m-1);
    end
end

S_p = cell(len_p,1);
S_m = cell(len_p,1);
S_x = cell(len_p,1);
for i = 1:len_p
    S_p{i} = zeros(N_s(i)+1);
    S_m{i} = zeros(N_s(i)+1);
    for m = 1:N_s(i)
        S_p{i}(m,m+1) = sqrt(S(i)*(S(i)+1)-S_z{i}(m+1)*(S_z{i}(m+1)+1));
        S_m{i}(m+1,m) = sqrt(S(i)*(S(i)+1)-S_z{i}(m)*(S_z{i}(m)-1));
    end
    
    S_p{i} = sparse(S_p{i});
    S_m{i} = sparse(S_m{i});
    
    S_x{i} = (S_p{i} + S_m{i})/2;
end

% construction of Hamiltonian
H1 = zeros(NN,1);
for i = 1:len_p
    temp = 1;
    for k = 1:i-1
        temp = kron(temp,ones(N_s(k)+1,1));
    end
    temp = kron(temp,S_z{i}.^2);
    for k = i+1:len_p
        temp = kron(temp,ones(N_s(k)+1,1));
    end
    H1 = H1 + 4*temp;
    
    for j = i+1:len_p
        if Jij(i,j) == 0
            continue
        end
        temp = ones(1);
        for k = 1:i-1
            temp = kron(temp,ones(N_s(k)+1,1));
        end
        temp = kron(temp,S_z{i});
        for k = i+1:j-1
            temp = kron(temp,ones(N_s(k)+1,1));
        end
        temp = kron(temp,S_z{j});
        for k = j+1:len_p
            temp = kron(temp,ones(N_s(k)+1,1));
        end
        H1 = H1 + Jij(i,j)*temp;
    end
end
H1 = -J*(H1)/(2*N);
H1 = sparse(H1);
H1 = diag(H1);

H2 = sparse(NN,NN);
for i = 1:len_p
    temp = sparse(1);
    for k = 1:i-1
        temp = kron(temp,diag(sparse(ones(N_s(k)+1,1))));
    end
    temp = kron(temp,S_x{i});
    for k = i+1:len_p
        temp = kron(temp,diag(sparse(ones(N_s(k)+1,1))));
    end
    H2 = H2 + h*temp;
end
H = H1 + H2;

spin_s = zeros(NN,p);
spin_sa = zeros(NN,p);
spin_s(1,1) = 1;
spin_sa(end,1) = 1;

for i = 2:p
    temp1 = ones(1);
    temp1a = ones(1);
    for k = 1:len_p
        temp2 = zeros(N_s(k)+1,1);
        temp3 = zeros(N_s(k)+1,1);
        if N_div(i,k) > 0
            temp2(1) = 1;
            temp3(end) = 1;
        else
            temp2(end) = 1;
            temp3(1) = 1;
        end
        temp1 = kron(temp1,temp2);
        temp1a = kron(temp1a,temp3);
    end
    spin_s(:,i) = temp1;
    spin_sa(:,i) = temp1a;
end

pro = zeros(p+1,nt);
m = zeros(p,nt);
M = zeros(NN,p);
for i = 1:p
    temp = zeros(1);
    for k = 1:len_p
        temp = kron_p(temp,S_z{k}*N_div(i,k));
    end
    M(:,i) = temp;
end

% time revolution
spin_x = 2*rand(NN,1)-1;
spin_x = spin_x./sqrt(sum(abs(spin_x).^2));
E = spin_x'*H*spin_x;
% spin = spin_s(:,1);
spin = spin_x;
H = -1i*H;
pro(1:p,1) = spin_s'*spin;
pro(p+1,1) = spin_x'*spin;
m(:,1) = M'*spin;

for i = 2:nt
    c1 = H*spin;
    c2 = H*(spin+c1.*(dt/2));
    c3 = H*(spin+c2.*(dt/2));
    c4 = H*(spin+c3.*dt);
    spin = spin + dt*(c1+2*c2+2*c3+c4)/6;
    
    pro(1:p,i) = abs(spin_s'*spin);
    pro(p+1,i) = abs(spin_x'*spin);
    
    spin_sq = abs(spin.^2);
    m(:,i) = M'*spin_sq;
end
m = m/N;
m_mean = mean(m);

% Fourier analysis
W = 2 * pi / dt;
dw = 2 * pi / ((nt - 1) * dt);
w = 0:dw:2 * pi / dt;

fpro = abs(fft(pro,nt,2));

% delete the peak at zero
fpro(:,1) = 0;

fm = abs(fft(m,nt,2));

% delete the peak at zero
fm(:,1) = 0;

le1 = cell(p,1);
le2 = cell(p,1);
le3 = cell(p,1);
le4 = cell(p,1);
for i = 1:p
    le1{i} = strcat('pro',num2str(i));
    le2{i} = strcat('m',num2str(i));
    le3{i} = strcat('fpro',num2str(i));
    le4{i} = strcat('fm',num2str(i));
end

figure;
set(gcf, 'position', [250 70 1400 900]);

subplot(2,2,1)
plot(t,pro)
xlabel('t')
ylabel('pro')
legend(le1)

subplot(2,2,2)
plot(t,m)
xlabel('t')
ylabel('m')
legend(le2)

subplot(2,2,3)
plot(w(1:cut),fpro(:,1:cut))
xlabel('w')
ylabel('fpro')
legend(le3)

subplot(2,2,4)
plot(w(1:cut),fm(:,1:cut))
xlabel('w')
ylabel('fm')
legend(le4)

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