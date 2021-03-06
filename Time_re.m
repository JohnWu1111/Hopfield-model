clear;
% close all;
clc;
format long
tic;

% Definition of parameters
N = 12; %size
No = 5;
dt = 1;
T = 1000;
t = 0:dt:T;
nt = length(t);

p = 2;

h = 0.3;


load(strcat('mem_N',num2str(N),'_p',num2str(p),'_No',num2str(No),'.mat'));
load(strcat('eig_N',num2str(N),'_p',num2str(p),'_No',num2str(No),'.mat'));
overlap = mean(mem_con{1}.*mem_con{2})

toc;
tic;

sigma_x = sparse([0,1;1,0]);
sigma_z = sparse([1,0;0,-1]);
I2 = sparse(eye(2));

temp1 = mem_con{1};
temp1a = round((-temp1+1)/2)+1;
temp1 = round((temp1+1)/2)+1;
spin1 = sigma_z(:,temp1(1));
spin1a = sigma_z(:,temp1a(1));

temp2 = mem_con{2};
temp2a = round((-temp2+1)/2)+1;
temp2 = round((temp2+1)/2)+1;
spin2 = sigma_z(:,temp2(1));
spin2a = sigma_z(:,temp2a(1));
for i = 2:N
    spin1 = kron(spin1,sigma_z(:,temp1(i)));
    spin1a = kron(spin1a,sigma_z(:,temp1a(i)));
    spin2 = kron(spin2,sigma_z(:,temp2(i)));
    spin2a = kron(spin2a,sigma_z(:,temp2a(i)));
end

% m1 = zeros(1,nt);
% m1a = zeros(1,nt);
% m1(1) = 1;
% m1a(1) = spin1a'*spin1;
% m2 = zeros(1,nt);
% m2a = zeros(1,nt);
% m2(1) = spin2'*spin1;
% m2a(1) = spin2a'*spin1;

% spin0 = V\spin1;
spin0 = V'*spin1;
% spin = spin0;
% trans = exp(-1i*e*dt);

toc;
tic;

trans = exp(-1i*e*t);
spin = spin0.*trans;
spint = V*spin;

m1 = abs(spin1'*spint);
m1a = abs(spin1a'*spint);
m2 = abs(spin2'*spint);
m2a = abs(spin2a'*spint);

figure;
plot(t,m1,t,m1a,t,m2,t,m2a);
legend('m1','m1a','m2','m2a');
% set(gcf,'position',[250 300 1200 500]);
xlabel('t');
ylabel('m');

% saveas(gcf,strcat('m_N',num2str(N),'_p',num2str(p),'_No',num2str(No),'.fig'));
% saveas(gcf,strcat('m_N',num2str(N),'_p',num2str(p),'_No',num2str(No),'.png'));
% saveas(gcf,strcat('m_N',num2str(N),'_p',num2str(p),'_No',num2str(No),'.pdf'));
% 
% save(strcat('result_N',num2str(N),'_p',num2str(p),'_No',num2str(No),'.mat'),'m1','m1a','m2','m2a','-v7.3');
% 
toc;