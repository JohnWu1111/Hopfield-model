clear;
% close all;
clc;
format long
tic;

% Definition of parameters
n = 20; %size
tmax = 3e7; % eventnumber
tstep = 5e5;
nt = round(tmax/tstep);
t=1:nt;

T = 0.1;

skip = 1e7;

% Cv = zeros(1,round((Tmax-Tmin)/Tstep+1));

% parameters of memeory
N = 3;
mem_con = cell(2,1);
for i = 1:N
    mem_con{i} = round(rand(n))*2-1;
end

spin = mem_con{1};

m2t = zeros(N,nt);
Et = zeros(1,nt);
figure;

subplot(2,2,1)
gM1 = plot(t,m2t(1,:));

subplot(2,2,2)
gM2 = plot(t,m2t(2,:));

subplot(2,2,3)
gM3 = plot(t,m2t(3,:));

subplot(2,2,4)
gEt = plot(t,Et);

Ett = 0;
ii = 1;

% Activation
for s = 1:skip
    r = ceil(rand(1,2)*n);
    su = spin(mod(r(1)-2,n)+1,r(2));
    sd = spin(mod(r(1),n)+1,r(2));
    sl = spin(r(1),mod(r(2)-2,n)+1);
    sr = spin(r(1),mod(r(2),n)+1);
    it = spin(r(1),r(2));
    dEt = 0;
    for k = 1:N
        temp = mem_con{k};
        mu = temp(mod(r(1)-2,n)+1,r(2));
        md = temp(mod(r(1),n)+1,r(2));
        ml = temp(r(1),mod(r(2)-2,n)+1);
        mr = temp(r(1),mod(r(2),n)+1);
        mm = temp(r(1),r(2));
        dEt = dEt + 2*it*su*mm*mu;
        dEt = dEt + 2*it*sd*mm*md;
        dEt = dEt + 2*it*sl*mm*ml;
        dEt = dEt + 2*it*sr*mm*mr;
    end
    dEt = dEt/N;
    %     dEt = spin(r(1),r(2))*(su+sd+sl+sr)*2;
    if rand < exp(-dEt/T)
        spin(r(1),r(2)) = -spin(r(1),r(2));
        Ett = Ett + dEt;
    end
end

spintotal = zeros(N,tmax);
Etotal = zeros(1,tmax);

for i = 1:tmax %模拟次数循环
    r = ceil(rand(1,2)*n);
    su = spin(mod(r(1)-2,n)+1,r(2));
    sd = spin(mod(r(1),n)+1,r(2));
    sl = spin(r(1),mod(r(2)-2,n)+1);
    sr = spin(r(1),mod(r(2),n)+1);
    %         dEt = spin(r(1),r(2))*(su+sd+sl+sr)*2;
    it = spin(r(1),r(2));
    dEt = 0;
    for k = 1:N
        temp = mem_con{k};
        mu = temp(mod(r(1)-2,n)+1,r(2));
        md = temp(mod(r(1),n)+1,r(2));
        ml = temp(r(1),mod(r(2)-2,n)+1);
        mr = temp(r(1),mod(r(2),n)+1);
        mm = temp(r(1),r(2));
        dEt = dEt + 2*it*su*mm*mu;
        dEt = dEt + 2*it*sd*mm*md;
        dEt = dEt + 2*it*sl*mm*ml;
        dEt = dEt + 2*it*sr*mm*mr;
    end
    dEt = dEt/N;
    
    if rand < exp(-dEt/T) %判断是否作出状态改变
        spin(r(1),r(2)) = -spin(r(1),r(2));
        Ett = Ett + dEt;
    end
    
    Etotal(i) = Ett;
    for j = 1:N
        temp = mem_con{j};
        spintotal(j,i) = mean(spin.*temp,'all')^2;
    end
    
    if mod(i-1,tstep)==0
        m2t(:,ii) = sqrt(spintotal(:,i));
        Et(ii) = Ett;
        ii = ii +1;
        set(gM1,'Ydata',double(m2t(1,:)));
        set(gM2,'Ydata',double(m2t(2,:)));
        set(gM3,'Ydata',double(m2t(3,:)));
        set(gEt,'Ydata',double(Et));
        drawnow
    end
end

m2 = sqrt(mean(spintotal,2))

% for k = 1:N
%     figure;
%     temp = mem_con{k};
%     pcolor(double(temp));
%     colormap(gray(2));
% end
% figure;
% pcolor(double(spin));
% colormap(gray(2));

Cv = var(Etotal)/T


toc;