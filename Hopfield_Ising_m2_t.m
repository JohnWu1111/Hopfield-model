clear;
% close all;
clc;
format long
tic;

% Definition of parameters
n = 20; %size
tmax = 1e8; % eventnumber

Tmax = 0.8; %temperature
Tmin = 0.05;
Tstep = 0.05;
T = Tmin:Tstep:Tmax;
nT = length(T);

skip = 1e7;

spin = round(rand(n))*2-1;

m2_std = zeros(1,round((Tmax-Tmin)/Tstep+1));

% parameters of memeory
N = 5;
mem_con = cell(2,1);
for i = 1:N
    mem_con{i} = round(rand(n))*2-1;
end

spintotal = zeros(N,tmax);
% Temperature loop
for j = 1:nT
    T(j)
    Ett = 0;
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
        if rand < exp(-dEt/T(j))
            spin(r(1),r(2)) = -spin(r(1),r(2));
        end
    end
    
%     spintotal = zeros(1,tmax);
    Etotal = zeros(1,tmax);
    
    for  i = 1:tmax %模拟次数循环
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
        
        if rand < exp(-dEt/T(j)) %判断是否作出状态改变
            spin(r(1),r(2)) = -spin(r(1),r(2));
            Ett = Ett + dEt;
        end
        
        for k = 1:N
            temp = mem_con{k};
            spintotal(k,i) = mean(spin.*temp,'all')^2;
        end
        
%         Etotal(i) = Ett;
%         spintotal(i) = mean(mean(spin));
    end
    m2 = sqrt(mean(spintotal,2));
    m2_std(j) = std(m2);
    
%     Cv(j) = var(Etotal)/T(j)^2;  
    
end

figure;
plot(T,m2_std);
xlabel('T')
ylabel('m2-std')

toc;