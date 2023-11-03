clear;
% close all;
clc;
format long
tic;

% Definition of parameters
N = 20;%size
p = 4;
num = 1;

ol = zeros(p,p);
len_p = 2^(p-1);

N_div = zeros(p,len_p);
for i = p:-1:1
    a = 2^(i-1);
    for j = 1:len_p
        temp1 = floor((j-1)/a);
        temp2 = mod(temp1,2);
        N_div(p-i+1,j) = -2*temp2+1;
    end
end

% trans = zeros(len_p,len_p);
% count = 1;
% for i = 1:p
%     for j = i+1:p
%         trans(count,:) = N_div(i,:).*N_div(j,:);
%         count = count + 1;
%     end
% end
% trans(end,:) = ones(1,len_p);

N_s = zeros(len_p,1);

% parameters of memeory
mem_con = zeros(p,N);
mem_con(1,:) = ones(1,N);
mem_con(2:p,:) = round(rand(p-1,N))*2-1;

for i = 1:p
    for j = i+1:p
        ol(i,j) = sum(mem_con(i,:).*mem_con(j,:));
    end
end

for i = 1:N
    for j = 1:len_p
        temp = mean(mem_con(:,i).*N_div(:,j));
        if temp == 1
            N_s(j) = N_s(j) + 1;
            break
        end
    end
end
% NN = sum(N_s);
NN = 1;
for i = 1:len_p
    NN = NN*(N_s(i)+1);
end

Jij = zeros(len_p,len_p);
for i = 1:len_p
    for j = i+1:len_p
        Jij(i,j) = 2*sum(N_div(:,i).*N_div(:,j));
    end
end
Jij = Jij + eye(len_p)*p;


fname = ['mem_N',num2str(N),'_p',num2str(p),'_No2.mat'];
save(fname,'mem_con','ol','N_s','N_div','Jij','len_p','-v7.3');

toc;