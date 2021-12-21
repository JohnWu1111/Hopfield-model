% clear;
% close all;
clc;
format long
tic;

N = 40; %size
No = 102;
p = 2;

% load(strcat('H\H_N',num2str(N),'_p',num2str(p),'_No',num2str(No),'.mat'));
% 
% H = full(H);
% [V,D] = eig(H);
% e = diag(D);

len = length(e);
e = sort(e);
s = zeros(len-1,1);
r = zeros(len-2,1);
for i = 1:len-1
    s(i) = e(i+1) - e(i);
end

s = sort(s);

% for i = 1:len-2
%     r(i) = s(i+1)/s(i);
% end
% 
% q = min([r 1./r],[],2);
% 
% q = sort(q);

figure;
histogram(s,1000,'Normalization','pdf','DisplayStyle','stairs');

% saveas(gcf,strcat('Energy_diff_N',num2str(N),'_p',num2str(p),'_No',num2str(No),'.fig'));
% print(strcat('Energy_diff_N',num2str(N),'_p',num2str(p),'_No',num2str(No)),'-dpng','-r0');
% save(strcat('Energy_diff_N',num2str(N),'_p',num2str(p),'_No',num2str(No),'.mat'),'q','s','e','-v7.3');

toc;