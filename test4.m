N = 26;
p = 4;
No = 1;
NN = length(e);
h = '0_1';
figure;
plot(1:NN, mm(1,:)');
xlabel('\epsilon')
ylabel('order parameter')
print(strcat('m2_N', num2str(N), '_p', num2str(p), '_h', num2str(h), '_No', num2str(No)), '-dpng', '-r0');

figure;
plot(1:NN, corre');
xlabel('\epsilon')
ylabel('correlation between largest spins')
print(strcat('corre_N', num2str(N), '_p', num2str(p), '_h', num2str(h), '_No', num2str(No)), '-dpng', '-r0');

figure;
plot(1:NN, corre_in');
xlabel('\epsilon')
ylabel('internal correlation')
print(strcat('corre_in_N', num2str(N), '_p', num2str(p), '_h', num2str(h), '_No', num2str(No)), '-dpng', '-r0');


figure;
plot(1:NN, corre_ext');
xlabel('\epsilon')
ylabel('external correlation')
print(strcat('corre_ext_N', num2str(N), '_p', num2str(p), '_h', num2str(h), '_No', num2str(No)), '-dpng', '-r0');

figure;
plot(1:NN, phi4');
xlabel('\epsilon')
ylabel('\phi^4')
print(strcat('phi4_N', num2str(N), '_p', num2str(p), '_h', num2str(h), '_No', num2str(No)), '-dpng', '-r0');

figure;
histogram(q,200,'Normalization','pdf','DisplayStyle','stairs');
print(strcat('En_diff_distri_N', num2str(N), '_p', num2str(p), '_h', num2str(h), '_No', num2str(No)), '-dpng', '-r0');

figure;
plot(1:NN, sus');
xlabel('\epsilon')
ylabel('\chi')
print(strcat('SG_order_N', num2str(N), '_p', num2str(p), '_h', num2str(h), '_No', num2str(No)), '-dpng', '-r0');