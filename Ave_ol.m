clear;
% close all;
clc;
format long
tic;

% global dt H

% Definition of parameters
N = 2; %size
p = 2;
NN = 2^N;
sum = 0;

% parameters of memeory

mem_con2 = ones(1,N);
sum = sum + 1;

for i = 1:NN-1
    mem_con2 = add(mem_con2,1);
    sum = sum + abs(mean(mem_con2));
end
sum = sum/NN

toc;

% load('mem_N12_p2_No5.mat')

function y=add(x,pos)
    if x(pos) == 1
        x(pos) = -1;
        y = x;
    else
        x(pos) = 1;
        y = add(x,pos+1);
    end    
end