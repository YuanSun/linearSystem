% Assignment 5
% Question 2
% Yuan Sun
% 2015-11-25

clear
clc

%% Contruct the system

A = [0 0 0; 1 0 -16; 0 1 -10];
B = [10; 0; 0];
C = [0 0 1];

sys = ss(A, B, C, 0);

%% Desired Characteristic polynomial of A-BK

poles = [-1.42 -1.04+2.14*1i -1.04-2.14*1i];

K = place(A, B, poles);