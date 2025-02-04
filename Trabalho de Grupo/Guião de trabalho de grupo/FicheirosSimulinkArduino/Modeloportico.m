clear all; close all; clc
n = 3;
Ts = 0.002;
Kv = 4.0816e3;
m = 0.93;
K = Kv*[2 -1 0;-1 2 -1;0 -1 1];
M = m*eye(n);
vP = [-0.3401 -0.6128 -0.7642;0.7642 0.3401 -0.6128;-0.6128 0.7642 -0.3401]';

freq = [4.6924; 13.1479; 18.9992];

Aol = [zeros(n) eye(n);-inv(M)*K zeros(n)];
Bol = [0 1/m 0 0 0 0]';
Col = [0 0 1 0 0 0];
Dol = [0];
%%
Co = ctrb(Aol,Bol);
Ob = obsv(Aol,Col);

rank (Co)
rank(Ob)
eePortico = ss(Aol, Bol, Col, Dol);
esPortico_disc = c2d(eePortico,Ts);
eig(Aol) %% tem um valor próprio positivo (sistema instavel
Kr=  [16.624 110.3 -55.716 0.51505 -0.43561 0.10403]

%%
filName = CreateCodeSSArduino('matrizA.txt',esPortico_disc,Kr,6,1)