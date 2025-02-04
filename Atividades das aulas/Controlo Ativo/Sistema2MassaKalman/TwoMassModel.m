%%
clear all; close all; clc;
M1 = 1; M2 = 2; K1=10; K2 = 20; C1 = 0.1; C2 = 0.2; %% constantes
x10 = 0.618; x20 = 1; % Posições iniciais das massas
Aol = [0 1 0 0; -(K1+K2)/M1  -(C1+C2)/M1 K2/M1 C2/M1; 0 0 0 1; K2/M2 C2/M2 -K2/M2 -C2/M2]; % matriz A
Bol = [0; 0; 0; 1/M2]; % matrixz B - somente força aplicada em m1
Col = eye(4); %% matriz C - Identidade: todos os estados são sídas
Dol  = [0;0;0;0]; % matriz D
x0 = [x10; 0; x20;0]; % vetor de estado inicial
sys = ss(Aol,Bol,Col,Dol); % espaço de estados
initial(sys,x0,5); % simulação só com condições iniciais durante 5 segundos

x10 = 1.618; x20 = -1;
%K1=K2; M1=M2; C1=0; C2=0;


%% Controlo por retorno de estado
P = [-2+4*i; -2-4*i; -1; -4];  %polos desejados para o controlo por retorno de estado
K = place(Aol,Bol,P); % posicionamento de polos com a função place
Acl = Aol-Bol*K; % matriz do sistema em malha fechada
sysControl = ss(Acl,Bol,Col,Dol); % espaço de estados
hold on;   % para sobreposição dos resultados de simulação
initial(sysControl,x0,5); % simulação do sistema em malha fechada só com condições iniciais duurante 5 segundos

%% Controlo por retorno de estado com observador
C1 = [1 0 0 0]; % o systema tem só a saída Y = x1
D1 = 0;
AT = Aol'; CT = C1'; %Matrizes transpostas para cálculo da matris de ganho L com a função "place"
Po = [-6+2*i; -6-2*i; -5; -8];  % posição desejada dos polos do observador, para cálculo matriz L (utilizar a matriz transposta de A de L e de C
LT = place(AT, CT,Po); % matriz L transposta
L = LT'; % matriz de ganho L
Ao = Aol-L*C1; % matriz A do observador
Co=eye(4); % matriz de saída do observador: monitoriza todos os estados estimados
Bo = [Bol L]; % matriz de entrada do observador: o observador tem a entrada u e y
Do = [0 0 0 0; 0 0 0 0]'; % matriz de alimentação direta do observador 
%  espaço de estados para o obesrvador "convencional"
esObs = ss(Ao,Bo,Co,Do);
%%
rank(ctrb(Aol,Bol))
rank(obsv(Aol,C1))
Q = [1 0 0 0; 0 0.2 0 0; 0 0 10 0; 0 0 0 0.2]*100;
R = 0.5;
% projeto do controlaad
KL = lqr(Aol,Bol,Q,R);
eig(Aol-Bol*KL)
% x dx teta dteta
x0 = [0.2 0 0.5 0]'

% ganho para ajuste de erro nulo à resposta ao degrau
N1 = -inv(C1*inv(Aol-Bol*KL)*Bol)
%%Observador ótimo de Kalman

% sistema com ruído na saída:
CK = [1 0 0 0];
DK= zeros (size(Col,1),size(Bol, 2));
rank(ctrb(Aol,Bol))
rank(obsv(Aol,C1))
%% Augment system with disturbances and noise
Vd = 0.1*eye(4); % disturbance variance
Vn = 1; % noise variance
BF = [Bol Vd 0*Bol];
% augment inputs with disturbance and noise
sysC= ss (Aol, BF, C1, [0 0 0 0 0 Vn]); % build state space system 
sysFullOutput = ss (Aol, BF, eye (4), zeros (4, size (BF, 2))); % system with full state output, disturbance, no noise
%%Build Kalman filter
% design Kalman filter with lqe function
[Kf, P,E] = lqe (Aol, Vd, C1, Vd, Vn);
%Kf = (lqr (Aol',CK', Vd, Vn))'; %% possible to design with "LQR" code
sysKF = ss(Aol-Kf*C1,[Bol Kf], eye(4), 0*[Bol Kf]); % Kalman filter

out = sim('TwoMassComObservadorKalman.slx');
x1 = out.Data.signals(1).values;
dx1 = out.Data.signals(2).values;
x2 = out.Data.signals(3).values;
dx2 = out.Data.signals(4).values;
u_noise = out.Data.signals(5).values;
y_noise =out.Data.signals(6).values;
t = out.Data.time;
figure
plot(t,x1(:,1),'r',t,x1(:,2),'-g',t,x1(:,3),'--b','LineWidth',2);


