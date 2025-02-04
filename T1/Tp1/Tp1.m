%% Trabalho Elaborado por:

% Magner Gusse (110180)
% Pedro Mendes (108063)

%% Ex. 1 b
clear all; clc;

% Variávies M (Massa), C (amortecimento) e K (rigidez)
M1=1; M2=1;M3=1; 
K1=10; K2=10; K3=10; K4=10; K5=10;
C1=1; C2=1; C3=1; C4=1; C5=1;

% Posições iniciais:
x10 = 0; x20=0.3; x30=0.3;
x0=[x10;0;x20;0;x30;0];

% Variáveis de estado em open loop
A=[0 1 0 0 0 0;-(K1+K2+K3)/M1 -(C1+C2+C3)/M1 K2/M1 C2/M1 K3/M1 C3/M1;...
    0 0 0 1 0 0; K2/M2 C2/M2 -(K2+K4+K5)/M2 -(C2+C4+C5)/M2 K4/M2 C4/M2;...
    0 0 0 0 0 1;K5/M3 C3/M3 K4/M3 C4/M3 -(K3+K4)/M3 -(C3+C4)/M3];
B=[0; 0; 0 ; 1/M2; 0; 0 ];
C=[eye(6)];  % Ver todas variáveis de estado à saída
D=[0;0;0;0;0;0];

sys=ss(A,B,C,D); % Define o modelo de espaço de estados com as matrizes A, 
% B, C e D

% Simulação da resposta do sistema descrito por 'sys' a partir da condição 
% inicial 'x0' por 15 segundos:
[Y,T,X]=initial(sys, x0,15);

% Plot da resposta do sistema:
figure(1);
plot(T,Y(:,1),'r',T,Y(:,3),'k',T,Y(:,5),'b','LineWidth',1);
legend('x1','x2','x3');
title('Posições das Massas - Ex.1');
xlabel('Tempo (s)');
ylabel('Amplitude (m)');
xlim([0 15]);
grid on;

%% 1. c

% Como para uma respota limitada, o sistema tende para um valor fixo,
% limite, 0, o sistema é estável. Para confirmar calcularam-se os valores
% próprios da matriz:

eig(A)

% Polos (malha aberta):
% -0.2929 + 2.4025i
% -0.2929 - 2.4025i
% -2.0000 + 6.0000i
% -2.0000 - 6.0000i
% -1.7071 + 5.5882i
% -1.7071 - 5.5882i


% Confirma-se a estabilidade, já que todos os polos encontram-se no 
% semi plano esquerdo.

%% Ex. 1 d

% Forma 1: Uso direto da matriz C para selecionar a posição da massa 3
C_x3 = [0 0 0 0 1 0]; % Matriz de saída para observar apenas m3
sys_x3 = ss(A, B, C_x3, 0); % Modelo de espaço de estados para observar x3
figure(2);
bode(sys_x3); % Plot do diagrama de bode
grid on;
title('Diagrama de Bode para a Massa 3');

% Forma 2: Conversão do sistema para uma função de transferência
f_x3 = zpk(tf(sys_x3));
figure(3);
bode(f_x3);
grid on;
title('Diagrama de Bode para a Massa 3');

%% Ex. 2 Perceber se o sistema é observável e controlável
n = size(A, 1); % Determinação do número de estados (dimensão da matriz A)
Cc = ctrb(A,B); % Cálculo da matriz de controlabilidade
if rank(Cc) == n
    fprintf('O sistema é controlável.\n');
else
    fprintf('O sistema não é controlável.\n');
end

Mo = obsv(A,C); % Cálculo da matriz de observabilidade
if rank(Mo) == n
    fprintf('O sistema é observável.\n');
else
    fprintf('O sistema não é observável.\n');
end

% R: O sistema é controlável.
%    O sistema é observável.

%% Ex. 2 Escolha dos Polos

% Definição dos polos desejados do sistema fechado:
P = [-3.2-2.1j,-3.2+2.1j,-2.1+1.05j,-2.1-1.05j,-2.05,-3.05];

K=place(A,B,P);% Cálculo da matriz de ganho K para realimentação de estados

Acl=A-B*K; % Cálculo da nova matriz A do sistema em malha fechada
syscl = ss(Acl,B,C,D); % Criação do sistema em malha fechada com 
% realimentação de estados

[YCL,TCL,XCL]=initial(syscl,x0,10); % 'initial' simula a resposta a x0 
% durante 10s

% Plot da resposta do estado x2 (massa 2) do sistema original e do sistema 
% em malha fechada

% Foi utilizada a posição x2, já que era nesta massa que se verificava
% maiores oscilações

figure(4);
plot(T,Y(:,3),'r',TCL,YCL(:,3),'k','LineWidth',2);
grid on
xlabel('Tempo (s)');
ylabel('Amplitude (m)');
xlim([0 15]);


%% Ex.2 Com Função Sim

% Simulação (Simulink) deste sistema, e obtenção de valores após controlo
% do sistema em malha fechada.
out2 = sim('Ex2.slx');
x1 = out2.Data2.signals(1).values;
dx1 = out2.Data2.signals(2).values;
x2 = out2.Data2.signals(3).values;
dx2 = out2.Data2.signals(4).values;
x3 = out2.Data2.signals(5).values;
dx3 = out2.Data2.signals(6).values;
t = out2.Data2.time;


% Obtenção de dados Open Loop (Sem controlo)
x10 = out2.Data2s.signals(1).values;
dx10 = out2.Data2s.signals(2).values;
x20 = out2.Data2s.signals(3).values;
dx20 = out2.Data2s.signals(4).values;
x30 = out2.Data2s.signals(5).values;
dx30 = out2.Data2s.signals(6).values;

% Plot das Posições das 3 Massas, com e sem controlo.
figure(5)
subplot(3,1,1)
plot(t,x1(:,1),'r',t,x10(:,1),'b','LineWidth',1);
legend('x1 c/controlo','x1 s/controlo','Location','best')
title('Massa 1');
grid on
xlabel('Tempo (s)');
ylabel('Amplitude (m)');
xlim([0 10]);


subplot(3,1,2)
plot(t,x2(:,1),'r',t,x20(:,1),'b','LineWidth',1);
legend('x2 c/controlo','x2 s/controlo','Location','best')
title('Massa 2');
grid on
xlabel('Tempo (s)');
ylabel('Amplitude (m)');
xlim([0 10]);


subplot(3,1,3)
plot(t,x3(:,1),'r',t,x30(:,1),'b','LineWidth',1);
legend('x3 c/controlo','x3 s/controlo','Location','best')
title('Massa 3');
grid on
xlabel('Tempo (s)');
ylabel('Amplitude (m)');
xlim([0 10]);


%% Ex. 3 

% Verificação da observabilidade e controlabilidade de X3 (enunciado):
C_x3 = [0 0 0 0 1 0]; % O sistema tem só a saída Y = x3
Mo = obsv(A,C_x3); % Cálculo da matriz de observabilidade

if rank(Mo) == n
    fprintf('O sistema com x3 é observável.\n');
else
    fprintf('O sistema com x3 não é observável.\n');
end

% R: O sistema não é observável, logo deve-se optar por outra variável de
% estado, como por exemplo x2, e voltar a verificar a observabilidade

C_x2 = [0 0 1 0 0 0];
Mo = obsv(A,C_x2);
if rank(Mo) == n
    fprintf('O sistema com x2 é observável.\n');
else
    fprintf('O sistema com x2 não é observável.\n');
end

% A controlabilidade do sistema foi assegurada anteriormente, já que apenas
% implica as matrizes A e B. 


% R: O sistema é observável, pelo que se pode observar x2.
AT = A'; CT = C_x2';

% Definição da posição desejada dos polos do observador
Pobs = [-9.6 - 6.3j, -9.6 + 6.3j, -6.3 + 3.15j, -6.3 - 3.15j, -6.15, -9.15];

% Cálculo da matriz de ganhos do observador (L transposta)
LT = place(AT, CT,Pobs); %'place' calcula ganhos para colocar polos em Pobs
L = LT'; % Transposição para obter a matriz de ganho L do observador

% Definição das matrizes do observador
Ao = A-L*C_x2; % Matriz A ajustada para o observador
Co = eye(6); % Matriz de saída do observador: todos os estados à saída
Bo = [B L]; % Matriz de entrada: o observador usa como entradas u e y
Do = [0 0 0 0 0 0; 0 0 0 0 0 0]'; % Matriz de alimentação direta

% Criação do sistema em espaço de estados para o observador
sysobs = ss(Ao,Bo,Co,Do);

% Simulação da resposta do observador com condição inicial
% 'initial' calcula a resposta temporal para as condições iniciais (x0)
[Y0,T0,X0]= initial(sysobs,x0,10);

% Simulação (Simulink) deste sistema, e obtenção de valores após controlo
% do sistema em malha fechada, com observador
out3 = sim('Ex3.slx');
x1_3 = out3.Data3.signals(1).values;
x2_3 = out3.Data3.signals(2).values;
x3_3 = out3.Data3.signals(3).values;
t3 = out3.Data3.time;

% Obtenção de dados Open Loop (Sem controlo), com observador
x10_3 = out3.Data3s.signals(1).values;
x20_3 = out3.Data3s.signals(2).values;
x30_3 = out3.Data3s.signals(3).values;

% Plot das Posições das 3 Massas, com e sem controlo.
figure(6)
subplot(3,1,1)
plot(t3,x1_3(:,1),'r',t3,x10_3(:,1),'b','LineWidth',1);
legend('x1 c/controlo', 'x1 s/controlo','Location','best')
title(' Ex.3 - Massa 1');
grid on
xlabel('Tempo (s)');
ylabel('Amplitude (m)');
xlim([0 10]);

subplot(3,1,2)
plot(t3,x2_3(:,1),'r',t3,x20_3(:,1),'b','LineWidth',1);
legend('x2 c/controlo', 'x2 s/controlo','Location','best')
title(' Ex.3 - Massa 2');
grid on
xlabel('Tempo (s)');
ylabel('Amplitude (m)');
xlim([0 10]);

subplot(3,1,3)
plot(t3,x3_3(:,1),'r',t3,x30_3(:,1),'b','LineWidth',1);
legend('x3 c/controlo', 'x3 s/controlo','Location','best')
title('Ex.3 - Massa 3');
grid on
xlabel('Tempo (s)');
ylabel('Amplitude (m)');
xlim([0 10]);

%% Ex. 4

% Configuração do controlador LQR
Q = diag([1,1,1,1,100,20]); % Matriz de ponderação dos estados
R=0.5; % Ponderação para o esforço de controlo
Kr=lqr(A,B,Q,R);  % Cálculo do ganho do controlador LQR
Ar=A-B*Kr; % Matriz dinâmica do sistema controlado (A ajustada)

% Simulação do sistema controlado (usando Simulink)
out4 = sim('Ex4.slx');

% Extração dos dados das variáveis de estado do sistema controlado
x1_4 = out4.Data4.signals(1).values;
x2_4 = out4.Data4.signals(3).values;
x3_4 = out4.Data4.signals(5).values;

dx1_4 = out4.Data4.signals(2).values;
dx2_4 = out4.Data4.signals(4).values;
dx3_4 = out4.Data4.signals(6).values;

t4 = out4.Data4.time; % Vetor de tempo da simulação com controlo


% Extração dos dados das variáveis de estado sem controlo
x10_4 = out4.Data4s.signals(1).values;
x20_4 = out4.Data4s.signals(3).values;
x30_4 = out4.Data4s.signals(5).values;

dx10_4 = out4.Data4s.signals(2).values;
dx20_4 = out4.Data4s.signals(4).values;
dx30_4 = out4.Data4s.signals(6).values;


% Plot das Posições das 3 Massas, com e sem controlo:
figure(7)
subplot(3,1,1)
plot(t4,x1_4(:,1),'r',t4,x10_4(:,1),'b','LineWidth',1); 
legend('x1 c/controlo', 'x1 s/controlo','Location','best')
title(' Ex.4 - Posição Massa 1');

grid on
xlabel('Tempo (s)');
ylabel('Amplitude (m)');
xlim([0 10]);

subplot(3,1,2)
plot(t4,x2_4(:,1),'r',t4,x20_4(:,1),'b','LineWidth',1); 
legend('x2 c/controlo', 'x2 s/controlo','Location','best')
title(' Ex.4 - Posição Massa 2');

grid on
xlabel('Tempo (s)');
ylabel('Amplitude (m)');
xlim([0 10]);

subplot(3,1,3)
plot(t4,x3_4(:,1),'r',t4,x30_4(:,1),'b','LineWidth',1); 
legend('x3 c/controlo', 'x3 s/controlo','Location','best')
title('Ex.4 - Posição Massa 3');

grid on
xlabel('Tempo (s)');
ylabel('Amplitude (m)');
xlim([0 10]);

% Conclusão sobre a resposta do sistema:
% Observa-se que a massa 3 (m3) estabiliza mais rapidamente que as outras
% variáveis de estado, especialmente em comparação com a massa 2 (m2).

% Plot das Velocidades das 3 Massas, com e sem controlo:
figure(8)
subplot(3,1,1)
plot(t4,dx1_4(:,1),'r',t4,dx10_4(:,1),'b','LineWidth',1); 
legend('x1 c/controlo', 'x1 s/controlo','Location','best')
title(' Ex.4 - Velocidade Massa 1');
grid on
xlabel('Tempo (s)');
ylabel('Velocidade (m/s)');
xlim([0 10]);

subplot(3,1,2)
plot(t4,dx2_4(:,1),'r',t4,dx20_4(:,1),'b','LineWidth',1); 
legend('x2 c/controlo', 'x2 s/controlo','Location','best')
title(' Ex.4 - Velocidade Massa 2');
grid on
xlabel('Tempo (s)');
ylabel('Velocidade (m/s)');
xlim([0 10]);

subplot(3,1,3)
plot(t4,dx3_4(:,1),'r',t4,dx30_4(:,1),'b','LineWidth',1); 
legend('x3 c/controlo', 'x3 s/controlo','Location','best')
title('Ex.4 - Velocidade Massa 3');
grid on
xlabel('Tempo (s)');
ylabel('Velocidade (m/s)');
xlim([0 10]);

% Simulação e obtenção do esforço da força F2
f2 = out4.forcaf2.signals(1).values;
t4 = out4.forcaf2.time;

figure(9)
plot(t4, f2,'r','LineWidth',2);
legend('Força F2','Location', 'best');
title('Esforço de F2');
grid on;
xlabel('t (s)');
ylabel('Força (N)');
xlim([0 10]);

% com este controlador, o esforço de f2 é bastante inferior a 50
% f2_max=6.16

% Justificação de Valores: Optou-se por atribuir um peso de 20 a dx3, já 
% que existe um bom compromisso entre esforço f2 e estabilização (pouca 
% oscilação). De facto x1 e x3 estabilizam aproximadamente ao mesmo tempo, 
% apesar das diferentes tentativas de valores de R e Q, pelo que o menor 
% esforço é o ideal. 

%% 5

% Cálculo do fator de pré-alimentação (N) para eliminar o erro em regime 
% permanente:

N = -1 / (C_x3 * ((A - B * Kr) \ B));

% N é calculado com base na relação entre as matrizes do sistema (A, B, C) 
% e o ganho de controlo (Kr).
% O objetivo é ajustar o sistema para garantir que o sinal de referência 
% esteja em concordância com o step unitário definido.

% Simulação do sistema com controlo (usando Simulink)
out5 = sim('Ex5.slx');

% Extração dos dados das variáveis de estado do sistema controlado
x1_5 = out5.Data5.signals(1).values;
x2_5 = out5.Data5.signals(3).values;
x3_5 = out5.Data5.signals(5).values;

dx1_5 = out5.Data5.signals(2).values;
dx2_5 = out5.Data5.signals(4).values;
dx3_5 = out5.Data5.signals(6).values;
t5 = out5.Data5.time; % Vetor de tempo da simulação com controlo

% Extração dos dados das variáveis de estado do sistema não controlado
x10_5 = out5.Data5s.signals(1).values;
x20_5 = out5.Data5s.signals(3).values;
x30_5 = out5.Data5s.signals(5).values;

dx10_5 = out5.Data5s.signals(2).values;
dx20_5 = out5.Data5s.signals(4).values;
dx30_5 = out5.Data5s.signals(6).values;

% Plot das Posições das 3 Massas, com e sem controlo:
figure(10)
subplot(3,1,1)
plot(t5,x1_5(:,1),'r',t5,x10_5(:,1),'b','LineWidth',1); 
legend('x1 c/controlo', 'x1 s/controlo','Location','best')
title(' Ex.5 - Posição Massa 1');
grid on
xlabel('Tempo (s)');
ylabel('Amplitude (m)');
xlim([0 10]);

subplot(3,1,2)
plot(t5,x2_5(:,1),'r',t5,x20_5(:,1),'b','LineWidth',1); 
legend('x2 c/controlo', 'x2 s/controlo','Location','best')
title(' Ex.5 - Posição Massa 2');
grid on
xlabel('Tempo (s)');
ylabel('Amplitude (m)');
xlim([0 10]);

subplot(3,1,3)
plot(t5,x3_5(:,1),'r',t5,x30_5(:,1),'b','LineWidth',1); 
legend('x3 c/controlo', 'x3 s/controlo','Location','best')
title('Ex.5 - Posição Massa 3');
grid on
xlabel('Tempo (s)');
ylabel('Amplitude (m)');
xlim([0 10]);

% Gráfico da Evolução do Erro:
erro = out5.erro.signals(1).values;
t5 = out5.erro.time;

figure(11)
plot(t5,erro(:,1),'r','LineWidth',1); 
title('Ex.5 - Erro');
grid on
xlabel('Tempo (s)');
ylabel('Amplitude (m)');
xlim([0 10]);