%%%% exercício Pórtico 2DOF %%%%
%%       _________
%       |___m2____|    m1=m2=1.2kg
%       |         |    K1=K2=4xKviga
%       | K2      |    Kviga=12EI/L^3
%       |_________|    E=70GPa
%       |___m1____|    I=b*h^3/12
%       |         |    b=10mm
%       | K1      |    h=2mm
%     __|_________|__
%     ///////////////
%
%%
m1=1.2; m2=m1;
b=10e-3;h=2e-3;L=0.150;
E=70e9;
I=b*h^3/12
Kviga=12*E*I/L^3
K1=4*Kviga
K2=K1;

M=[m1 0;0 m2];
K=[K1+K2 -K2; -K2 K2];

[a,b]=eig(K,M)

Phi=a
Wn=sqrt(diag(b))
Fn=Wn/2/pi

%normalização massa modal unitária

Phi'*M*Phi
Phi'*K*Phi


%% exercícios



% sistema 2dof com variantes

m1=1;m2=1;K1=1000;K2=100;
M=[m1 0;0 m2];
K=[K1+K2 -K2; -K2 K2];
[a,b]=eig(K,M);
Phi=a
Wn=sqrt(diag(b));Fn=Wn/2/pi

%
m1=1;m2=2;K1=1000;K2=1000;
M=[m1 0;0 m2];
K=[K1+K2 -K2; -K2 K2];
[a,b]=eig(K,M);
Phi=a
Wn=sqrt(diag(b));Fn=Wn/2/pi

%
m1=1;m2=1;K1=0;K2=1000; %sistema semidefinido
M=[m1 0;0 m2];
K=[K1+K2 -K2; -K2 K2];
[a,b]=eig(K,M);
Phi=a
Wn=sqrt(diag(b));Fn=Wn/2/pi

% SISTEMA 3DOF A
m1=1;m2=1;m3=1;
K1=1000;K2=1000;K3=1000;
M=[m1 0 0; 0 m2 0; 0 0 m3];
K=[K1+K2  -K2  0; -K2 K2+K3 -K3; 0 -K3 K3];

[a,b]=eig(K,M);
Phi=a
Wn=sqrt(diag(b));Fn=Wn/2/pi

%
m1=1;m2=1;m3=1;
K1=1000;K2=100000;K3=1000;
M=[m1 0 0; 0 m2 0; 0 0 m3];
K=[K1+K2  -K2  0; -K2 K2+K3 -K3; 0 -K3 K3];
[a,b]=eig(K,M);
Phi=a
Wn=sqrt(diag(b));Fn=Wn/2/pi


m1=1;m2=1;m3=1;
K1=1000;K2=1000;K3=100000;
M=[m1 0 0; 0 m2 0; 0 0 m3];
K=[K1+K2  -K2  0; -K2 K2+K3 -K3; 0 -K3 K3];
[a,b]=eig(K,M);
Phi=a
Wn=sqrt(diag(b));Fn=Wn/2/pi

% 3DOF B
m1=1;m2=1;m3=1;
K1=1000;K2=1000;K3=1000;K4=1000;
M=[m1 0 0; 0 m2 0; 0 0 m3];
K=[K1+K2+K3  -K2  -K3; -K2 K2+K4 -K4; -K3 -K4 K3+K4];
[a,b]=eig(K,M);
Phi=a
Wn=sqrt(diag(b));Fn=Wn/2/pi

%4dof xy
m1=1;m2=1;
K1x=1000;K1y=1000;K2x=1000;K2y=1000;
M=[m1 0 0 0; 0 m1 0 0; 0 0 m2 0; 0 0 0 m2];
K=[K1x 0 0 0; 0 K1y+K2y 0 -K2y;  0 0 K2x 0; 0 -K2y 0 K2y];
[a,b]=eig(K,M);
Phi=a
Wn=sqrt(diag(b));Fn=Wn/2/pi

m1=1;m2=1;
K1x=1000;K1y=10;K2x=1000;K2y=10;
M=[m1 0 0 0; 0 m1 0 0; 0 0 m2 0; 0 0 0 m2];
K=[K1x 0 0 0; 0 K1y+K2y 0 -K2y;  0 0 K2x 0; 0 -K2y 0 K2y];
[a,b]=eig(K,M);
Phi=a
Wn=sqrt(diag(b));Fn=Wn/2/pi

m1=1;m2=1;
K1x=10;K1y=1000;K2x=10;K2y=1000;
M=[m1 0 0 0; 0 m1 0 0; 0 0 m2 0; 0 0 0 m2];
K=[K1x 0 0 0; 0 K1y+K2y 0 -K2y;  0 0 K2x 0; 0 -K2y 0 K2y];
[a,b]=eig(K,M);
Phi=a
Wn=sqrt(diag(b));Fn=Wn/2/pi