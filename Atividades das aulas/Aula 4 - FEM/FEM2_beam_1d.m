%%beam 1D

% k=EI/L^3*[12 6L -12 6L;6L 4L^2 -6L 2L^2;-12 -6L 12 -6L;6L 2L^2 -6L 4L^2];
% beam with 2 nodes
% each node has 2 dofs (degrees-of-freedom): 
%      transverse displacement v
%      rotation theta
%
%     ^           ^
%     |v1,theta1  |v2,theta2
%     o-----------o
%
%
%

%% Example

% Beam with 4 elements fixed at one end and 
% simply supported at the other end 2
% with vertical load at midspan

%% FEM model parameters
clear; clc;close all;
nelem=4; %n�mero de elementos de viga que constitui a viga
nnode=5; %n�mero de n�s
ndof=nnode*2; %considerando 2 graus de liberdade por n� (v, theta)

%% Par�metros materiais e geom�tricos
E=210e9; %Pa (a�o)
Es=E*ones(nelem,1); %vetor de m�dulo de young das vigas

I=0.04*0.08^3/12; %Momento de in�rcia de sec��o retangular (40x80mm)
Is=I*ones(nelem,1); %vetor de momentos de in�rcia para as vigas


L=1; % comprimento de cada viga (1 metro)
Ls=L*ones(nelem,1);


%localiza��es dos n�s e liga��es nodais
xyU=zeros(nnode,2);
lnodal={nelem};
for ielem=1:nelem
    xyU(ielem,1)=L*(ielem-1);
    xyU(ielem+1,1)=L*ielem;
    
    
    lnodal{ielem}=[ielem, ielem+1];
end

%% Condi��es de fronteira
% graus de liberdade restritos: n�1_v, n�1_theta, n�_end_v
Rdof=[1,2];

% carregamento vertical de 1000N no n� 3 (corresponde ao dof 5)
Fdof=[10]; %vetor de dofs carregados
VF=[1000]; %Vetor de cargas

%% C�lculo da Matriz de rigidez

K=zeros(ndof,ndof); %inicializa��o da matriz de rigidez global

for ielem=1:nelem
% matriz de rigidez do elemento de viga (2dof por n�=4dof por elemento)
Kbeam= Es(ielem)*Is(ielem)/Ls(ielem)^3;

Kelem=Kbeam*[12 6*L -12 6*L;6*L 4*L^2 -6*L 2*L^2;-12 -6*L 12 -6*L;6*L 2*L^2 -6*L 4*L^2];


%montagem da matriz de rigidez do elemento na matriz de rigidez global
dof1=2*lnodal{ielem}(1)-1; %1� dof do 1� n� do elemento ielem
dof2=2*lnodal{ielem}(2)-1; %1� dof do 2� n� do elemento ielem

K(dof1:dof1+1,dof1:dof1+1)=K(dof1:dof1+1,dof1:dof1+1)+Kelem(1:2,1:2);
K(dof1:dof1+1,dof2:dof2+1)=K(dof1:dof1+1,dof2:dof2+1)+Kelem(1:2,3:4);
K(dof2:dof2+1,dof1:dof1+1)=K(dof2:dof2+1,dof1:dof1+1)+Kelem(3:4,1:2);
K(dof2:dof2+1,dof2:dof2+1)=K(dof2:dof2+1,dof2:dof2+1)+Kelem(3:4,3:4);
end

%% Carregamento e condi��es de fronteira
%Determina��o do vetor de carregamento
F=zeros(ndof,1);
F(Fdof)=VF;

%Aplica��o das condi��es de fronteira
free_dofs=setxor(1:ndof,Rdof); %graus de liberdade n�o restritos
Kp=K(free_dofs,free_dofs);
Fp=F(free_dofs,1);

%% Resolu��o do sistema de equa��es
Up=Kp\Fp;
U=zeros(ndof,1);
U(free_dofs)=Up;

%% output de resultados

% listagem de deslocamentos
for inode=1:nnode
    v=U(inode*2-1);
    theta=U(inode*2);
    fprintf('d%iv = %7.3e\n', inode,v);
    fprintf('d%itheta = %7.3e\n', inode,theta);
end
  
%representa��o da deformada


aF =10* (max(U)-min(U))/(max(max(xyU))-min(min(xyU))); %factor de amplifica��o da deformada


for ielem=1:nelem
   no1=lnodal{ielem}(1);%1� n� do elemento ielem
   no2=lnodal{ielem}(2); %2� n� do elemento ielem
   dof1=2*lnodal{ielem}(1)-1; %1� dof do n� 1 do elemento ielem
   dof2=2*lnodal{ielem}(2)-1; %1� dof do n� 2 do elemento ielem 
   LinesU(ielem,1:4)=[xyU(no1,1) xyU(no2,1) xyU(no1,2)  xyU(no2,2)];
   LinesUD(ielem,1:4)=LinesU(ielem,1:4)+1/aF*[0 0 U(dof1) U(dof2)];
   
end

figure(1); hold on; axis equal
for i=1:size(LinesU,1);
    undef(i)=line(LinesU(i,1:2),LinesU(i,3:4)); set(undef,'LineStyle',':','color','b','Linewidth',.5);
    def(i)=line(LinesUD(i,1:2),LinesUD(i,3:4)); set(def,'LineStyle','-','color','b','Linewidth',2);
end
hold off

% anima��o

%animdef=[NaN,NaN];
max_anim_steps=100;time_step=0.01;
for anim_step=1:max_anim_steps
    pause(time_step);
    for ielem=1:nelem
        no1=lnodal{ielem}(1);%1� n� do elemento ielem
        no2=lnodal{ielem}(2); %2� n� do elemento ielem
        dof1=2*lnodal{ielem}(1)-1; %1� dof do n� 1 do elemento ielem
        dof2=2*lnodal{ielem}(2)-1; %1� dof do n� 2 do elemento ielem
        LinesUD(ielem,1:4)=LinesU(ielem,1:4)+anim_step/max_anim_steps*...
            1/aF*[0 0 U(dof1) U(dof2)];
        set(def(ielem),'XData',LinesUD(ielem,1:2),'YData',LinesUD(ielem,3:4));
    end
   % set(def(:),'XData',LinesUD(:,1:2),'YData',LinesUD(:,3:4));
    
end

