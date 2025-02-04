%% BEAM 2D + BAR/TRUSS2D 
% 
% this finite element overlaps the stiffness for a BEAM in a 2D space with
% a BAR or Truss element in a 2D space
%
% % k=EI/L^3* (beam component)
%             [ 12s^2 | -12cs | -6Ls  || -12s^2 |  12cs  | -6Ls ]
%             [-12cs  | 12c^2 |  6Lc  ||  12cs  | -12c^2 |  6Lc ]
%             [ -6Ls  |   6Lc | 4L^2  ||  +6Ls  |  -6Lc  | 2L^2 ]
%              --------------------------------------------------
%             [-12s^2 |  12cs |  6Ls  ||  12s^2 | -12cs  |  6Ls ]
%             [ 12cs  |-12c^2 | -6Lc  || -12cs  |  12c^2 | -6Lc ]
%             [ -6Ls  |   6Lc | 2L^2  ||  +6Ls  |  -6Lc  | 4L^2 ]%
%  +
%     EA/L* (truss/bar component
%             [  c^2  |   cs  |   0   ||  -c^2  |  -cs   |   0  ]
%             [  cs   |  s^2  |   0   ||  -cs   |  -s^2  |   0  ]
%             [   0   |   0   |   0   ||   0    |    0   |   0  ]
%              -------------------------------------------------
%             [ -c^2  |  -cs  |   0   ||   c^2  |   cs   |   0  ]
%             [  -cs  | -s^2  |   0   ||   cs   |   s^2  |   0  ]
%             [   0   |   0   |   0   ||   0    |    0   |   0  ]

% beam with 2 nodes
% each node has 3 dofs (degrees-of-freedom): 
%      longitudinal displacement u
%      transverse displacement v  
%      rotation Rx
%
%     ^           ^
%     |u1,v1,Rx1  |u2,v2,Rx2
%     o-----------o
%
%
%

%% Example
%
%    20kN
% ^  <--- o----------o 12kN.m(moment)
% |     2 |          | 3
% 3       |          |
% |       |          | 
% v     1 |          | 4
%        FFF  <-4-> FFF 
%%



%% Parâmetros do modelo de FEM
clear; clc;close all;
nelem=3; %número de elementos (vigas) que constitui a frame
nnode=4; %número de nós
ndof=nnode*3; %considerando 3 graus de liberdade por nó

%% Parâmetros materiais e geométricos
E=210e9;
Es=E*ones(nelem,1); %vetor de módulo de young das barras

A=2e-2; %área da secção das barras (tubo de secção circular com diâmetro exterior 50mm e interior 40mm
As=A*ones(nelem,1); %vetor de áreas para a totalidade de barras que constituem a treliça

I=5e-5;
Is=I*ones(nelem,1);

Ls=zeros(nelem,1);Ls=[3;4;3];

theta=[pi/2 0 -pi/2]; %orientação de cada barra
%localizações dos nós
xyU=[0 0;0 3; 4 3; 4 0 ];

lnodal={[1 2],[2 3],[3 4]}; %pares de nós que definem cada elemento (barra)

%% Condições de fronteira
% graus de liberdade restritos: 
% nó1 e nó 4: direção X(u), Y(v) e rotação(Rx)
Rdof=[1,2,3,10,11,12];

% força horizontal de 20kN no nó 2 (corresponde ao dof 4)
% momento de 12kN.m no nó 3 (correspondente ao dof 9)
Fdof=[4 9]; %vetor de dofs carregados
VF=[-20e3 12e3]; %Vetor de cargas

%% Cálculo da Matriz de rigidez

K=zeros(ndof,ndof); %inicialização da matriz de rigidez global

for ielem=1:nelem
    
% matriz de rigidez do elemento de barra (2dof por nó=4dof por elemento)
Kbar= As(ielem)*Es(ielem)/Ls(ielem);

%matriz de transformação referencial barra para referencial global
c=cos(theta(ielem));s=sin(theta(ielem));c2=c^2;s2=s^2;cs=c*s;
T=[c2 cs 0 -c2 -cs 0;cs s2 0 -cs -s2 0;zeros(1,6);-c2 -cs 0 c2 cs 0;-cs -s2 0 cs s2 0;zeros(1,6)];
Kelem1=Kbar*T;


Kbeam= Es(ielem)*Is(ielem)/Ls(ielem)^3;

L=Ls(ielem);
Kelem2=Kbeam*[ 12*s2  -12*cs  -6*L*s  -12*s2  12*cs -6*L*s;
   -12*cs   12*c2   6*L*c   12*cs -12*c2  6*L*c;
   -6*L*s   6*L*c   4*L^2   6*L*s  -6*L*c  2*L^2;
   -12*s2  12*c*s   6*L*s  12*s2  -12*cs  6*L*s;
    12*cs  -12*c2  -6*L*c  -12*c*s  12*c2 -6*L*c;
   -6*L*s   6*L*c  2*L^2   6*L*s  -6*L*c  4*L^2 ];

Kelem=Kelem1+Kelem2;%

%montagem da matriz de rigidez do elemento na matriz de rigidez global
dof1=2*lnodal{ielem}(1)-1; %1º dof do 1º nó do elemento ielem
dof2=2*lnodal{ielem}(2)-1; %1º dof do 2º nó do elemento ielem

K(dof1:dof1+2,dof1:dof1+2)=K(dof1:dof1+2,dof1:dof1+2)+Kelem(1:3,1:3);
K(dof1:dof1+2,dof2:dof2+2)=K(dof1:dof1+2,dof2:dof2+2)+Kelem(1:3,4:6);
K(dof2:dof2+2,dof1:dof1+2)=K(dof2:dof2+2,dof1:dof1+2)+Kelem(4:6,1:3);
K(dof2:dof2+2,dof2:dof2+2)=K(dof2:dof2+2,dof2:dof2+2)+Kelem(4:6,4:6);
end


%Determinação do vetor de carregamento
F=zeros(ndof,1);
F(Fdof)=VF;

%Aplicação das condições de fronteira
free_dofs=setxor(1:ndof,Rdof); %graus de liberdade não restritos
Kp=K(free_dofs,free_dofs);
Fp=F(free_dofs,1);

% Resolução do sistema de equações
Up=Kp\Fp;
U=zeros(ndof,1);
U(free_dofs)=Up;

fvector=K*U;


% stress=zeros(nelem,1);
% 
% for ielem=1:nelem
%    no1=lnodal{ielem}(1);%1º nó do elemento ielem
%    no2=lnodal{ielem}(2); %2º nó do elemento ielem
%    dof1=2*lnodal{ielem}(1)-1; %1º dof do nó 1 do elemento ielem
%    dof2=2*lnodal{ielem}(2)-1; %1º dof do nó 2 do elemento ielem 
% 
%     T=[cos(theta(ielem)) sin(theta(ielem)) 0 0; 0 0 cos(theta(ielem)) sin(theta(ielem))];
%     stress(ielem,1)=Es(ielem)/Ls(ielem)*[-1,1]*T*[U(dof1:dof1+1);U(dof2:dof2+1)];
% end

%% output de resultados

% listagem de deslocamentos
for inode=1:nnode
    dx=U(inode*3-2);
    dy=U(inode*3-1);
    Rx=U(inode*3);
    fprintf('d%ix = %7.3e\n', inode,dx);
    fprintf('d%iy = %7.3e\n', inode,dy);
    fprintf('d%iRx = %7.3e\n', inode,Rx);
end
%     % listagem de forças
% for inode=1:nnode
%     fx=fvector(inode*2-1);
%     fy=fvector(inode*2);
%     fprintf('F%ix = %7.3e\n', inode,fx);
%     fprintf('F%iy = %7.3e\n', inode,fy);
% end

% listagem de tensões normais

%here

%representação da deformada

%aF =10* (max(U)-min(U))/(max(max(xyU))-min(min(xyU))); %factor de amplificação da deformada
aF=1e-4;
% cm=colormap('jet');
% deltaS=max(stress)-min(stress);

for ielem=1:nelem
   no1=lnodal{ielem}(1);%1º nó do elemento ielem
   no2=lnodal{ielem}(2); %2º nó do elemento ielem
   dof1=(no1-1)*3+1; %1º dof do nó 1 do elemento ielem
   dof2=(no2-1)*3+1; %1º dof do nó 2 do elemento ielem 
   LinesU(ielem,1:4)=[xyU(no1,1) xyU(no2,1) xyU(no1,2)  xyU(no2,2)];
   LinesUD(ielem,1:4)=LinesU(ielem,1:4)+1/aF*[U(dof1) U(dof2) U(dof1+1) U(dof2+1)];
   %stress_color(ielem,1:3)=cm(floor((stress(ielem,1)-min(stress))*63/deltaS)+1,1:3);
end

figure(1); hold on; axis equal
for i=1:size(LinesU,1);
    undef=line(LinesU(i,1:2),LinesU(i,3:4)); set(undef,'LineStyle',':','color','b','Linewidth',.5);
    def(i)=line(LinesU(i,1:2),LinesU(i,3:4)); set(def,'LineStyle','-','color','b','Linewidth',2);
end%
%cbh=colorbar;
%cbh.Ticks=linspace(0,1,10);
%cbh.TickLabels=num2cell(linspace(min(stress),max(stress),10));
hold off

% animação

%animdef=[NaN,NaN];
max_anim_steps=100;time_step=0.01;
for anim_step=1:max_anim_steps
    pause(time_step);
    for ielem=1:nelem
        no1=lnodal{ielem}(1);%1º nó do elemento ielem
        no2=lnodal{ielem}(2); %2º nó do elemento ielem
        dof1=(no1-1)*3+1; %1º dof do nó 1 do elemento ielem
        dof2=(no2-1)*3+1; %1º dof do nó 2 do elemento ielem 
        LinesUD(ielem,1:4)=LinesU(ielem,1:4)+anim_step/max_anim_steps*...
            1/aF*[U(dof1) U(dof2) U(dof1+1) U(dof2+1)];
       
        set(def(ielem),'XData',LinesUD(ielem,1:2),'YData',LinesUD(ielem,3:4));
    end
   % set(def(:),'XData',LinesUD(:,1:2),'YData',LinesUD(:,3:4));
    
end


        

