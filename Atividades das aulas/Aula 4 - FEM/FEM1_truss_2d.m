%% treli�a 2D
%  5 ______6_____7 
%   /\    /\    /\
%  /  \  /  \  /  \
% /____\/____\/____\
%1     2      3     4

%% Par�metros do modelo de FEM
clear; clc;close all;
nelem=11; %n�mero de elementos (barras) que constitui a treli�a
nnode=7; %n�mero de n�s
ndof=nnode*2; %considerando 2 graus de liberdade por n�

%% Par�metros materiais e geom�tricos
E=70e9; %Pa (alum�nio)
Es=E*ones(nelem,1); %vetor de m�dulo de young das barras

A=(0.05^2-0.04^2)*pi; %�rea da sec��o das barras (tubo de sec��o circular com di�metro exterior 50mm e interior 40mm
As=A*ones(nelem,1); %vetor de �reas para a totalidade de barras que constituem a treli�a

L=1; % comprimento de cada barra (1 metro)
Ls=L*ones(nelem,1);

theta=[0 0 0 pi/3 2*pi/3 pi/3 2*pi/3 pi/3 2*pi/3 0 0 0 pi/2]; %orienta��o de cada barra
%localiza��es dos n�s
xyU=[0 0;L 0; 2*L 0; 3*L 0; L*cos(pi/3) L*sin(pi/3); L*cos(pi/3)+L L*sin(pi/3);L*cos(pi/3)+2*L L*sin(pi/3) ];

lnodal={[1 2],[2 3],[3 4],[1 5],[2 5],[2 6],[3 6],[3 7],[4 7], [5 6], [6 7]}; %pares de n�s que definem cada elemento (barra)

%% Condi��es de fronteira
% graus de liberdade restritos: n�1 dire��o X, n�1 dire��o y, n�3 dire��o y
Rdof=[1,2,8];

% carregamento vertical de 1000N no n� 6 (corresponde ao dof 12)
Fdof=[7]; %vetor de dofs carregados
VF=[-1000]; %Vetor de cargas

%% C�lculo da Matriz de rigidez

K=zeros(ndof,ndof); %inicializa��o da matriz de rigidez global

for ielem=1:nelem
% matriz de rigidez do elemento de barra (2dof por n�=4dof por elemento)
Kbar= As(ielem)*Es(ielem)/Ls(ielem);

%matriz de transforma��o referencial barra para referencial global
c=cos(theta(ielem));s=sin(theta(ielem));c2=c^2;s2=s^2;cs=c*s;
T=[c2 cs -c2 -cs;cs s2 -cs -s2;-c2 -cs c2 cs;-cs -s2 cs s2];
Kelem=Kbar*T;

%montagem da matriz de rigidez do elemento na matriz de rigidez global
dof1=2*lnodal{ielem}(1)-1; %1� dof do 1� n� do elemento ielem
dof2=2*lnodal{ielem}(2)-1; %1� dof do 2� n� do elemento ielem

K(dof1:dof1+1,dof1:dof1+1)=K(dof1:dof1+1,dof1:dof1+1)+Kelem(1:2,1:2);
K(dof1:dof1+1,dof2:dof2+1)=K(dof1:dof1+1,dof2:dof2+1)+Kelem(1:2,3:4);
K(dof2:dof2+1,dof1:dof1+1)=K(dof2:dof2+1,dof1:dof1+1)+Kelem(3:4,1:2);
K(dof2:dof2+1,dof2:dof2+1)=K(dof2:dof2+1,dof2:dof2+1)+Kelem(3:4,3:4);
end


%Determina��o do vetor de carregamento
F=zeros(ndof,1);
F(Fdof)=VF;

%Aplica��o das condi��es de fronteira
free_dofs=setxor(1:ndof,Rdof); %graus de liberdade n�o restritos
Kp=K(free_dofs,free_dofs);
Fp=F(free_dofs,1);

% Resolu��o do sistema de equa��es
Up=Kp\Fp;
U=zeros(ndof,1);
U(free_dofs)=Up;

fvector=K*U;


stress=zeros(nelem,1);

for ielem=1:nelem
   no1=lnodal{ielem}(1);%1� n� do elemento ielem
   no2=lnodal{ielem}(2); %2� n� do elemento ielem
   dof1=2*lnodal{ielem}(1)-1; %1� dof do n� 1 do elemento ielem
   dof2=2*lnodal{ielem}(2)-1; %1� dof do n� 2 do elemento ielem 

    T=[cos(theta(ielem)) sin(theta(ielem)) 0 0; 0 0 cos(theta(ielem)) sin(theta(ielem))];
    stress(ielem,1)=Es(ielem)/Ls(ielem)*[-1,1]*T*[U(dof1:dof1+1);U(dof2:dof2+1)];
end

%% output de resultados

% listagem de deslocamentos
for inode=1:nnode
    dx=U(inode*2-1);
    dy=U(inode*2);
    fprintf('d%ix = %7.3e\n', inode,dx);
    fprintf('d%iy = %7.3e\n', inode,dy);
end
    % listagem de for�as
for inode=1:nnode
    fx=fvector(inode*2-1);
    fy=fvector(inode*2);
    fprintf('F%ix = %7.3e\n', inode,fx);
    fprintf('F%iy = %7.3e\n', inode,fy);
end

% listagem de tens�es normais

%here

%representa��o da deformada


aF =10* (max(U)-min(U))/(max(max(xyU))-min(min(xyU))); %factor de amplifica��o da deformada
cm=colormap('jet');
deltaS=max(stress)-min(stress);

for ielem=1:nelem
   no1=lnodal{ielem}(1);%1� n� do elemento ielem
   no2=lnodal{ielem}(2); %2� n� do elemento ielem
   dof1=2*lnodal{ielem}(1)-1; %1� dof do n� 1 do elemento ielem
   dof2=2*lnodal{ielem}(2)-1; %1� dof do n� 2 do elemento ielem 
   LinesU(ielem,1:4)=[xyU(no1,1) xyU(no2,1) xyU(no1,2)  xyU(no2,2)];
   LinesUD(ielem,1:4)=LinesU(ielem,1:4)+1/aF*[U(dof1) U(dof2) U(dof1+1) U(dof2+1)];
   stress_color(ielem,1:3)=cm(floor((stress(ielem,1)-min(stress))*255/deltaS)+1,1:3);
end

figure(1); hold on; axis equal
for i=1:size(LinesU,1);
    undef=line(LinesU(i,1:2),LinesU(i,3:4)); set(undef,'LineStyle',':','color','b','Linewidth',.5);
    def=line(LinesUD(i,1:2),LinesUD(i,3:4)); set(def,'LineStyle','-','color',stress_color(i,1:3),'Linewidth',2);
end
cbh=colorbar;
cbh.Ticks=linspace(0,1,10);
cbh.TickLabels=num2cell(linspace(min(stress),max(stress),10));
hold off



        

