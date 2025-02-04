%% Treliça 2D
%
%   5 ______6_____7
%    /\    /\    /\
%   /  \  /  \  /  \
%  /____\/____\/____\
% 1     2      3     4
%
%% Parâmetros do modelo de FEM
clear; clc;close all;
nelem=11; %número de elementos (barras) que constitui a treliça
nnode=7; %número de nós
ndof=nnode*2; %considerando 2 graus de liberdade por nó

%% Parâmetros materiais e geométricos
E=70e9; %Pa (alumínio)
Es=E*ones(nelem,1); %vetor de módulo de young das barras

A=(0.05^2-0.04^2)*pi; %área da secção das barras (tubo de secção circular com diâmetro exterior 50mm e interior 40mm
As=A*ones(nelem,1); %vetor de áreas para a totalidade de barras que constituem a treliça

L=1; % comprimento de cada barra (1 metro)
Ls=L*ones(nelem,1);

rho=2700; % massa volúmina do material
rhos=rho*ones(nelem,1);

theta=[0 0 0 pi/3 2*pi/3 pi/3 2*pi/3 pi/3 2*pi/3 0 0 0 pi/2]; %orientação de cada barra
%localizações dos nós
xyU=[0 0;L 0; 2*L 0; 3*L 0; L*cos(pi/3) L*sin(pi/3); L*cos(pi/3)+L L*sin(pi/3);L*cos(pi/3)+2*L L*sin(pi/3)];

lnodal={[1 2],[2 3],[3 4],[1 5],[2 5],[2 6],[3 6],[3 7],[4 7], [5 6], [6 7]}; %pares de nós que definem cada elemento (barra)

%% Condições de fronteira
% graus de liberdade restritos: nó1 direção X, nó1 direção y, nó4 direção y
Rdof=[1, 2, 7, 8];


%% Cálculo da Matriz de rigidez

M=zeros(ndof,ndof); %inicialização da matriz de rigidez global

K=zeros(ndof,ndof); %inicialização da matriz de rigidez global

for ielem=1:nelem

    % matriz de massa do elemento de barra
    Mbar= As(ielem)*Ls(ielem)*rhos(ielem);
    Melem=Mbar*0.5*eye(ndof);
    % matriz de rigidez do elemento de barra (2dof por nó=4dof por elemento)
    Kbar= As(ielem)*Es(ielem)/Ls(ielem);

    %matriz de transformação referencial barra para referencial global
    c=cos(theta(ielem));s=sin(theta(ielem));c2=c^2;s2=s^2;cs=c*s;
    T=[c2 cs -c2 -cs;cs s2 -cs -s2;-c2 -cs c2 cs;-cs -s2 cs s2];
    Kelem=Kbar*T;

    %montagem da matriz de massa e de rigidez do elemento nas matrizes globais
    no1=2*lnodal{ielem}(1)-1; %1º dof do 1º nó do elemento ielem
    no2=2*lnodal{ielem}(2)-1; %1º dof do 2º nó do elemento ielem

    M(no1:no1+1,no1:no1+1)=M(no1:no1+1,no1:no1+1)+Melem(1:2,1:2);
    M(no1:no1+1,no2:no2+1)=M(no1:no1+1,no2:no2+1)+Melem(1:2,3:4);
    M(no2:no2+1,no1:no1+1)=M(no2:no2+1,no1:no1+1)+Melem(3:4,1:2);
    M(no2:no2+1,no2:no2+1)=M(no2:no2+1,no2:no2+1)+Melem(3:4,3:4);

    K(no1:no1+1,no1:no1+1)=K(no1:no1+1,no1:no1+1)+Kelem(1:2,1:2);
    K(no1:no1+1,no2:no2+1)=K(no1:no1+1,no2:no2+1)+Kelem(1:2,3:4);
    K(no2:no2+1,no1:no1+1)=K(no2:no2+1,no1:no1+1)+Kelem(3:4,1:2);
    K(no2:no2+1,no2:no2+1)=K(no2:no2+1,no2:no2+1)+Kelem(3:4,3:4);
end




%Aplicação das condições de fronteira
free_dofs=setxor(1:ndof,Rdof); %graus de liberdade não restritos
Mp=M(free_dofs,free_dofs);
Kp=K(free_dofs,free_dofs);

%% Resolução do problema de valores e vetores proprios

[a,b]=eig(Kp,Mp);
Wn=sqrt(diag(b));Fn=Wn/2/pi();
Phi=a;



%% output de resultados

for imode=1:length(Wn)
    figure(1);clf;
    U=zeros(ndof,1);
    U(free_dofs)=Phi(:,imode);
    % listagem de deslocamentos
    for inode=1:nnode
        dx=U(inode*2-1);
        dy=U(inode*2);
        fprintf('d%ix = %7.3e\n', inode,dx);
        fprintf('d%iy = %7.3e\n', inode,dy);
    end


    %representação da deformada
    aF =10* (max(U)-min(U))/(max(max(xyU))-min(min(xyU))); %factor de amplificação da deformada


    for ielem=1:nelem
        no1=lnodal{ielem}(1);%1º nó do elemento ielem
        no2=lnodal{ielem}(2); %2º nó do elemento ielem
        dof1=2*lnodal{ielem}(1)-1; %1º dof do nó 1 do elemento ielem
        dof2=2*lnodal{ielem}(2)-1; %1º dof do nó 2 do elemento ielem
        LinesU(ielem,1:4)=[xyU(no1,1) xyU(no2,1) xyU(no1,2)  xyU(no2,2)];
        LinesUD(ielem,1:4)=LinesU(ielem,1:4)+1/aF*[U(dof1) U(dof2) U(dof1+1) U(dof2+1)];
        defmodal(ielem,1:4)=1/aF*[U(dof1) U(dof2) U(dof1+1) U(dof2+1)];
        LinesUD2(ielem,1:4)=LinesU(ielem,1:4)-1/aF*[U(dof1) U(dof2) U(dof1+1) U(dof2+1)];
    end

    figure(1); hold on; axis equal
    for i=1:size(LinesU,1);
        undef=line(LinesU(i,1:2),LinesU(i,3:4)); set(undef,'LineStyle',':','color','b','Linewidth',.5);
        def=line(LinesUD(i,1:2),LinesUD(i,3:4)); set(def,'LineStyle','-','color','k','Linewidth',2);
        def2=line(LinesUD2(i,1:2),LinesUD2(i,3:4)); set(def,'LineStyle','--','color','r','Linewidth',1);
    end


    hold off
    stringtitle=sprintf('modo: %i  @ [ %7.1f Hz ]',imode,round(Fn(imode)));
    title(stringtitle);
    pause
end


%% animacao dos modos
for imode=1:length(Wn)

figure(2);clf;
pause
U=zeros(ndof,1);
U(free_dofs)=Phi(:,imode);

aF =10* (max(U)-min(U))/(max(max(xyU))-min(min(xyU))); %factor de amplificação da deformada

for ielem=1:nelem
    no1=lnodal{ielem}(1);%1º nó do elemento ielem
    no2=lnodal{ielem}(2); %2º nó do elemento ielem
    dof1=2*lnodal{ielem}(1)-1; %1º dof do nó 1 do elemento ielem
    dof2=2*lnodal{ielem}(2)-1; %1º dof do nó 2 do elemento ielem
    LinesU(ielem,1:4)=[xyU(no1,1) xyU(no2,1) xyU(no1,2)  xyU(no2,2)];
end
hold on; axis equal
for i=1:size(LinesU,1);
    undef=line(LinesU(i,1:2),LinesU(i,3:4)); set(undef,'LineStyle',':','color','b','Linewidth',.5);
    def(i)=line(LinesU(i,1:2),LinesU(i,3:4)); set(def(i),'LineStyle','-','color','k','Linewidth',2);
end
stringtitle=sprintf('modo: %i  @ [ %7.1f Hz ]',imode,round(Fn(imode)));
    title(stringtitle);
    
cicle_num=4;
anim_steps=100;time_step=0.001;cicle_sin=sin(linspace(0,2*pi,anim_steps));
for icicle=1:cicle_num
    for anim_step=1:anim_steps
        for ielem=1:nelem
            no1=lnodal{ielem}(1);%1º nó do elemento ielem
            no2=lnodal{ielem}(2); %2º nó do elemento ielem
            dof1=2*lnodal{ielem}(1)-1; %1º dof do nó 1 do elemento ielem
            dof2=2*lnodal{ielem}(2)-1; %1º dof do nó 2 do elemento ielem

            LinesUD(ielem,1:4)=LinesU(ielem,1:4)+cicle_sin(anim_step)*...
                1/aF*[U(dof1) U(dof2) U(dof1+1) U(dof2+1)];

            set(def(ielem),'XData',LinesUD(ielem,1:2),'YData',LinesUD(ielem,3:4));
            pause(time_step);
            
        end
    end
end
end


