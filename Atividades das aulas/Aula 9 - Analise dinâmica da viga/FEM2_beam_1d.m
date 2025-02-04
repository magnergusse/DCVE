%%beam 1D - with mass formulation

% k=EI/L^3*[12 6L -12 6L;6L 4L^2 -6L 2L^2;-12 -6L 12 -6L;6L 2L^2 -6L 4L^2];
% m=rho*A*L/420*[156 22*L 54 -13*L; 22*L 4*L^2 13*L -3*L^2; ...
% 54 13*L 156 -22*L; -13*L -3*L^2 -22*L 4*L^2];
%
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
clear; clc;close all;
%% Examples

% analisisT='static';
  analisisT='dynamic';mode=3;


%% FEM model parameters

nelem=5; %n�mero de elementos de viga que constitui a viga
nnode=nelem+1; %n�mero de n�s
ndof=nnode*2; %considerando 2 graus de liberdade por n� (v, theta)

%% Par�metros materiais e geom�tricos
E=70e9;   % [Pa]     M�dulo de Young (alum�nio=70GPa, A�o=210GPa...)
rho=2700; % [kg/m^3] Massa vol�mica do material

Lt=0.5;  % [m] Comprimento da viga
h=0.003; % [m] espessura da viga
b=0.030; % [m] largura da viga
%------------------------------------------------------------------

I=b*h^3/12; % [m^4] Momento de in�rcia de sec��o retangular
A=b*h;      % [m^2] �rea da sec��o da viga

L=Lt/nelem; % [m]   comprimento de cada elemento de viga

%%------ defini��o da malha de elementos finitos --------------------

% localiza��es dos n�s e liga��es nodais
xyU=zeros(nnode,2); % matriz com as coordenadas dos n�s [x_n�1 y_n�1; 
%                                                        x_n�2 y_n�2;..]
lnodal={nelem};     % string com as liga��es nodais {[1 2], [2 3],...}

for ielem=1:nelem
    xyU(ielem,1)=L*(ielem-1);
    xyU(ielem+1,1)=L*ielem;   
    lnodal{ielem}=[ielem, ielem+1];
end

%% Condi��es de fronteira

Rdof=[1,2]; %FL - encastrada % graus de liberdade restritos: n�1_v, n�1_theta, n�_end_v
%Rdof=[];   %LL - livre-livre
%Rdof=[1 (nelem+1)*2-1]; % SS simplesmente apoiada

%% Carregamento est�tico

switch analisisT
    case 'static'
     % carregamento vertical de 1000N no n� 3 (corresponde ao dof 5)
     Fdof=[5]; %vetor de dofs carregados
     VF=[1000]; %Vetor de cargas
    case 'dynamic'
        %sem carregamento est�tico
end
%% C�lculo da Matriz de rigidez e massa

K=zeros(ndof,ndof); %inicializa��o da matriz de rigidez global
M=zeros(ndof,ndof); %inicializa��o da matriz de massa global

for ielem=1:nelem
% matriz de rigidez e de massa do elemento de viga (2dof por n�=4dof por elemento)
Kbeam= E*I/L^3;%  
Kelem=Kbeam*[12 6*L -12 6*L;6*L 4*L^2 -6*L 2*L^2;...
    -12 -6*L 12 -6*L;6*L 2*L^2 -6*L 4*L^2];%%
Mbeam=rho*A*L;
Melem=Mbeam/420*[156 22*L 54 -13*L; 22*L 4*L^2 13*L -3*L^2; 54 13*L 156 -22*L;...
    -13*L -3*L^2 -22*L 4*L^2];

% diagonaliza��o da matrix consistente para obter a massa diagonal 
% (m�todo HRZ - E.Hinton, T. Rock, O.C. Zienkiewicz.
% A note on mass lumping and related processes in the finite element method. 
% Earthquake Eng.&Struct Dyn, 4(3):245-249,1976.)

Melem=diag(diag(Melem)*420/312); %comentar se se pretender a matriz de
%massa consistente


%montagem da matriz de rigidez/massa do elemento na matriz de rigidez global
dof1=2*lnodal{ielem}(1)-1; %1� dof do 1� n� do elemento ielem
dof2=2*lnodal{ielem}(2)-1; %1� dof do 2� n� do elemento ielem

K(dof1:dof1+1,dof1:dof1+1)=K(dof1:dof1+1,dof1:dof1+1)+Kelem(1:2,1:2);
K(dof1:dof1+1,dof2:dof2+1)=K(dof1:dof1+1,dof2:dof2+1)+Kelem(1:2,3:4);
K(dof2:dof2+1,dof1:dof1+1)=K(dof2:dof2+1,dof1:dof1+1)+Kelem(3:4,1:2);
K(dof2:dof2+1,dof2:dof2+1)=K(dof2:dof2+1,dof2:dof2+1)+Kelem(3:4,3:4);

M(dof1:dof1+1,dof1:dof1+1)=M(dof1:dof1+1,dof1:dof1+1)+Melem(1:2,1:2);
M(dof1:dof1+1,dof2:dof2+1)=M(dof1:dof1+1,dof2:dof2+1)+Melem(1:2,3:4);
M(dof2:dof2+1,dof1:dof1+1)=M(dof2:dof2+1,dof1:dof1+1)+Melem(3:4,1:2);
M(dof2:dof2+1,dof2:dof2+1)=M(dof2:dof2+1,dof2:dof2+1)+Melem(3:4,3:4);
end

%% Carregamento e condi��es de fronteira

%Determina��o do vetor de carregamento
%F=zeros(ndof,1);
%F(Fdof)=VF;

%Aplica��o das condi��es de fronteira
free_dofs=setxor(1:ndof,Rdof); %graus de liberdade n�o restritos
Kp=K(free_dofs,free_dofs);
Mp=M(free_dofs,free_dofs);

switch analisisT
    case 'static'
      %Determina��o do vetor de carregamento
      F=zeros(ndof,1);
      F(Fdof)=VF;
      Fp=F(free_dofs,1);
    case 'dynamic'
        %nothing to do
end

%% Resolu��o do sistema de equa��es
U=zeros(ndof,1);

switch analisisT
    case 'static'
        Up=Kp\Fp;
        U(free_dofs)=Up;
    case 'dynamic'
        [a,b]=eig(Kp,Mp);
        freqs=sqrt(diag(b))/2/pi; %[Hz]
        U(free_dofs)=a(:,mode);
end

%% output de resultados

switch analisisT
    case 'static'
        % listagem de deslocamentos
        for inode=1:nnode
            v=U(inode*2-1);
            theta=U(inode*2);
            fprintf('d%iv = %7.3e\n', inode,v);
            fprintf('d%itheta = %7.3e\n', inode,theta);
        end
    case 'dynamic'
        for ifreq=1:length(freqs)
            fprintf('M%i = %12.6e [Hz] \n',ifreq,freqs(ifreq))
        end
end

  

aF =(max(U)-min(U))/(max(max(xyU))-min(min(xyU))); %factor de amplifica��o da deformada


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


max_anim_steps=100;time_step=0.01;
switch analisisT
    case 'static'
        for anim_step=1:max_anim_steps
            pause(time_step);
            for ielem=1:nelem
                no1=lnodal{ielem}(1); %1� n� do elemento ielem
                no2=lnodal{ielem}(2); %2� n� do elemento ielem
                dof1=2*lnodal{ielem}(1)-1; %1� dof do n� 1 do elemento ielem
                dof2=2*lnodal{ielem}(2)-1; %1� dof do n� 2 do elemento ielem
                LinesUD(ielem,1:4)=LinesU(ielem,1:4)+anim_step/max_anim_steps*...
                    1/aF*[0 0 U(dof1) U(dof2)];
                set(def(ielem),'XData',LinesUD(ielem,1:2),'YData',LinesUD(ielem,3:4));
            end
        end

    case 'dynamic'
        max_anim_steps=6;
        limY=get(gca,'ylim');MY=max(abs(limY));limY=[-MY,+MY];
        set(gca,'Ylim',limY,'Xlim',[0,Lt]);

        for anim_step=1:max_anim_steps
            sinus=sin(linspace(0,2*pi,180));
            for sin_step=1:length(sinus)
                pause(time_step);
                for ielem=1:nelem
                    no1=lnodal{ielem}(1); %1� n� do elemento ielem
                    no2=lnodal{ielem}(2); %2� n� do elemento ielem
                    dof1=2*lnodal{ielem}(1)-1; %1� dof do n� 1 do elemento ielem
                    dof2=2*lnodal{ielem}(2)-1; %1� dof do n� 2 do elemento ielem
                    LinesUD(ielem,1:4)=LinesU(ielem,1:4)+sinus(sin_step)*...
                        1/aF*[0 0 U(dof1) U(dof2)];
                    set(def(ielem),'XData',LinesUD(ielem,1:2),'YData',LinesUD(ielem,3:4));
                end
            end
        end
end



