%% Beam 1D

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
%% Example

% Beam with 4 elements fixed at one end and 
% simply supported at the other end 2
% with vertical load at midspan

%% FEM model parameters
clear; clc;close all;
nelem=18; %número de elementos de viga que constitui a viga
nnode=nelem+1; % Número de nós
ndof=nnode*2; % Considerando 2 graus de liberdade por nó (v, theta)

%% Parâmetros materiais e geométricos
E=70e9; %Pa (material)
Es=E*ones(nelem,1); % Vetor de módulo de young das vigas

larg=30e-3;
h=3e-3;

I=larg*h^3/12; % Momento de inércia de secção retangular (40x80mm)
Is=I*ones(nelem,1); %vetor de momentos de inércia para as vigas

L=500e-3; % Comprimento total da viga (0.5 ou 0.4 metros)
Lt=L/nelem;
Ls=L*ones(nelem,1);

rho=2700;
rhos=rho*ones(nelem,1);
Area=larg*h;
Areas=Area*ones(nelem,1);

% Localizações dos nós e ligações nodais
xyU=zeros(nnode,2);
lnodal={nelem};
for ielem=1:nelem
    xyU(ielem,1)=Lt*(ielem-1);
    xyU(ielem+1,1)=Lt*ielem;
    lnodal{ielem}=[ielem, ielem+1];
end

%% Condições de fronteira
% graus de liberdade restritos: nó1_v, nó1_theta, nó_end_v

%Rdof=[1,2]; %FL - encastrada % graus de liberdade restritos: nó1_v, nó1_theta, nó_end_v
Rdof=[];   %LL - livre-livre
% Rdof=[1 (nelem+1)*2-1]; % SS simplesmente apoiada

%% Cálculo da Matriz de rigidez
%alpha=1/50;

K=zeros(ndof,ndof); %inicialização da matriz de rigidez global
M=zeros(ndof,ndof);

for ielem=1:nelem
% matriz de rigidez do elemento de viga (2dof por nó=4dof por elemento)
Kbeam= E*I/Lt^3;

Kelem=Kbeam*[12 6*L -12 6*L;6*L 4*L^2 -6*L 2*L^2;-12 -6*L 12 -6*L;6*L...
    2*L^2 -6*L 4*L^2];

%montagem da matriz de rigidez do elemento na matriz de rigidez global
dof1=2*lnodal{ielem}(1)-1; %1º dof do 1º nó do elemento ielem
dof2=2*lnodal{ielem}(2)-1; %1º dof do 2º nó do elemento ielem

K(dof1:dof1+1,dof1:dof1+1)=K(dof1:dof1+1,dof1:dof1+1)+Kelem(1:2,1:2);
K(dof1:dof1+1,dof2:dof2+1)=K(dof1:dof1+1,dof2:dof2+1)+Kelem(1:2,3:4);
K(dof2:dof2+1,dof1:dof1+1)=K(dof2:dof2+1,dof1:dof1+1)+Kelem(3:4,1:2);
K(dof2:dof2+1,dof2:dof2+1)=K(dof2:dof2+1,dof2:dof2+1)+Kelem(3:4,3:4);

Mbeam= rho*Area*Lt/420;


Melem=Mbeam*[156 22*L 54 -13*L; 22*L 4*L^2 13*L -3*L^2; 54 13*L 156 -22*L; -13*L -3*L^2 -22*L 4*L^2];


M(dof1:dof1+1,dof1:dof1+1)=M(dof1:dof1+1,dof1:dof1+1)+Melem(1:2,1:2);
M(dof1:dof1+1,dof2:dof2+1)=M(dof1:dof1+1,dof2:dof2+1)+Melem(1:2,3:4);
M(dof2:dof2+1,dof1:dof1+1)=M(dof2:dof2+1,dof1:dof1+1)+Melem(3:4,1:2);
M(dof2:dof2+1,dof2:dof2+1)=M(dof2:dof2+1,dof2:dof2+1)+Melem(3:4,3:4);
end

%% Carregamento e condições de fronteira
%Determinação do vetor de carregamento
% F=zeros(ndof,1);
% F(Fdof)=VF;

%Aplicação das condições de fronteira
free_dofs=setxor(1:ndof,Rdof); %graus de liberdade não restritos
Kp=K(free_dofs,free_dofs);
% Fp=F(free_dofs,1);
Mp=M(free_dofs,free_dofs);
%% Resolução do sistema de equações
U=zeros(ndof,1);
[a,b]=eig(Kp, Mp);
freqs=sqrt(diag(b))/2/pi; %[Hz]

for ifreq=1:10
    fprintf('F%i = %12.6e [Hz] \n',ifreq,freqs(ifreq));
end


%% output de resultados: Movimento até parar

%representação da deformada
% U(free_dofs)=a(:,6);
% aF =(max(U)-min(U))/(max(max(xyU))-min(min(xyU))); %factor de amplificação da deformada
% 
% for ielem=1:nelem
%    no1=lnodal{ielem}(1);%1º nó do elemento ielem
%    no2=lnodal{ielem}(2); %2º nó do elemento ielem
%    dof1=2*lnodal{ielem}(1)-1; %1º dof do nó 1 do elemento ielem
%    dof2=2*lnodal{ielem}(2)-1; %1º dof do nó 2 do elemento ielem 
%    LinesU(ielem,1:4)=[xyU(no1,1) xyU(no2,1) xyU(no1,2)  xyU(no2,2)];
%    LinesUD(ielem,1:4)=LinesU(ielem,1:4)+1/aF*[0 0 U(dof1) U(dof2)];
% 
% end
% grid on
% figure(1); hold on; axis equal
% for i=1:size(LinesU,1)
%     undef(i)=line(LinesU(i,1:2),LinesU(i,3:4)); set(undef,'LineStyle',':','color','b','Linewidth',.5);
%     def(i)=line(LinesUD(i,1:2),LinesUD(i,3:4)); set(def,'LineStyle','-','color','b','Linewidth',2);
% end
% hold off
% 
% % animação
% 
% xlabel('Secções da viga, em comprimento, em x (m)');
% ylabel('Amplitude em y (m)')
% max_anim_steps=100;time_step=0.01;
% for anim_step=1:max_anim_steps
%     pause(time_step);
%     for ielem=1:nelem
%         no1=lnodal{ielem}(1);%1º nó do elemento ielem
%         no2=lnodal{ielem}(2); %2º nó do elemento ielem
%         dof1=2*lnodal{ielem}(1)-1; %1º dof do nó 1 do elemento ielem
%         dof2=2*lnodal{ielem}(2)-1; %1º dof do nó 2 do elemento ielem
%         LinesUD(ielem,1:4)=LinesU(ielem,1:4)+anim_step/max_anim_steps*...
%             1/aF*[0 0 U(dof1) U(dof2)];
%         set(def(ielem),'XData',LinesUD(ielem,1:2),'YData',LinesUD(ielem,3:4));
%     end

% end



%% Animacao Inf

U(free_dofs)=a(:,3);
aF =(max(U)-min(U))/(max(max(xyU))-min(min(xyU))); %factor de amplificação da deformada


for ielem=1:nelem
   no1=lnodal{ielem}(1);%1º nó do elemento ielem
   no2=lnodal{ielem}(2); %2º nó do elemento ielem
   dof1=2*lnodal{ielem}(1)-1; %1º dof do nó 1 do elemento ielem
   dof2=2*lnodal{ielem}(2)-1; %1º dof do nó 2 do elemento ielem 
   LinesU(ielem,1:4)=[xyU(no1,1) xyU(no2,1) xyU(no1,2)  xyU(no2,2)];
   LinesUD(ielem,1:4)=LinesU(ielem,1:4)+1/aF*[0 0 U(dof1) U(dof2)];

end

figure(1); hold on; axis equal
for i=1:size(LinesU,1);
    undef(i)=line(LinesU(i,1:2),LinesU(i,3:4)); set(undef,'LineStyle',':','color','b','Linewidth',.5);
    def(i)=line(LinesUD(i,1:2),LinesUD(i,3:4)); set(def,'LineStyle','-','color','b','Linewidth',2);
end
hold off

time_step=0.01;
max_anim_steps=6;
limY=get(gca,'ylim');MY=max(abs(limY));limY=[-MY,+MY];
set(gca,'Ylim',limY,'Xlim',[0,L]);

for anim_step=1:max_anim_steps
    sinus=sin(linspace(0,2*pi,180));
    for sin_step=1:length(sinus)
        pause(time_step);
        for ielem=1:nelem
            no1=lnodal{ielem}(1); %1º nó do elemento ielem
            no2=lnodal{ielem}(2); %2º nó do elemento ielem
            dof1=2*lnodal{ielem}(1)-1; %1º dof do nó 1 do elemento ielem
            dof2=2*lnodal{ielem}(2)-1; %1º dof do nó 2 do elemento ielem
            LinesUD(ielem,1:4)=LinesU(ielem,1:4)+sinus(sin_step)*...
            1/aF*[0 0 U(dof1) U(dof2)];
            set(def(ielem),'XData',LinesUD(ielem,1:2),'YData',LinesUD(ielem,3:4));
        end
   end
end