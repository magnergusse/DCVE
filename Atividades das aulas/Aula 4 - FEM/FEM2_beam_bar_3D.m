clearvars
close all


%% modelo de vigas simples
%xyzNode={[0 4 0], [3 4 0],[0 4 -3],[0 0 0]};
%lnodal={[1 2],[4 1],[1 3]};
%% portico cubo
xyzNode={[0 0 0],[0 1 0],[0 2 0],[0 3 0],[0 4 0],[0 0 4],[0 1 4], [0 2 4], ...
    [0 3 4],[0 4 4],[4 0 4],[4 1 4],[4 2 4],[4 3 4], [4 4 4],...
    [4 0 0],[4 1 0],[4 2 0], [4 3 0],[4 4 0],[0 4 1],[0 4 2],[0 4 3],...
    [1 4 4],[2 4 4],[3 4 4],[4 4 3],[4 4 2],[4 4 1],[3 4 0],[2 4 0],[1 4 0]};
lnodal={[1 2],[2 3],[3 4],[4 5],[6 7],[7 8],[8 9],[9 10],...
    [11 12],[12 13],[13 14],[14 15],[16 17],[17 18],[18 19],[19 20],...
    [5 21],[21 22],[22 23],[23 10],[10 24],[24 25],[25 26],[26 15],...
    [15 27],[27 28],[28 29],[29 20],[20 30],[30 31],[31 32],[32 5]};

[~,nnode]=size(xyzNode);
[~,nelem]=size(lnodal);
ndof=nnode*6;
%% parametros de geometria e material
%Lel=1;                   % comprimento de cada elemento
%b = 0.01;                % Largura secção da viga
%h = 0.01;                % Altura secção da viga
r = 0.1;                  % raio da secção da viga
A = 2e-2;                 % Área da secção da viga
E = 210e9;                % Módulo de Young do material da viga
Iz = 20e-5;               % Momento de Inércia da área da viga em z
Iy =10e-5;                % Momento de Inércia da área da viga em y
J = 5e-5;                 % Momento Polar de Inércia da seccão da viga
G = 84e9;                 % Shear Modulus
%% condições de fronteira e carregamentos
F=zeros(ndof,1);                       % inicialização do vetor de carregamento
 
fixed_dofs= [1 2 3 4 5 6 31 32 33 34 35 36 61 62 63 64 65 66 91 92 93 94 95 96]; %  Nós fixos 1 6 11 16

free_dofs=setxor(1:ndof,fixed_dofs);   % Todos os outros são livres

force_dofs=[121 127 133];      % Aplicação dos carregamentos nós 21 22 e 23 Fx
force_val=[10e3 10e3 10e3];               % Valor dos carregamentos

F(force_dofs)=force_val;  


%% geração da matriz de rigidez (montagem da matriz)
K=zeros(ndof,ndof); %inicialização da matriz de rigidez

for ielem=1:nelem
    
    % coordenadas do no 1
    x1 = xyzNode{lnodal{ielem}(1)}(1);
    y1 = xyzNode{lnodal{ielem}(1)}(2);
    z1 = xyzNode{lnodal{ielem}(1)}(3);
    
    % coordenadas do no 2
    x2 = xyzNode{lnodal{ielem}(2)}(1);
    y2 = xyzNode{lnodal{ielem}(2)}(2);
    z2 = xyzNode{lnodal{ielem}(2)}(3);
    
    % Ponto 3, para definir a torsão ao longo da barra
    x3 = (x2-x1)/2 + x1 + 1;
    y3 = (y2-y1)/2 + y1 + 1;    
    z3 = (z2-z1)/2 + z1;
    
    V1 = [x1 y1 z1];
    V2 = [x2 y2 z2];
    V3 = [x3 y3 z3];
    
    % comprimento da barra
    L = sqrt((x2-x1)^2+(y2-y1)^2+(z2-z1)^2);
    
    % cossenos diretores
    lx = (x2-x1)/L;
    mx = (y2-y1)/L;
    nx = (z2-z1)/L;
      
    aux = cross((V2-V1),(V3-V1))/norm(cross((V2-V1),(V3-V1)));
    lz = aux(1);
    mz = aux(2);
    nz = aux(3);
    
    aux = cross(aux,[lx,mx,nx]);
    ly = aux(1);
    my = aux(2);
    ny = aux(3);
    
    % matriz de transformação
    T_aux = [lx mx nx
        ly my ny
        lz mz nz]
    
    T = [T_aux zeros(3) zeros(3) zeros(3)
        zeros(3) T_aux zeros(3) zeros(3)
        zeros(3) zeros(3) T_aux zeros(3)
        zeros(3) zeros(3) zeros(3) T_aux];

    % Matriz de Rigidez ( sem transformação de variáveis )
    S1 = [E*A/L 0 0 0 0 0
        0 12*E*Iz/L^3 0 0 0 6*E*Iz/L^2
        0 0 12*E*Iy/L^3 0 -6*E*Iy/L^2 0
        0 0 0 G*J/L 0 0
        0 0 -6*E*Iy/L^2 0 4*E*Iy/L 0
        0 6*E*Iz/L^2 0 0 0 4*E*Iz/L];
    
    S2 = [-E*A/L 0 0 0 0 0
        0 -12*E*Iz/L^3 0 0 0 6*E*Iz/L^2
        0 0 -12*E*Iy/L^3 0 -6*E*Iy/L^2 0
        0 0 0 -G*J/L 0 0
        0 0 6*E*Iy/L^2 0 2*E*Iy/L 0
        0 -6*E*Iz/L^2 0 0 0 2*E*Iz/L];
    
    S3 = [-E*A/L 0 0 0 0 0
        0 -12*E*Iz/L^3 0 0 0 -6*E*Iz/L^2
        0 0 -12*E*Iy/L^3 0 6*E*Iy/L^2 0
        0 0 0 -G*J/L 0 0
        0 0 -6*E*Iy/L^2 0 2*E*Iy/L 0
        0 6*E*Iz/L^2 0 0 0 2*E*Iz/L];
    
    S4 = [E*A/L 0 0 0 0 0
        0 12*E*Iz/L^3 0 0 0 -6*E*Iz/L^2
        0 0 12*E*Iy/L^3 0 6*E*Iy/L^2 0
        0 0 0 G*J/L 0 0
        0 0 6*E*Iy/L^2 0 4*E*Iy/L 0
        0 -6*E*Iz/L^2 0 0 0 4*E*Iz/L];
    
    Kelem = [S1 S2 ; S3 S4];    % Antes da transformação de variáveis
    Kelem = T'*Kelem*T;         % Após a transformação de variáveis
    
    % Montagem da matriz de rigidez do elemento na matriz de rigidez global
    dof1=6*lnodal{ielem}(1)-5; %1º dof do 1º nó do elemento ielem
    dof2=6*lnodal{ielem}(2)-5; %1º dof do 2º nó do elemento ielem
    
    indexB = [dof1:dof1+5,dof2:dof2+5];
    
    K(indexB,indexB) = K(indexB,indexB)+Kelem;
    
end

%% aplicação das condições de fronteira
Kp=K(free_dofs, free_dofs);
Fp=F(free_dofs,1);

%% Resolução do sistema de equações
Up=Kp\Fp;
U=zeros(ndof,1);
U(free_dofs)=Up;

%% output de resultados

% listagem de deslocamentos
for inode=1:nnode
    dx(inode)=U(inode*6-5);
    dy(inode)=U(inode*6-4);
    dz(inode)=U(inode*6-3);
    Rx(inode)=U(inode*6-2);
    Ry(inode)=U(inode*6-1);
    Rz(inode)=U(inode*6);
    nodes(inode)=inode; 
end
    VarNames = {'Node', 'Tx', 'Ty', 'Tz', 'Rx', 'Ry', 'Rz'};
    T = table(nodes',dx',dy',dz',Rx',Ry',Rz','VariableNames',VarNames)
    
% listagem de forças


% listagem de tensões normais

%here

%representação da deformada

aF=0.002;


for ielem=1:nelem
    
   no1=lnodal{ielem}(1); % 1º nó do elemento ielem
   no2=lnodal{ielem}(2); % 2º nó do elemento ielem
   dof1=no1*6-5; %1º dof do nó 1 do elemento ielem
   dof2=no2*6-5; %1º dof do nó 2 do elemento ielem 
   
   x1=xyzNode{no1}(1); y1=xyzNode{no1}(2); z1=xyzNode{no1}(3);
   x2=xyzNode{no2}(1); y2=xyzNode{no2}(2); z2=xyzNode{no2}(3);
 
   LinesU(ielem,1:6)=[x1 x2 y1 y2 z1 z2];
      
   LinesUD(ielem,1:6)=LinesU(ielem,1:6)+1/aF*[U(dof1) U(dof2) U(dof1+1) U(dof2+1) U(dof1+2) U(dof2+2)];
   %stress_color(ielem,1:3)=cm(floor((stress(ielem,1)-min(stress))*63/deltaS)+1,1:3);
   
end

figure(1); clf; hold on; axis equal

xlabel('X');ylabel('Y');zlabel('Z');

for i=1:size(LinesU,1);
    undef=line(LinesU(i,1:2),LinesU(i,3:4),LinesU(i,5:6)); set(undef,'LineStyle',':','color','b','Linewidth',.5);
    def(i)=line(LinesU(i,1:2),LinesU(i,3:4),LinesU(i,5:6)); set(def(i),'LineStyle','-','color','b','Linewidth',2);
end%

hold off

view(3); grid on
view(50,30);camup([0 1 0])

% animação
grid off;
%xbound=get(gca,'xlim');ybound=get(gca,'ylim');
set(gca,'xlim',[-.5,5.5],'YLim',[-.5 4.5],'ZLim',[-.5, 4.5]);


%animdef=[NaN,NaN];
max_anim_steps=100;time_step=0.03;
for anim_step=1:max_anim_steps
    pause(time_step);
    for ielem=1:nelem
        no1=lnodal{ielem}(1);%1º nó do elemento ielem
        no2=lnodal{ielem}(2); %2º nó do elemento ielem
        dof1=no1*6-5; %1º dof do nó 1 do elemento ielem
        dof2=no2*6-5; %1º dof do nó 2 do elemento ielem 
      
        LinesUD(ielem,1:6)=LinesU(ielem,1:6)+anim_step/max_anim_steps*...
            1/aF*[U(dof1) U(dof2) U(dof1+1) U(dof2+1) U(dof1+2) U(dof2+2)];
       
        set(def(ielem),'XData',LinesUD(ielem,1:2),'YData',LinesUD(ielem,3:4),...
            'ZData',LinesUD(ielem,5:6));
    end

    
end

