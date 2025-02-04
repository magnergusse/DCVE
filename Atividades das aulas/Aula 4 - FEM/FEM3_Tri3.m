clearvars
close all

% -------------------------
% \\ Elemento Triangular //
% -------------------------

%% parametros de geometria e material

E=70e8;     % Módulo de Young do material da placa
nu = 0.3;   % Poisson
t = 0.1;    % espessura  

lnodal={[1 2 9] [2 10 9] [2 11 10] [2 3 11] [3 4 11] [4 12 11] [4 13 12]...
    [4 5 13] [5 6 13] [6 14 13] [6 15 14] [6 7 15] [7 8 15] [8 16 15] [8 9 16] [8 1 9]};
% pares de nós que definem cada elemento 

c=cos(pi()/4);
xyzNode = {[-2 -2], [0 -2],[2 -2], [2 0] [2 2 ] [0 2] [-2 2] [-2 0] ...
    [-c -c] [0 -1] [c -c] [1 0] [c c ] [0 1] [-c c] [-1 0]};    


% localização de cada nó (x,y,z)
[~ , nelem] = size(lnodal);         % numero de elementos
[~ , nnode]=size(xyzNode);        % numero de nós
ndof = nnode*2;
%% condições de fronteira e carregamentos
F=zeros(ndof,1); %inicialização do vetor de carregamento

fixed_dofs= [1:16];                 % Nós fixos (1 ao 8)---dofs 1->16
free_dofs=setxor(1:ndof,fixed_dofs);    % Todos os outros são livres

force_dofs=[17:ndof];                % Aplicação do carregamento
f=1000;
force_val=f*[-c -c 0 -1 c -c 1 0 c c 0 1 -c c -1 0 ];            % valor do carregamento

F(force_dofs)=force_val;  



%% geração da matriz de rigidez (montagem da matriz)
K=zeros(ndof,ndof); %inicialização da matriz de rigidez

% Matrizes constante
D = E/(1-nu^2)*[1 nu 0 ; nu 1 0 ; 0 0 (1-nu)/2];
for ielem=1:nelem
 
    x1 = xyzNode{lnodal{ielem}(1)}(1);
    y1 = xyzNode{lnodal{ielem}(1)}(2);
    
    x2 = xyzNode{lnodal{ielem}(2)}(1);
    y2 = xyzNode{lnodal{ielem}(2)}(2);
    
    x3 = xyzNode{lnodal{ielem}(3)}(1);
    y3 = xyzNode{lnodal{ielem}(3)}(2);
    
    A = 1/2*det([1 x1 y1 ; 1 x2 y2 ; 1 x3 y3]); % Área
    
    B = 1/(2*A)*[y2-y1 0 y3-y1 0 y1-y2 0
        0 x3-x2 0 x1-x3 0 x2-x1
        x3-x2 y2-y3 x1-x3 y3-y1 x2-x1 y1-y2];
    
    Kelem = B'*D*B*t;
    
    %montagem da matriz de rigidez do elemento na matriz de rigidez global
    no1=2*lnodal{ielem}(1)-1; %1º dof do 1º nó do elemento ielem
    no2=2*lnodal{ielem}(2)-1; %1º dof do 2º nó do elemento ielem
    no3=2*lnodal{ielem}(3)-1; %1º dof do 3º nó do elemento ielem
    
    indexB = [no1:no1+1,no2:no2+1,no3:no3+1];
    
    K(indexB,indexB) = K(indexB,indexB)+Kelem;
    
end

%% aplicação das condições de fronteira
Kp=K(free_dofs, free_dofs);
Fp=F(free_dofs,1);

%% Resolução do sistema de equações
Up=Kp\Fp;
U=zeros(ndof,1);
U(free_dofs)=Up;


%% Pós-Processamento


% listagem de deslocamentos
for inode=1:nnode
    dx(inode)=U(inode*2-1);
    dy(inode)=U(inode*2);
    nodes(inode)=inode; 
end
    VarNames = {'Node', 'Tx', 'Ty'};
    T = table(nodes',dx',dy','VariableNames',VarNames)
%representação da deformada

MU=max(max(abs(dx)),max(abs(dy)));

aF=0.2/MU; %fator de amplificação (que vai multiplicar pelo vetor de deslocamentos para adicionar ao vetor da geometria para representar a deformada
 

for ielem=1:nelem
    
    no1=lnodal{ielem}(1); %1º nó do elemento ielem
    no2=lnodal{ielem}(2); %2º nó do elemento ielem
    no3=lnodal{ielem}(3); %3º nó do elemento ielem
    
    x1 = xyzNode{no1}(1);
    y1 = xyzNode{no1}(2);
    
    x2 = xyzNode{no2}(1);
    y2 = xyzNode{no2}(2);
    
    x3 = xyzNode{no3}(1);
    y3 = xyzNode{no3}(2);
    
     
    
    dof1=no1*2-1; %1º dof do nó 1 do elemento ielem
    dof2=no2*2-1; %1º dof do nó 2 do elemento ielem
    dof3=no3*2-1; %1º dof do nó 3 do elemento ielem
  
    LinesU(ielem,1:8)=[x1 x2 x3 x1 y1 y2 y3 y1];
    LinesUD(ielem,1:8)=LinesU(ielem,1:8)+aF*[U(dof1) U(dof2) U(dof3)...
        U(dof1) U(dof1+1) U(dof2+1) U(dof3+1) U(dof1+1)];
    %stress_color(ielem,1:3)=cm(floor((stress(ielem,1)-min(stress))*63/deltaS)+1,1:3);
end

 

figure(1); hold on; 
axis([min(min(LinesUD(:,1:3))),max(max(LinesUD(:,1:3))),...
    min(min(LinesUD(:,5:7))),max(max(LinesUD(:,5:7)))])
axis equal
for i=1:size(LinesU,1)
    undef=line(LinesU(i,1:4),LinesU(i,5:8)); set(undef,'LineStyle',':','color','b','Linewidth',.5);
    def(i)=line(LinesU(i,1:4),LinesU(i,5:8)); set(def(i),'LineStyle','-','color','b','Linewidth',2);
end%
%cbh=colorbar;
%cbh.Ticks=linspace(0,1,10);
%cbh.TickLabels=num2cell(linspace(min(stress),max(stress),10));
hold off

 

% animação


max_anim_steps=100;time_step=0.01;
for anim_step=1:max_anim_steps
    pause(time_step);
    for ielem=1:nelem
            no1=lnodal{ielem}(1); %1º nó do elemento ielem
            no2=lnodal{ielem}(2); %2º nó do elemento ielem
            no3=lnodal{ielem}(3); %3º nó do elemento ielem
          
            dof1=no1*2-1; %1º dof do nó 1 do elemento ielem
            dof2=no2*2-1; %1º dof do nó 2 do elemento ielem
            dof3=no3*2-1; %1º dof do nó 3 do elemento ielem
         
            
        LinesUD(ielem,1:8)=LinesU(ielem,1:8)+anim_step/max_anim_steps*...
        aF*[U(dof1) U(dof2) U(dof3) U(dof1) U(dof1+1) U(dof2+1) U(dof3+1)  U(dof1+1)];

        set(def(ielem),'XData',LinesUD(ielem,1:4),'YData',LinesUD(ielem,5:8));
        
              
        
    end
    

end


%% new - colormap

nf=figure(2);

 hold on; 
 axis([min(min(LinesUD(:,1:3))) , max(max(LinesUD(:,1:3))), min(min(LinesUD(:,5:7))) , max(max(LinesUD(:,5:7)))])
axis equal
for i=1:size(LinesU,1)
    undef=line(LinesU(i,1:4),LinesU(i,5:8)); set(undef,'LineStyle',':','color','b','Linewidth',.5);
    def(i)=line(LinesU(i,1:4),LinesU(i,5:8)); set(def,'LineStyle','-','color','b','Linewidth',2);
end%
for ielem=1:nelem
            no1=lnodal{ielem}(1); %1º nó do elemento ielem
            no2=lnodal{ielem}(2); %2º nó do elemento ielem
            no3=lnodal{ielem}(3); %3º nó do elemento ielem
           

            dof1=no1*2-1; %1º dof do nó 1 do elemento ielem
            dof2=no2*2-1; %1º dof do nó 2 do elemento ielem
            dof3=no3*2-1; %1º dof do nó 3 do elemento ielem
           
            
        LinesUD(ielem,1:8)=LinesU(ielem,1:8)+...
        aF*[U(dof1) U(dof2) U(dof3) U(dof1) U(dof1+1) U(dof2+1) U(dof3+1)  U(dof1+1)];

        set(def(ielem),'XData',LinesUD(ielem,1:4),'YData',LinesUD(ielem,5:8));
        
        data=[(abs(U(dof1))+abs(U(dof1+1)))/2, (abs(U(dof2))+abs(U(dof2+1)))/2,...
            (abs(U(dof3))+abs(U(dof3+1)))/2, ...
            (abs(U(dof1))+abs(U(dof1+1)))/2];
        
        fill(LinesUD(ielem,1:4),LinesUD(ielem,5:8),data);
        colormap(jet)
        
    end