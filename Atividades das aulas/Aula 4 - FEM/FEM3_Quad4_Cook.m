clearvars
close all

% ---------------------------
% \\ Elemento Quadrilátero //
% ---------------------------

%% parametros de geometria e material

E=1;        % Módulo de Young do material da placa
nu = 0.333; % Poisson
t = 1;      % espessura da placa 
% 
%% exemplo simples com 4 elementos
% xyzNode = {[0 0],[24 22],[48 44],[0 22], [24 37], [48 52], [0 44], [24 52], [48 60]};
% lnodal ={[1 2 5 4], [2 3 6 5], [4 5 8 7], [5 6 9 8]};

%% exemplo genérico com nxm elementos

n=10; % #elementos na vertical
m=8; % #elementos na horizontal

% geração de nós
k=0;
for i=0:n
    for j=0:m
        k=k+1;
        xyzNode{k}=[48/m*j,44/n*i+44/m*j-28/(m*n)*i*j];
    end
end
% geração dos elementos
k=0;
for i=1:n
    for j=1:m
        k=k+1;
        lnodal{k}=[(i-1)*(m+1)+j,(i-1)*(m+1)+j+1,i*(m+1)+j+1,i*(m+1)+j ];
    end
end

[~ , nelem] = size(lnodal);         % numero de elementos
[~ , n_nodes]=size(xyzNode);        % numero de nós

ndof = n_nodes*2;
%% condições de fronteira e carregamentos
F=zeros(ndof,1); %inicialização do vetor de carregamento

% Fixar os nós em x=0
fixed_dofs=[];
for i=1:n_nodes
    if xyzNode{i}(1)==0
        fixed_dofs=[fixed_dofs,2*i-1:2*i];
    end
end


free_dofs=setxor(1:ndof,fixed_dofs);    % Todos os outros são livres

force_dofs=[ndof/2];           % Aplicação do carregamento
force_val=[-1];            % valor do carregamento

F(force_dofs)=force_val;  


%% geração da matriz de rigidez (montagem da matriz)
K=zeros(ndof,ndof); %inicialização da matriz de rigidez

% Matriz constante
D = E/(1-nu^2)*[1 nu 0 ; nu 1 0 ; 0 0 (1-nu)/2];
% quadratura
qd=1/sqrt(3);
Q = [-qd -qd ; qd -qd ; qd qd ; -qd qd];
W = [1;1;1;1];

for ielem=1:nelem
 
    x1 = xyzNode{lnodal{ielem}(1)}(1);
    y1 = xyzNode{lnodal{ielem}(1)}(2);
    
    x2 = xyzNode{lnodal{ielem}(2)}(1);
    y2 = xyzNode{lnodal{ielem}(2)}(2);
    
    x3 = xyzNode{lnodal{ielem}(3)}(1);
    y3 = xyzNode{lnodal{ielem}(3)}(2);

    x4 = xyzNode{lnodal{ielem}(4)}(1);
    y4 = xyzNode{lnodal{ielem}(4)}(2);


    
    Kelem=zeros(8,8);
    for q=1:size(W,1)
        xi = Q(q,1);
        eta= Q(q,2);
       
       % Funções de forma
        N = 1/4*[ (1-xi)*(1-eta);(1+xi)*(1-eta);
           (1+xi)*(1+eta);(1-xi)*(1+eta)];
       
        dNdxi=1/4*[-(1-eta) , -(1-xi)
                     1-eta  , -(1+xi)
                     1+eta  ,   1+xi
                   -(1+eta) ,   1-xi];
       
        J = [x1 x2 x3 x4 ; y1 y2 y3 y4 ]*dNdxi;
        J = J';
        dNdx = dNdxi/J;
        
        B=zeros(3,8);
        B(1, [1 3 5 7])= dNdx(:,1)';
        B(2, [2 4 6 8])= dNdx(:,2)';
        B(3, [1 3 5 7])= dNdx(:,2)'; B(3, [2 4 6 8])= dNdx(:,1)';
        
        Kelem = Kelem+B'*D*B*W(q)*det(J)*t;
    end
        %montagem da matriz de rigidez do elemento na matriz de rigidez global
        no1=2*lnodal{ielem}(1)-1; %1º dof do 1º nó do elemento ielem
        no2=2*lnodal{ielem}(2)-1; %1º dof do 2º nó do elemento ielem
        no3=2*lnodal{ielem}(3)-1; %1º dof do 3º nó do elemento ielem
        no4=2*lnodal{ielem}(4)-1; %1º dof do 4º nó do elemento ielem

        indexB = [no1:no1+1,no2:no2+1,no3:no3+1,no4:no4+1];

        K(indexB,indexB) = K(indexB,indexB)+Kelem;
 
end

%% aplicação das condições de fronteira
Kp=K(free_dofs, free_dofs);
Fp=F(free_dofs,1);

%% Resolução do sistema de equações
Up=Kp\Fp;
U=zeros(ndof,1);
U(free_dofs)=Up;
U_C = U(ndof);

%% Pós-Processamento

%representação da deformada

aF=10; %fator de amplificação (que vai multiplicar pelo vetor de deslocamentos para adicionar ao vetor da geometria para representar a deformada
 

for ielem=1:nelem
    
    x1 = xyzNode{lnodal{ielem}(1)}(1);
    y1 = xyzNode{lnodal{ielem}(1)}(2);
    
    x2 = xyzNode{lnodal{ielem}(2)}(1);
    y2 = xyzNode{lnodal{ielem}(2)}(2);
    
    x3 = xyzNode{lnodal{ielem}(3)}(1);
    y3 = xyzNode{lnodal{ielem}(3)}(2);
    
    x4 = xyzNode{lnodal{ielem}(4)}(1);
    y4 = xyzNode{lnodal{ielem}(4)}(2);
    
    no1=lnodal{ielem}(1); %1º nó do elemento ielem
    no2=lnodal{ielem}(2); %2º nó do elemento ielem
    no3=lnodal{ielem}(3); %3º nó do elemento ielem
    no4=lnodal{ielem}(4); %4º nó do elemento ielem
    
    dof1=no1*2-1; %1º dof do nó 1 do elemento ielem
    dof2=no2*2-1; %1º dof do nó 2 do elemento ielem
    dof3=no3*2-1; %1º dof do nó 3 do elemento ielem
    dof4=no4*2-1; %1º dof do nó 4 do elemento ielem
    
    
    LinesU(ielem,1:10)=[x1 x2 x3 x4 x1 y1 y2 y3 y4 y1];
    LinesUD(ielem,1:10)=LinesU(ielem,1:10)+1/aF*[U(dof1) U(dof2) U(dof3) U(dof4) U(dof1) U(dof1+1) U(dof2+1) U(dof3+1) U(dof4+1) U(dof1+1)];
    %stress_color(ielem,1:3)=cm(floor((stress(ielem,1)-min(stress))*63/deltaS)+1,1:3);
end

 

figure(1); hold on; 
axis equal;
boundsD=[min(min(LinesUD(:,1:4))) , max(max(LinesUD(:,1:4))), min(min(LinesUD(:,6:9))) , max(max(LinesUD(:,6:9)))];
boundsU=[min(min(LinesU(:,1:4))) , max(max(LinesU(:,1:4))), min(min(LinesU(:,6:9))) , max(max(LinesU(:,6:9)))];
bounds=[min(boundsD(:,1),boundsU(:,1)),max(boundsD(:,2),boundsU(:,2)), min(boundsD(:,3),boundsU(:,3)),max(boundsD(:,4),boundsU(:,4))];
axis(bounds);
for i=1:size(LinesU,1)
    undef=line(LinesU(i,1:5),LinesU(i,6:10)); set(undef,'LineStyle',':','color','b','Linewidth',.5);
    def(i)=line(LinesU(i,1:5),LinesU(i,6:10)); set(def,'LineStyle','-','color','b','Linewidth',2);
end%

hold off

 

% animação

%animdef=[NaN,NaN];
max_anim_steps=100;time_step=0.01;
for anim_step=1:max_anim_steps
    pause(time_step);
    for ielem=1:nelem
            no1=lnodal{ielem}(1); %1º nó do elemento ielem
            no2=lnodal{ielem}(2); %2º nó do elemento ielem
            no3=lnodal{ielem}(3); %3º nó do elemento ielem
            no4=lnodal{ielem}(4); %4º nó do elemento ielem

            dof1=no1*2-1; %1º dof do nó 1 do elemento ielem
            dof2=no2*2-1; %1º dof do nó 2 do elemento ielem
            dof3=no3*2-1; %1º dof do nó 3 do elemento ielem
            dof4=no4*2-1; %1º dof do nó 4 do elemento ielem
            
        LinesUD(ielem,1:10)=LinesU(ielem,1:10)+anim_step/max_anim_steps*...
        1/aF*[U(dof1) U(dof2) U(dof3) U(dof4) U(dof1) U(dof1+1) U(dof2+1) U(dof3+1) U(dof4+1) U(dof1+1)];

        set(def(ielem),'XData',LinesUD(ielem,1:5),'YData',LinesUD(ielem,6:10));
    end
    % set(def(:),'XData',LinesUD(:,1:2),'YData',LinesUD(:,3:4));

end


%% new - colormap

nf=figure(2);

 hold on; axis(bounds);
axis equal
for i=1:size(LinesU,1)
    undef=line(LinesU(i,1:5),LinesU(i,6:10)); set(undef,'LineStyle',':','color','b','Linewidth',.5);
    def(i)=line(LinesU(i,1:5),LinesU(i,6:10)); set(def,'LineStyle','-','color','b','Linewidth',2);
end%
for ielem=1:nelem
            no1=lnodal{ielem}(1); %1º nó do elemento ielem
            no2=lnodal{ielem}(2); %2º nó do elemento ielem
            no3=lnodal{ielem}(3); %3º nó do elemento ielem
            no4=lnodal{ielem}(4); %4º nó do elemento ielem

            dof1=no1*2-1; %1º dof do nó 1 do elemento ielem
            dof2=no2*2-1; %1º dof do nó 2 do elemento ielem
            dof3=no3*2-1; %1º dof do nó 3 do elemento ielem
            dof4=no4*2-1; %1º dof do nó 4 do elemento ielem
            
        LinesUD(ielem,1:10)=LinesU(ielem,1:10)+...
        1/aF*[U(dof1) U(dof2) U(dof3) U(dof4) U(dof1) U(dof1+1) U(dof2+1) U(dof3+1) U(dof4+1) U(dof1+1)];

        set(def(ielem),'XData',LinesUD(ielem,1:5),'YData',LinesUD(ielem,6:10));
        
        data=[(abs(U(dof1))+abs(U(dof1+1)))/2, (abs(U(dof2))+abs(U(dof2+1)))/2,...
            (abs(U(dof3))+abs(U(dof3+1)))/2, (abs(U(dof4))+abs(U(dof4+1)))/2,...
            (abs(U(dof1))+abs(U(dof1+1)))/2];
        
        fill(LinesUD(ielem,1:5),LinesUD(ielem,6:10),data);
        colormap(jet)
        
    end