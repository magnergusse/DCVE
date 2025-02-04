%% exemplo elemento de barra Bar2dof
clearvars;
close all;

%% parametros de geometria e material
L=10;            % Comprimento da barra [m]
A=0.02*0.02;    % �rea da sec��o da barra [m x m]
E=210e9;        % M�dulo de Young do material da barra [Pa]
rho=7800;       % massa vol�mica do material da barra [kg/m3]

%% inicializa��o de vari�veis
nelem=20;       % n�mero de elementos
nnode=nelem+1;  % n�mero de n�s;
ndof=nnode;     % n�mero de graus de liberdade

Lel=L/nelem;    %comprimento de cada elemento de barra

%% condi��es de fronteira e carregamentos

fixed_dofs= [1 ]; %identifica��o dos gdl fixos
free_dofs=setxor(1:nnode,fixed_dofs); %


%% matriz de rigidez e de massa do elemento de barra

Kel=A*E/Lel*[1 -1; -1 1];
Mel=A*rho*L*[0.5 0;0 0.5];

%% gera��o da matriz de rigidez e de massa (montagem da matriz)
K=zeros(ndof,ndof); %inicializa��o da matriz de rigidez
M=zeros(ndof,ndof); %inicializa��o da matriz de massa

for ielem=1:nelem
    range=[ielem ielem+1];
    K(range,range)=K(range,range)+Kel;
    M(range,range)=M(range,range)+Mel;
end

%% aplica��o das condi��es de fronteira
Kp=K(free_dofs, free_dofs);
Mp=M(free_dofs, free_dofs);

%% Resolu��o do problema de valores e vetores proprios

[a,b]=eig(Kp,Mp);
Wn=sqrt(diag(b));Fn=Wn/2/pi();
Phi=a;



%% Pos-processamento de resultados
X=zeros(ndof,1);Y=zeros(ndof,1);
for ielem=1:nelem
    X(ielem+1)=X(ielem) + Lel;
end

for imode=1:length(Wn)
    U=zeros(ndof,1); %inicializa��o do vetor de resultados de deslocamentos
    phi=Phi(:,imode);U(free_dofs)=phi;
    wn=Wn(imode);
    
Afactor=0.1*L/max(abs(U));
Xdef=X+Afactor*U;
Ydef=Y;

%% 

figure(1);clf(1);
Xdef(end+1)=NaN;
Ydef(end+1)=NaN;
c=U;
c(end+1)=0;


%% 
%figure(2)
xP=zeros(2*nnode,1);
yP=zeros(2*nnode,1);
cP=zeros(2*nnode,1);
h=L/20;
for j=1:2:2*nnode
    xU(j)=X((j+1)/2);%
    xU(j+1)=xU(j);%
    yU(j)=0+3*h;%
    yU(j+1)=0+h;%
    
    xP(j)=Xdef((j+1)/2);
    xP(j+1)=xP(j);
    yP(j)=0+h;
    yP(j+1)=0-h;
    cP(j)=U((j+1)/2);
    cP(j+1)=cP(j);
end
xUF=zeros(4,nelem);
yUF=zeros(4,nelem);

xF=zeros(4,nelem);
yF=zeros(4,nelem);
cF=zeros(4,nelem);
for j=1:nelem
    xUcur=xU((2*j)-1:(2*j)+2); %
    yUcur=yU((2*j)-1:(2*j)+2); %
    xcur=xP((2*j)-1:(2*j)+2);
    ycur=yP((2*j)-1:(2*j)+2);
    ccur=cP((2*j)-1:(2*j)+2);
    
    xUF(:,j)=xUcur([1 3 4 2]);%
    yUF(:,j)=yUcur([1 3 4 2]);%
    xF(:,j)=xcur([1 3 4 2]);
    yF(:,j)=ycur([1 3 4 2]);
    cF(:,j)=ccur([1 3 4 2]);
    
end
cUF=zeros(size(cF));
undefH=patch(xUF,yUF,cUF,'Linewidth',1,'Facecolor','none');
defH=patch(xF,yF,cF);

%colormap jet;
cb=colorbar;
%t=get(cb,'Limits');
%set(cb,'Ticks',linspace(t(1),t(2),10));
stringcb=sprintf('modo: %i  @ [ %7.1f Hz ]',imode,round(Fn(imode)));
cb.Label.String=stringcb;

axis equal;
pause;
end