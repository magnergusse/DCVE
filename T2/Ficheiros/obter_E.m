%% Calculo de E ótimo analiticamente

num=500; % Numero de iterações

% Valores obtidos experimentalmente:
Exp=[93.75,256.25,503.75,837.5,1247.5,1738.75];
nfreq=length(Exp);

% Matrizes com as frequências [Hz] nas linhas e iterações nas colunas

SS=zeros(nfreq,num); 
fSS=zeros(nfreq,num); % SS simplesmente apoiada
LL=zeros(nfreq,num);
fLL=zeros(nfreq,num); % LL - livre-livre
FL=zeros(nfreq,num);
fFL=zeros(nfreq,num); % FL - encastrada

E=zeros(1,num); % Vetor de Módulos de Young
E(2)=40E9; % [Pa] Módulo de Young Inicial
h=0.003; % [m] Espessura da viga
b=0.030; % [m] Largura da viga
L=0.4; % [m] Comprimento da viga
rho=2571; % [kg/m^3] Massa volúmica do material
m=0.23139; % [kg/m] Massa por comprimento
I=b*h^3/12; % [m^4] Momento de inércia de secção retangular

difer=zeros(nfreq,num);
adifer=zeros(nfreq,num);
media=zeros(1,num);

k=2; % Para cálculo para vários módulos de Young

% Inicializações apenas para garantir funcionamento do algoritmo
media(1,1)=10000000;
media(1,2)=999999;


while(media(1,k)<media(1,k-1) && media(1,k+1)<media(1,k) && k<num)
    k=k+1;
    E(1,k)=E(1,k-1)+0.1E9;
    for i=1:nfreq
        SS(i,k)=i*pi;
        fSS(i,k)=(1/(2*pi))*(SS(i,k)^2/L^2)*sqrt((E(1,k)*I)/m);
        LL(i,k)=(2*i+1)*pi/2;
        fLL(i,k)=(1/(2*pi))*(LL(i,k)^2/L^2)*sqrt((E(1,k)*I)/m);
        FL(i,k)=(2*i-1)*pi/2;
        fFL(i,k)=(1/(2*pi))*(FL(i,k)^2/L^2)*sqrt((E(1,k)*I)/m);
        difer(i,k)=fLL(i,k)-Exp(i);
        adifer(i,k)=abs(difer(i,k));
    end
    media(1,k)=sum(adifer(:,k))/nfreq;
end

% Para um valor a seguir ao valor corrigo
fprintf(['Analiticamente: Módulo de Young Corrigido: %.1f GPa, com ...' ...
    'variação %.3f \n'],E(1,k-1)/1E9,media(1,k-1));

%% Calculo de E ótimo numericamente

num=500; % Numero de iterações


E=zeros(1,num); % Vetor de Módulos de Young
E(2)=55e9;   % [Pa] Módulo de Young Inicial
rho=2571; % [kg/m^3] Massa volúmica do material
Lt=0.4;  % [m] Comprimento da viga
h=0.003; % [m] espessura da viga
b=0.030; % [m] largura da viga

I=b*h^3/12; % [m^4] Momento de inércia de secção retangular
A=b*h; % [m^2] Área da secção da viga

% Valores obtidos experimentalmente:
Exp=[93.75,256.25,503.75,837.5,1247.5,1738.75];
nfreq=length(Exp);

difer=zeros(nfreq,num);
adifer=zeros(nfreq,num);
media=zeros(1,num);

k=2;  % Para cálculo para vários módulos de Young

% Inicializações apenas para garantir funcionamento do algoritmo
media(1,1)=10000000;
media(1,2)=999999;

% Para dar apenas display de 6 frequências
seisfreqs=zeros(nfreq,num);

while(media(1,k)<media(1,k-1) && media(1,k+1)<media(1,k) && k<num)
    k=k+1;
    E(1,k)=E(1,k-1)+0.1E9;
    freqs = beam_solver(E(1,k),I, rho, A, Lt);

    seisfreqs(:,k)= freqs(3:nfreq+2); % começa em 3 devido a modos de corpo
                                      % rígido
    difer(:,k)=seisfreqs(:,k)-Exp(:);
    adifer(:,k)= abs(difer(:,k));
    media(1,k)=sum(adifer(:,k))/nfreq;
end

% Para um valor a seguir ao valor corrigo
fprintf(['Numericamente: Módulo de Young Corrigido: %.1f GPa, com ...' ...
    'variação %.3f \n'],E(1,k-1)/1E9,media(1,k-1));


%% Cálculo numérico das frequências

function [freqs] = beam_solver(E, I, rho, A, Lt)

    % Número de Elementos
    nelem=18;
    L=Lt/nelem;
    Rdof=[];

    % Inicialização de parâmetros
    nnode = nelem + 1;  
    ndof = nnode * 2;        
    xyU = zeros(nnode, 2);  
    lnodal = cell(nelem, 1);   
    K = zeros(ndof, ndof);    
    M = zeros(ndof, ndof);     

    for ielem = 1:nelem
        xyU(ielem, 1) = L * (ielem - 1);
        xyU(ielem + 1, 1) = L * ielem;
        lnodal{ielem} = [ielem, ielem + 1];
    end

    % Matrizes de rigidez e massa 
    for ielem = 1:nelem

        Kbeam = E * I / L^3;
        Kelem = Kbeam * [12, 6 * Lt, -12, 6 * Lt;
                         6 * Lt, 4 * Lt^2, -6 * Lt, 2 * Lt^2;
                         -12, -6 * Lt, 12, -6 * Lt;
                         6 * Lt, 2 * Lt^2, -6 * Lt, 4 * Lt^2];
                     

        Mbeam = rho * A * L;
        Melem = Mbeam / 420 * [156, 22 * Lt, 54, -13 * Lt;
                               22 * Lt, 4 * Lt^2, 13 * Lt, -3 * Lt^2;
                               54, 13 * Lt, 156, -22 * Lt;
                               -13 * Lt, -3 * Lt^2, -22 * Lt, 4 * Lt^2];

        dof1 = 2 * lnodal{ielem}(1) - 1;
        dof2 = 2 * lnodal{ielem}(2) - 1; 

        K(dof1:dof1 + 1, dof1:dof1 + 1) = K(dof1:dof1 + 1, dof1:dof1 + 1) + Kelem(1:2, 1:2);
        K(dof1:dof1 + 1, dof2:dof2 + 1) = K(dof1:dof1 + 1, dof2:dof2 + 1) + Kelem(1:2, 3:4);
        K(dof2:dof2 + 1, dof1:dof1 + 1) = K(dof2:dof2 + 1, dof1:dof1 + 1) + Kelem(3:4, 1:2);
        K(dof2:dof2 + 1, dof2:dof2 + 1) = K(dof2:dof2 + 1, dof2:dof2 + 1) + Kelem(3:4, 3:4);

        M(dof1:dof1 + 1, dof1:dof1 + 1) = M(dof1:dof1 + 1, dof1:dof1 + 1) + Melem(1:2, 1:2);
        M(dof1:dof1 + 1, dof2:dof2 + 1) = M(dof1:dof1 + 1, dof2:dof2 + 1) + Melem(1:2, 3:4);
        M(dof2:dof2 + 1, dof1:dof1 + 1) = M(dof2:dof2 + 1, dof1:dof1 + 1) + Melem(3:4, 1:2);
        M(dof2:dof2 + 1, dof2:dof2 + 1) = M(dof2:dof2 + 1, dof2:dof2 + 1) + Melem(3:4, 3:4);
    end

    free_dofs = setxor(1:ndof, Rdof);
    Kp = K(free_dofs, free_dofs);
    Mp = M(free_dofs, free_dofs);


    [~, b] = eig(Kp, Mp);
    freqs=sqrt(diag(b))/2/pi; %[Hz]
end