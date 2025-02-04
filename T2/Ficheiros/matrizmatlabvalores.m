%% Display de matrizes no formato fLL | fFL | fSS (Analiticamente), com 
% correção do Módulo de Young e Massa Volúmica

nfreqs=8; % Numero de Frequências a obter

% Matrizes com as frequências [Hz] nas linhas e iterações nas colunas
SS=zeros(8,1);
fSS=zeros(8,1);
LL=zeros(8,1);
fLL=zeros(8,1);
FL=zeros(8,1);
fFL=zeros(8,1);

E=60.3E9; % [Pa] Módulo de Young 
b=3E-2; % [m] largura da viga
h=3E-3; % [m] espessura da viga
L=0.4; % [m] Comprimento da viga
rho=2571; % [kg/m^3] Massa volúmica do material
m=0.23139; % [kg/m] Massa por comprimento
I=b*h^3/12; % [m^4] Momento de inércia de secção retangular

matrizfreqs = zeros(nfreqs,3);

for i=1:nfreqs
    SS(i)=i*pi;
    fSS(i)=(1/(2*pi))*(SS(i)^2/L^2)*sqrt((E*I)/m);
    LL(i)=(2*i+1)*pi/2;
    fLL(i)=(1/(2*pi))*(LL(i)^2/L^2)*sqrt((E*I)/m);
    FL(i)=(2*i-1)*pi/2;
    fFL(i)=(1/(2*pi))*(FL(i)^2/L^2)*sqrt((E*I)/m);
end

matrizfreqs(:,1)=fLL;
matrizfreqs(:,2)=fFL;
matrizfreqs(:,3)=fSS;


%% Display de matrizes no formato fLL | fFL | fSS (Numericamente), com 
% correção do Módulo de Young e Massa Volúmica

num=500;

E=60.4E9; % [Pa] Módulo de Young 
rho=2571; % [kg/m^3] Massa volúmica do material
Lt=0.5; % [m] Comprimento da viga
h=0.003; % [m] espessura da viga
b=0.030; % [m] largura da viga

I=b*h^3/12; % [m^4] Momento de inércia de secção retangular
A=b*h; % [m^2] Área da secção da viga
nfreq=8; % Numero de Frequências a obter

matrizfreq = zeros(nfreq,3); % Matriz de Resultados

for i = 1:3
    freqs = beam_solver_LL(18, E, I, rho, A, Lt, i);
    if i == 1
        matrizfreq(:, i) = freqs(3:nfreq+2); % LL - livre-livre
    elseif i == 2
        matrizfreq(:, i) = freqs(1:nfreq); % FL - encastrada
    else
        matrizfreq(:, i) = freqs(1:nfreq); % SS - simplesmente apoiada
    end
end


%% Cálculo numérico das frequências

function [freqs] = beam_solver_LL(nelem,E, I, rho, A, Lt,rdofs)

    L=Lt/nelem;

    switch rdofs
        case 1
            Rdof = []; % LL - livre-livre
        case 2
            Rdof = [1, 2]; % FL - encastrada
        case 3
            Rdof = [1, (nelem + 1) * 2 - 1]; % SS - simplesmente apoiada
    end

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