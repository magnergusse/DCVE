


%% Modelo 3dof Modelos  de representação
%
% Script sobre modelos de representação 
% 
% Modelo espacial: [M], [C] , [K]
% Modelo modal: {omega_r}, [Phi]
% Modelo de resposta: H_jk(omega)
%
%%
% Modelo espacial
m1=1; m2=1;m3=1;
K1=1000;K2=1000;K3=1000;K4=1000;K5=2000;

M=[m1 0 0;
   0 m2 0;
   0 0 m3];

K=[K1+K2+K3    -K2      -K3; 
    -K2      K2+K3+K5   -K4;
    -K3        -K4     K3+K4];

alpha=1e-2; %factor de amortecimento
C=K*alpha; %matriz de amortecimento [C]=alpha*[K] (proporcional)

%% Modelo modal (real- considera C=0)


[a,b]=eig(K,M);     % Resolução do problema de valores e vetores próprios

Phi=a;              % matriz de vetores modais [adimensional relativa]
Wn=sqrt(diag(b));   % vetor de frequências naturais [rad/s]
Fn=Wn/2/pi;         % vetor de frequências naturais [Hz]

%output de resultados

disp(['Frequência natural fundamental w1 (f1): ',num2str(Wn(1)),' rad/s (',num2str(Fn(1)),' Hz))']);
disp(['Vetor modal para w1: ',num2str(Phi(1,1)),' [dof1] e ', num2str(Phi(2,1)),' [dof2] e ',num2str(Phi(3,1)),' [dof3].']);

disp(['2a Frequência natural w2 (f2): ',num2str(Wn(2)),' rad/s (',num2str(Fn(2)),' Hz))']);
disp(['Vetor modal para w2: ',num2str(Phi(1,2)),' [dof1] e ', num2str(Phi(2,2)),' [dof2] e ', num2str(Phi(3,2)),' [dof3].']);

disp(['3a Frequência natural w3 (f3): ',num2str(Wn(3)),' rad/s (',num2str(Fn(3)),' Hz))']);
disp(['Vetor modal para w3: ',num2str(Phi(1,3)),' [dof1] e ', num2str(Phi(2,3)),' [dof2] e ', num2str(Phi(3,3)),' [dof3].']);


%% Modelo de resposta gerado pelo modelo espacial

% [M],[C],[K] --> H_jk(omega)

disp('>> Modelo de resposta gerado pelo modelo espacial (press key)');
pause

f=linspace(0,2*Fn(3),10000); %vetor de frequências 0 2*Fn3, 10000 valores
w=2*pi()*f;w2=w.^2; 

for i=1:length(f)
    wi=w(i);wi2=w2(i);
    Zw=(-wi2*M+1i*wi*C+K); % matriz da inversa da função de transferência
    iZw=Zw^-1;             % matriz da função de transferência
    H(1:3,1:3,i)=iZw;      % matriz das funções de resposta em frequência
end

%output de resultados
for j=1:3
    for k=1:3
        figure(3*(j-1)+k);
        frf(:)=H(j,k,:);
        sf1=subplot(2,1,1);semilogy(f,abs(frf));
        sf2=subplot(2,1,2);plot(f,angle(frf));
        set(sf1,'Position',[0.13,0.3, 0.775, 0.6]);
        set(sf2,'Position',[0.13,0.1, 0.775, 0.2]);
        frfname=['H_',num2str(j),'_',num2str(k),' [m]'];
        set(get(sf2,'XLabel'),'String','Frequência [Hz]');
        set(get(sf1,'YLabel'),'String',frfname);
        set(get(sf2,'YLabel'),'String','Fase [rad]');
        fh=get(get(sf2,'Ylabel'),'Fontsize');
        set(get(sf1,'YLabel'),'Fontsize',fh);
        set(sf1,'XTicklabel',{});
    end;
end;

% representação 3D da Função de Resposta em Frequência

figure(10);frf1(:)=H(1,1,:);
plot3(f,real(frf1),imag(frf1)); ylabel('REAL'),zlabel('IMAG'),xlabel('Freq')

figure(11);frf2(:)=H(1,2,:);frf3(:)=H(1,3,:);
plot3(f,real(frf1),imag(frf1),f,real(frf2),imag(frf2),f,real(frf3),imag(frf3)); ylabel('REAL'),zlabel('IMAG'),xlabel('Freq')


figure(12); plot(real(frf1),imag(frf1),real(frf2),imag(frf2),real(frf3),imag(frf3)); xlabel('REAL');ylabel('IMAG');axis equal;

pause

%% Modelo de resposta gerado pelo modelo modal real

% {omega_r},[Phi]--> Alpha_jk(omega)

disp('>> Modelo de resposta gerado pelo modelo modal real (press key)');
pause

f=linspace(0,2*Fn(3),10000); %vetor de frequências 0 2*Fn3, 10000 valores
w=2*pi()*f;w2=w.^2;

for j=1:3
    for k=1:3
        frf=0;
        for r=1:3
            wr=Wn(r);
            frf=frf+ Phi(j,r)*Phi(r,k)./(wr^2-w.^2);
        end
        Alpha(j,k,:)=frf;
    end
end

%output de resultados
for j=1:3
    for k=1:3
        figure(12+3*(j-1)+k);
        frf(:)=Alpha(j,k,:);
        sf1=subplot(2,1,1);semilogy(f,abs(frf));
        sf2=subplot(2,1,2);plot(f,angle(frf));
        set(sf1,'Position',[0.13,0.3, 0.775, 0.6]);
        set(sf2,'Position',[0.13,0.1, 0.775, 0.2]);
        frfname=['H_',num2str(j),'_',num2str(k),' [m]'];
        set(get(sf2,'XLabel'),'String','Frequência [Hz]');
        set(get(sf1,'YLabel'),'String',frfname);
        set(get(sf2,'YLabel'),'String','Fase [rad]');
        fh=get(get(sf2,'Ylabel'),'Fontsize');
        set(get(sf1,'YLabel'),'Fontsize',fh);
        set(sf1,'XTicklabel',{});
    end;
end;


