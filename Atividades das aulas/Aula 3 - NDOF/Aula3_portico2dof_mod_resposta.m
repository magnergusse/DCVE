%% Pórtico Modelo de resposta
%
%%%%% exercício Pórtico 2DOF %%%%
%%       _________
%       |___m2____|    m1=m2=1.2kg
%       |         |    K1=K2=4xKviga
%       | K2      |    Kviga=12EI/L^3
%       |_________|    E=70GPa
%       |___m1____|    I=b*h^3/12
%       |         |    b=10mm
%       | K1      |    h=2mm
%     __|_________|__
%     ///////////////
%
%%
m1=1.2; m2=m1;
b=10e-3;h=2e-3;L=0.150;
E=70e9;
I=b*h^3/12;
Kviga=12*E*I/L^3;
K1=4*Kviga;
K2=K1;

M=[m1 0;0 m2];
K=[K1+K2 -K2; -K2 K2];

[a,b]=eig(K,M)

Phi=a;
Wn=sqrt(diag(b));
Fn=Wn/2/pi;


disp(['Frequência natural fundamental w1 (f1): ',num2str(Wn(1)),' rad/s (',num2str(Fn(1)),' Hz))']);
disp(['Vetor modal para w1: ',num2str(Phi(1,1)),' dof1 e ', num2str(Phi(2,1)),' dof 2.']);

disp(['2a Frequência natural w2 (f2): ',num2str(Wn(2)),' rad/s (',num2str(Fn(2)),' Hz))']);
disp(['Vetor modal para w2: ',num2str(Phi(1,2)),' dof1 e ', num2str(Phi(2,2)),' dof 2.']);


alpha=6e-3; %amortecimento
C=K*alpha; %matriz de amortecimento [C]=alpha*[K]
% [M]a(t)+[C]v(t)+[K]x(t)=f(t)

F=1;

w=linspace(0,2*Wn(2),1000); %vetor de frequências 0 2*Wn2, 1000 valores
f=w/2/pi; %Hz
for i=1:length(w)
    Z11(i)=-w(i)^2*M(1,1)+1i*w(i)*C(1,1)+K(1,1); %-w^2*m11+jw*c11+k11
    Z12(i)=1i*w(i)*C(1,2)+K(1,2); %-w^2*m12+jw*c12+k12
    Z21(i)=1i*w(i)*C(2,1)+K(2,1); %-w^2*m21+jw*c21+k21
    Z22(i)=-w(i)^2*M(2,2)+1i*w(i)*C(2,2)+K(2,2); %-w^2*m22+jw*c22+k22
    Det(i)=Z11(i)*Z22(i)-Z12(i)*Z21(i); %determinante da matriz Z
    X1(i)=Z22(i)/Det(i)*F;
    X2(i)=-Z21(i)/Det(i)*F;
end


figure(1);
sf1=subplot(2,1,1);semilogy(f,abs(X1));
sf2=subplot(2,1,2);plot(f,angle(X1));
set(sf1,'Position',[0.13,0.3, 0.775, 0.6]); 
set(sf2,'Position',[0.13,0.1, 0.775, 0.2]); 
set(get(sf2,'XLabel'),'String','Frequência [Hz]');
set(get(sf1,'YLabel'),'String','|X_2| [m]');
set(get(sf2,'YLabel'),'String','Fase [rad]');
fh=get(get(sf2,'Ylabel'),'Fontsize');
set(get(sf1,'YLabel'),'Fontsize',fh);
set(sf1,'XTicklabel',{});

figure(2);
sf3=subplot(2,1,1);semilogy(f,abs(X2));
sf4=subplot(2,1,2);plot(f,angle(X2));
set(sf3,'Position',[0.13,0.3, 0.775, 0.6]); 
set(sf4,'Position',[0.13,0.1, 0.775, 0.2]); 
set(get(sf4,'XLabel'),'String','Frequência [Hz]');
set(get(sf3,'YLabel'),'String','|X_2| [m]');
set(get(sf4,'YLabel'),'String','Fase [rad]');
fh=get(get(sf4,'Ylabel'),'Fontsize');
set(get(sf3,'YLabel'),'Fontsize',fh);
set(sf3,'XTicklabel',{});



    
    
    
