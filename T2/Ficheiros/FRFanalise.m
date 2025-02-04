%% Obter Frequencias Experimentais para viga de 400mm
load('viga400.mat');

data_viga400 = load('viga400.mat');

frequencias_viga400 = data_viga400.f;  % Frequencias associadas a viga 400

n_tentativas = 20;  % Numero de tentativas experimentais

% Matriz para armazenar as frequencias naturais as quais ocorrem picos
frequencias_naturais_viga400 = zeros(n_tentativas, 8);

tol = 1;

hold on 
for i = 1:n_tentativas
    % Obtem cada FRF
    FRF_i = data_viga400.(['FRF_6_' num2str(i)]);

    % Encontrar os picos na componente imaginaria
    [picos, idx] = findpeaks(abs(imag(FRF_i)));

    j=1;
    for k=1:length(picos)
        if picos(k) > tol % Verifica quais os picos relevantes
           frequencias_naturais_viga400(i,j) = frequencias_viga400(idx(k));
           j=j+1;
        end
    end
    % figure(1)
    % % Plot do gráfico do modulo
    % plot(frequencias_viga400,log(abs(FRF_i)));
    % title('FRF Viga 400 mm');
    % xlabel('Frequência (Hz)');
    % ylabel('Mag(Acelerância) [m s^{-2}/N]');
    % legend('show', 'Location', 'bestoutside');  
    % grid on;

    % Plot do gráfico da parte imaginária
    % plot(frequencias_viga400,imag(FRF_i));
    % title('FRF Viga 400 mm');
    % xlabel('Frequência (Hz)');
    % ylabel('Amplitude (Parte Imaginária)');
    % legend('show', 'Location', 'bestoutside');
    % grid on;

    % Plot do gráfico da parte real
    % plot(frequencias_viga400,real(FRF_i));
    % title('FRF Viga 400 mm');
    % xlabel('Frequência (Hz)');
    % ylabel('Amplitude (Parte Real)');
    % legend('show', 'Location', 'bestoutside');
    % grid on;
end

% Obter todos os valores unicos de todas as tentativas, ignorando os zeros
valores_unicos = unique(frequencias_naturais_viga400(frequencias_naturais_viga400 ~= 0));
% Exibir os valores com a mensagem personalizada
disp('Valores unicos de frequencias naturais para l=400mm:');

% Iterar sobre os valores unicos
modosnormais_400 = [];
r=0; % Numero de Modos de Corpo Rigido
for i = 1:length(valores_unicos)
    if valores_unicos(i) < 10
        % Exibir como "Modos de Corpo Rigido"
        fprintf('Modos de Corpo Rigido: %.2f\n', valores_unicos(i));
        r=r+1;
    else
        % Exibir como "Modo X"
        modosnormais_400 = [modosnormais_400; valores_unicos(i)];
    end
end


guardar=0; % Flag que indica se foi detetao um valor semelhante
modosnormais_final_400 = []; % Valores Finais dos Modos Normais Detetados
semelhantes = []; % Vetor que guarda valores semelhantes

p=1;

for i = 1:length(modosnormais_400)-1

    % Valores de freqs distintas:
    if modosnormais_400(i+1)-modosnormais_400(i)>100
        if guardar == 1
            if i == length(modosnormais_400)-1
                semelhantes=[semelhantes;modosnormais_400(i)];
                modosnormais_final_400(p)=mean(semelhantes);
                modosnormais_final_400(p+1) = modosnormais_400(i+1);
            else 
                semelhantes=[semelhantes;modosnormais_400(i)];
                modosnormais_final_400(p)=mean(semelhantes);
                p=p+1;
                semelhantes = [];
                guardar=0;
            end
        else % guardar == 0
            if i == length(modosnormais_400)-1
                modosnormais_final_400(p) = modosnormais_400(i);
                modosnormais_final_400(p+1) = modosnormais_400(i+1);
            else 
                modosnormais_final_400(p) = modosnormais_400(i);
                p=p+1;
                guardar=0;
            end
        end

    % Valores de Freq Semelhantes:
    elseif modosnormais_400(i+1)-modosnormais_400(i)<100

        if i == length(modosnormais_400)-1 % Penúltimo Elemento
            semelhantes=[semelhantes;modosnormais_400(i)];
            modosnormais_final_400(p)=mean(semelhantes);

        else 
            semelhantes=[semelhantes;modosnormais_400(i)];
            guardar=1;
        end
    end
end


for i=1:length(modosnormais_final_400)
    fprintf('Modo %d: %.2f\n', i, modosnormais_final_400(i));
end



% Apenas para perceber os diferentes casos visualmente
% figure(1)
% plot(frequencias_viga400,log(abs(FRF_6_1)))
% title('FRF Viga 400 mm');
% xlabel('Frequência (Hz)');
% ylabel('Mag(Acelerância) [m s^{-2}/N]');
% grid on;

figure(2)
plot(frequencias_viga400,abs(imag(FRF_6_1)))
title('FRF Viga 400 mm');
xlabel('Frequência (Hz)');
ylabel('Amplitude (Componente Imaginária)');
grid on;

% figure(3)
% plot(frequencias_viga400,real(FRF_6_1))
% title('FRF Viga 400 mm');
% xlabel('Frequência (Hz)');
% ylabel('Amplitude (Componente Real)');
grid on;
hold off

%% Ver Modos Normais Viga 400 mm num gráfico 3D

% frequencias (frequencias_viga400)
% FRF_real (parte real de FRF)
% FRF_imag (parte imaginária de FRF)

FRF_real = real(FRF_6_1);  % Parte real da FRF
FRF_imag = imag(FRF_6_1);  % Parte imaginária da FRF
frequencias = frequencias_viga400;  % Frequências associadas

figure;
plot3( FRF_imag,FRF_real,frequencias, 'b', 'LineWidth', 1.5);
grid on;

xlabel('Parte Real da FRF');
ylabel('Parte Imaginária da FRF');
zlabel('Frequência (Hz)');
title('Gráfico 3D da FRF para Modos Normais');
view(45, 30);




%% Obter Frequencias Experimentais para viga de 500mm

data_viga500 = load('viga500.mat');

frequencias_viga500 = data_viga500.f;  % Frequencias associadas a viga 500

n_tentativas = 20;  % Numero de tentativas experimentais

% Matriz para armazenar as frequencias naturais as quais ocorrem picos
frequencias_naturais_viga500 = zeros(n_tentativas, 11);

tol=0.9;

hold on
for i = 1:n_tentativas
    % Obtem cada FRF
    FRF_i = data_viga500.(['FRF_7_' num2str(i)]);

    % Encontrar os picos na componente imaginaria
    [picos, idx] = findpeaks(abs(imag(FRF_i)));

    j=1;
    
    for k=1:length(picos)
        if picos(k) > tol % Verifica quais os picos relevantes
           frequencias_naturais_viga500(i,j) = frequencias_viga500(idx(k));
           j=j+1;
        end
    end

    % Plot do gráfico do modulo
    % figure(1)
    % plot(frequencias_viga500,log(abs(FRF_i)));
    % title('FRF Viga 500 mm');
    % xlabel('Frequência (Hz)');
    % ylabel('Mag(Acelerância) [m s^{-2}/N]');
    % legend('show', 'Location', 'bestoutside');  
    % grid on;

    % Plot do gráfico da parte imaginária
    % plot(frequencias_viga500,imag(FRF_i));
    % title('FRF Viga 500 mm');
    % xlabel('Frequência (Hz)');
    % ylabel('Amplitude Real');
    % legend('show', 'Location', 'bestoutside');
    % grid on;

    % Plot do gráfico da parte real
    % plot(frequencias_viga500,real(FRF_i));
    % title('FRF Viga 500 mm');
    % xlabel('Frequência (Hz)');
    % ylabel('Amplitude Imaginária');
    % legend('show', 'Location', 'bestoutside');
    % grid on;

end

title('FRFs - Parte Real');
xlabel('Frequência (Hz)');
ylabel('Amplitude');
legend('show', 'Location', 'Best');  % Mostra a legenda com os nomes das variáveis
grid on;
% Obter todos os valores unicos de todas as tentativas, ignorando os zeros
valores_unicos = unique(frequencias_naturais_viga500(frequencias_naturais_viga500 ~= 0));

% Exibir valores unicos
disp('Valores unicos de frequencias naturais para l=500mm:');

% Iterar sobre os valores unicos
modosnormais_500 = [];
r=0; % Numero de Modos de Corpo Rigido
for i = 1:length(valores_unicos)
    if valores_unicos(i) < 10
        % Exibir como "Modos de Corpo Rigido"
        fprintf('Modos de Corpo Rigido: %.2f\n', valores_unicos(i));
        r=r+1;
    else
        % Exibir como "Modo X"
        % fprintf('Modo %d: %.2f\n', i-r, valores_unicos(i));
        modosnormais_500 = [modosnormais_500; valores_unicos(i)];
    end
end

% Algoritmo de deteção e correção de valores semelhantes
guardar=0;
modosnormais_final_500 = [];
semelhantes = [];
p=1;

for i = 1:length(modosnormais_500)-1

    if modosnormais_500(i+1)-modosnormais_500(i)>100
        if guardar == 1
            if i == length(modosnormais_500)-1
                semelhantes=[semelhantes;modosnormais_500(i)];
                modosnormais_final_500(p)=mean(semelhantes);
                modosnormais_final_500(p+1) = modosnormais_500(i+1);
            else 
                semelhantes=[semelhantes;modosnormais_500(i)];
                modosnormais_final_500(p)=mean(semelhantes);
                p=p+1;
                semelhantes = [];
                guardar=0;
            end
            
        else
            if i == length(modosnormais_500)-1

                modosnormais_final_500(p) = modosnormais_500(i);
                modosnormais_final_500(p+1) = modosnormais_500(i+1);
            else 
                modosnormais_final_500(p) = modosnormais_500(i);
                p=p+1;
                guardar=0;
            end
        end

    elseif modosnormais_500(i+1)-modosnormais_500(i)<100

        if i == length(modosnormais_500)-1
            semelhantes=[semelhantes;modosnormais_500(i)];
            modosnormais_final_500(p)=mean(semelhantes);

        else 
            semelhantes=[semelhantes;modosnormais_500(i)];
            guardar=1;
        end
    end
end


for i=1:length(modosnormais_final_500)
    fprintf('Modo %d: %.2f\n', i, modosnormais_final_500(i));
end

% Apenas para perceber os diferentes casos visualmente
% figure(1)
% plot(frequencias_viga500,log(abs(FRF_7_1)))
% title('FRF Viga 500 mm');
% xlabel('Frequência (Hz)');
% ylabel('Mag(Acelerância) [m s^{-2}/N]');
% grid on;

figure(2)
plot(frequencias_viga500,abs(imag(FRF_7_1)))
title('FRF Viga 500 mm');
xlabel('Frequência (Hz)');
ylabel('Amplitude (Componente Imaginária)');
grid on;

% figure(3)
% plot(frequencias_viga500,real(FRF_7_1))
% title('FRF Viga 500 mm');
% xlabel('Frequência (Hz)');
% ylabel('Amplitude (Componente Real)');
% grid on;
hold off