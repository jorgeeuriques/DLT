clc
clear all
% MONORRESTITUIÇÃO PELA DIRECT LINEAR TRANSFORMATION

%% ESPAÇO IMAGEM 
%Leitura de coordenadas dos pontos na imagem. Estas medidas (observações)
%devem ser efetuadas através de algum software (ex. multispec)
%[ID PONTO,Coluna,Linha]
pts_digital = [865  1048 %02 AP
               129  317  %57 AP
               2205 456  %123 AP
               2370 1575 %129 AP
               1816 1432 %133 AP
               280 1530  %136 AP
               2686 1900 %84 VR
               1338 1477 %94 VR
               1075 89   %108 VR
               3057 503  %127 VR
               85   1111 %138 VR
               1692 1199]; %145 VR


%% PONTOS DE APOIO
%[ID PONTO,Coluna,Linha,X,Y,Z]

pts_ap=[57  129  316  444.191 669.765 21.42;%57
        84  2686 1900 675.587 374.946 1.432;%84
        94  1338 1477 525.366 492.427 7.269;%94
        108 1075 89   554.357 662.660 12.16;%108
        127 3057 503  770.173 537.567 8.069;%127
        138 85   1111 395.067 588.307 8.878;];%138


%% PONTOS DE VERIFICACAO
%[ID PONTO,Coluna,Linha,X,Y,Z]

pts_ver =[02  865  1049 488.402 563.665 8.267;%02
          123 2207 457  669.969 577.359 9.741;%123
          129 2373 1577 648.682 432.527 2.117;%129
          133 1816 1433 585.134 476.260 3.493;%133
          136 279  1532 395.148 530.693 6.973;%136
          145 1692 1199 579.864 510.927 5.448];%145


%Medidas dos pontos de apoio
xp_ap=pts_ap(:, 2);
yp_ap=pts_ap(:, 3);


%Medidas dos pontos de verificação
xp_ver=pts_ver(:, 2);
yp_ver=pts_ver(:, 3);


%% TRANSFORMACAO LINEAR DIRETA
%  AJUSTAMENTO
XYZ = [pts_ap(:,4) pts_ap(:,5)  pts_ap(:,6)];    

xy =[xp_ap,yp_ap];


 %Numero de pontos de apoio
  n=6;
    
 %Construir vetor das observacoes (Lb)
  Lb = zeros(n*2,1);
  for i=1:n
   Lb(i*2-1,1)= xp_ap(i);
   Lb(i*2,1) = yp_ap(i);
  end
  
%Parâmetros Iniciais (Aproximados)
%X0=zeros(11); %L1,L2,L3,L4,L5,L6,L7,L8,L9,L10,L11
  dx = zeros(11,11);
  b = zeros(11,1);
  X0 = zeros(11,1);
  X = XYZ(:,1);
  Y = XYZ(:,2);
  Z = XYZ(:,3);
  
 for i = 1:n
    
    % DERIVADA DA FUNCAO x
    dx(i*2-1,1) = X(i);
    dx(i*2-1,2) = Y(i);
    dx(i*2-1,3) = Z(i);
    dx(i*2-1,4) = 1;
    dx(i*2-1,9) = xy(i,1)*X(i);
    dx(i*2-1,10) = xy(i,1)*Y(i); 
    dx(i*2-1,11) = xy(i,1)*Z(i);
    
    % DERIVADA DA FUNCAO Y
    dx(i*2,5) = X(i);
    dx(i*2,6) = Y(i);
    dx(i*2,7) = Z(i);
    dx(i*2,8) = 1;
    dx(i*2,9) = xy(i,2)*X(i);
    dx(i*2,10) = xy(i,2)*Y(i);
    dx(i*2,11) = xy(i,2)*Z(i);
    
    b(i*2-1,1) = xy(i,1);
    b(i*2,1) = xy(i,2);
    b(11,1) = xy(6,1);
end  
    
X0 = dx\b

%precisao de um pixel 
%Matriz dos pesos
 %P=eye(n*2)/(tp)^2; (Usar esta se processar com coordenadas em milímetros)
 
 P=eye(n*2)/1^2; %(Usar este se as coordenadas estiverem no referencial digital)
  %======= ITERACOES ========%
  fim=false;
  iteracoes=-1;
  while ~fim %Enquanto o fim nao for true
     iteracoes=iteracoes+1;
     fprintf('Iteracao: %d\n', iteracoes);

     %Definir a matriz A e Vetor L0
    
      X = XYZ(:,1);
      Y = XYZ(:,2);
      Z = XYZ(:,3);

      A=zeros(n*2,11);
      L0=zeros(n*2,1);
  
      %MODELO FUNCIONAL: EQUAÇÕES DLT direta
      %x=(L1*X+L2*Y+L3*Z+L4)/(L9*X+L10*Y+L11*Z+1)
      %Y=(L5*X+L6*Y+L7*Z+L8)/(L9*X+L10*Y+L11*Z+1)
      
      for i = 1:n
      % VETOR L0 (L0=f(x0))
      L0(i*2-1) = (X0(1)*X(i)+X0(2)*Y(i)+X0(3)*Z(i)+X0(4))/(X0(9)*X(i)+X0(10)*Y(i)+X0(11)*Z(i)+1);
      L0(i*2) = (X0(5)*X(i)+X0(6)*Y(i)+X0(7)*Z(i)+X0(8))/(X0(9)*X(i)+X0(10)*Y(i)+X0(11)*Z(i)+1);
    
       % DERIVADA DA FUNCAO x em relação aos parâmetros
       A(i*2-1,1) = X(i)/(X0(9)*X(i)+X0(10)*Y(i)+X0(11)*Z(i)+1);%dx/dL1
       A(i*2-1,2) = Y(i)/(X0(9)*X(i)+X0(10)*Y(i)+X0(11)*Z(i)+1);%dx/dL2
       A(i*2-1,3) = Z(i)/(X0(9)*X(i)+X0(10)*Y(i)+X0(11)*Z(i)+1);%dx/dL3
       A(i*2-1,4) = 1/(X0(9)*X(i)+X0(10)*Y(i)+X0(11)*Z(i)+1);   %dx/dL4
       A(i*2-1,9) = -X(i)*(X0(1)*X(i)+X0(2)*Y(i)+X0(3)*Z(i)+X0(4))/(X0(9)*X(i)+X0(10)*Y(i)+X0(11)*Z(i)+1)^2; %dx/dL9
       A(i*2-1,10) = -Y(i)*(X0(1)*X(i)+X0(2)*Y(i)+X0(3)*Z(i)+X0(4))/(X0(9)*X(i)+X0(10)*Y(i)+X0(11)*Z(i)+1)^2;%dx/dL10
       A(i*2-1,11) = -Z(i)*(X0(1)*X(i)+X0(2)*Y(i)+X0(3)*Z(i)+X0(4))/(X0(9)*X(i)+X0(10)*Y(i)+X0(11)*Z(i)+1)^2;%dx/dL11

       % DERIVADA DA FUNCAO Y
       A(i*2,5) = X(i)/(X0(9)*X(i)+X0(10)*Y(i)+X0(11)*Z(i)+1);%dy/dL5
       A(i*2,6) = Y(i)/(X0(9)*X(i)+X0(10)*Y(i)+X0(11)*Z(i)+1);%dy/dL6
       A(i*2,7) = Z(i)/(X0(9)*X(i)+X0(10)*Y(i)+X0(11)*Z(i)+1);%dy/dL7
       A(i*2,8) = 1/(X0(9)*X(i)+X0(10)*Y(i)+X0(11)*Z(i)+1);   %dy/dL8
       A(i*2,9) = -X(i)*(X0(5)*X(i)+X0(6)*Y(i)+X0(7)*Z(i)+X0(8))/(X0(9)*X(i)+X0(10)*Y(i)+X0(11)*Z(i)+1)^2; %dy/dL9
       A(i*2,10) = -Y(i)*(X0(5)*X(i)+X0(6)*Y(i)+X0(7)*Z(i)+X0(8))/(X0(9)*X(i)+X0(10)*Y(i)+X0(11)*Z(i)+1)^2;%dy/dL10
       A(i*2,11) = -Z(i)*(X0(5)*X(i)+X0(6)*Y(i)+X0(7)*Z(i)+X0(8))/(X0(9)*X(i)+X0(10)*Y(i)+X0(11)*Z(i)+1)^2;%dy/dL11
  
    end  

     %Matriz L
     L=L0-Lb;
     %Equacoes normais
     N=A'*P*A;
     U=A'*P*L;
     X=-inv(N)*U;
     Xa=X0+X; %Parâmetros ajustados
     fim=(max(abs(X))<10^-6) || (iteracoes>50);
     X0=Xa;
  end
  V=A*X+L; %RESIDUOS
  
%Avaliação da qualidade do ajustamento
%Graus de Liberdade (número de equações - número de parâmetros)
GL=12-11;
%Variância a posteriori
varpost=(V'*P*V)/GL %VARIANCIA POSTERIORI
varXa=varpost*inv(N);
%Precisão dos parâmetros
dpXa=sqrt(diag(varXa));
  
%Teste do Chi^2
qui_cal = varpost*(GL);
% Chi^2 tabelado (Verificado na tabela em função do GL e nível de confiança)
qui1 = 0.001; qui2 = 5.024;

OBS_AJUST=[Xa dpXa];

disp('===================================================================')
disp('RESULTADO APÓS ITERAÇÕES')
disp('===================================================================')

fprintf('Parâmetro L1=  %f +- %f \n', OBS_AJUST(1,1), OBS_AJUST(1,2));
fprintf('Parâmetro L2=  %f +- %f \n', OBS_AJUST(2,1), OBS_AJUST(2,2));
fprintf('Parâmetro L3=  %f +- %f \n', OBS_AJUST(3,1), OBS_AJUST(3,2));
fprintf('Parâmetro L4=  %f +- %f \n', OBS_AJUST(4,1), OBS_AJUST(4,2));
fprintf('Parâmetro L5=  %f +- %f \n', OBS_AJUST(5,1), OBS_AJUST(5,2));
fprintf('Parâmetro L6=  %f +- %f \n', OBS_AJUST(6,1), OBS_AJUST(6,2));
fprintf('Parâmetro L7=  %f +- %f \n', OBS_AJUST(7,1), OBS_AJUST(7,2));
fprintf('Parâmetro L8=  %f +- %f \n', OBS_AJUST(8,1), OBS_AJUST(8,2));
fprintf('Parâmetro L9=  %f +- %f \n', OBS_AJUST(9,1), OBS_AJUST(9,2));
fprintf('Parâmetro L10=  %f +- %f \n', OBS_AJUST(10,1), OBS_AJUST(10,2));
fprintf('Parâmetro L11=  %f +- %f \n', OBS_AJUST(11,1), OBS_AJUST(11,2));
   

fprintf('%f<=%f<=%f \n',qui1,qui_cal,qui2);
disp('xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx')
disp('Hipótese aceita no teste do Chi^2');
disp('xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx')
  
%% Calcular X,Y pontos de verificação
%O cálculo das coordenadas pela DLT se dá através de sua forma inversa
%[X;Y]= inv (A)*C, em que A=[L1-x*L9

xy = [xp_ver,yp_ver];
Z = pts_ver(:,6);

%Definição das matrizes
  n = size(xy);
  A = zeros(2,2);
  C = zeros(2,1);
  XY = zeros(n);
  
  for i = 1:n
    A = [Xa(1)-xy(i,1)*Xa(9), Xa(2)-xy(i,1)*Xa(10);
         Xa(5)-xy(i,2)*Xa(9), Xa(6)-xy(i,2)*Xa(10)];
     
    C = [-Z(i)*(Xa(3)-xy(i,1)*Xa(11))-Xa(4)+xy(i,1);
         -Z(i)*(Xa(7)-xy(i,2)*Xa(11))-Xa(8)+xy(i,2)];
    
    XY(i,:) = [inv(A)*C]';
  end

%% COMPARACAO DAS COORDENADAS DOS PONTOS DE CONTROLE

dx= XY(:,1)- pts_ver(:,4);
dy= XY(:,2) - pts_ver(:,5);


%%
fprintf('Variância a posteriori');
varpost 

fprintf('Discrepância entre as coordenadas X dos pontos calculados e dos pontos de verificação');
dx

fprintf('Discrepância entre as coordenadas Y dos pontos calculados e dos pontos de verificação');
dy

%Média das discrepâncias na componente X
fprintf('Média dos resíduos na componente X');
mdx=mean(dx)

fprintf('Média dos resíduos na componente Y');
mdy=mean(dy)

%% Avaliação de tendências na modelagem pela DLT
%TESTE T Student
n=6;
%Hipósete 0= modelagem livre de tendência (discrepâncias com media nula a um determinado nível de significância)

fprintf('tx para avaliação da tendência na componente X');
tx=(mdx/(std(dx)))*sqrt(n)
%tx (em módulo)deve ser menor que o t tabelado, dado em função do intervalo de confiança e do tamanho da amostra.
fprintf('ty para avaliação da tendência na componente Y');
ty=(mdy/(std(dy)))*sqrt(n)
