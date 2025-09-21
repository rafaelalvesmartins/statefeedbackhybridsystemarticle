function saida = estHInfSintPolyLMILab(poly,h,delta,tol)
    % estHInfSintPolyLMILab
    %
    % Sintase Ganho K com custo HInf garantido de um sistema híbrido
    % intervalar
    %
    % input:  poly -> Matrizes do sistema
    %         h -> período de amostragem
    %         delta -> relacionada a aproximação da derivada
    %         tol -> tolerância do solver
    %
    % output: saida.S{i},saida.Z{i},saida.L{i},saida.W -> Matrizes Lyapunov
    %         saida.K -> Ganho de alimentação de estado
    %         saida.gamma -> Custo garantido HInf
    % Criado data: 12 Set 2020
    % Revisão data: 12 Set 2020
    % Autor: rafaelalves@dca.fee.unicamp.br;msouza@fee.unicamp.br
    
     % Matrizes Incertas
    A = poly.APoly;
    B = poly.BPoly;
    E = poly.EPoly;
    C = poly.CPoly;
    D = poly.DPoly;
%     poly.D1Poly = D1Poly;
    
%     Tamanho das matrizes
    nx = length(A{1});
    nu = size(B{1},2);
    nw = size(E{1},2);
    ny = size(C{1},1);
    N = length(A);
    
    % Inicialização das variáveis LMIs
    %Inicia o sistema de LMI
    setlmis([]);
	%Configura contador LMI
    c = 1;
    
    % Cria as matrizes S_i,L_i,Z_i     
    quantInt = round(h/delta);
    S = cell(N,quantInt);
    sS = cell(N,quantInt);
    
    L = cell(N,quantInt);
    sL = cell(N,quantInt);
    
    Z = cell(N,quantInt);
    sZ = cell(N,quantInt);
    for i = 1:N
        for j=1:quantInt
            [S{i,j},nVar,sS{i,j}] = lmivar(1,[nx 1]);
            [L{i,j},nVar,sL{i,j}] = lmivar(2,[nx nu]);
            [Z{i,j},nVar,sZ{i,j}] = lmivar(1,[nu 1]);
        end
    end
    
    % Cria as matrizes Yi
    X.X1_11 = lmivar(2,[nx nx]);
    X.X1_12 = lmivar(2,[nx nu]);
    X.X1_21 = lmivar(2,[nu nx]);
    X.X1_22 = lmivar(2,[nu nu]);
    
    X.X2_11 = lmivar(2,[nw nx]);
    X.X2_12 = lmivar(2,[nw nu]);
  
    
    X.X3_11 = lmivar(2,[ny nx]);
    X.X3_12 = lmivar(2,[ny nu]);
    
    X.X4_11 = lmivar(2,[nx nx]);
    X.X4_12 = lmivar(2,[nx nu]);
    X.X4_21 = lmivar(2,[nu nx]);
    X.X4_22 = lmivar(2,[nu nu]);
    
    
    % Função objetiva
	mu = lmivar(1,[1 0]);
    lmiterm([-c 1 1 mu],1,1); % mu > 0
    c = c + 1;
    obj = mu;
    
    % Montar as restrições LMIs
    % Primeira restrição
    for i = 1:N
        for j=1:(quantInt-1)

            % Para j
            % Derivadas
            lmiterm([c 1 1 S{i,j+1}],-1/delta,1);
            lmiterm([c 1 1 S{i,j}],1/delta,1);
            lmiterm([c 2 1 -L{i,j+1}],-1/delta,1);
            lmiterm([c 2 1 -L{i,j}],1/delta,1);
            lmiterm([c 2 2 Z{i,j+1}],-1/delta,1);
            lmiterm([c 2 2 Z{i,j}],1/delta,1);
                                 
            gerarSinLMIAlternadaDerivadaPolyLMILab(A{i},B{i},E{i},C{i},D{i},S{i,j},Z{i,j},L{i,j},X,mu,c);
             
            c = c + 1;

            % Para i+1

            % Primeira restrição
            % Derivadas
            lmiterm([c 1 1 S{i,j+1}],-1/delta,1);
            lmiterm([c 1 1 S{i,j}],1/delta,1);
            lmiterm([c 2 1 -L{i,j+1}],-1/delta,1);
            lmiterm([c 2 1 -L{i,j}],1/delta,1);
            lmiterm([c 2 2 Z{i,j+1}],-1/delta,1);
            lmiterm([c 2 2 Z{i,j}],1/delta,1);

             gerarSinLMIAlternadaDerivadaPolyLMILab(A{i},B{i},E{i},C{i},D{i},S{i,j+1},Z{i,j+1},L{i,j+1},X,mu,c)
            c = c + 1;

        end
    end
    
    % Segunda restrição
    % Restrição para garantir que a Lyapunov sempre decai
    W = lmivar(1,[nx 1]);
    KChap = lmivar(2,[nu nx]);
    for i = 1:N
        lmiterm([-c 1 1 S{i,quantInt}],1,1);
        lmiterm([-c 2 1 -L{i,quantInt}],1,1);
        lmiterm([-c 2 2 Z{i,quantInt}],1,1);
        lmiterm([-c 3 1 S{i,quantInt}],1,1);
        lmiterm([-c 3 2 L{i,quantInt}],1,1);
        lmiterm([-c 3 3 W],1,1);
        c = c + 1;
        lmiterm([-c 1 1 W],1,1);
        lmiterm([-c 2 1 W],1,1);
        lmiterm([-c 2 2 S{i,1}],1,1);
        lmiterm([-c 3 1 KChap],1,1);
        lmiterm([-c 3 2 -L{i,1}],1,1);
        lmiterm([-c 3 3 Z{i,1}],1,1);
        c = c + 1;
    end
    % Terceira restrição
    for i = 1:N
        for j=1:quantInt
            lmiterm([-c 1 1 S{i,j}],1,1); % S{i} > 0
            lmiterm([-c 2 1 -L{i,j}],1,1); 
            lmiterm([-c 2 2 Z{i,j}],1,1); % Z{i} > 0
            c = c + 1;
        end
    end
        
    LMIs = getlmis;

    c = zeros(decnbr(LMIs),1);
    c(diag(decinfo(LMIs,obj))) = 1;
    options = [tol,2000,0,200,1];
    [copt,xopt] = mincx(LMIs,c,options);

    
    if isempty(xopt)
       saida.feas = 0;
    else
        saida.factivel = 1;
         for i = 1:N
            for j=1:quantInt
                saida.S{i,j} = dec2mat(LMIs,xopt,S{i,j});
                saida.Z{i,j} = dec2mat(LMIs,xopt,Z{i,j});
                saida.L{i,j} = dec2mat(LMIs,xopt,L{i,j});
            end
        end
        saida.mu = dec2mat(LMIs,xopt,mu);
        saida.gamma = sqrt(saida.mu);
        saida.W = dec2mat(LMIs,xopt,W);
        saida.KChap = dec2mat(LMIs,xopt,KChap);
        saida.K = saida.KChap*(saida.W^-1);
    end
    

end