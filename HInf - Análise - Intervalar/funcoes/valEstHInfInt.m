function saida = valEstHInfInt(ACal,ECal,CCal,KCal,h,delta,tol)
    % validarLMIHInf
    %
    % Avaliar norma HInf de um sistema híbrido
    %
    % input:  ACal,ECal,CCal,KCal -> Matrizes aumentadas
    %         h -> período de amostragem
    %         delta -> relacionada a aproximação da derivada
    %         tol -> tolerância do solver
    %
    % output: saida.X -> Matriz Lyapunov
    %
    % Criado date: 15 Fev 2020
    % Revisão date: 15 Fev 2020
    % Autor: rafaelmartinsalves@gmail.com
    
    % Matrizes Incertas
    ACal0 = mid(ACal);
    deltaACal = rad(ACal);
    
    ECal0 = mid(ECal);
    deltaECal = rad(ECal);
    
    CCal0 = mid(CCal);
    deltaCCal = rad(CCal);
    
    % Dimensões das matrizes de entrada
    nXi = length(ACal.inf);
    nW = size(ECal.inf, 2);
    nY = size(CCal.inf, 1);
    
    % Inicialização das variáveis LMIs
    saida.cpuSeg_m = clock;
    %LMI contador de linhas
    saida.linhas = 0;
    %Inicia o sistema de LMI
    LMIs = set([]);
    
    % Cria as matrizes X_i     
    quantInt = round(h/delta);
    X = cell(quantInt);
    for i=1:quantInt
        X{i} = sdpvar(nXi,nXi,'symmetric');
        LMIs = LMIs + (X{i} >= 0);
        saida.linhas = saida.linhas + size(X{i},1);
    end
    
    % Função objetiva
    mu = sdpvar(1,1);
    obj = mu;
    
    %   Alphas 1
    for i=1:nXi
        for j=1:nXi
            alphas(i,j) = sdpvar(1,1); % Scalar
        end
    end
    
    %    Gammas
    for i=1:nY
        for j=1:nXi
            gammas(i,j) = sdpvar(1,1); % Scalar
        end
    end
    
    %    Epsilons
    for i=1:nXi
        for j=1:nW
            epsilons(i,j) = sdpvar(1,1); % Scalar
        end
    end
    
    
    %   Calcula S
    for i=1:nXi  
       for j=1:nXi       
            ej = generatesEi(j,nXi); 
            productOfDeltaAAlphaij = (abs(deltaACal(i,j))^2)*alphas(i,j)*(ej*ej');
            
            % Primeira vez i = j = 1
            if (i==1 && j == 1) 
                S = productOfDeltaAAlphaij;
            else
                S = productOfDeltaAAlphaij + S;
            end
            
       end
    end
    
    
     %   Calcula T
    for i=1:nXi %
       for j=1:nW       
            fj = generatesEi(j,nW); 
            productOfDeltaEEpsilonij = (abs(deltaECal(i,j))^2)*epsilons(i,j)*(fj*fj');
            
            % Primeir vez i = j = 1
            if (i==1 && j == 1) 
                T = productOfDeltaEEpsilonij;
            else
                T = productOfDeltaEEpsilonij + T;
            end
            
       end
    end
    
    
    %   Calcula U
    for i=1:nY %
       for j=1:nXi       
            gi = generatesEi(i,nY); 
            productOfDeltaCGammasnij = (abs(deltaCCal(i,j))^2)*gammas(i,j)*(gi*gi');
            
            % Primeir vez i = j = 1
            if (i==1 && j == 1) 
                U = productOfDeltaCGammasnij;
            else
                U = productOfDeltaCGammasnij + U;
            end
            
       end
    end
    
    
    Alphas = diag(reshape(alphas',[nXi*nXi, 1]));
    Gammas = diag(reshape(gammas',[nY*nXi, 1]));
    Epsilons = diag(reshape(epsilons',[nW*nXi, 1]));
    
  
    
    % Montar as restrições LMIs
    for i=1:(quantInt-1)
        
        
        der = (X{i+1}-X{i})/delta;
        T_i = gerarLMIAlternadaDerivadaInt(ACal0,ECal0,CCal0,der,X,S,T,U,Alphas,Gammas,Epsilons,mu,nXi,nW,nY,i);
        T_iMais1 = gerarLMIAlternadaDerivadaInt(ACal0,ECal0,CCal0,der,X,S,T,U,Alphas,Gammas,Epsilons,mu,nXi,nW,nY,i+1);
        
        LMIs = LMIs + (T_i <= 0);
        LMIs = LMIs + (T_iMais1 <= 0);
        
        % Conta quantidade de linhas LMIs
        saida.linhas = saida.linhas + size(T_i,1);
        saida.linhas = saida.linhas + size(T_iMais1,1);
        
    end
    
    
    % Restrição para garantir que a Lyapunov sempre decai
    U11 = X{quantInt};
    U21 = X{1}*KCal;
    U22 = X{1};
   
    
    T =  [U11 U21';
          U21 U22];
    
    LMIs = LMIs + (T >= 0);
    
    % Conta qtd de linhas LMIs, tempo de montagem e qtd de variáveis
    saida.linhas = saida.linhas + size(T,1);
    saida.cpuSeg_m = etime(clock,saida.cpuSeg_m);
    saida.var = size(getvariables(LMIs),2);
    
    % Resolve as LMIs
    sol = optimize(LMIs,obj,sdpsettings('verbose',0,'solver','mosek'));
    saida.cpusec = sol.solvertime;
    
    % Checar a solução
    p = min(checkset(LMIs));
    saida.delta = p;
    saida.factivel = 0;
    
    if(sol.problem ~= 0 ) % Factível
       saida.factivel = 0;
    else
        saida.factivel = 1;
        for i=1:quantInt
            saida.X{i} = double(X{i});
        end
        saida.mu = double(mu);
        saida.gamma = sqrt(saida.mu);
    end
    

end