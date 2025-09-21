function saida = estHInfSint(A,B,E,C,D,h,delta,tol)
    % validarLMIHInf
    %
    % Síntase norma HInf de um sistema híbrido
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
    A0 = mid(A);
    deltaA = rad(A);
    
    B0 = mid(B);
    deltaB = rad(B);
    
    E0 = mid(E);
    deltaE = rad(E);
    
    C0 = mid(C);
    deltaC = rad(C);
    
    D0 = mid(D);
    deltaD = rad(D);
    
    %     Tamanho das matrizes
    nx = length(A.inf);
    nu = size(B.inf,2);
    nw = size(E.inf,2);
    ny = size(C.inf,1);
    
    
   
    
    % Inicialização das variáveis LMIs
    saida.cpuSeg_m = clock;
    %LMI contador de linhas
    saida.linhas = 0;
    %Inicia o sistema de LMI
    LMIs = set([]);
    
    % Cria as matrizes X_i     
    quantInt = round(h/delta);
    S = cell(quantInt,1);
    L = cell(quantInt,1);
    Z = cell(quantInt,1);
    for i=1:quantInt
        S{i} = sdpvar(nx,nx,'symmetric');
        L{i} = sdpvar(nx,nu);
        Z{i} = sdpvar(nu,nu,'symmetric');
    end
    
    % Função objetiva
    mu = sdpvar(1,1);
    obj = mu;
    
    %   Alphas 
    for i=1:nx
        for j=1:nx
            alphas(i,j) = sdpvar(1,1); % Scalar
        end
    end
    
    %   betas 
    for i=1:nx
        for j=1:nu
            betas(i,j) = sdpvar(1,1); % Scalar
        end
    end
    
    %   gammas 
    for i=1:ny
        for j=1:nx
            gammas(i,j) = sdpvar(1,1); % Scalar
        end
    end
    
    %   deltas 
    for i=1:ny
        for j=1:nu
            deltas(i,j) = sdpvar(1,1); % Scalar
        end
    end
    
    %   epsilons 
    for i=1:nx
        for j=1:nw
            epsilons(i,j) = sdpvar(1,1); % Scalar
        end
    end
    
     %   Calcula S
    for i=1:nx  
       for j=1:nx       
            ei = generatesEi(i,nx); 
            aux = 2*(abs(deltaA(i,j))^2)*alphas(i,j)*(ei*ei');
            
            % Primeira vez i = j = 1
            if (i==1 && j == 1) 
                S_ = aux;
            else
                S_ = aux + S_;
            end
            
       end
    end
    
     for i=1:nx  
       for j=1:nu       
            ei = generatesEi(i,nx); 
            aux = 2*(abs(deltaB(i,j))^2)*betas(i,j)*(ei*ei');
            
            S_ = aux + S_;
            
       end
     end
    
     
    %   Calcula U
    for i=1:ny  
       for j=1:nx       
            gi = generatesEi(i,ny); 
            aux = 2*(abs(deltaC(i,j))^2)*gammas(i,j)*(gi*gi');
            
            % Primeira vez i = j = 1
            if (i==1 && j == 1) 
                U_ = aux;
            else
                U_ = aux + U_;
            end
            
       end
    end
    
     for i=1:ny  
       for j=1:nu       
            gi = generatesEi(i,ny); 
            aux = 2*(abs(deltaD(i,j))^2)*deltas(i,j)*(gi*gi');
            
            U_ = aux + U_;
            
       end
     end
    
    %   Calcula T
    for i=1:nx  
       for j=1:nw       
            hj = generatesEi(j,nw); 
            aux = (abs(deltaE(i,j))^2)*epsilons(i,j)*(hj*hj');
            
            % Primeira vez i = j = 1
            if (i==1 && j == 1) 
                T_ = aux;
            else
                T_ = aux + T_;
            end
            
       end
    end
    
    Alphas = diag(reshape(alphas',[nx*nx, 1]));
    Betas = diag(reshape(betas',[nx*nu, 1]));
    Gammas = diag(reshape(gammas',[ny*nx, 1]));
    Deltas = diag(reshape(deltas',[ny*nu, 1]));
    Epsilons = diag(reshape(epsilons',[nx*nw, 1]));
   
    
    % Montar as restrições LMIs
    for i=1:(quantInt-1)
        % Para i
        % Derivadas
        derS = (S{i+1}-S{i})/delta;
        derL = (L{i+1}-L{i})/delta;
        derZ = (Z{i+1}-Z{i})/delta;
    
        T_i = gerarSinLMIAlternadaDerivada(A0,B0,E0,C0,D0,Alphas,Betas,Gammas,Deltas,Epsilons,S_,T_,U_,nx,nw,ny,S,Z,L,derS,derL,derZ,mu,i);
        T_iMais1 = gerarSinLMIAlternadaDerivada(A0,B0,E0,C0,D0,Alphas,Betas,Gammas,Deltas,Epsilons,S_,T_,U_,nx,nw,ny,S,Z,L,derS,derL,derZ,mu,i+1);
        
        LMIs = LMIs + (T_i <= 0);
        LMIs = LMIs + (T_iMais1 <= 0);
        
        % Conta quantidade de linhas LMIs
        saida.linhas = saida.linhas + size(T_i,1);
        saida.linhas = saida.linhas + size(T_iMais1,1);
    end
    
    % Segunda restrição
    % Restrição para garantir que a Lyapunov sempre decai
     W = sdpvar(nx,nx,'symmetric');
     KChap = sdpvar(nu,nx);
     U11 = S{quantInt};
     U21 = L{quantInt}';
     U22 = Z{quantInt};
     U31 = S{quantInt};
     U32 = L{quantInt};
     U33 = W;
    
     T =  [U11 U21' U31';
           U21 U22  U32';
           U31 U32  U33];
    LMIs = LMIs + (T >= 0);
    saida.linhas = saida.linhas + size(T,1);
    
    U11 = W;
    U21 = W;
    U22 = S{1};
    U31 = KChap;
    U32 = L{1}';
    U33 = Z{1};

    T =  [U11 U21' U31';
          U21 U22  U32';
          U31 U32  U33];
   LMIs = LMIs + (T >= 0);
   saida.linhas = saida.linhas + size(T,1);
    
   % Terceira restrição
   for i=1:quantInt
       U11 = S{i};
       U21 = L{i}';
       U22 = Z{i};
        T =  [U11 U21';
              U21 U22];
      LMIs = LMIs + (T >= 0);  
      saida.linhas = saida.linhas + size(T,1);
   end
  
    
    % Conta qtd de linhas LMIs, tempo de montagem e qtd de variáveis
  
    saida.cpuSeg_m = etime(clock,saida.cpuSeg_m);
    saida.var = size(getvariables(LMIs),2);
    
    % Resolve as LMIs  
    sol = optimize(LMIs,obj,sdpsettings('verbose',0,'solver','sedumi'));
    saida.cpusec = sol.solvertime;
    
    % Checar a solução
    p = min(checkset(LMIs));
    saida.delta = p;
    saida.factivel = 0;
    
    if(p > -tol ) % Factível
        saida.factivel = 1;
        for i=1:quantInt
            saida.S{i} = double(S{i});
            saida.L{i} = double(L{i});
            saida.Z{i} = double(Z{i});
        end
        saida.W = double(W);
        saida.KChap = double(KChap);
        saida.K = saida.KChap*inv(saida.W);
        saida.mu = double(mu);
        saida.gamma = sqrt(saida.mu);
    end
    

end