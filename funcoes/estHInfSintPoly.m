function saida = estHInfSintPoly(poly,h,delta,tol)
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
    saida.cpuSeg_m = clock;
    %LMI contador de linhas
    saida.linhas = 0;
    %Inicia o sistema de LMI
    LMIs = set([]);
    
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
            S{i,j} = sdpvar(nx,nx,'symmetric');
            L{i,j} = sdpvar(nx,nu);
            Z{i,j} = sdpvar(nu,nu,'symmetric');
        end
    end
    
    % Cria as matrizes Yi
    X.X1_11 = sdpvar(nx,nx);
    X.X1_12 = sdpvar(nx,nu);
    X.X1_21 = sdpvar(nu,nx);
    X.X1_22 = sdpvar(nu,nu);
    
    X.X2_11 = sdpvar(nw,nx);
    X.X2_12 = sdpvar(nw,nu);
  
    
    X.X3_11 = sdpvar(ny,nx);
    X.X3_12 = sdpvar(ny,nu);
    
    X.X4_11 = sdpvar(nx,nx);
    X.X4_12 = sdpvar(nx,nu);
    X.X4_21 = sdpvar(nu,nx);
    X.X4_22 = sdpvar(nu,nu);
    
    
    % Função objetiva
	mu = sdpvar(1,1);
    obj = mu;
    
    % Montar as restrições LMIs
    % Primeira restrição
    for i = 1:N
        for j=1:(quantInt-1)

            % Para j
            % Derivadas
            derS = (S{i,j+1}-S{i,j})/delta;
            derL = (L{i,j+1}-L{i,j})/delta;
            derZ = (Z{i,j+1}-Z{i,j})/delta;
                                 
            T_i = gerarSinLMIAlternadaDerivadaPoly(A{i},B{i},E{i},C{i},D{i},S{i,j},Z{i,j},L{i,j},X,derS,derL,derZ,mu);
            T_iMais1 = gerarSinLMIAlternadaDerivadaPoly(A{i},B{i},E{i},C{i},D{i},S{i,j+1},Z{i,j+1},L{i,j+1},X,derS,derL,derZ,mu);
             
            LMIs = LMIs + (T_i <= 0);
            LMIs = LMIs + (T_iMais1 <= 0);

            % Conta quantidade de linhas LMIs
            saida.linhas = saida.linhas + size(T_i,1);
            saida.linhas = saida.linhas + size(T_iMais1,1);
        end
    end
    
  
     % Segunda restrição
    % Restrição para garantir que a Lyapunov sempre decai
     W = sdpvar(nx,nx,'symmetric');
     KChap = sdpvar(nu,nx);
     for i = 1:N
         U11 = S{i,quantInt};
         U21 = L{i,quantInt}';
         U22 = Z{i,quantInt};
         U31 = S{i,quantInt};
         U32 = L{i,quantInt};
         U33 = W;
    
         T =  [U11 U21' U31';
               U21 U22  U32';
               U31 U32  U33];
        LMIs = LMIs + (T >= 0);
        
        saida.linhas = saida.linhas + size(T,1);
       
    
        U11 = W;
        U21 = W;
        U22 = S{i,1};
        U31 = KChap;
        U32 = L{i,1}';
        U33 = Z{i,1};

        T =  [U11 U21' U31';
              U21 U22  U32';
              U31 U32  U33];
       LMIs = LMIs + (T >= 0);
       
       saida.linhas = saida.linhas + size(T,1);
   
     end
    
   % Terceira restrição
    for i = 1:N
        for j=1:quantInt
           U11 = S{i,j};
           U21 = L{i,j}';
           U22 = Z{i,j};
            T =  [U11 U21';
                  U21 U22];
          LMIs = LMIs + (T >= 0);  
          
          saida.linhas = saida.linhas + size(T,1);
        end
   end
  
        
    % Conta qtd de linhas LMIs, tempo de montagem e qtd de variáveis
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
        for i = 1:N
            for j=1:quantInt
				saida.S{i,j} = double(S{i,j});
				saida.L{i,j} = double(L{i,j});
				saida.Z{i,j} = double(Z{i,j});
			end
        end
        saida.W = double(W);
        saida.KChap = double(KChap);
        saida.K = saida.KChap*inv(saida.W);
        saida.mu = double(mu);
        saida.gamma = sqrt(saida.mu);
     end
    

end