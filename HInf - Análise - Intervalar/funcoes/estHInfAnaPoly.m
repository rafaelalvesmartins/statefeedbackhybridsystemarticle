function saida = estHInfAnaPoly(poly,KCal,h,delta,tol)
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
    APoly = poly.APoly;
    EPoly = poly.EPoly;
    CPoly = poly.CPoly;
%     poly.D1Poly = D1Poly;
    
%     Tamanho das matrizes
    nXi = length(APoly{1});
    nwi = size(EPoly{1},2);
    nyi = size(CPoly{1},1);
    N = length(APoly);
    
    % Inicialização das variáveis LMIs
    saida.cpuSeg_m = clock;
    %LMI contador de linhas
    saida.linhas = 0;
    %Inicia o sistema de LMI
    LMIs = set([]);
    
    % Cria a matriz X 
    quantInt = round(h/delta);
    Y = cell(N,quantInt);
    
     for i=1:N
        for j = 1:quantInt
            Y{i,j} = sdpvar(nXi,nXi,'symmetric');
            saida.linhas = saida.linhas + size(Y{i,j},1);
        end
     end

    % Função objetiva
    mu = sdpvar(1,1);
    obj = mu;
	
    % Cria as matrizes Yi
    X = cell(4,1);
    X{1} = sdpvar(nXi,nXi);
	X{2} = sdpvar(nwi,nXi);
	X{3} = sdpvar(nyi,nXi);
	X{4} = sdpvar(nXi,nXi);
    
  
    % Montar as restrições LMIs
    % Primeira restrição
    for i = 1:N
        for j=1:(quantInt-1)
            % Para j
            % Derivadas
			der = (Y{i,j+1}-Y{i,j})/delta;
            T_j = gerarAnaLMIAlternadaDerivadaPoly(APoly{i},EPoly{i},CPoly{i},der,Y{i,j},X,mu);
			T_jMais1 = gerarAnaLMIAlternadaDerivadaPoly(APoly{i},EPoly{i},CPoly{i},der,Y{i,j+1},X,mu);
            
			LMIs = LMIs + (T_j <= 0);
			LMIs = LMIs + (T_jMais1 <= 0);
			
			% Conta quantidade de linhas LMIs
			saida.linhas = saida.linhas + size(T_j,1);
			saida.linhas = saida.linhas + size(T_jMais1,1);

        end
    end
    
    for i = 1:N
         % Restrição para garantir que a Lyapunov sempre decai
        U11 = Y{i,quantInt};
        U21 = Y{i,1}*KCal;
        U22 = Y{i,1};
        T =  [U11 U21';
          U21 U22];
    
        LMIs = LMIs + (T >= 0);
        
        % Conta quantidade de linhas LMIs
        saida.linhas = saida.linhas + size(T,1);
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
        for i=1:quantInt
            saida.Y{i} = double(Y{i});
        end
        
        for i=1:N
             for j = 1:quantInt
                saida.Y{i,j} = double(Y{i,j});
             end
        end
        
		for i=1:4
            saida.X{i} =  double(X{i});
        end
		
        saida.mu = double(mu);
        saida.gamma = sqrt(saida.mu);
    end
end