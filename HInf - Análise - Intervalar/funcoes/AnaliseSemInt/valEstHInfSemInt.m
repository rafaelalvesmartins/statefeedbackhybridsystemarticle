function saida = valEstHInfSemInt(ACal,ECal,CCal,KCal,h,delta,tol)
    % validarLMIHInf
    %
    % Avaliar norma HInf de um sistema híbrido da tese do Matheus
    %
    % input:  ACal,ECal,CCal,KCal -> Matrizes aumentadas
    %         h -> período de amostragem
    %         delta -> relacionada a aproximação da derivada
    %         tol -> tolerância do solver
    %
    % output: saida.X -> Matriz Lyapunov
    %
    % Criado date: 21 Mar 2020
    % Revisão date: 21 Mar 2020
    % Autor: rafaelmartinsalves@gmail.com

    % Dimensões das matrizes de entrada
    nXi = length(ACal);
    nW = size(ECal, 2);
    
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
    end
    
    % Função objetiva
    mu = sdpvar(1,1);
    obj = mu;
    
    % Montar as restrições LMIs
    for i=1:(quantInt-1)
        der = (X{i+1}-X{i})/delta;
   
        T_i = gerarLMIAlternadaDerivadaSemInt(ACal,ECal,CCal,der,X,mu,nW,i);
        T_iMais1 = gerarLMIAlternadaDerivadaSemInt(ACal,ECal,CCal,der,X,mu,nW,i+1);
       
            
        
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
    sol = optimize(LMIs,obj,sdpsettings('verbose',0,'solver','sedumi'));
    saida.cpusec = sol.solvertime;
    
    % Checar a solução
    p = min(checkset(LMIs));
    saida.delta = p;
    saida.factivel = 0;
    
    if(p > -tol ) % Factível
        saida.factivel = 1;
        for i=1:quantInt
            saida.X{i} = double(X{i});
        end
        saida.mu = double(mu);
        saida.gamma = sqrt(saida.mu);
    end
    

end