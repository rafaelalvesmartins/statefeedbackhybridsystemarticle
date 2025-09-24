function saida = valEstHInfLMILabSemInt(ACal,ECal,CCal,KCal,h,delta,tol)
    % validarLMIHInfLab
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

    % Dimensões das matrizes de entrada
    nXi = length(ACal);
 
	% Inicialização das variáveis LMIs
    %Inicia o sistema de LMI
    setlmis([]);
	%Configura contador LMI
    c = 1;
	
    % Cria as matrizes X_i     
    quantInt = round(h/delta);
	X = cell(quantInt,1);
    sX = cell(quantInt,1);
    for i=1:quantInt
		[X{i},nVar,sX{i}] = lmivar(1,[nXi 1]);
		lmiterm([-c 1 1 X{i}],1,1); % X{i} > 0
		c = c + 1;
    end
	
	 % Função objetiva
	mu = lmivar(1,[1 0]);
    obj = mu;
    
    % Montar as restrições LMIs
    for i=1:(quantInt-1)
	
		lmiterm([c 1 1 X{i+1}],1/delta,1);
        lmiterm([c 1 1 X{i}],-1/delta,1);
        gerarLMIAlternadaDerivadaLMILabSemInt(ACal,ECal,CCal,X,mu,i,c);
		c = c + 1;
		
		lmiterm([c 1 1 X{i+1}],1/delta,1);
        lmiterm([c 1 1 X{i}],-1/delta,1);
        gerarLMIAlternadaDerivadaLMILabSemInt(ACal,ECal,CCal,X,mu,i+1,c);
		c = c + 1;
        
    end
    
    % Restrição para garantir que a Lyapunov sempre decai
	lmiterm([-c 1 1 X{quantInt}],1,1);
	lmiterm([-c 2 1 X{1}],1,KCal);
	lmiterm([-c 2 2 X{1}],1,1);
	c = c + 1;
    
	% Resolve LMI
    LMIs = getlmis;

    c = zeros(decnbr(LMIs),1);
    c(diag(decinfo(LMIs,obj))) = 1;
    options = [tol,2000,0,200,1];
    [copt,xopt] = mincx(LMIs,c,options);

    
    if isempty(xopt)
       saida.feas = 0;
    else
        saida.factivel = 1;
        for i=1:quantInt
            saida.X{i} = dec2mat(LMIs,xopt,X{i});
        end
        saida.mu = dec2mat(LMIs,xopt,mu);
        saida.gamma = sqrt(saida.mu);
    end
	

end