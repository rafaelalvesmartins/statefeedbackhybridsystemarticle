function saida = estHInfAnaPolyLMILabOriginal(poly,KCal,h,delta,tol)
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
    nxi = length(APoly{1});
    nwi = size(EPoly{1},2);
    nyi = size(CPoly{1},1);
    N = length(APoly);
    
    % Inicialização das variáveis LMIs
    %Inicia o sistema de LMI
    setlmis([]);
	%Configura contador LMI
    c = 1;
    
    % Cria a matriz X 
    quantInt = round(h/delta);
    Y = cell(quantInt,1);
    sX = cell(quantInt,1);
    
    for i=1:quantInt
		[Y{i},nVar,sX{i}] = lmivar(1,[nxi 1]);
    end

  
    
    
    
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
            lmiterm([c 1 1 Y{j+1}],-1/delta,1);
            lmiterm([c 1 1 Y{j}],1/delta,1);
            gerarAnaLMIAlternadaDerivadaPolyLMILabOriginal(APoly{i},EPoly{i},CPoly{i},Y{j},mu,c);
            c = c + 1;

            % Para j+1
            lmiterm([c 1 1 Y{j+1}],-1/delta,1);
            lmiterm([c 1 1 Y{j}],1/delta,1);
            gerarAnaLMIAlternadaDerivadaPolyLMILabOriginal(APoly{i},EPoly{i},CPoly{i},Y{j+1},mu,c);
            c = c + 1;
        end
    end
    
    % Segunda restrição
    % Restrição para garantir que a Lyapunov sempre decai
    lmiterm([-c 1 1 Y{quantInt}],1,1);
	lmiterm([-c 2 1 Y{quantInt}],KCal,1);
	lmiterm([-c 2 2 Y{1}],1,1);
    
    % Resolve LMI
    LMIs = getlmis;
    
           
  
    Nm = decnbr(LMIs);
    Cz = zeros(Nm,1);
    for j=1:Nm
        Miz = defcx(LMIs,j,obj);
        Cz(j) = Cz(j)+ Miz;
    end
%     c = zeros(decnbr(LMIs),1);
%     c(diag(decinfo(LMIs,obj))) = 1;
%     options = [tol,2000,0,200,1];
%     [copt,xopt] = mincx(LMIs,c,options);
    options = [tol,2000,1e9,10,1];
    [copt,xopt] = mincx(LMIs,Cz,options);
    
    if isempty(xopt)
       saida.feas = 0;
    else
        saida.factivel = 1;
        for i=1:quantInt
            saida.Y{i} = dec2mat(LMIs,xopt,Y{i});
        end
       
        saida.mu = dec2mat(LMIs,xopt,mu);
        saida.gamma = sqrt(saida.mu);
    end
end