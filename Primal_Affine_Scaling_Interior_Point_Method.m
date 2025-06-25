%_________________________________________________________________________%
%                                                                         %
%                                                                         %
%                 PRIMAL AFFINE SCALING INTERIOR-POINT METHOD             %
%                     (with or without Central Trajectory)                %
%                                                                         %
%                                          Developed by:                  % 
%                                                 Joao Augusto Silva Ledo %
%_________________________________________________________________________%

function result = PrimalAfim()
    clear all;
    clc;
    format long;      
    data = loadInputData(); % Load the input data
    result = solvePrimalAfim(data.A, data.C, data.pontoInicial, data.Erro, data.e, data.answer, data.Mi); % Cria a Instancia que resolve juntamento com os dados
end

%_________________________________________________________________________%
 function result = loadInputData() % Load the input data function
         resultado.Erro = 10^-3;
         
%         resultado.answer = 1; % Without Central Trajectory
%         resultado.C{1} = [-2; 1; 0; 0];
%         resultado.A{1} = [1, -1, 1, 0; 1, 1, 0, 1];     % <------------- Problem 1
%         resultado.pontoInicial = [2,1,1,3];
%         resultado.Mi = 0;

%         resultado.answer = 2; % With Central Trajectory
%         resultado.C{1} = [-2; 1; 0; 0];
%         resultado.A{1} = [1, -1, 1, 0; 1, 1, 0, 1];     % <------------- Problem 1
%         resultado.pontoInicial = [2,1,1,3];
%         resultado.Mi = 1;

         resultado.answer = 1; % Without Central Trajectory
         resultado.C{1} = [-2; 1; 0; 0];
         resultado.A{1} = [1, -1, 1, 0; 0, 1, 0, 1];   % <------------- Problem 2
         resultado.pontoInicial = [8; 2; 9; 13];
         resultado.Mi = 0;
         
%          resultado.answer = 2; % With Central Trajectory
%          resultado.C{1} = [-2; 1; 0; 0];
%          resultado.A{1} = [1, -1, 1, 0; 0, 1, 0, 1];   % <------------- Problem 2
%          resultado.pontoInicial = [8; 2; 9; 13];
%          resultado.Mi = 1;

        resultado.e = controiE(resultado.pontoInicial);
     result = resultado;
 end
%_________________________________________________________________________% 

 function result = solvePrimalAfim(A, C, pontoInicial, Erro, e, answer, Mi_ini)
   k = 1;
   flag = false;
   Mi{k} = Mi_ini;
   y{k} = e;
   x{k} = pontoInicial;
   Alpha = {};
   while(flag == false)      
        X{k} = constroiX(A{1},x{k});
        W{k} = EstimativaDual(A, X{k}, C, answer, Mi{k}, e);
        r{k} = CustoRelativo(C, A, X{k}, W{k}, answer, Mi{k}, e);
        erro{k} = calculaErro(e, X{k}, r{k});
        FO{k} = calculoFO(x{k}, C{1});
           if((verificanegatividade(r{k}) == false) && (calculaErro(e, X{k}, r{k}) <= Erro) )
                   resultado.Nome = VerificaNome(answer);
                   resultado.Alpha = Alpha;
                   resultado.W = W;
                   resultado.R = r;
                   resultado.Erro = erro;
                   resultado.Direcao = Dy;
                   resultado.MatrizX = X;
                   resultado.FO = FO;
                   resultado.x = x;
                   resultado.Mi = Mi;
                   resultado.y = y;
                   resultado.X = X;
                   resultado.Iteracoe = k-1;
                   flag = true;
           else
                Dy{k} = Direcao(X{k},r{k});
                if(verificanegatividade(Dy{k}) == false)
                    resultado.Nome = VerificaNome(answer);
                    resultado.x = 'Unlimited Problem';
                    resultado.Mi = Mi;
                    resultado.Iteracoe = k-1;
                    flag = true;
                else if(sum(abs(Dy{k})) <= Erro)
                        resultado.Nome = VerificaNome(answer);
                        resultado.Alpha = Alpha;
                        resultado.W = W;
                        resultado.R = r;
                        resultado.Erro = erro;
                        resultado.Direcao = Dy;
                        resultado.MatrizX = X;
                        resultado.FO = FO;
                        resultado.x = x;
                        resultado.Mi = Mi;
                        resultado.y = y;
                        resultado.X = X;
                        resultado.Iteracoes = k-1;
                        flag = true;
                    else
                        Alpha{k} = passo(e, Dy{k});                   
                        Mi{k+1} = ProximaBarreira(Mi{k}, answer);                    
                        y{k+1} = NextY(Alpha{k},Dy{k},e);
                        x{k+1} = NextX(X{k},y{k+1});
                        k = k + 1;
                    end   
                end    
           end      
       end
   result = resultado;
 end
 
 function result = VerificaNome(answer)
    if(answer == 1)
        resposta = 'Primal Affine Scaling Method';
    else if (answer == 2)
            resposta = 'Primal Affine Scaling Method with Central Trajectory';
        end
    end
    result = resposta;
 end
 
 function result = calculoFO(X, C)
 for i = 1 : length(X)
    fo(i) = X(i)*C(i);
 end
    result = sum(fo);
 end
 
function result = controiE(x)
    for i = 1 : length(x)
       aux(i) = 1; 
    end
    result = aux';
end
 
function result = calculaErro(e, X, r)
    result = e' * X * r;
end
 
 function result = passo(e, Dy)   
     for i = 1 : length(e)
        if Dy(i) < 0 
            alpha(i) = 0.9995*e(i)/-Dy(i);
        else
            alpha(i) = 0;
        end
     end    
     alpha = nonzeros(alpha)';
    result = min(alpha);
 end
 
 function result = Direcao(X,r)
    result = - X*r;
 end
 
 function result = EstimativaDual(A, X, C, answer, Mi, e)
    if(answer == 1)
        resposta = (inv(A{1}*power(X,2)*A{1}')*A{1}*power(X,2))*C{1};
    else if (answer == 2)
            resposta = (inv(A{1}*power(X,2)*A{1}')*A{1}*power(X,2))*(C{1}-Mi*inv(X)*e);
        end
    end
    result =  resposta;
 end
 
 function result = CustoRelativo(C, A, X, W, answer, Mi, e)
   if(answer == 1)
        resposta = C{1}-(A{1}'* W);
    else if (answer == 2)
            resposta = (C{1}-Mi*inv(X)*e) - (A{1}'* W);
        end
    end
    result = resposta;
 end
 
 function resultado = NextY(alfa,dy,e)
   resultado = e+(alfa*dy);
end

function resultado = NextX(X,y)
   resultado = X*y;
end
 
 function result = constroiX(z,x)
    for i = 1 : length(z)
        for j = 1 : length(z)
            if i == j
                DiagX(i,j) = x(i);
            end
        end
    end
     result = DiagX;
 end
 
 function result = verificanegatividade(S)
    flag = false;
    for i = 1 : length(S)
       if (S(i) < 0)
           flag = true;
       end
    end
    result = flag;
 end

 function result = ProximaBarreira(Mi, answer)
   if(answer == 1)
        resposta = Mi;
    else if (answer == 2)
            resposta = Mi/2;
        end
    end
    result = resposta;
end
