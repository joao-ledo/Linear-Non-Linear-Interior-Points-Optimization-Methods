%_________________________________________________________________________%
%                                                                         %
%                                                                         %
%                            DUAL AFIM ALGORITHM                          %
%                                                                         %
%                                                                         %
%                                          Developed by:                  % 
%                                                 Joao Augusto Silva Ledo %
%_________________________________________________________________________%

function result = DualAfim()
    clear all;
    clc;
    format long;      
    data = loadInputData(); % Carregar os dados de entrada.
    result = solveDualAfim(data.A, data.C, data.pontoInicial, data.Erro, data.e, data.Sinicial, data.b); % Cria a Instancia que resolve juntamento com os dados
end

%_________________________________________________________________________%
 function result = loadInputData() % Local onde se carregam as informacoes
        resultado.Erro = 10^-3;
        
        resultado.C{1} = [-2; -1];
        resultado.C{2} = ['x1';'x2'];
        resultado.A{1} = [-1, 1; 1, 0; 0, 1; -1, 0; 0, -1];    % <------------- Problema 1
        resultado.A{2} = ['x1';'x2'];
        resultado.b = [1; 3; 2; 0; 0];
        resultado.pontoInicial = [1; 1];
        resultado.Sinicial = [1; 2; 1; 1; 1];       
        resultado.pontoInicial = [2; 1];
        resultado.Sinicial = [2; 1; 1; 2; 1];
      
%         resultado.C{1} = [-2; 1];
%         resultado.C{2} = ['x1'; 'x2'];
%         resultado.A{1} = [1, -1; 0, 1; -1, 0; 0, -1];   % <------------- Problema 2
%         resultado.A{2} = ['x1';'x2'];
%         resultado.b = [15; 15; 0; 0];
%         resultado.pontoInicial = [8; 2];
%         resultado.Sinicial = [9; 13; 8; 2];
             
        resultado.e = controiE(resultado.Sinicial);
     result = resultado;
 end
%_________________________________________________________________________% 

 function result = solveDualAfim(A, C, pontoInicial, Erro, e, SInicial, b)
   k = 1;
   flag = false;
   x{k} = pontoInicial;
   s{k} = SInicial;
   Beta = {};
   while (flag == false)
           S{k} = constroiX(s{k});
           Dx{k} = DirecaoX(A{1}, S{k}, C{1});
           Ds{k} = DirecaoS(A{1}, Dx{k});
           FO{k} = calculoFO(x{k}, C{1});
           if((sum(abs(Ds{k})) <= Erro))
                   resultado.Nome = 'Algoritmo Dual Afim';
                   resultado.Beta = Beta;
                   resultado.W = W;
                   resultado.Matriz_S = S;
                   resultado.Dx = Dx;
                   resultado.Ds = Ds;
                   resultado.S = s;
                   resultado.FO = FO;
                   resultado.x = x;
                   resultado.Iteracoe = k-1;
                   flag = true;
           else
                if(verificanegatividade(Ds{k}) == false)
                    resultado.Nome = 'Algoritmo Dual Afim';
                    resultado.x = 'Problema Ilimitado';
                    resultado.Iteracoe = k-1;
                    flag = true;
                else
                    W{k} = EstimativaDual(S{k}, Ds{k});
                    if((verificanegatividade(W{k}) == false) && (calculaErroDual(e, S{k}, W{k}) <= Erro))
                        resultado.Nome = 'Algoritmo Dual Afim';
                        resultado.Beta = Beta;
                        resultado.W = W;
                        resultado.Matriz_S = S;
                        resultado.Dx = Dx;
                        resultado.Ds = Ds;
                        resultado.S = s;
                        resultado.FO = FO;
                        resultado.x = x;
                        resultado.Iteracoes = k-1;
                        flag = true;
                    else if((verificanegatividade(W{k}) == false) && (calculaErroPrimal(x{k}, b, W{k}) <= Erro))
                            resultado.Nome = 'Algoritmo Dual Afim';
                            resultado.Beta = Beta;
                            resultado.W = W;
                            resultado.Matriz_S = S;
                            resultado.Dx = Dx;
                            resultado.Ds = Ds;
                            resultado.S = s;
                            resultado.FO = FO;
                            resultado.x = x;
                            resultado.Iteracoes = k-1;
                            flag = true;                            
                        else if((verificanegatividade(W{k}) == false) && (VerificaErro(S{k}, Erro) == false))
                                resultado.Nome = 'Algoritmo Dual Afim';
                                resultado.Beta = Beta;
                                resultado.W = W;
                                resultado.Matriz_S = S;
                                resultado.Dx = Dx;
                                resultado.Ds = Ds;
                                resultado.S = s;
                                resultado.FO = FO;
                                resultado.x = x;
                                resultado.Iteracoes = k-1;
                                flag = true;                                       
                            else
                                Beta{k} = passo(s{k}, Ds{k});
                                x{k+1} = NextValue(x{k}, Beta{k}, Dx{k});
                                s{k+1} = NextValue(s{k}, Beta{k}, Ds{k});
                                k = k + 1;
                            end
                        end
                    end   
                end    
           end      
       end
   result = resultado;
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

function result = calculaErroDual(e, S, w)
    result = e' * S * w;
end

function result = calculaErroPrimal(x, b, w)
    result = x - (b' * w);
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
 
 function result = DirecaoX(A, S, c)
    result = - inv(A' * S^(-2) * A) * c;
 end
 
  function result = DirecaoS(A, Dx)
    result = - A * Dx;
 end
 
 function result = EstimativaDual(S, Ds)        
    result = S^(-2) * Ds;
 end
 
 function result = CustoRelativo(C, A, W)
    result = C{1}-(A{1}'* W);
 end
 
 function resultado = NextValue(valoratual, passo, direcao)
   resultado = valoratual + passo * direcao;
end
 
 function result = constroiX(x)
    for i = 1 : length(x)
        for j = 1 : length(x)
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

 function result = VerificaErro(s, Erro)
    flag = false;
    [linha, coluna] = size(s);
    for i = 1 : linha
        for j = 1 : coluna
            if(s(i,j) > Erro)
                flag = true;
            end
        end
    end
    result = flag;
 end