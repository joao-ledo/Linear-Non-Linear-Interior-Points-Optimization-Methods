%_________________________________________________________________________%
%                                                                         %
%                                                                         %
%                 DUAL AFFINE SCALING INTERIOR-POINTS METHOD              %
%         (polynomial quadratic/linear order with Central Trajectory)     % 
%                                                                         %
%                                                                         %
%                                                                         %
%                                          Developed by:                  % 
%                                                 Joao Augusto Silva Ledo %
%_________________________________________________________________________%

function result = DualAfimPPQTragetoriaCentral()
    clear all;
    clc;
    format long;      
    data = loadInputData(); % Load the input data
    result = solvePrimalAfimPPQ(data.A, data.C, data.C_quadratic, data.pontoInicial, data.Erro, data.b, data.Q, data.Sinicial, data.e, data.Mi, data.answer); % Creates the instances and solve it by loading the input data
end

%_________________________________________________________________________%
 function result = loadInputData() % Load the input data function
    resultado.Erro = 10^-1;
    resultado.Mi = 1;
    
   resultado.answer = 2;
    resultado.C = [4; -6];
    resultado.C_quadratic = [2; 3];
    resultado.A = [-1, 1; -1, -1;1, 0;0, 1; -1, 0; 0, -1]; %<------------- Problem 2
    resultado.b = [1; -1; 3; 2; 0; 0]; 
    resultado.Q = [4, 0; 0, 6];
    resultado.pontoInicial{1} = inv(-resultado.Q) * resultado.C;
    resultado.pontoInicial{2} = [2.5; 1];% [1.5; 1.5]; %[1; 1]; %%               [2.5; 1];
    resultado.Sinicial = [2.5; 2.5; 0.5; 1; 2.5; 1];%[1; 2; 1.5; 0.5; 1.5; 1.5];%[1; 1; 2; 1; 1; 1]; %           [2.5; 2.5; 0.5; 1; 2.5; 1];
    resultado.e = controiE(resultado.Sinicial); 

%     resultado.answer = 3;
%     resultado.C = [-2; -1];
%     resultado.C_quadratic = [2; 3];
%     resultado.A = [-1, 1; 1, 0; 0, 1; -1, 0; 0, -1]; %<------------- Problem 2
%     resultado.b = [1; 3; 2; 0; 0]; 
%     resultado.Q = [0, 0; 0, 0];
%     resultado.pontoInicial{1} = inv(-resultado.Q) * resultado.C;
%     resultado.pontoInicial{2} = [1; 1];% [1.5; 1.5]; %[1; 1]; %%               [2.5; 1];
%     resultado.Sinicial = [1; 2; 1; 1; 1];%[1; 2; 1.5; 0.5; 1.5; 1.5];%[1; 1; 2; 1; 1; 1]; %           [2.5; 2.5; 0.5; 1; 2.5; 1];
%     resultado.e = controiE(resultado.Sinicial); 
    
    result = resultado;
 end
%_________________________________________________________________________% 

 function result = solvePrimalAfimPPQ(A, C, C_quadratic, pontoInicial, Erro, b, Q_ini, Sinicial, e, Mi_ini, answer)
   k = 1;
   flag = false;
   x{k} = pontoInicial{1};
   s{k} = Sinicial;
   w{k} = {};
%    answer = input('Type 1 if wants to solve using polynomial quadratic order or 2 for using N-dimensional polynomial quadratic order method: \n');
%    answer = VerificaOpcao(answer);
%    VerificaOpcao(answer);
   Fo{k} = calculoFO(x{k}, C, C_quadratic);
   if((VerificaEquacao(A, x{k}, b) == true) && (verificanegatividade(x{k}) == false))
       resultado.Nome = VerificaNome(answer);
       resultado.Fo = Fo;
       resultado.x = x;
       resultado.Iteracoe = k-1;
       flag = true;
   else
       x{k} = pontoInicial{2};
   end
   Mi{k} = Mi_ini;
   while(flag == false)        
        Fo{k} = calculoFO(x{k}, C, C_quadratic);
        Q{k} = Q_ini;%Hessian(x{k});%ChooseMethod(Q_ini, x{k});
        X{k} = diag(x{k});
        S{k} = diag(s{k});
        Dx{k} = DirecaoX(A, S{k}, Q{k}, x{k}, C, Mi{k}, answer, e);
        Ds{k} = DirecaoS(A, Dx{k});
        if(VerificaMenorErro(Ds{k}, Erro) == true)
            resultado.Nome = VerificaNome(answer);
            resultado.Fo = Fo;
            resultado.Q = Q;
            resultado.S = S;
            resultado.X = X;
            resultado.x = x;
            resultado.s = s;
            resultado.Dx = Dx;
            resultado.Ds = Ds;
            resultado.w = w;
            resultado.Mi = Mi;
            resultado.Iteracoe = k-1;
            flag = true;
        else         
            if((VerificaPositividade(Ds{k}) == true))
                resultado.Nome = VerificaNome(answer);
                resultado.x = 'Problema Ilimitado';
                resultado.Iteracoe = k-1;
                flag = true;
            else
                w{k} = Calcula_W(S{k}, Ds{k});
                if((verificanegatividade(w{k}) == false) || (VerificaErro(e, S{k}, w{k}, Erro) == true))
                    resultado.Nome = VerificaNome(answer);
                    resultado.Fo = Fo;
                    resultado.Q = Q;
                    resultado.S = S;
                    resultado.x = x;
                    resultado.X = X;
                    resultado.Matriz_X = X;
                    resultado.s = s;
                    resultado.Dx = Dx;
                    resultado.Ds = Ds;
                    resultado.w = w;
                    resultado.Mi = Mi;
                    resultado.Iteracoe = k-1;
                    flag = true;
                else
                    beta{k} = passo(x{k}, Dx{k}, Q{k}, C, s{k}, Ds{k});
                    Mi{k+1} = ProximaBarreira(Mi{k});
                    x{k+1} = NextValue(x{k}, beta{k}, Dx{k});
                    s{k+1} = NextValue(s{k}, beta{k}, Ds{k});
                    k = k + 1;
                end
            end
        end
   end    
   result = resultado;
 end
 
 function result = VerificaNome(answer)
    if(answer == 1)
        resposta = 'Algoritmo Dual Afim (PPQ)';
    else if (answer == 2)
            resposta = 'Algoritmo Dual Afim (PPQ) Trajetoria Central';
        else
            if(answer == 3)
                resposta = 'Algoritmo Dual Afim (PPL) Trajetoria Central';
            end
        end
    end
    result = resposta;
end
 
 function result = VerificaOpcao(answer)
    while((answer ~= 1) && (answer ~= 2))
       answer = input('Favor Escolher corretamente o metodo pelo qual deseja Executar o Algoritmo \nDigite 1 se deseja resolver pelo metodo de polinomio ordem quadratica ou 2 por metodo de polinomio de ordem N-dimensional: \n');
    end
    result = answer;
 end
 
 function result = ChooseMethod(answer, Q, x)
    if(answer == 1)
        resposta = Q;
    else if (answer == 2)
            resposta = Hessian(x);
        end
    end
    result = resposta;
 end
 
function result = ProximaBarreira(Mi)
    result = Mi/2;
end
 
 function result = Hessian(x)
    result = [12*x(1)^2 - 96*x(1)+194, -8; -8, 32];
 end
 
 function result = Grad(x)
    result = [4*(x(1)^3) - 48*(x(1)^2) + 194*x(1)-8*x(2)-256; 32*x(2)-8*x(1)];
 end
 
function result = controiE(x)
    for i = 1 : length(x)
       aux(i) = 1; 
    end
    result = aux';
end
 
 function result = Calcula_W(S, Ds) 
    result = inv(S^2) * Ds;
 end
 
 function result = DirecaoX(A, S, Q, x, C, Mi, answer, e)
    if(answer == 1)
        resposta = - inv(A' * inv(S^2) * A) * (Q * x + C);
    else if (answer == 2)
            resposta = - inv(A' * inv(S^2) * A) * (Q * x + C + (Mi * A'*inv(S)*e));
        else
            if(answer == 3)
                resposta = - inv(A' * inv(S^2) * A) * (C + (Mi * A'*inv(S)*e));
            end
        end
    end
    result = resposta;
 end
 
 function result = DirecaoS(A, Dx)
    result = -A * Dx;
 end
 
 function result = calculoFO(X, C, C_quadratic)
     for i = 1 : length(X)
        fo(i) = X(i)*C(i) + (X(i)^2)*C_quadratic(i);
     end
    result = sum(fo);
 end
 
%  function result = calculoFO(x)
%     result = ((x(1)-4)^4) + (x(1)-4*x(2))^2;
%  end
 
 function result = passo(x, Dx, Q, c, s, Ds)
     for i = 1 : length(s)
        if Ds(i) < 0 
            alpha(i) = -0.9995*s(i)/Ds(i);
        else
            alpha(i) = 0;
        end
     end    
     beta1 = min(nonzeros(alpha));    
     beta2 = - (Dx' * (Q*x + c))/(Dx' * Q  * Dx);
     if(sum(alpha) == 0)
         resposta = beta2;
     else
         resposta = 0.9995*min(beta1, beta2);
     end
    result = resposta;
 end
  
 function resultado = NextValue(valoratual, passo, direcao)
   resultado = valoratual + passo * direcao;
 end
 
 function result = verificanegatividade(valor)
    flag = false;
    for i = 1 : length(valor)
       if (valor(i) < 0)
           flag = true;
       end
    end
    result = flag;
 end
 
 function result = VerificaPositividade(valor)
    flag = true;
    for i = 1 : length(valor)
       if(valor(i) <= 0) 
           flag = false;
       end
    end
    result = flag;
 end

 function result = VerificaEquacao(A, x, b)
     aux1 = A*x;
     flag = true;
    for i = 1 : length(aux1)
       if(aux1(i) ~= b(i))
           flag = false;
       end
    end
    result = flag;
 end
 
 function result = VerificaErro(e, S, w, Erro)     
    flag = true;
    if(abs(e'*S*w) > Erro)
        flag = false;
    end
    result = flag;
 end
 
 function result = VerificaMenorErro(Dx, Erro)
    flag = true;
    for i = 1 : length(Dx)
       if(Dx(i) > Erro)
           flag = false;
       end
    end
    result = flag;
 end
