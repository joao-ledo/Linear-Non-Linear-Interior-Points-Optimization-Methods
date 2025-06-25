%_________________________________________________________________________%
%                                                                         %
%                                                                         %
%                PRIMAL AFFINE SCALING INTERIOR-POINTS METHOD             %
%              (polynomial quadratic order - Central Trajectory)          % 
%                     (With or without Rosenbrock Search)                 % 
%                     (For the Economic Dispatch Problem)                 %  
%                                                                         %
%                                          Developed by:                  % 
%                                                 Joao Augusto Silva Ledo %
%_________________________________________________________________________%

function result = PrimalAfimPPQTragetoriaCentral()
    clear all;
    clc;
    format long;      
    data = loadInputData(); % Load the input data
    result = SolvePrimalAfimPPQTragetoriaCentral(data.A, data.C, data.C_quadratic, data.pontoInicial, data.Erro, data.b, data.Q, data.Mi, data.e, data.answer); % Creates the instances and solve it by loading the input data
end

%_________________________________________________________________________%
 function result = loadInputData() % Load input data function
    resultado.Erro{1} = 10^-4;
    resultado.Erro{2} = 10^-5;
    resultado.Erro{3} = 10^-3;
    resultado.Mi = 1;
      
%     resultado.answer = 1; % Without Rosenbrock
%     resultado.C = [1; 2; -3];
%     resultado.C_quadratic = [2; 3; 5];
%     resultado.A = [1, 1, 0; 0, 1, 1];
%     resultado.b = [5; 10];
%     resultado.Q = [4, 0, 0; 0, 6, 0; 0, 0, 10];
%     resultado.pontoInicial{1} = inv(-resultado.Q) * resultado.C;
%     resultado.pontoInicial{2} = [3; 2; 8];
    
    resultado.answer = 1; % Without Rosenbrock
    resultado.C = [7.92; 7.97; 7.85; 0; 0; 0; 0; 0; 0];
    resultado.C_quadratic = [0.001562; 0.00482; 0.001940; 0; 0; 0; 0; 0; 0];
    resultado.A = [1, 1, 1, 0, 0, 0, 0, 0, 0;                                   %<------------- Economic Dispatch Problem with 3 generators
                   -1, 0, 0, 1, 0, 0, 0, 0, 0;
                   1, 0, 0, 0, 1, 0, 0, 0, 0;
                   0, -1, 0, 0, 0, 1, 0, 0, 0;
                   0, 1, 0, 0, 0, 0, 1, 0, 0;
                   0, 0, -1, 0, 0, 0, 0, 1, 0;
                   0, 0, 1, 0, 0, 0, 0, 0, 1];
    resultado.b = [850; -100; 600; -50; 200; -100; 400];
    resultado.Q = [0.003124, 0, 0, 0, 0, 0, 0, 0, 0; 
                   0, 0.00964, 0, 0, 0, 0, 0, 0, 0;
                   0, 0, 0.00388, 0, 0, 0, 0, 0, 0;
                   0, 0, 0, 0, 0, 0, 0, 0, 0;
                   0, 0, 0, 0, 0, 0, 0, 0, 0;     %13X18
                   0, 0, 0, 0, 0, 0, 0, 0, 0;
                   0, 0, 0, 0, 0, 0, 0, 0, 0;
                   0, 0, 0, 0, 0, 0, 0, 0, 0;
                   0, 0, 0, 0, 0, 0, 0, 0, 0];
    resultado.pontoInicial{1} = inv(-resultado.Q) * resultado.C;
    resultado.pontoInicial{2} = [500; 150; 200; 400; 100; 100; 50; 100; 200];
    
%     resultado.answer = 2; % With Rosenbrock
%     resultado.C = [-2; 0; 0; 0; 0];
%     resultado.C_quadratic = [0; 0; 0];
%     resultado.A = [1, -1, 1, 0, 0; 0, -1, 0, 1, 0; 0, 1, 0, 0, 1];  
%     resultado.b = [15; -2; 15];
%     resultado.pontoInicial{2} = [2; 4; 17; 2; 11];%[10; 4; 9; 2; 11];%[10; 3; 8; 1; 12];
%     resultado.Q = Hessian(resultado.pontoInicial{2});
%     resultado.pontoInicial{1} = inv(-resultado.Q) * resultado.C;

    resultado.e = constroeE(resultado.pontoInicial{2});
    result = resultado;
 end
%_________________________________________________________________________% 

function result = SolvePrimalAfimPPQTragetoriaCentral(A, C, C_quadratic, pontoInicial, Error, b, Q_ini, Mi_ini, e, answer)
   Erro = Error{1};
   ErroPrimal = Error{2};
   ErroDual = Error{3};
   k = 1;
   flag = false;
   alpha{k} = {};
   x{k} = pontoInicial{1};
   Mi{k} = Mi_ini;
   Fo{k} = calculoFO(x{k}, C, C_quadratic, answer);
   if((VerificaEquacao(A, x{k}, b) == true) && (verificanegatividade(x{k}) == false))
       resultado.Nome = VerificaNome(answer);
       resultado.Fo = Fo;
       resultado.x = x;
       resultado.Iteracoe = k-1;
       flag = true;
   else
       x{k} = pontoInicial{2};
   end
   while(flag == false)
        Fo{k} = calculoFO(x{k}, C, C_quadratic, answer);
        X{k} = diag(x{k});
        Q{k} = CalculaQ(Q_ini, x{k}, answer);
        H{k} = Calcula_H(Q{k}, X{k}, Mi{k});
        w{k} = Calcula_W(A, H{k}, Q{k}, x{k}, C, Mi{k}, X{k}, e, answer);
        s{k} = Calcula_S(Q{k}, x{k}, C, A, w{k}, Mi{k}, X{k}, e, answer);
        if((VerificaErro(x{k}, s{k}, Erro) == true)||(VerificaPositividade(x{k})== true) && (VerificaPositividade(s{k})== true) && ((factbilidadePrimal(A, x{k}, b, ErroPrimal) == true) && (factibilidadeDual(s{k}, Q{k}, x{k}, C, ErroDual) == true)))
            resultado.Nome = VerificaNome(answer);
            resultado.Fo = Fo;
            resultado.Q = Q;
            resultado.x = x;
            resultado.Matriz_X = X;
            resultado.H = H;
            resultado.w = w;
            resultado.s = s;
            resultado.Alpha = alpha;
            resultado.Mi = Mi;
            resultado.Iteracoe = k-1;
            flag = true;
        else
            Dx{k} = Direcao(H{k}, s{k});
            if((VerificaPositividade(Dx{k}) == true))
                resultado.Nome = VerificaNome(answer);
                resultado.x = 'Problema Ilimitado';
                resultado.Iteracoe = k-1;
                flag = true;
            else
                if((VerificaMenorErro(Dx{k}, Erro) == true) && ((VerificaPositividade(x{k}) == true) || (VerificaPositividade(s{k}) == true)))
                    resultado.Nome = VerificaNome(answer);
                    resultado.Fo = Fo;
                    resultado.Q = Q;
                    resultado.x = x;
                    resultado.Matriz_X = X;
                    resultado.H = H;
                    resultado.w = w;
                    resultado.s = s;
                    resultado.Dx = Dx;
                    resultado.Alpha = alpha;
                    resultado.Mi = Mi;
                    resultado.Iteracoe = k-1;
                    flag = true;
                else
                    alpha{k} = passo(x{k}, Dx{k}, Q{k}, C, answer);
                    Mi{k+1} = ProximaBarreira(Mi{k});
                    x{k+1} = NextValue(x{k}, alpha{k}, Dx{k});
                    k = k + 1;
                end
            end
        end
   end
    result = resultado;
end

function result = VerificaNome(answer)
    if(answer == 1)
        resposta = 'Primal Affine Scalling polynomial quadratic order method with Central Trajectory';
    else if (answer == 2)
            resposta = 'Primal Affine Scalling polynomial quadratic order method with Rosenbrock';
        end
    end
    result = resposta;
end

function result = constroeE(x)
    for i = 1 : length(x)
       aux(i) = 1; 
    end
    result = aux';
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

function result = Calcula_H(Q, X, Mi)
   result = inv(Q + (Mi + 1)*inv(X^2));
end

function result = Calcula_W(A, H, Q, x, c, Mi, X, e, answer)
    if(answer == 1)
        resposta = inv(A*H*A')*A*H*(Q*x+c-Mi*inv(X)*e);
    else if (answer == 2)
            resposta = inv(A*H*A')*A*H*(Grad(x)-Mi*inv(X)*e);
        end
    end
   result = resposta;  
end

function result = Calcula_S(Q, x, c, A, w, Mi, X, e, answer)
    if(answer == 1)
        resposta = (Q*x + c - Mi*inv(X)*e) - A'*w;
    else if (answer == 2)
            resposta = (Grad(x) - Mi*inv(X)*e) - A'*w;
        end
    end
    result = resposta;
end

function result = VerificaErro(x, s, Erro)     
    flag = true;
    if(x'*s > Erro)
        flag = false;
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

function result = Direcao(H, s)
    result = -H*s;
end

function result = CalculaQ(Q, x, answer)
    if(answer == 1)
        resposta = Q;
    else if (answer == 2)
            resposta = Hessian(x);
        end
    end
    result = resposta;
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

function result = passo(x, Dx, Q, c, answer)
     for i = 1 : length(x)
        if Dx(i) < 0 
            alpha(i) = -0.9995*x(i)/Dx(i);
        else
            alpha(i) = 0;
        end
     end    
     alpha1 = min(nonzeros(alpha));
    if(answer == 1)
        alpha2 = - (Dx' * (Q*x + c))/(Dx' * Q  * Dx);
    else if (answer == 2)
            alpha2 = - (Dx' * (Grad(x)))/(Dx' * Q  * Dx);
        end
    end    
    if((sum(alpha) == 0))
        resposta = alpha2;
    else
        resposta = 0.9995*min(alpha1, alpha2);
    end
    result = resposta;
end

function result = ProximaBarreira(Mi)
    result = Mi/2;
end

function resultado = NextValue(valoratual, passo, direcao)
   resultado = valoratual + passo * direcao;
end

function result = factbilidadePrimal(A, x, b, Erro) 
    if((VerificaErroMatriz(abs(A*x - b)/(abs(b)+1), Erro) == true) && (verificanegatividade(x) == false))
        flag = true;
    else
        flag = false;
    end   
    result = flag;
end

function result = factibilidadeDual(s, Q, x, c, Erro)   
    if((VerificaErroMatriz(abs(s)/(abs(Q*x + c) + 1), Erro) == true) && (verificanegatividade(s) == false))
        flag = true;
    else
        flag = false;
    end
    result = flag;
end

function result = calculoFO(X, C, C_quadratic, answer)
    if(answer == 1)
        for i = 1 : length(X)
            fo(i) = X(i)*C(i) + (X(i)^2)*C_quadratic(i);
        end
        resposta = sum(fo);
    else if (answer == 2)
            resposta = 100*(X(1)- X(2)^2)^2 + (1- X(1))^2;
        end
    end

    result = resposta;
end

function result = Hessian(x)
    result = [202, -400*x(2), 0, 0, 0; -400*x(2), 1200*x(2)^2-400*x(1), 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0];
end

function result = Grad(x)
    result =[202*x(1)-200*(x(2)^2)-2; 400*(x(2)^3)-400*x(1)*x(2);0; 0; 0];
end
