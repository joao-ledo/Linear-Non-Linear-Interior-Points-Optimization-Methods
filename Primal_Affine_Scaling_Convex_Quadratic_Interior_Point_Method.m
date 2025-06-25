%_________________________________________________________________________%
%                                                                         %
%                                                                         %
%        PRIMAL AFFINE SCALING CONVEX QUADRATIC INTERIOR-POINT METHOD     %
%                                                                         %
%                                                                         %
%                                          Developed by:                  % 
%                                                 Joao Augusto Silva Ledo %
%_________________________________________________________________________%

function result = PrimalAfimPPQ()
    clear all;
    clc;
    format long;      
    data = loadInputData(); % Carregar os dados de entrada.
    result = solvePrimalAfimPPQ(data.A, data.C, data.C_quadratic, data.pontoInicial, data.Erro, data.b, data.Q); % Cria a Instancia que resolve juntamento com os dados
end

%_________________________________________________________________________%
 function result = loadInputData() % Local onde se carregam as informacoes
    resultado.Erro{1} = 10^-4;
    resultado.Erro{2} = 10^-5;
    resultado.Erro{3} = 10^-3;
    
    resultado.C = [-2; 0; 0; 0; 0];%[1; 2; -3];
    resultado.C_quadratic = [2; 3; 5];
    resultado.A = [1, -1, 1, 0, 0; 0, -1, 0, 1, 0; 0, 1, 0, 0, 1];   % <------------- Problema 2
    resultado.b = [15; -2; 15];%[5; 10];

    resultado.pontoInicial{2} = [2; 4; 17; 2; 11];%[10; 4; 9; 2; 11];%[10; 3; 8; 1; 12];%[3; 2; 8];
    resultado.Q = Hessian(resultado.pontoInicial{2});
    resultado.pontoInicial{1} = inv(-resultado.Q) * resultado.C;
        
     result = resultado;
 end
%_________________________________________________________________________% 

 function result = solvePrimalAfimPPQ(A, C, C_quadratic, pontoInicial, Error, b, Q_ini)
   Erro = Error{1};
   ErroPrimal = Error{2};
   ErroDual = Error{3};
   k = 1;
   flag = false;
   alpha{k} = {};
   x{k} = pontoInicial{1};
   Fo{k} = FuncaoObjetivo(x{k});%calculoFO(x{k}, C, C_quadratic);
   if((VerificaEquacao(A, x{k}, b) == true) && (verificanegatividade(x{k}) == false))
       resultado.Nome = 'Algoritmo Primal Afim (PPQ)';
       resultado.Fo = Fo;
       resultado.x = x;
       resultado.Iteracoe = k-1;
       flag = true;
   else
       x{k} = pontoInicial{2};
   end
%    answer = input('Digite 1 se deseja resolver pelo metodo de polinomio ordem quadratica ou 2 por metodo de polinomio de ordem N-dimensional: \n');
%    answer = VerificaOpcao(answer);
%    VerificaOpcao(answer);
   while(flag == false)
        Fo{k} = FuncaoObjetivo(x{k});%calculoFO(x{k}, C, C_quadratic);
        Q{k} = Hessian(x{k});%ChooseMethod(answer, Q_ini, x{k});
        X{k} = diag(x{k});
        H{k} = Calcula_H(Q{k}, X{k});
        w{k} = Calcula_W(A, H{k}, Q{k}, x{k}, C);
        s{k} = Calcula_S(Q{k}, x{k}, C, A, w{k});
        if((VerificaErro(x{k}, s{k}, Erro) == true) && (VerificaPositividade(x{k})== true) && (VerificaPositividade(s{k})== true) && ((factbilidadePrimal(A, x{k}, b, ErroPrimal) == true) && (factibilidadeDual(s{k}, Q{k}, x{k}, C, ErroDual) == true)))
            resultado.Nome = 'Algoritmo Primal Afim (PPQ)';
            resultado.Fo = Fo;
            resultado.Q = Q;
            resultado.x = x;
            resultado.Matriz_X = X;
            resultado.H = H;
            resultado.w = w;
            resultado.s = s;
            resultado.Alpha = alpha;
            resultado.Iteracoe = k-1;
            flag = true;
        else
            Dx{k} = Direcao(H{k}, s{k});
            if((VerificaPositividade(Dx{k}) == true))
                resultado.Nome = 'Algoritmo Primal Afim (PPQ)';
                resultado.x = 'Problema Ilimitado';
                resultado.Iteracoe = k-1;
                flag = true;
            else
                if((VerificaMenorErro(Dx{k}, Erro) == true) && ((VerificaPositividade(x{k}) == true) || (VerificaPositividade(s{k}) == true)))
                    resultado.Nome = 'Algoritmo Primal Afim (PPQ)';
                    resultado.Fo = Fo;
                    resultado.Q = Q;
                    resultado.x = x;
                    resultado.Matriz_X = X;
                    resultado.H = H;
                    resultado.w = w;
                    resultado.s = s;
                    resultado.Dx = Dx;
                    resultado.Alpha = alpha;
                    resultado.Iteracoe = k-1;
                    flag = true;
                else
                    alpha{k} = passo(x{k}, Dx{k}, Q{k}, C);
                    x{k+1} = NextValue(x{k}, alpha{k}, Dx{k});
                    k = k + 1;
                end
            end
        end
   end
   result = resultado;
 end
 
 function result = FuncaoObjetivo(x)   
    result = 100*(x(1)- x(2)^2)^2 + (1- x(1))^2;
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
 
 function result = Calcula_H(Q, X)
    result = inv(Q + inv(X^2));
 end
 
 function result = Hessian(x)
    result = [202, -400*x(2), 0, 0, 0; -400*x(2), 1200*x(2)^2-400*x(1), 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0];
 end
 
 function result = Calcula_W(A, H, Q, x, c) 
    result = inv(A*H*A')*A*H*(Grad(x));%inv(A*H*A')*A*H*(Q*x+c);  
 end
 
 function result = Calcula_S(Q, x, c, A, w)
    result = (Grad(x)) - A'*w; %(Q*x + c) - A'*w;
 end
 
 
 function result = Direcao(H, s)
    result = -H*s;
 end
 
 function result = calculoFO(X, C, C_quadratic)
     for i = 1 : length(X)
        fo(i) = X(i)*C(i) + (X(i)^2)*C_quadratic(i);
     end
    result = sum(fo);
 end
 
 function result = passo(x, Dx, Q, c)
     for i = 1 : length(x)
        if Dx(i) < 0 
            alpha(i) = -0.9995*x(i)/Dx(i);
        else
            alpha(i) = 0;
        end
     end    
     alpha1 = min(nonzeros(alpha));    
     alpha2 = - (Dx' * (Grad(x)))/(Dx' * Q  * Dx);
     %alpha2 = - (Dx' * (Q*x + c))/(Dx' * Q  * Dx);
     if((sum(alpha) == 0))
         resposta = alpha2;
     else
         resposta = 0.9995*min(alpha1, alpha2);
     end
    result = resposta;
 end
 
 function result = Grad(x)
    result =[202*x(1)-200*(x(2)^2)-2; 400*(x(2)^3)-400*x(1)*x(2);0; 0; 0];
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
 
 function result = VerificaErro(x, s, Erro)     
    flag = true;
    if(x'*s > Erro)
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
 
  function result = VerificaErroMatriz(s, Erro)
    flag = true;
    [linha, coluna] = size(s);
    for i = 1 : linha
        for j = 1 : coluna
            if(s(i,j) > Erro)
                flag = false;
            end
        end
    end
    result = flag;
 end
