%_________________________________________________________________________%
%                                                                         %
%                                                                         %
%        PRIMAL DUAL AFFINE PREDICTOR CORRECTOR INTERIOR-POINTS METHOD    %
%                         (polynomial linear order)                       % 
%                                                                         %
%                                                                         %
%                                                                         %
%                                          Developed by:                  % 
%                                                 Joao Augusto Silva Ledo %
%_________________________________________________________________________%

function result = PrimalDualPPL()
    clear all;
    clc;
    format long;      
    data = loadInputData(); % Carregar os dados de entrada.
    result = solveDualAfim(data.A, data.C, data.pontoInicial, data.Erro, data.e, data.Sinicial, data.b, data.Mi, data.w); % Cria a Instancia que resolve juntamente com os dados
end

%_________________________________________________________________________%
 function result = loadInputData() % Local onde se carregam as informacoes
 
        resultado.Erro{1} = 10^-2;
        resultado.Erro{2} = 10^-2;
        resultado.Erro{3} = 10^-1;   
        
%         resultado.C = [-2; 1; 0; 0];      %<---------- Problem 1
%         resultado.A = [1, -1, 1, 0; 0, 1, 0, 1];
%         resultado.b = [15; 15];
%         resultado.pontoInicial = [0.5; 1; 1; 0.5];
%         resultado.Sinicial = [1; 1; 1; 1];
%         resultado.w = [0;0];

        resultado.C = [-1; -2; 0; 0];      %<---------- Problem 2 (Mehrotra)
        resultado.A = [-1, 1, 1, 0; 1, 1, 0, 1];
        resultado.b = [1; 5];
        resultado.pontoInicial = [0.5; 1; 1; 0.5];
        resultado.Sinicial = [1; 1; 1; 1];
        resultado.w = [0;0];

        resultado.Mi = 1;
        resultado.e = constroiE(resultado.Sinicial);
     result = resultado;
 end
%_________________________________________________________________________% 

 function result = solveDualAfim(A, C, pontoInicial, Erro, e, SInicial, b, mi_ini, w_ini)
   resposta = ChoseSolutionMethod();
   k = 1;
   flag = false;
   x{k} = pontoInicial;
   s{k} = SInicial;
   w{k} = w_ini;    
   BetaP = {};
   BetaD = {};
   Dw = {};
   Dx = {};
   Ds = {};
   mi{k} = mi_ini;
   while (flag == false)           
       t{k} = constroiT(b, A, x{k});
       u{k} = constroiU(C, A, w{k}, s{k}); 
       X{k} = diag(x{k});
       S{k} = diag(s{k});
       v{k} = constroiV(mi{k}, e, X{k}, S{k});
       p{k} = constroiP(X{k}, v{k});
       teta{k} = constroiTETA(S{k}, X{k});
       FO{k} = calculoFO(x{k}, C);
       erro2{k} = Factbilidade(t{k}, b, Erro{2});
       erro3{k} = Factbilidade(u{k}, C, Erro{3});
       if((mi{k} < Erro{1})  && (erro2{k}.flag == true) && (erro3{k}.flag == true) || (isfinite(sum(x{k})) == false) || (k==130))
               resultado.Nome = VerificaNome(resposta);
               resultado.FO = FO;
               resultado.x = x;
               resultado.s = s;
               resultado.w = w;
               resultado.Dx = Dx;
               resultado.Ds = Ds;
               resultado.Dw = Dw; 
               resultado.mi = mi;
               resultado.BetaP = BetaP;
               resultado.BetaD = BetaD;
               resultado.Diagonal_S = S;
               resultado.Diagonal_X = X;
               resultado.t = t;
               resultado.u = u; 
               resultado.v = v;
               resultado.p = p;
               resultado.teta = teta;
               resultado.erro2 = erro2;
               resultado.erro3 = erro3;
               if(resposta.answer1 == 2)
                resultado.PrevisorCorretor = PrevisorCorretor;
               end
               resultado.Iteracoe = k-1;
               flag = true;
       else
            Escolha{k} = EscolheMetodoResolucao(A, e, x{k}, w{k}, s{k}, t{k}, u{k}, v{k}, p{k}, teta{k}, X{k}, S{k}, resposta, k, mi{k});
            Dw{k} = Escolha{k}.Dw;
            Dx{k} = Escolha{k}.Dx;
            Ds{k} = Escolha{k}.Ds;
            if(resposta.answer1 == 2)
                PrevisorCorretor{k} = Escolha{k}.PrevCorre;
            end
            if((FactbilidadePrimal(t{k}, Dx{k}, C, 10^-12) == true) || (FactbilidadeDual(u{k}, Ds{k}, b, Dw{k}, 10^-12) == true))
                resultado.Nome = VerificaNome(resposta);
                resultado.x = 'Problema Ilimitado';
                resultado.Iteracoe = k-1;
                flag = true;
            else              
                BetaP{k} = passo(x{k}, Dx{k});
                BetaD{k} = passo(s{k}, Ds{k});
                x{k+1} = NextValue(x{k}, BetaP{k}, Dx{k});
                w{k+1} = NextValue(w{k}, BetaD{k}, Dw{k});
                s{k+1} = NextValue(s{k}, BetaD{k}, Ds{k});
                mi{k+1} = ProximaBarreira(mi{k},x{k+1}, s{k+1});
                k = k + 1;
            end    
       end      
   end
   result = resultado;
 end
 
function result = VerificaNome(resposta)
    if(resposta.answer1 == 1)
        resultado = 'Primal Dual Affine Logarithmic Barrier Method for Polynomial Linear Order Problems';
    else
        if (resposta.answer1 == 2)
            if(resposta.answer2 == 1)
                resultado = 'Primal Dual Affine Mehrotra Predictor and Corrector Method for Polynomial Linear Order Problems';
            else
                if(resposta.answer2 == 2)
                     resultado = 'Primal Dual Affine Tanabe-Todd-Ye Predictor and Corrector Method for Polynomial Linear Order Problems';
                else
                    if(resposta.answer2 == 3)
                        resultado = 'Primal Dual Affine Gondzio Predictor and Corrector Method for Polynomial Linear Order Problems';
                    else
                        if(resposta.answer2 == 4)
                            resultado = 'Primal Dual Affine Mixed Predictor and Corrector Method for Polynomial Linear Order Problems';
                        end
                    end
                end
            end         
        end
    end
    result = resultado;
end

 function result = FactbilidadePrimal(t, Dx, c, epslon)
    flag = false;
    if((max(abs(t)) <= epslon) && (verificanegatividade(Dx) == false) && (c'*Dx < 0))
        flag = true;
    end
    result = flag;
 end

  function result = FactbilidadeDual(u, Ds, b, Dw, epslon)
    flag = false;
    if((max(abs(u)) <= epslon) && (verificanegatividade(Ds) == false) && (b'*Dw > 0))
        flag = true;
    end
    result = flag;
 end
 
 function result = Factbilidade(parametro1, parametro2, Erro)
    resultado.flag = false;
    resultado.valor = (max(abs(parametro1))/(max(abs(parametro2))+1));
    if((max(abs(parametro1))/(max(abs(parametro2))+1)) < Erro)
        resultado.flag = true;       
    end
    result = resultado;
 end
 
 function result = constroiT(b, A, x)
    result = b - A*x;
 end
 
 function result = constroiU(c, A, w, s)
    result = c - A'*w - s;
 end
 
 function result = constroiV(mi, e, X, S)
    result = mi*e - X*S*e;
 end
 
 function result = constroiP(X, v)
    result = inv(X)*v;
 end
 
 function result = constroiTETA(S, X)
    result = inv(S)*X;
 end
 
 function result = calculoFO(X, C)
     for i = 1 : length(X)
        fo(i) = X(i)*C(i);
     end
  result = sum(fo);
 end
 
function result = constroiE(x)
    for i = 1 : length(x)
       aux(i) = 1; 
    end
    result = aux';
end

 function result = passo(e, Dy)   
     for i = 1 : length(e)
        if Dy(i) < 0 
            alpha(i) = -0.9995*e(i)/Dy(i);
        else
            alpha(i) = 0;
        end
     end
     if(sum(abs(alpha)) == 0)
         alpha = 1;
     else
        alpha = nonzeros(alpha);
     end
    result = min(alpha);
 end
 
 function result = Direcaow(A, teta, t, u, p)
    result = inv(A*teta*A')*(t + A*teta*(u - p));
 end
 
  function result = DirecaoX(teta, A, Dw, u, p) 
    result = teta*(A'*Dw - u + p);
 end
 
 function result = DirecaoS(X,p, S, Dx) 
    result = p - inv(X)*S*Dx;
 end
 
 function resultado = NextValue(valoratual, passo, direcao)
   resultado = valoratual + passo * direcao;
 end
 
 function result = verificanegatividade(S)
    flag = false;
    for i = 1 : length(S)
       if (S(i) <= 0)
           flag = true;
       end
    end
    result = flag;
 end
 
function result = ProximaBarreira(mi, x, s)
    result = mi/2;
end

function result = ProximaBarreiraPrevCorr(mi, x, s)
    result = ((x'*s)/length(x));
end
 
function result = Previsor(A, e, x, w, s, t, u, v, teta, X, S, resposta, mi)
    if(resposta == 1) % Mehrotra
        resultado.Dw = inv(A*teta*A')*(t + A*teta*u);
        resultado.Dx = teta*(A'*resultado.Dw - u);
        resultado.Ds = -inv(X)*S*resultado.Dx;
        resultado.BetaP = passo(x, resultado.Dx);
        resultado.BetaD = passo(s, resultado.Ds);
        resultado.x = NextValue(x, resultado.BetaP, resultado.Dx);
        resultado.w = NextValue(w, resultado.BetaD, resultado.Dw);
        resultado.s = NextValue(s, resultado.BetaD, resultado.Ds);
        resultado.mi = ProximaBarreiraPrevCorr(mi, resultado.x, resultado.s);
        resultado.v = constroiV(resultado.mi, e, X, S) - diag(resultado.Dx) * diag(resultado.Ds) * e;
    end
    if(resposta == 2) % Tanabe-Todd-Ye
        resultado.Dw = inv(A*teta*A')*(t + A*teta*u);
        resultado.Dx = teta*(A'*resultado.Dw - u);
        resultado.Ds = -inv(X)*S*resultado.Dx;
        resultado.BetaP = passo(x, resultado.Dx);
        resultado.BetaD = passo(s, resultado.Ds);
        resultado.x = NextValue(x, resultado.BetaP, resultado.Dx);
        resultado.w = NextValue(w, resultado.BetaD, resultado.Dw);
        resultado.s = NextValue(s, resultado.BetaD, resultado.Ds);
        resultado.mi = ProximaBarreiraPrevCorr(mi, resultado.x, resultado.s);
        resultado.v = constroiV(resultado.mi, e, X, S);
    end
    if(resposta == 3) % Gondzio
        resultado.Dw = inv(A*teta*A')*(t + A*teta*u);
        resultado.Dx = teta*(A'*resultado.Dw - u);
        resultado.Ds = -inv(X)*S*resultado.Dx;
        resultado.BetaP = passo(x, resultado.Dx);
        resultado.BetaD = passo(s, resultado.Ds);
        resultado.x = NextValue(x, resultado.BetaP, resultado.Dx);
        resultado.w = NextValue(w, resultado.BetaD, resultado.Dw);
        resultado.s = NextValue(s, resultado.BetaD, resultado.Ds);
        resultado.mi = resultado.x'*resultado.s;
        resultado.v = constroiV(resultado.mi, e, X, S) - diag(resultado.Dx) * diag(resultado.Ds) * e;
    end
    if(resposta == 4) %Misto
        resultado.Dw = inv(A*teta*A')*(t + A*teta*(u - inv(X)*v));
        resultado.Dx = teta*(A'*resultado.Dw - u + inv(X)*v);
        resultado.Ds = -inv(X)*(v - S*resultado.Dx);
        resultado.BetaP = passo(x, resultado.Dx);
        resultado.BetaD = passo(s, resultado.Ds);
        resultado.x = NextValue(x, resultado.BetaP, resultado.Dx);
        resultado.w = NextValue(w, resultado.BetaD, resultado.Dw);
        resultado.s = NextValue(s, resultado.BetaD, resultado.Ds);
        resultado.mi = ProximaBarreiraPrevCorr(mi, resultado.x, resultado.s);
        resultado.v = constroiV(resultado.mi, e, X, S) - diag(resultado.Dx) * diag(resultado.Ds) * e;
    end
    result = resultado;
end

function result = Corretor(A, e, x, w, s, t, u, v, teta, X, S, resposta, k, p, mi)
    if(resposta == 1) %Mehotra
        previsor = Previsor(A, e, x, w, s, t, u, v, teta, X, S, resposta, mi);
        resultado.Dw = inv(A*teta*A')*(t + A*teta*(u - inv(X)*previsor.v));
        resultado.Dx = teta*(A'*resultado.Dw-u+inv(X)*previsor.v);
        resultado.Ds = inv(X)*(previsor.v-S*resultado.Dx);
        resultado.Prev = previsor;
    end
    if(resposta == 2)
        previsor = Previsor(A, e, x, w, s, t, u, v, teta, X, S, resposta, mi);
        if(rem(k,2) ~= 0) %Tanabi       
            resultado.Dw = inv(A*teta*A')*(t + A*teta*(u - inv(X)*previsor.v));
            resultado.Dx = teta*(A'*resultado.Dw-u+inv(X)*previsor.v);
            resultado.Ds = inv(X)*(previsor.v-S*resultado.Dx);
            resultado.Prev = previsor;
        else
            resultado.Dw = previsor.Dw;
            resultado.Dx = previsor.Dx;
            resultado.Ds = previsor.Ds;
            resultado.Prev = previsor;
        end
        
    end
    if(resposta == 3) %Gondzio
        previsor = Previsor(A, e, x, w, s, t, u, v, teta, X, S, resposta, mi);
        correction.Dw = inv(A*teta*A')*(t + A*teta*(u - inv(X)*previsor.v));
        correction.Dx = teta*(A'*correction.Dw-u+inv(X)*previsor.v);
        correction.Ds = inv(X)*(previsor.v-S*correction.Dx);      
        correction.BetaP = passo(x, correction.Dx);
        correction.BetaD = passo(s, correction.Ds);
        correction.x = NextValue(x, correction.BetaP, correction.Dx);
        correction.w = NextValue(w, correction.BetaD, correction.Dw);
        correction.s = NextValue(s, correction.BetaD, correction.Ds);
        correction.mi = correction.x'*correction.s;
        if(previsor.mi < correction.mi)
            resultado.Dw = previsor.Dw;
            resultado.Dx = previsor.Dx;
            resultado.Ds = previsor.Ds;
            resultado.Prev = previsor;
        else
                resultado.Dw = correction.Dw;
                resultado.Dx = correction.Dx;
                resultado.Ds = correction.Ds;
                resultado.Prev = previsor;
        end
    end
    if(resposta == 4) %Mixed
        previsor = Previsor(A, e, x, w, s, t, u, v, teta, X, S, resposta, mi);
        resultado.Dw = inv(A*teta*A')*(t + A*teta*(u - inv(X)*previsor.v));
        resultado.Dx = teta*(A'*resultado.Dw-u+inv(X)*previsor.v);
        resultado.Ds = inv(X)*(previsor.v-S*resultado.Dx);
        resultado.Prev = previsor;
    end
    result = resultado;
end
 
function result = ChoseSolutionMethod()
    resultado.answer1 = input('Select the Solution Method: \n 1-Primal Dual Method \n 2-Primal Dual with Predictor and Corrector Method \n');
    resultado.answer1 = VerificaOpcao(resultado.answer1);
    if(resultado.answer1 == 2)
        resultado.answer2 = input('Select the Predictor and Corrector Method: \n 1-Mehrotra \n 2-Tanabe-Todd-Ye \n 3-Gondzio \n 4-Mixed \n');
        resultado.answer2 = VerificaOpcaoMetodo(resultado.answer2); 
    end
    result = resultado;
end

function result = VerificaOpcao(answer)
    while((answer ~= 1) && (answer ~= 2))
       answer = input('Favor Escolher corretamente o metodo pelo qual deseja Executar the Solution Method: \n 1-Primal-Dual \n 2-Primal Dual com Previsor-Corretor \n');
    end
    result = answer;
end

function result = VerificaOpcaoMetodo(answer)
    while((answer ~= 1) && (answer ~= 2) && (answer ~= 3) && (answer ~= 4))
       answer = input('Favor Escolher corretamente o metodo pelo qual deseja Executar the Predictor and Corrector Method: \n 1-Mehrotra \n 2-Tanabe-Todd-Ye \n 3-Gondzio \n 4-Mixed \n');
    end
    result = answer;
end 
 
function result = EscolheMetodoResolucao(A, e, x, w, s, t, u, v, p, teta, X, S, resposta, k, mi)
    if(resposta.answer1 == 1)
        resultado.Dw = Direcaow(A, teta, t, u, p);
        resultado.Dx = DirecaoX(teta, A, resultado.Dw, u, p);
        resultado.Ds = DirecaoS(X, p, S, resultado.Dx);       
    else
        if(resposta.answer1 == 2)
           correct = Corretor(A, e, x, w, s, t, u, v, teta, X, S, resposta.answer2, k, p, mi);
           resultado.Dw = correct.Dw;
           resultado.Dx = correct.Dx;
           resultado.Ds = correct.Ds;
           resultado.PrevCorre = correct;
        end
    end
    result = resultado;
end
