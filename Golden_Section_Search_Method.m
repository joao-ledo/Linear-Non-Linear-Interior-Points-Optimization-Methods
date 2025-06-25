%_________________________________________________________________________%
%                       GOLDEN SECTION SEARCH METHOD                      %
%__*Developed by Joao Augusto Silva Ledo*_________________________________%

function result = Busca_Seguimento_Aureo()
  clear all;
  clc;
  % =======================================================================
  %                         LOADING INPUT DATA
  % =======================================================================
  alpha = (sqrt(5.0)-1)/2;
  
%   a=1;
%   b=3;
%   
  a=0;
  b=1;

%   a=-3.0;
%   b=6.0;
  
  l=0.2;
  k = 1;
  
  % =======================================================================
  %               STARTING THE GOLDEN SECTION SEARCH METHOD
  % =======================================================================
  segmento_aureo.nome = 'Golden Section Search Method';
  segmento_aureo.a(k) = a;
  segmento_aureo.b(k) = b;  
     segmento_aureo.lambda(k) = lambda(segmento_aureo.a(k),segmento_aureo.b(k),alpha);
      segmento_aureo.mi(k) = mi(segmento_aureo.a(k),segmento_aureo.b(k),alpha);
  while abs(segmento_aureo.b(k)-segmento_aureo.a(k))>=l 
      segmento_aureo.numero_de_iteracoes = k+1;
      segmento_aureo.lambda(k+1) = lambda(segmento_aureo.a(k),segmento_aureo.b(k),alpha);
      segmento_aureo.mi(k+1) = mi(segmento_aureo.a(k),segmento_aureo.b(k),alpha);
      segmento_aureo.funcao_em_lambda(k) = objetivo(segmento_aureo.lambda(k));
      segmento_aureo.funcao_em_mi(k) = objetivo(segmento_aureo.mi(k));
      if segmento_aureo.funcao_em_lambda(k) <= segmento_aureo.funcao_em_mi(k)
          segmento_aureo.a(k+1) = segmento_aureo.a(k);
          segmento_aureo.b(k+1) = segmento_aureo.mi(k);
          segmento_aureo.mi(k+1) = segmento_aureo.lambda(k);
          segmento_aureo.funcao_em_mi(k+1) = segmento_aureo.funcao_em_lambda(k);
          segmento_aureo.funcao_em_lambda(k+1) = lambda(segmento_aureo.a(k+1),segmento_aureo.b(k+1),alpha);
          segmento_aureo.funcao_em_lambda(k+1) = objetivo(segmento_aureo.funcao_em_lambda(k+1));
      else
          segmento_aureo.a(k+1) = segmento_aureo.lambda(k);
          segmento_aureo.b(k+1) = segmento_aureo.b(k);
          segmento_aureo.lambda(k+1) = segmento_aureo.mi(k);
          segmento_aureo.funcao_em_lambda(k+1) = segmento_aureo.funcao_em_mi(k);
          segmento_aureo.funcao_em_mi(k+1) = mi(segmento_aureo.a(k+1),segmento_aureo.b(k+1),alpha);
          segmento_aureo.funcao_em_mi(k+1) = objetivo(segmento_aureo.funcao_em_mi(k+1));
      end     

      segmento_aureo.objetivo = min(segmento_aureo.funcao_em_mi,segmento_aureo.funcao_em_lambda);
      k = k + 1;  
  end
  result = segmento_aureo;
end

  % =======================================================================
  %                           CALCULATING LAMBDA
  % =======================================================================
    function result = lambda(a,b,alpha)
    result = a+(1-alpha)*(b-a);
    end

  % =======================================================================
  %                           CALCULATING MU
  % =======================================================================
    function result = mi(a,b,alpha)
    result = a+alpha*(b-a);
    end

  % =======================================================================
  %                  FINDING THE OBJECTIVE FUNCTION VALUES
  % =======================================================================
    function result = objetivo(x)
     % result =  x^2-5*x+6;
    
     result =  4*x^2+3*x;
    
     % result = x^2+2*x;
    end 
