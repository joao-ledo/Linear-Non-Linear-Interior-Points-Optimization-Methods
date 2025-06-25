% ======================================================================= %
%                                                                         %
%                           Busca Dicotomica                              %
%                                                                         %
%                                   Developed by Jo?o Augusto Silva L?do  %
% ======================================================================= %
function result = Busca_Dicotomica()
  clear all;
  clc;
  % =======================================================================
  %                  Inicializacao dos dados de entrada
  % =======================================================================
  E = 0.01;
  
  a=1;
  b=3;
%   
%    a=0;
%    b=1;

%  a=-3.0;
%  b=6.0;
  
  l=0.2;
  k = 1;
  
  % =======================================================================
  %                         Inicio da Busca Dicotomica
  % =======================================================================
  dicotomica.nome = 'Busca Dicotomica';
  dicotomica.a(k) = a;
  dicotomica.b(k) = b;
  while abs(dicotomica.b(k)-dicotomica.a(k))>=l 
      dicotomica.numero_de_iteracoes = k+1;
      dicotomica.lambda(k) = lambda(dicotomica.a(k),dicotomica.b(k),E);
      dicotomica.mi(k) = mi(dicotomica.a(k),dicotomica.b(k),E); 
      dicotomica.funcao_em_lambda(k) = objetivo(dicotomica.lambda(k));
      dicotomica.funcao_em_mi(k) = objetivo(dicotomica.mi(k));
      if dicotomica.funcao_em_lambda(k) <= dicotomica.funcao_em_mi(k)
          dicotomica.a(k+1) = dicotomica.a(k);
          dicotomica.b(k+1) = dicotomica.mi(k);
      else
          dicotomica.a(k+1) = dicotomica.lambda(k);
          dicotomica.b(k+1) = dicotomica.b(k);
      end     
          dicotomica.objetivo = min(dicotomica.funcao_em_mi,dicotomica.funcao_em_lambda);
      k = k + 1;  
  end
  result = dicotomica;
end

  % =======================================================================
  %                            Calculo do Lambda
  % =======================================================================
    function result = lambda(a,b,E)
    result = ((a+b)/2)-E;
    end

  % =======================================================================
  %                          Calculo do Mi
  % =======================================================================
    function result = mi(a,b,E)
    result = ((a+b)/2)+E;
    end

  % =======================================================================
  %                        Calculo da Funcao Objetivo
  % =======================================================================
    function result = objetivo(x)
     result =  x^2-5*x+6;
    
    % result =  4*x^2+3*x;
    
     % result = x^2+2*x;
    end 