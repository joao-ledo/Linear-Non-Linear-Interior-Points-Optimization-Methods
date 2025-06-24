%_________________________________________________________________________%
%                           FAKE POSITION METHOD                          %
%__*Developed by Joao Augusto Silva Ledo*_________________________________%

function result = Falsa_posicao()
  clear all;
  clc;
  
%   
  a=1;
  b=3;
%   
%   a=0;
%   b=1;
  
%   a=-3.0;
%   b=6.0;
  
  E=0.2;
  k = 1;
  metodo.nome = 'Fake Position Method';
  metodo.x(k)=a; 
  metodo.x(k+1) = b;
  metodo.derivadas(k) = deriva_objetivo(metodo.x(k));
   metodo.objetivo(k) = abs(metodo.derivadas(k));
  while(abs(metodo.derivadas(k))<E) || (abs(metodo.x(k+1)-metodo.x(k))/abs(metodo.x(k+1))>=E)
      k = k+1;
      metodo.derivadas(k+1) = deriva_objetivo(metodo.x(k));
      metodo.objetivo(k+1) = metodo.derivadas(k);
      metodo.x(k+1) = metodo.x(k) - deriva_objetivo(metodo.x(k))*((metodo.x(k-1)-metodo.x(k))/(deriva_objetivo(metodo.x(k-1)-deriva_objetivo(metodo.x(k)))));   
  end
  metodo.numero_iteracoes = k;
  result = metodo;
end

function result = deriva_objetivo(x)
   % syms y;
   % f = y^2+2*y;
   % f_linha = diff(f,y);
   % f_2linhas = diff(f_linha,y);
   % y = x;
   % result = [subs(f_linha), subs(f_2linhas)];
   
   result = 2*x-5;
   
   %result = 8*x+3;
   
   %result = 2*x+2;
end
