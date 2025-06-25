%_________________________________________________________________________%
%                         CUBIC ADJUSTMENT METHOD                         %
%__*Developed by Joao Augusto Silva Ledo*_________________________________%

function result = Ajuste_Cubico()
  clear all;
  clc;
  
%   a=1;
%   b=3;
%   
%   a=0;
%   b=1;
  
  a=-3.0;
  b=6.0;
  
  E=0.2;
  k = 1;
  metodo.nome = 'Cubic Adjustment Method';
  metodo.x(k)=a; 
  metodo.x(k+1)=b;
  metodo.derivadas(k) = deriva_objetivo(metodo.x(k));
   metodo.objetivo(k) = abs(metodo.derivadas(k));
  while((abs(metodo.x(k+1)-metodo.x(k))/abs(metodo.x(k+1))>=E) || (abs(metodo.derivadas(k))<E))
      k = k+1;    
      metodo.derivadas(k+1) = deriva_objetivo(metodo.x(k));
      metodo.objetivo(k+1) = metodo.derivadas(k);
      metodo.x(k+1) = metodo.x(k) - (metodo.x(k)-metodo.x(k-1))*((deriva_objetivo(metodo.x(k)) + u2(metodo.x,k) - u1(metodo.x,k))/(deriva_objetivo(metodo.x(k))-deriva_objetivo(metodo.x(k-1))+2*u2(metodo.x,k)));     
  end
  metodo.numero_iteracoes = k;

  result = metodo;
end

function result = u1(x,k)

    result = deriva_objetivo(x(k-1)) + deriva_objetivo(x(k))-3* ( (objetivo(x(k-1))-objetivo(x(k)))/(x(k-1)-x(k)));
end

function result = u2(x,k)
    result = (u1(x,k)^2 - deriva_objetivo(x(k-1))* deriva_objetivo(x(k)))^(1/2);
end

function result = objetivo(x)
   % result =  x^2-5*x+6;
    
   % result =  4*x^2+3*x;
    
    result = x^2+2*x;
end

function result = deriva_objetivo(x)
   % syms y;
   % f = y^2+2*y;
   % f_linha = diff(f,y);
   % f_2linhas = diff(f_linha,y);
   % y = x;
   % result = [subs(f_linha), subs(f_2linhas)];
  % result = 2*x-5;
   
 %  result = 8*x+3;
   
   result = 2*x+2;
end
