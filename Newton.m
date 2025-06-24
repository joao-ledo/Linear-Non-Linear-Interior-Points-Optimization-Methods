%_________________________________________________________________________%
%                              NEWTON METHOD                              %
%_________________________________________________________________________%

function result = Newton()
  clear all;
  clc;
     
%   a=1;
%   b=3;
%   
 % a=0;
 % b=1;
  
   a=-3.0;
   b=6.0;
  
  
  E=0.01;
  k = 1;
  metodo.nome = 'Newton';
  metodo.x(k)=a; 
  metodo.derivadas{k} = deriva_objetivo(metodo.x(k));
  metodo.x(k+1) = b;%metodo.x(k)-(metodo.derivadas{k}(1)/metodo.derivadas{k}(2));
 % metodo.derivadas{k}(1) <---- First order derivative
  % metodo.derivadas{k}(2) <---- Second order derivative
  while(abs(metodo.derivadas{k}(1))>E) && (abs(metodo.x(k+1)-metodo.x(k))/abs(metodo.x(k+1))>=E) 
      k = k+1;
      metodo.derivadas{k} = deriva_objetivo(metodo.x(k));
      metodo.x(k+1) = metodo.x(k)-(metodo.derivadas{k}(1)/metodo.derivadas{k}(2));           
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
   
  %result = [2*x-5,2];
   
  %result = [8*x+3,8];
   
   result = [2*x+2, 2];
end
