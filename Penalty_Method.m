%_________________________________________________________________________%
%                              PENALTY METHOD                            %
%__*Developed by Joao Augusto Silva Ledo*_________________________________%


function result = penalidade()
    clear all;
    clc; 
    k = 1;
    x{k} = [2;1];
    mi(k) = 0.1;
    beta(k) = 10;
    epslon = 1*10^-3;

  grad{k} = gradiente(x{k},mi(k));
  compara = max(abs(grad{k}));
  while(compara >= epslon)
      compara = max(abs(grad{k}));
      % d{k}=-grad{k};
      d{k}=-inv(hessiana(x{k},mi(k)))*grad{k};
      alpha(k) = fibonacci(d{k}, x{k}, mi(k));
      x{k + 1} = x{k} + alpha(k)*d{k}; 
      mi(k+1) = atualizami(mi(k),beta);
      grad{k+1} = gradiente(x{k+1},mi(k+1));
      alpha(k+1) = fibonacci(-grad{k+1}, x{k+1}, mi(k+1));
      
      if (mi(k)*alpha(k+1) < epslon)%Passo essencial que determina quando parar pois encontrou o ponto que satisfa?a a restricao.
          compara = 1*10^-4;
      end
      k = k + 1;
      
  end 
  valor = [x{k}(1),x{k}(2)];
  resultado.nome = 'Metodo da Penalidade';
  resultado.valor = valor;
  resultado.alpha = alpha(k);
  resultado.mi = mi(k);
  resultado.iteracoes = k;
    
    
    result = resultado;
end

function result = penalty(x)
    p = max(0,x(1)^2-x(2))^2;
    result = p;
end

function result = funcao_auxiliar_penalisada(x,mi)
    FAP = (x(1)-2)^4 + (x(1)-2*x(2))^2 + mi*penalty(x);
    result = FAP;
end

function result = hessiana(x,mi)
    x1 = x(1);
    x2 = x(2);
    result = [ 12*(x1 - 2)^2 + 8*mi*x1^2 - 4*mi*(- x1^2 + x2) + 2, - 4*mi*x1 - 4 ;    - 4*mi*x1 - 4,      2*mi + 8];
end

function result = fibonacci(d,x,mi)

% This portion aims to find alpha | alpha is in between a and b points
% Also the goal in this unidimensional search is finding function f(xk + alpha*dk) minimum
% by feeding dk as the main function gradient
  a=-1.0;
  b=1.0;
  l=0.2;  
  k = 1;
  metodo.nome = 'Fibonacci';
  metodo.a(k) = a;
  metodo.b(k) = b;
  metodo.fibonacci = determinar_N(metodo.a(k),metodo.b(k),l); % the tolerance gap is in calculating n
  metodo.n = length(metodo.fibonacci);
  metodo.lambda(k) = lambd(metodo.a(k),metodo.b(k),metodo.fibonacci,metodo.n,k);
  metodo.mi(k) = calculo_mi(metodo.a(k),metodo.b(k),metodo.fibonacci,metodo.n,k);
  
  metodo.f_lambda(k) = objetivo(metodo.lambda(k),d,x,mi);
  
  metodo.f_mi(k) = objetivo(metodo.mi(k),d,x,mi);
  
  while (k<metodo.n-1)      
      if (objetivo(metodo.lambda(k),d,x,mi) > objetivo(metodo.mi(k),d,x,mi))
        metodo.a(k+1) = metodo.lambda(k);
        metodo.b(k+1) = metodo.b(k);
        metodo.lambda(k+1) = metodo.mi(k);
        metodo.mi(k+1) = calculo_mi(metodo.a(k+1),metodo.b(k+1),metodo.fibonacci,metodo.n,k);
        
        metodo.f_mi(k+1) = objetivo(metodo.mi(k+1),d,x,mi);
        
        metodo.f_lambda(k+1) =  metodo.f_mi(k+1);
      else
          metodo.a(k+1) = metodo.a(k);
          metodo.b(k+1) = metodo.mi(k);
          metodo.mi(k+1) = metodo.lambda(k);
          metodo.lambda(k+1) = lambd(metodo.a(k+1),metodo.b(k+1),metodo.fibonacci,metodo.n,k);
          
          metodo.f_lambda(k+1) = objetivo(metodo.lambda(k+1),d,x,mi);
          
          metodo.f_mi(k+1) =  metodo.f_lambda(k+1);
      end
    k = k+1;
  end 
  metodo.lambda(k) = metodo.lambda(k-1);
  metodo.mi(k) = metodo.mi(k-1)+l;
  if (objetivo(metodo.lambda(k),d,x,mi) < objetivo(metodo.mi(k),d,x,mi))
      metodo.a(k+1) = metodo.a(k);
      metodo.b(k+1) = metodo.mi(k);
  else
      metodo.a(k+1) = metodo.lambda(k);
      metodo.b(k+1) = metodo.b(k);
  end
  metodo.numero_iteracoes = k;
  if objetivo(metodo.a(k),d,x,mi) < objetivo(metodo.b(k),d,x,mi)
      saida = metodo.a(k);
  else
      saida = metodo.b(k);
  end
  result = saida;
end

function result = objetivo(alpha, variavel, z, mi)
    x1 = z(1) + alpha*variavel(1);
    x2 = z(2) + alpha*variavel(2);
    x = [x1,x2];
    result = funcao_auxiliar_penalisada(x,mi);
end

function result = determinar_N(a,b,l)
    fn = (b-a)/l;
    n=1;
    fibonacci = sequencia_fibonacci(10);
    while (fibonacci(n)<fn)
        n = n + 1;
    end
    result = sequencia_fibonacci(n);
end

function result = sequencia_fibonacci(n)
    f(1)=1;
    f(2)=1;
    for i=3:n
        f(i)=f(i-2)+f(i-1);
    end
    result = f;
end

function result = lambd(a,b,fibonacci,n,k)
        lambda = a+(fibonacci(n-1-k)/fibonacci(n-k+1))*(b-a);
    result = lambda;
end

function result = gradiente(variavel,mi)
    x1=variavel(1);
    x2=variavel(2);
    result = [ 2*x1 - 4*x2 + 4*(x1 - 2)^3+0.1*(-4*x1*(- x1^2 + x2));8*x2 - 4*x1+0.1*(- 2*x1^2 + 2*x2)] ;
end

function result = calculo_mi(a,b,fibonacci,n,k)
        mi = a+(fibonacci(n-k)/fibonacci(n-k+1))*(b-a);
    result = mi;
end

function result = atualizami(mi,beta)
    result = mi*beta;
end
