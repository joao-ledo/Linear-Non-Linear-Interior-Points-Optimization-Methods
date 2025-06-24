%_________________________________________________________________________%
%                       AUGMENTED LAGRANGIAN METHOD                       %
%__*Developed by Joao Augusto Silva Ledo*_________________________________%

function result = LagrangianaAumentada()
    clear all;
    clc; 
    k = 1;
    i = 1;
    j = 1;
    lambda(i) = 0;
    c(j) = 1;
    x{k} = [0,0];
    beta = 10;
    viola(k) = 100000;
    epslon = 0.0001;
    while viola(k) > epslon
        x{k+1} = Newton(x{k},c(j),lambda(i));
        viola(k+1) = abs(residuo(x{k+1}));
        if viola(k+1) <= (1/4)*viola(k)        
           lambda(i+1) = lambda(i) + c(j)*restrictions(x{k+1});
           i = i + 1;
        else
           c(j+1) = beta * c(j);
           j = j+1;
        end
        k = k + 1;
    end
    metodo.nome = 'Augmented Lagrangiana Method';
    metodo.valor = x{k};
    metodo.iteracoes = k;
    %viola(k) = residuo(x(k));
    result = metodo;
end

function result = residuo(x)
    x1 = x(1);
    x2 = x(2);
    g(1) = x1 + x2 - 1; % This array keeps the found value of each constraint to the X of such iteration
    result = max(g);
end

function result = restrictions(x)
   result = x(1) + x(2) - 1; 
end

function result = hessiana(x,lambda,c)
    result = [ 2*c + 2,     2*c ; 2*c, 2*c + 2];
end

function result = Newton(x,c,lambda)
  grad = gradiente(x,lambda,c);
 % while(max(abs(grad{contador})) >= epslon)
      d=-inv(hessiana(x,lambda,c))*grad;
      alpha = fibonacci(d, x, lambda,c);
      x = x + alpha*d';    
      %grad{contador+1} = gradiente(x{contador+1},lambda(contador+1));    
 % end 

result = x;
end

function result = gradiente(variavel,lambda,c)
    x1=variavel(1);
    x2=variavel(2);
    result = [lambda + 2*x1 + c*(2*x1 + 2*x2 - 2);lambda + 2*x2 + c*(2*x1 + 2*x2 - 2)] ;
end

function result = fibonacci(d,x,lambda,c)

% o intuito dessa parte ? descobrir o alpha, a e b s?o pontos que entre
% eles exista um alpha, o objetivo dessa parte da busca unidimensional ?
% encontrar o minimo da fun??o f(xk + alpha*dk), alimentando o dk como
% ogradiente da fun??o principal.
  a=-1.0;
  b=1.0;
  l=0.2;  
  k = 1;
  metodo.nome = 'Fibonacci';
  metodo.a(k) = a;
  metodo.b(k) = b;
  metodo.fibonacci = determinar_N(metodo.a(k),metodo.b(k),l); % o erro esta no calculo de n
  metodo.n = length(metodo.fibonacci);
  metodo.lambda(k) = lambd(metodo.a(k),metodo.b(k),metodo.fibonacci,metodo.n,k);
  metodo.mi(k) = mi(metodo.a(k),metodo.b(k),metodo.fibonacci,metodo.n,k);
  
  metodo.f_lambda(k) = objetivo(metodo.lambda(k),d,x,lambda,c);
  
  metodo.f_mi(k) = objetivo(metodo.mi(k),d,x,lambda,c);
  
  while (k<metodo.n-1)      
      if (objetivo(metodo.lambda(k),d,x,lambda,c) > objetivo(metodo.mi(k),d,x,lambda,c))
        metodo.a(k+1) = metodo.lambda(k);
        metodo.b(k+1) = metodo.b(k);
        metodo.lambda(k+1) = metodo.mi(k);
        metodo.mi(k+1) = mi(metodo.a(k+1),metodo.b(k+1),metodo.fibonacci,metodo.n,k);
        
        metodo.f_mi(k+1) = objetivo(metodo.mi(k+1),d,x,lambda,c);
        
        metodo.f_lambda(k+1) =  metodo.f_mi(k+1);
      else
          metodo.a(k+1) = metodo.a(k);
          metodo.b(k+1) = metodo.mi(k);
          metodo.mi(k+1) = metodo.lambda(k);
          metodo.lambda(k+1) = lambd(metodo.a(k+1),metodo.b(k+1),metodo.fibonacci,metodo.n,k);
          
          metodo.f_lambda(k+1) = objetivo(metodo.lambda(k+1),d,x,lambda,c);
          
          metodo.f_mi(k+1) =  metodo.f_lambda(k+1);
      end
    k = k+1;
  end 
  metodo.lambda(k) = metodo.lambda(k-1);
  metodo.mi(k) = metodo.mi(k-1)+l;
  if (objetivo(metodo.lambda(k),d,x,lambda,c) < objetivo(metodo.mi(k),d,x,lambda,c))
      metodo.a(k+1) = metodo.a(k);
      metodo.b(k+1) = metodo.mi(k);
  else
      metodo.a(k+1) = metodo.lambda(k);
      metodo.b(k+1) = metodo.b(k);
  end
  metodo.numero_iteracoes = k;
  if objetivo(metodo.a(k),d,x,lambda,c) < objetivo(metodo.b(k),d,x,lambda,c)
      saida = metodo.a(k);
  else
      saida = metodo.b(k);
  end
  result = saida;
end

function result = objetivo(alpha, variavel, z, lambda, c) % recebe alpha, recebe direcao, recebe z = [x,y]
    x1 = z(1) + alpha*variavel(1);
    x2 = z(2) + alpha*variavel(2);
    result = x1^2 + x2^2 + lambda*(x1 + x2 - 1)+c*(x1 + x2 - 1)^2;
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

function result = mi(a,b,fibonacci,n,k)
        mi = a+(fibonacci(n-k)/fibonacci(n-k+1))*(b-a);
    result = mi;
end
