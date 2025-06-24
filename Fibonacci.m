function result = Fibonacci()
  clear all;
  clc;
  
%   a=1;
%   b=3;
%   
%   a=0;
%   b=1;

  a=-3.0;
  b=6.0;
  
  l=0.2;
  k = 1;
  metodo.nome = 'Fibonacci';
  metodo.a(k) = a;
  metodo.b(k) = b;
  metodo.fibonacci = determinar_N(metodo.a(k),metodo.b(k),l); % o erro esta no calculo de n
  metodo.n = length(metodo.fibonacci);
  metodo.lambda(k) = lambda(metodo.a(k),metodo.b(k),metodo.fibonacci,metodo.n,k);
  metodo.mi(k) = mi(metodo.a(k),metodo.b(k),metodo.fibonacci,metodo.n,k);
  metodo.f_lambda(k) = objetivo(metodo.lambda(k));
  metodo.f_mi(k) = objetivo(metodo.mi(k));
  while (k<metodo.n-1)      
      if (objetivo(metodo.lambda(k)) > objetivo(metodo.mi(k)))
        metodo.a(k+1) = metodo.lambda(k);
        metodo.b(k+1) = metodo.b(k);
        metodo.lambda(k+1) = metodo.mi(k);
        metodo.mi(k+1) = mi(metodo.a(k+1),metodo.b(k+1),metodo.fibonacci,metodo.n,k);
        metodo.f_mi(k+1) = objetivo(metodo.mi(k+1));
        metodo.f_lambda(k+1) =  metodo.f_mi(k+1);
      else
          metodo.a(k+1) = metodo.a(k);
          metodo.b(k+1) = metodo.mi(k);
          metodo.mi(k+1) = metodo.lambda(k);
          metodo.lambda(k+1) = lambda(metodo.a(k+1),metodo.b(k+1),metodo.fibonacci,metodo.n,k);
          metodo.f_lambda(k+1) = objetivo(metodo.lambda(k+1));
          metodo.f_mi(k+1) =  metodo.f_lambda(k+1);
      end
    k = k+1;
  end 
  metodo.lambda(k) = metodo.lambda(k-1);
  metodo.mi(k) = metodo.mi(k-1)+l;
  if (objetivo(metodo.lambda(k)) < objetivo(metodo.mi(k)))
      metodo.a(k+1) = metodo.a(k);
      metodo.b(k+1) = metodo.mi(k);
  else
      metodo.a(k+1) = metodo.lambda(k);
      metodo.b(k+1) = metodo.b(k);
  end
  metodo.numero_iteracoes = k;
  result = metodo;
end

function result = objetivo(x)
    %result =  x^2-5*x+6;
    
    %result =  4*x^2+3*x;
    
    result = x^2+2*x;
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

function result = lambda(a,b,fibonacci,n,k)
    result = a+(fibonacci(n-1-k)/fibonacci(n-k+1))*(b-a);
end

function result = mi(a,b,fibonacci,n,k)
    result = a+(fibonacci(n-k)/fibonacci(n-k+1))*(b-a);
end