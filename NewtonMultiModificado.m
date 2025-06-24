%_________________________________________________________________________%
%                 MODIFIED MULTIDIMENSIONAL NEWTON METHOD                 %
%__*Developed by Joao Augusto Silva Ledo*_________________________________%

function result = NewtonMultiModificado()
  clear all;
  clc; 
  epslon = 0.01;
  contador = 1;
  %x{contador} = [0,0];  
  %x{contador} = [8,9];
  x{contador} = [-1/2,1];
  grad{contador} = gradiente(x{contador});
  hessi{contador} = hessiana(x{contador});
  while(max(abs(grad{contador})) >= epslon)
      d{contador}= - inv(hessi{contador}) * grad{contador};
      alpha = fibonacci(d{contador}, x{contador}); % Busca Unidimensional!
      x{contador + 1} = x{contador} + alpha*d{contador}';    
      grad{contador+1} = gradiente(x{contador+1});
      hessi{contador + 1} = hessiana(x{1}); % In the modified version of multidimensional newton method the hessian matrix is constant in X1
      contador = contador + 1;
  end 
resultado.nome = 'Modified Newton Multidimensional Method';
resultado.valor = x{contador};
result = resultado;
end

function result = fibonacci(d,x)

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
  metodo.fibonacci = determinar_N(metodo.a(k),metodo.b(k),l); % o erro esta no calculo de n
  metodo.n = length(metodo.fibonacci);
  metodo.lambda(k) = lambda(metodo.a(k),metodo.b(k),metodo.fibonacci,metodo.n,k);
  metodo.mi(k) = mi(metodo.a(k),metodo.b(k),metodo.fibonacci,metodo.n,k);
  
  metodo.f_lambda(k) = objetivo(metodo.lambda(k),d,x);
  
  metodo.f_mi(k) = objetivo(metodo.mi(k),d,x);
  
  while (k<metodo.n-1)      
      if (objetivo(metodo.lambda(k),d,x) > objetivo(metodo.mi(k),d,x))
        metodo.a(k+1) = metodo.lambda(k);
        metodo.b(k+1) = metodo.b(k);
        metodo.lambda(k+1) = metodo.mi(k);
        metodo.mi(k+1) = mi(metodo.a(k+1),metodo.b(k+1),metodo.fibonacci,metodo.n,k);
        
        metodo.f_mi(k+1) = objetivo(metodo.mi(k+1),d,x);
        
        metodo.f_lambda(k+1) =  metodo.f_mi(k+1);
      else
          metodo.a(k+1) = metodo.a(k);
          metodo.b(k+1) = metodo.mi(k);
          metodo.mi(k+1) = metodo.lambda(k);
          metodo.lambda(k+1) = lambda(metodo.a(k+1),metodo.b(k+1),metodo.fibonacci,metodo.n,k);
          
          metodo.f_lambda(k+1) = objetivo(metodo.lambda(k+1),d,x);
          
          metodo.f_mi(k+1) =  metodo.f_lambda(k+1);
      end
    k = k+1;
  end 
  metodo.lambda(k) = metodo.lambda(k-1);
  metodo.mi(k) = metodo.mi(k-1)+l;
  if (objetivo(metodo.lambda(k),d,x) < objetivo(metodo.mi(k),d,x))
      metodo.a(k+1) = metodo.a(k);
      metodo.b(k+1) = metodo.mi(k);
  else
      metodo.a(k+1) = metodo.lambda(k);
      metodo.b(k+1) = metodo.b(k);
  end
  metodo.numero_iteracoes = k;
  if objetivo(metodo.a(k),d,x) < objetivo(metodo.b(k),d,x)
      saida = metodo.a(k);
  else
      saida = metodo.b(k);
  end
  result = saida;
end

function result = objetivo(alpha, variavel, z)
    x = z(1) + alpha*variavel(1);
    y = z(2) + alpha*variavel(2);
    %result = x^2 + y^2 + x*y - 3*x;
    %result = 4*(x-5)^2+(y-6)^2;
    result = -12*y+4*x^2+4*y^2-4*x*y;
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
        lambda = a+(fibonacci(n-1-k)/fibonacci(n-k+1))*(b-a);
    result = lambda;
end

function result = hessiana(variavel)
    x=variavel(1);
    y=variavel(2);
    %result = [2, 1 ;1, 2];
    %result = [8,0;0,2];
    result = [8,-4;-4,8];
end

function result = gradiente(variavel)
    x=variavel(1);
    y=variavel(2);
    %result = [2*x+y-3;2*y+x] ;
    %result = [8*x-40;2*y-12] ;
    result = [8*x-4*y;8*y-4*x-12];
end

function result = mi(a,b,fibonacci,n,k)
        mi = a+(fibonacci(n-k)/fibonacci(n-k+1))*(b-a);
    result = mi;
end
