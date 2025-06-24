function result = NewtonDFP()
  clear all;
  clc; 
  epslon = 0.01;
  contador = 1;
  x{contador} = [0;0];  
  %x{contador} = [8;9];
  %x{contador} = [-1/2;1];
  s{contador} = [1,0;0,1];
  grad{contador} = gradiente(x{contador}');
  
  if (max(abs(grad{contador})) >= epslon)     
        d{contador}= - s{contador} * grad{contador};
        alpha(contador) = fibonacci(d{contador}, x{contador}');
        x{contador + 1} = x{contador} + alpha(contador)*d{contador};
  end

  while(max(abs(gradiente(x{contador}')) >= epslon) )      
%       if (((gradiente(x{contador+1}) - gradiente(x{contador}))' * ( (x{contador + 1}-x{contador}) - s{contador} * (gradiente(x{contador+1} - gradiente(x{contador}))))) <= 0)
%           s{contador + 1} = s{contador};
%           contador = contador + 1
%       else
          p{contador} = x{contador + 1}-x{contador};
          q{contador} = gradiente(x{contador+1}) - gradiente(x{contador});
          
          s{contador + 1} = s{contador} + ((p{contador}*p{contador}')/(p{contador}'*q{contador}))-((s{contador}*q{contador}*q{contador}'*s{contador})/(q{contador}'*s{contador}*q{contador}));
          contador = contador + 1;
          
%       end
      d{contador} = -s{contador}*gradiente(x{contador});

      alpha(contador) = fibonacci(d{contador}, x{contador}');
      
      x{contador + 1} = x{contador} + alpha(contador)*d{contador};      
      %grad{contador} = gradiente(x{contador}');
  end 
resultado.nome = 'Metodo de Davidon-Fletcher-Powell';
resultado.valor = x{contador}';
result = resultado;
end

function result = fibonacci(d,x)

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

function result = objetivo(alpha, variavel, z) % recebe alpha, recebe dire??o, recebe z = [x,y]
    x = z(1) + alpha*variavel(1);
    y = z(2) + alpha*variavel(2);
    result = x^2 + y^2 + x*y - 3*x;
    %result = 4*(x-5)^2+(y-6)^2;
    %result = -12*y+4*x^2+4*y^2-4*x*y;
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

function result = gradiente(variavel)
    x=variavel(1);
    y=variavel(2);
    result = [2*x+y-3;2*y+x] ;
    %result = [8*x-40;2*y-12] ;
    %result = [8*x-4*y;8*y-4*x-12];
end

function result = mi(a,b,fibonacci,n,k)
        mi = a+(fibonacci(n-k)/fibonacci(n-k+1))*(b-a);
    result = mi;
end