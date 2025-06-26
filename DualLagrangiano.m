function result = DualLagrangiano()
  clear all;
  clc; 
  epslon = 0.0001;
  contador = 1;
  x{contador} = [0,0];  
  lambda(contador) = 1;
  grad{contador} = gradiente(x{contador},lambda(contador));
  while(max(abs(grad{contador})) >= epslon)
      
      d{contador}=-inv(hessiana(x{contador},lambda(contador)))*grad{contador};
      alpha = fibonacci(d{contador}, x{contador}, lambda(contador));
      x{contador + 1} = x{contador} + alpha*d{contador}';    
      restricao = calculoRestricaoIgualdade(x{contador + 1});
      lambda(contador+1) = atualizaLambda(lambda(contador),restricao);
      grad{contador+1} = gradiente(x{contador+1},lambda(contador+1));    
      contador = contador + 1;
  end 
  resultado.nome = 'Metodo Dual-Lagrangiano';
  resultado.valor = x{contador};
  resultado.iteracoes = contador;
result = resultado;
end

function result = hessiana(x,lambda)
    result = [ 2, 0 ; 0, 2];
end

function result = fibonacci(d,x,lambda)

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
  
  metodo.f_lambda(k) = objetivo(metodo.lambda(k),d,x,lambda);
  
  metodo.f_mi(k) = objetivo(metodo.mi(k),d,x,lambda);
  
  while (k<metodo.n-1)      
      if (objetivo(metodo.lambda(k),d,x,lambda) > objetivo(metodo.mi(k),d,x,lambda))
        metodo.a(k+1) = metodo.lambda(k);
        metodo.b(k+1) = metodo.b(k);
        metodo.lambda(k+1) = metodo.mi(k);
        metodo.mi(k+1) = mi(metodo.a(k+1),metodo.b(k+1),metodo.fibonacci,metodo.n,k);
        
        metodo.f_mi(k+1) = objetivo(metodo.mi(k+1),d,x,lambda);
        
        metodo.f_lambda(k+1) =  metodo.f_mi(k+1);
      else
          metodo.a(k+1) = metodo.a(k);
          metodo.b(k+1) = metodo.mi(k);
          metodo.mi(k+1) = metodo.lambda(k);
          metodo.lambda(k+1) = lambd(metodo.a(k+1),metodo.b(k+1),metodo.fibonacci,metodo.n,k);
          
          metodo.f_lambda(k+1) = objetivo(metodo.lambda(k+1),d,x,lambda);
          
          metodo.f_mi(k+1) =  metodo.f_lambda(k+1);
      end
    k = k+1;
  end 
  metodo.lambda(k) = metodo.lambda(k-1);
  metodo.mi(k) = metodo.mi(k-1)+l;
  if (objetivo(metodo.lambda(k),d,x,lambda) < objetivo(metodo.mi(k),d,x,lambda))
      metodo.a(k+1) = metodo.a(k);
      metodo.b(k+1) = metodo.mi(k);
  else
      metodo.a(k+1) = metodo.lambda(k);
      metodo.b(k+1) = metodo.b(k);
  end
  metodo.numero_iteracoes = k;
  if objetivo(metodo.a(k),d,x,lambda) < objetivo(metodo.b(k),d,x,lambda)
      saida = metodo.a(k);
  else
      saida = metodo.b(k);
  end
  result = saida;
end

function result = objetivo(alpha, variavel, z, lambda) % recebe alpha, recebe dire??o, recebe z = [x,y]
    x1 = z(1) + alpha*variavel(1);
    x2 = z(2) + alpha*variavel(2);
    result = x1^2 + x2^2 + lambda*(x1 + x2 - 1);
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

function result = gradiente(variavel,lambda)
    x1=variavel(1);
    x2=variavel(2);
    result = [2*x1+lambda;2*x2+lambda] ;
end

function result = mi(a,b,fibonacci,n,k)
        mi = a+(fibonacci(n-k)/fibonacci(n-k+1))*(b-a);
    result = mi;
end

function result = calculoRestricaoIgualdade(x)
    result =  x(1) + x(2) - 1;
end

function result = atualizaLambda(lambda,restricao)
    result = lambda+restricao;
end

