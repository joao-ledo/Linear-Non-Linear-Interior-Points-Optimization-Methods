%_________________________________________________________________________%
%                             BARRIER METHOD                              %
%__*Developed by Joao Augusto Silva Ledo*_________________________________%

function result = Barreira()
  clear all;
  clc; 
  
  epslon = 1*10^-3;
  contador = 1;
  x{contador} = [0,1];
  mi(contador) = 10;
  beta = 0.1;
  
  grad{contador} = gradiente(x{contador},mi(contador));
  hessi{contador} = hessiana(x{contador},mi(contador));
  verificador = max(abs(grad{contador}));
  while(verificador >= epslon)
      verificador = max(abs(grad{contador}));
      d{contador} = - inv(hessi{contador}) * grad{contador};
      alpha = fibonacci(d{contador}, x{contador},mi(contador)); % Unidimensional Search
      x{contador + 1} = x{contador} + alpha*d{contador}';
      %B = barreirainversa(x{contador+1});
      B = barreiralogaritimica(x{contador+1});   
      grad{contador+1} = gradiente(x{contador+1},mi(contador));
      hessi{contador + 1} = hessiana(x{contador+1},mi(contador));
      
      if (mi(contador)*B < epslon) % Determine when to stop based on finding the feasible point that satisfies the constraint
          verificador = 1*10^-14;
      end
      
      mi(contador+1) = beta*mi(contador);
      contador = contador + 1;
  end 
resultado.nome = 'Metodo de Barreira';
resultado.valor = x{contador};
resultado.mi = mi(contador);
resultado.iteracoes = contador;
result = resultado;
end

function result = calculo_mi(a,b,fibonacci,n,k)
        mi = a+(fibonacci(n-k)/fibonacci(n-k+1))*(b-a);
    result = mi;
end

function result = hessiana(variavel,mi)
    x1=variavel(1);
    x2=variavel(2);
    result = [ 12*(x1 - 2)^2 + 2+mi*(2/(x2 - x1^2)^2 + (8*x1^2)/(x2 - x1^2)^3), -4+mi*( -(4*x1)/(x2 - x1^2)^3) ; -4+mi*(-(4*x1)/(x2 - x1^2)^3),  8+mi*(2/(x2 - x1^2)^3)];
end

function result = gradiente(variavel,mi)
    x1=variavel(1);
    x2=variavel(2);
    result = [2*x1 - 4*x2 + 4*(x1 - 2)^3+mi*((2*x1)/(x2 - x1^2)^2);8*x2 - 4*x1+mi*(-1/(x2 - x1^2)^2)] ;
end

function result = barreirainversa(x)   
    b = -1/(x(1)^2-x(2));
    result = b;
end

function result = barreiralogaritimica(x)   
    b = -log(-(x(1)^2-x(2)));
    result = b;
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
  metodo.fibonacci = determinar_N(metodo.a(k),metodo.b(k),l); % o erro esta no calculo de n
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
 %   result = (x1-2)^4 + (x1-2*x2)^2 + mi*(barreirainversa(x));
    result = (x1-2)^4 + (x1-2*x2)^2 + mi*(barreiralogaritimica(x));    
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
