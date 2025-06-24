%_________________________________________________________________________%
%                        PROJECTED GRADIENT METHOD                        %
%__*Developed by Joao Augusto Silva Ledo*_________________________________%

function result = GradienteProjetado()
    clear all;
    clc; 
    k = 1;
    alphamax = [0,0];
    %v = [1,0;0,1];
    x{k} = [0;0];
    v{k} = [0,0,0,0];
    m{k} = funcaoAtiva(x{k},v{k});
    %I{k} = criaIdentidade(m{k});
    %P{k} = calculoP(I{k},m{k});%I{k}-m{k}'inv(m{k}*m{k}')*m{k};
    %gradiente{k} = grad(x{k},P{k});
    %d{k} = calculoD(P{k},gradiente{k});
    verifica = true;
    while verifica == true%(alphamax(1) <= alphamax(2)) || (verificaV(v) == 1) || (max(m{k}) == 0)

        m{k} = funcaoAtiva(x{k},v{k});
        
        identidade{k} = criaIdentidade(m{k});
        p{k} = calculoP(identidade{k},m{k});
        gradiente{k} = grad(x{k},p{k});
        d{k} = calculoD(p{k},gradiente{k});
        v{k} = -(inv(m{k}*m{k}'))*m{k}*gradiente{k};
        minimo = min(v{k});      
        if (minimo >= 0) || (alphamax(1) <= alphamax(2))
           verifica = false; 
        end
        
        if (d{k} == 0)  
           verifica = false;
        else
            alphamax = calculoAlpha(x{k},d{k});
            alpha(k) = verificaalpha(alphamax);
            for i = 1 : length(x{k})
               aux(i) = d{k}(i); 
            end
            d{k} = aux';
            x{k+1} = x{k} + alpha(k)*d{k};
        end
        k = k + 1;
    end

    metodo.nome = 'Projected Gradient';
    metodo.resultado = [x{k}(1),x{k}(2)];
    result = metodo;
end
function result = funcaoAtiva(x,v)

m = [0,0];
linhadeM(1).valor = [1,1];
linhadeM(1).checa = 0;
linhadeM(2).valor = [1,5];
linhadeM(2).checa = 0;
linhadeM(3).valor = [-1,0];
linhadeM(3).checa = 0;
linhadeM(4).valor =[0,-1];
linhadeM(4).checa = 0;
if x(1)+x(2)-2 <= 0
    linhadeM(1).checa = 1;
end
if x(1)+5*x(2)-5 <= 0
    linhadeM(2).checa = 1;
end
if -x(1) <= 0
    linhadeM(3).checa = 1;
end
if -x(2) <= 0
    linhadeM(4).checa = 1;
end

for q=1 : length(v)
    if (v(q) < 0)
        position(q) = q;
        position
    end
end




for i = 1 : length(linhadeM)
    for j = 1 : length(linhadeM(i).valor)
        if linhadeM(i).checa == 1
            m(i,j) = linhadeM(i).valor(j);
        end
    end
end
[l,c] = size(m);

if (l>c)
   for t = 1 : l
        for r=1 : c
            aux(t,r)=0;
        end
   end
  m = [m,aux]; 
end

if (c>l)
   for w = 1 : l
        for f=1 : c
            aux(w,f)=0;
        end
   end
   m = [m;aux];
end



    result = m;
end

function result = criaIdentidade(m)
    [l,c]=size(m);
    if l < c
        c=l;
    else
        l=c;
    end
    for i = 1 : l
        for j = 1 : c
            if i == j
               identidade(i,j) = 1; 
            else
               identidade(i,j) = 0;
            end
        end
    end
    result = identidade;
end

function result = calculoP(identidade,m)
    result = identidade - m'*inv(m * m')*m;
end

function result = grad(x,p)
gradiente = [4*x(1)-2*x(2)-4;4*x(2)-2*x(1)-6];
if(length(gradiente)<length(p))
    for i = 1 : (length(p)-length(gradiente))
        aux(i) = 0;
    end
    gradiente = [gradiente;aux'];
end
    result = gradiente;
end

function result = calculoD(p,grad)
    result = -p*grad;
end

function result = calculoAlpha(x,d)
    alpha = subs('alpha');   
    for i = 1 : length(x)
       aux(i) = d(i); 
    end
    d = aux';
    Xnew = x + alpha*d;
    x1 = Xnew(1);
    x2 = Xnew(2);
    g1 = x1+x2-2;
    g2 = x1+5*x2-5;
    g3 = -x1;
    g4 = -x2;
    alpha1 = solve(g1,alpha);
    alpha2 = solve(g2,alpha);
    alpha3 = solve(g3,alpha);
    alpha4 = solve(g4,alpha);
    alphamax1 = max([alpha1,alpha2,alpha3,alpha4]);
    alphamax2 = fibonacci(d,x);
    alphamax = [double(alphamax1),double(alphamax2)]; 
    %alpha =([2;5;0;0]\[1,1;1,5;-1,0;0,-1]-[x(1),x(2)])\[d(1),d(2)];
    result = alphamax;
end

function result = verificaV(v)
    verificador = 0;
    [l,c] = size(v);
    for i = 1 : l
        for j = 1 : c
            if (v(i,j) < 0)
                verificador = 1;
            end
        end
    end

    result = verificador;
end


function result = verificaalpha(alphamax)
    if alphamax(2) > alphamax(1)
        result = alphamax(1);
    else
        result = alphamax(2);
    end   
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
    x1 = z(1) + alpha*variavel(1);
    x2 = z(2) + alpha*variavel(2);
    result = 2*x1^2 + 2*x2^2 - 2*x1*x2 - 4*x1 - 6*x2;
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
