%_________________________________________________________________________%
%           SIMPLEX METHOD TO SOLVE LINEAR OPTIMIZATION PROBLEMS          %
%__*Developed by Joao Augusto Silva Ledo*_________________________________%

function result = Simplex()

    clear all;
    clc;
    format long;
    
    % Load the input data.
    data = loadInputData();
    z = data{1};
    C = data{2};
    answer = data{3};
    
    % Return the solution after loading the simplex method with the input data.
    resposta = solveSimplex(z,C,answer);
    result = resposta;
end

 function result = loadInputData()% Load input data
%_________________________________________________________________________%
%                               TEST 1
%         x1.valor = [7;3];
%         x1.nome = 'x1';
%         x2.valor = [-5;2];
%         x2.nome = 'x2';
%         x3.valor = [1;0];
%         x3.nome = 'x3';
%         x4.valor = [0;1];
%         x4.nome = 'x4';
% 
%         y1.valor = -5;
%         y1.nome = 'x1';
%         y2.valor = 1;
%         y2.nome = 'x2';
%         y3.valor = 0;
%         y3.nome = 'x3';
%         y4.valor = 0;
%         y4.nome = 'x4';
% 
%         answer = [13;17];
%         C = {y1,y2,y3,y4};
%         z = {x1,x2,x3,x4};
%_________________________________________________________________________%
%                               TEST 2
        x1.valor = [10;1;1;0];
        x1.nome = 'x1';
        
        x2.valor = [8;1;0;1];
        x2.nome = 'x2';
        
        x3.valor = [1;0;0;0];
        x3.nome = 'x3';
        x4.valor = [0;1;0;0];
        x4.nome = 'x4';
        x5.valor = [0;0;1;0];
        x5.nome = 'x5';
        x6.valor = [0;0;0;1];
        x6.nome = 'x6';
        
        y1.valor = -100;
        y1.nome = 'x1';
        y2.valor = -50;
        y2.nome = 'x2';
        y3.valor = 0;
        y3.nome = 'x3';
        y4.valor = 0;
        y4.nome = 'x4';
        y5.valor = 0;
        y5.nome = 'x5';
        y6.valor = 0;
        y6.nome = 'x6';
       
        answer = [25000;4500;1500;6000];
        C = {y1,y2,y3,y4,y5,y6};
        z = {x1,x2,x3,x4,x5,x6};
%_________________________________________________________________________%
%                                TEST 3
%         x1.valor = [1;1;0];
%         x1.nome = 'x1';
%         x2.valor = [1;0;1];
%         x2.nome = 'x2';
%         x3.valor = [1;0;0];
%         x3.nome = 'x3';
%         x4.valor = [0;1;0];
%         x4.nome = 'x4';
%         x5.valor = [0;0;1];
%         x5.nome = 'x5';
%         
%         y1.valor = -2;
%         y1.nome = 'x1';
%         y2.valor = -1;
%         y2.nome = 'x2';
%         y3.valor = 0;
%         y3.nome = 'x3';
%         y4.valor = 0;
%         y4.nome = 'x4';
%         y5.valor = 0;
%         y5.nome = 'x5';
%         
%         answer = [4;3;7/2];
%         C = {y1,y2,y3,y4,y5};
%         z = {x1,x2,x3,x4,x5};
%_________________________________________________________________________%    
     result = {z,C,answer};
 end

function result = solveSimplex(z,C,answer)
    B = constroiB1(z); 
    N = constroiN1(z,B);    
    CB = constroiCBeCN(C,B);
    CN = constroiCBeCN(C,N);
    lambda = calculaLambda(CB,B);
    CR = encontraCR(CN,lambda,N);
    sw = swap(B,N,CR,answer);
    interacoes = 1;
    while (verificanegatividade(CR) == true) && (sw.regiaoilimitada == false)          
          B = sw.B;
          N = sw.N;
          CB = constroiCBeCN(C,B); % Basic cost
          CN = constroiCBeCN(C,N); % Not basic cost
          lambda = calculaLambda(CB,B); % Simplex multiplier 
          CR = encontraCR(CN,lambda,N); % Relative Cost
          sw = swap(B,N,CR,answer); % Swap function between basic and not basic
          interacoes = interacoes + 1;
    end  
   if sw.regiaoilimitada == false 
        resolution = basicsolution(answer,B);
        resposta.Nome = resolution{2};
        for i = 1 : length(resolution{1})
            resposta.Valor(i) = resolution{1}(i);
        end
        resposta.Interacoes = interacoes;   
   else
       resposta = '- Infinito';
   end
    result = resposta;
end

function result = swap(b,n,cr,answer)
    outsideB = descobrequemsaideB(b,n,cr,answer);% Returns the output name from B
    entraEmB = descobrequementraemB(cr); % Returns the input name of B
    saiDeB = outsideB.variavel;    
    for i = 1 : length(b{2})
        if b{2}{i} == saiDeB % Searching in B the taking out values of B
            outposition = i;
        end
    end
    for j = 1 : length(n{2})
        if n{2}{j} == entraEmB
            inposition = j;
        end
    end
    aux = b{1}(:,outposition);
    b{1}(:,outposition) = n{1}(:,inposition);
    n{1}(:,inposition) = aux;   
    b{2}{outposition} = entraEmB;
    n{2}{inposition} = saiDeB;
    resultado.B = b;
    resultado.N = n;
    resultado.regiaoilimitada = outsideB.regiaoilimitada;
    result = resultado;
end

function result = direcaosimplex(b,n,cr)
    a = [];   
    for i = 1 : length(n{2})
        if n{2}{i} == descobrequementraemB(cr) % Serach for a name of a value taking out of N and taking in B to calculate y
            a = n{1}(:,i);
        end
    end
    y.valor = inv(b{1})*a;
    y.ilimitado = false;
   if (min(y.valor) <= 0) && (max(y.valor) <= 0) % If all y values are negative, this is a limited problem
       y.ilimitado = true;
   end    
    result = y;
end

function result = basicsolution(answer,b)
    resposta{1} = inv(b{1})*answer;
    resposta{2} = b{2};    
    result = resposta;
end
    
function result = descobrequemsaideB(b,n,cr,answer)
    solucao_basica = basicsolution(answer,b);
    y = direcaosimplex(b,n,cr);    
    for i = 1 : length(y.valor) % Step size
        if y.valor(i) < 0
            E.valor(i) = abs(solucao_basica{1}(i)/y.valor(i))*1000;
            E.nome{i} = solucao_basica{2}{i};
        else
            E.valor(i) = solucao_basica{1}(i)/y.valor(i);
            E.nome{i} = solucao_basica{2}{i};
        end
    end     
    [minimum,position] = min(E.valor);
    resposta.variavel = E.nome{position};
    resposta.regiaoilimitada = y.ilimitado;
    result = resposta;
end

function result = descobrequementraemB(cr)
    [minimum,position] = min(cr{1});
    resultado = cr{2}{position};
    result = resultado;
end

function result = verificanegatividade(cr)
    flag = false;
    for i = 1 : length(cr{1})
       if cr{1}(i) < 0
           flag = true;
       end
    end
result = flag;
end

 function result = constroiB1(z)
     Ident = identidade(z);
     i = 1;
         for (k = 1 : length(z))
            while (i <= length(z{1}.valor))
                for j = 1 : length(z{1}.valor)
                    B1{j} = buscaigual(Ident(:,i) , z);
                    i = i + 1;
                end
            end
         end
        B = encontra(B1);
     result = B;
 end
 
 function result = constroiN1(z,b)
    flag = 1;
    N1 = {};
      for i=1 : length(z)
          for j = 1 : length(b{2})
              if z{i}.nome == b{2}{j}
                  flag = flag + 1;
              end
          end
          if ((flag < ( length(z) - length(b{1})))) % flag ~= length(b{2}) && maybe use ~= 2
              N1 = [N1,z{i}];
          end
          flag = 1;
      end
    N = encontra(N1);
     result = N;
  end
  
  function result = constroiCBeCN(c,b)
    [line,column] = size(b{1});
    for i = 1 : length(c)
        for j = 1 : column
            if c{i}.nome == b{2}{j};
                CB1{j} = c{i};
            end
        end
    end
    CB = encontra(CB1);   
    result = CB;
  end
  
  function result = calculaLambda(cb,b)
    lambda{1} = cb{1}*inv(b{1});
    lambda{2} = cb{2};
    result = lambda;
  end
  
  function result = encontraCR(CN,lambda,N)
    for i = 1 : length(CN{1})
        CR1(i) = CN{1}(i)-lambda{1}*N{1}(:,i);
    end
    CR{1} = CR1;
    CR{2} = CN{2};
    result = CR; 
  end
  
 function result = identidade(z) % Creates the identity matrix with same size of the problem
    for i = 1 : length(z{1}.valor)
        for j = 1 : length(z{1}.valor)
            if i == j
                ident(i,j) = 1;
            end
        end
    end
     result = ident;
 end
 
 function result = buscaigual(a,b)
    for i = 1 : length(b)
        if a(:) == b{i}.valor;
            resultado = b{i};
        end
    end    
    result = resultado;
 end
 
 function result = encontra(x)  
    for j = 1 : length(x)
        nome{j} = x{j}.nome;
        for k = 1 : length(x{j}.valor)
            mat(k,j) = x{j}.valor(k);
        end 
    end
    result = {mat,nome};
 end
