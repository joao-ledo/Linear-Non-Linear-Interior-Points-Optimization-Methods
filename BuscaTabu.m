function resultado = BuscaTabu()
    clear all;
    clc;
    format long;    
    dados = carregadados();
    resultado = resolveBuscaTabu(dados.G, dados.Pgmin, dados.Pgmax, dados.Pd, dados.a, dados.b, dados.c, dados.e, dados.f, dados.penalizacao,dados.tamanhoListaTabu, dados.numInteracoes, dados.tamanhoPasso, dados.tamanhoVizinhanca);
end

function resultado = carregadados()
    dados.G = 3;
    dados.Pgmin = [100,50,100];
    dados.Pgmax = [600,200,400];
    dados.Pd = 850;
    dados.a = [0.001562,0.00482,0.00194];
    dados.b = [7.92,7.97,7.85];
    dados.c = [561,78,310];
    dados.e = [300,150,200];
    dados.f = [0.0315,0.063,0.042];
    dados.penalizacao = 1000;
    dados.tamanhoListaTabu = 25;
    dados.numInteracoes = 100;
    dados.tamanhoPasso = (dados.Pgmax-dados.Pgmin)*0.005; 
    dados.tamanhoVizinhanca = 100;
    resultado=dados;
end

function resultado = resolveBuscaTabu(G, Pgmin, Pgmax, Pd, a, b, c, e, f, penalizacao,tamanhoListaTabu, numInteracoes, passo, tamanhoVizinhanca)
    k=1;
    j=1;
    Sbest = constroiSolucaoInicial(G, Pgmin, Pgmax);
    ListaTabu{1} = Sbest;
    para = false;        
    while((para == false) && (k <= numInteracoes))%abs(Fo(a,b,c,e,f,Pgmin,Pgmax,Pd,Sbest,G,penalizacao) - FuncaoAvaliacao(a,b,c,e,f,Pgmin,Pgmax,Pd,Sbest,G,penalizacao)) >= 10^-1)%(para == false)%(condicaoParada()~= true) (para == false) && (k <= numInteracoes) && (sum(Sbest) <= Pd+1) || (sum(Sbest) >= Pd-1)
        ListaCandidatos = {};
        %auz = constroiSolucaoInicial(G, Pgmin, Pgmax);
        for s = 1 : tamanhoVizinhanca
            Svizinhos{s} = CalculoVizinho(Sbest, passo, Pd, Pgmin, Pgmax, G);%  constroivizinhanca(Sbest,G); %Sbest*rand();
        end
        for i=1:length(Svizinhos)
            if(VerificaTabu(Svizinhos{i},ListaTabu) == false)
                ListaCandidatos{j} = Svizinhos{i};
                j=j+1;
            end
        end
        if((length(ListaCandidatos) ~= 0) )
            ListaCandidatos = ListaCandidatos(~cellfun('isempty',ListaCandidatos));
            Scandidato = BuscaMelhor(ListaCandidatos,a,b,c,e,f,Pgmin,Pgmax,Pd,G,penalizacao);
            if(Scandidato.avaliacao <= FuncaoAvaliacao(a,b,c,e,f,Pgmin,Pgmax,Pd,Sbest,G,penalizacao))
                Sbest = Scandidato.candidato;
            end
            if(VerificaTabu(Sbest,ListaTabu) == false)
                if(k <= tamanhoListaTabu)
                    ListaTabu{k+1} = Sbest;
                else
                    ListaTabu{tamanhoListaTabu+1} = Sbest;
                end
            end
            Elite{k} = BuscaMelhor(ListaTabu,a,b,c,e,f,Pgmin,Pgmax,Pd,G,penalizacao);       
            if(length(ListaTabu) > tamanhoListaTabu)
%                 ListaTabu{1} = [];
%                 ListaTabu = ListaTabu(~cellfun('isempty',ListaTabu));
              for m=1:tamanhoListaTabu
                  auxiliar{m} = ListaTabu{m+1};                     
              end
              ListaTabu = auxiliar;
            end

        else
            para = true;
        end
       k=k+1
      % MelhorElite(Elite)
    end
    
    resposta = MelhorElite(Elite);%BuscaMelhor(ListaTabu,a,b,c,e,f,Pgmin,Pgmax,Pd,G,penalizacao);
    resposta.SomaPg = sum(Sbest);
    resposta.Pd = Pd;
    resposta.Pgmin = Pgmin;
    resposta.Pgmax = Pgmax;
    resposta.Tabu = ListaTabu;
    resposta.Candidatos = ListaCandidatos;
    resposta.Vizinhos = Svizinhos;
    resposta.Iteracoes = k;
    resposta.TamanhoTabu = tamanhoListaTabu;
    resposta.Penalidade = penalizacao;
    resposta.Elite = Elite;
    resultado = resposta;
end

function result = MelhorElite(Elite)
    for i = 1 : length(Elite)
        aux(i) = Elite{i}.avaliacao;
        aux2{i} = Elite{i}.candidato;
    end
    [valor, posicao] = min(aux);
    resposta.candidato = aux2{posicao};
    resposta.avaliacao = valor;
    result = resposta;
end

function resultado = constroiSolucaoInicial(G, Pgmin, Pgmax)
    for i = 1:G
         Pg(i) = Pgmin(i) + ((Pgmax(i)-Pgmin(i)) * rand());
    end
    resultado = Pg;
end

% function resultado = condicaoParada()
%     
% 
% end

function resultado = constroivizinhanca(Sbest,G)
   for i = 1:G
       SbestBin{i} = dec2bin(Sbest(i));
       for j = 1:length(SbestBin{i})
        
        Svizinhos{j}(i) = movimentos(str2num(SbestBin{i}),j);
       end
   end
   resultado = Svizinhos;
end

function resultado = movimentos(PosicaoSbestBin,j)

     for i = 1: length(PosicaoSbestBin)
        if(i==j)
            if(PosicaoSbestBin(i) == 1)
                PosicaoSbestBin(i) = 0;
            else
                PosicaoSbestBin(i) = 1;
            end
        end
     end
    
    resultado = bin2dec(num2str(PosicaoSbestBin));
end

function resultado = FuncaoAvaliacao(a,b,c,e,f,Pgmin,Pgmax,Pd,Pg,G,penalizacao)
    penalizaMin = 0;
    penalizaMax = 0;
    for i = 1:G
        if(Pg(i) < Pgmin(i))
           penalizaMin =  penalizacao;
           if(Pg(i) > Pgmax(i))
               penalizaMax = penalizacao;
           end
        end
       FuncaoAvaliacao(i) = a(i)*Pg(i)*Pg(i) + b(i)*Pg(i) + c(i) + abs(e(i)*sin(f(i)*(Pgmin(i)-Pg(i)))) + penalizaMin*abs((Pgmin(i)-Pg(i)))+ penalizaMax*abs((Pg(i)-Pgmax(i)));
    end    
    resultado = sum(FuncaoAvaliacao) + penalizacao*abs((sum(Pg)-Pd));
end

function resultado = Fo(a,b,c,e,f,Pgmin,Pgmax,Pd,Pg,G,penalizacao)
    for i = 1:G
        if(Pg(i) < Pgmin(i))
           penalizaMin =  penalizacao;
           if(Pg(i) > Pgmax(i))
               penalizaMax = penalizacao;
           end
        end
       FuncaoAvaliacao(i) = a(i)*Pg(i)*Pg(i) + b(i)*Pg(i) + c(i) + abs(e(i)*sin(f(i)*(Pgmin(i)-Pg(i))));
    end   
    resultado = sum(FuncaoAvaliacao);
end

function result = CalculoVizinho(Sbest, passo, pd, pgiMin, pgiMax, G) % passo dado para os dados vizinhos, o tamanho desse passo ? c?lculado pelo minmax, que ? um percentual m?dio tirado do valor que se deseja encontrar para o lado (vizinhos)
  position = size(Sbest);
  vetor = size(size(Sbest));
  cont = 0;
  totalpercent = 100.0;
   for i = 1: length(Sbest)
    if(i == length(pgiMin)) % Se i estiver apontando para a ultima posi??o fa?a:
        vetor(i) = pd * (totalpercent/100); % adicionar o percentual de restri??o de min?mos e m?ximos
    else
      minimo = max(pgiMin(i), Sbest(i)-passo(i));
      maximo = min(pgiMax(i), Sbest(i)+passo(i)); 
      aleatorio = rand_in_bounds((minimo*100.0)/pd, (maximo*100.0)/pd);  % Respons?vel por chutar valores convenientes 
      totalpercent = totalpercent - aleatorio ;
      vetor(i) = pd * (aleatorio/100);
    end
   end
        position = vetor;
  result = position;
end

function result = rand_in_bounds(min, max)
  result = min + ((max-min) * rand());
end

function resultado = BuscaMelhor(ListaCandidatos,a,b,c,e,f,Pgmin,Pgmax,Pd,G,penalizacao)
    for i=1:length(ListaCandidatos)
       avaliacao(i) = FuncaoAvaliacao(a,b,c,e,f,Pgmin,Pgmax,Pd,ListaCandidatos{i},G,penalizacao);
    end
    [melhorvalor,posicao] = min(avaliacao);
       
    Scandidato.candidato = ListaCandidatos{posicao}; %struct
    Scandidato.avaliacao = melhorvalor;
    resultado = Scandidato;
end

function resultado = VerificaTabu(Svizinhos,ListaTabu)
    controle = 0;
       for j = 1:length(ListaTabu)
            if(sum(ismember(Svizinhos , ListaTabu{j})) == length(Svizinhos))
                controle = controle+1;           
            end
        end
        if(controle ~= 0)
            resposta = true;
        else
            resposta = false;
        end
   
    resultado = resposta;
end