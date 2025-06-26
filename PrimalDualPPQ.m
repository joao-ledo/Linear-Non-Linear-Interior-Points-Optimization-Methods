%_________________________________________________________________________%
%                                                                         %
%                                                                         %
%                        PRIMAL DUAL AFIM ALGORITHM                       %
%                           (PPQ, PNL, Rosenbrock)                        %
%                                                                         %
%                                          Developed by:                  % 
%                                                 Joao Augusto Silva Ledo %
%_________________________________________________________________________%

function result = PrimalDualPPQ()
    clear all;
    clc;
    format long;      
    data = loadInputData(); % Carregar os dados de entrada.
    result = solveDualAfim(data.A, data.C, data.pontoInicial, data.Erro, data.e, data.Sinicial, data.b, data.Mi, data.w, data.Q, data.C_quadratic, data.constant); % Cria a Instancia que resolve juntamento com os dados
end

%_________________________________________________________________________%
 function result = loadInputData() % Local onde se carregam as informacoes
        
%         resultado.Erro{1} = 10^-2;
%         resultado.Erro{2} = 10^-2;
%         resultado.Erro{3} = 10^-2;  
%         resultado.C = [1; 2; -3];
%         resultado.C_quadratic = [2; 3; 5]; %<------------- TP8 Normal
%         resultado.constant = [0; 0; 0];
%         resultado.A = [1, 1, 0; 0, 1, 1];
%         resultado.b = [5; 10];
%         resultado.Q = [4, 0, 0; 0, 6, 0; 0, 0, 10];
%         resultado.pontoInicial{1} = -inv(resultado.Q) * resultado.C;
%         resultado.pontoInicial{2} = [1; 1; 1]; % o otimo ?: [0; 5; 5]
%         resultado.w = [-14; 40]; %[-15; 47]
%         resultado.Sinicial = [18; 1; 1];%[16; 0;0];
        
%         resultado.Erro{1} = 10^-2;
%         resultado.Erro{2} = 10^-2;
%         resultado.Erro{3} = 10^-2;  
%         resultado.C = [4; -6; 0; 0; 0; 0];
%         resultado.C_quadratic = [2; 3; 0; 0; 0; 0];
%         resultado.constant = [0; 0; 0; 0; 0; 0];
%         resultado.A =  [-1, 1, 1, 0, 0, 0; -1, -1, 0, 1, 0, 0;1, 0, 0, 0, 1, 0; 0, 1, 0, 0, 0, 1]; %<------------- PROVA %Mehotra Exercicio 1
%         resultado.b = [1; -1; 3; 2]; 
%         resultado.Q = [4, 0, 0, 0, 0, 0; 0, 6, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0];
%         resultado.pontoInicial{1} = -inv(resultado.Q) * resultado.C;
%         resultado.pontoInicial{2} = [2.5; 1; 2.5; 2.5; 0.5; 1];
%         resultado.w = [0; 0; 0; 0];
%         resultado.Sinicial = [1; 1; 1; 1; 1; 1];%[2.5; 2.5; 0.5; 1; 2.5; 1];%[1; 2; 1.5; 0.5; 1.5; 1.5];%[1; 1; 2; 1; 1; 1]; %           [2.5; 2.5; 0.5; 1; 2.5; 1];

        
%         resultado.Erro{1} = 10^-2;
%         resultado.Erro{2} = 10^-2;
%         resultado.Erro{3} = 10^-2;          
%         resultado.C = [7.92; 7.97; 7.85; 0; 0; 0; 0; 0; 0]; %->b
%         resultado.C_quadratic = [0.001562; 0.00482; 0.001940; 0; 0; 0; 0;0; 0]; %->a
%         resultado.constant = [561; 310; 78; 102; 51; 178; 0; 0; 0; 0; 0; 0];      %->c
%         resultado.A = [1, 1, 1, 0, 0, 0, 0, 0, 0;                                   %<-------------Problema de Despacho Economico com 3 geradores
%                        -1, 0, 0, 1, 0, 0, 0, 0, 0;
%                        1, 0, 0, 0, 1, 0, 0, 0, 0;
%                        0, -1, 0, 0, 0, 1, 0, 0, 0;
%                        0, 1, 0, 0, 0, 0, 1, 0, 0;
%                        0, 0, -1, 0, 0, 0, 0, 1, 0;
%                        0, 0, 1, 0, 0, 0, 0, 0, 1];
%         resultado.b = [850; -100; 600; -50; 200; -100; 400];
%         resultado.Q = [0.003124, 0, 0, 0, 0, 0, 0, 0, 0; 
%                        0, 0.00964, 0, 0, 0, 0, 0, 0, 0;
%                        0, 0, 0.00388, 0, 0, 0, 0, 0, 0;
%                        0, 0, 0, 0, 0, 0, 0, 0, 0;
%                        0, 0, 0, 0, 0, 0, 0, 0, 0;
%                        0, 0, 0, 0, 0, 0, 0, 0, 0;
%                        0, 0, 0, 0, 0, 0, 0, 0, 0;
%                        0, 0, 0, 0, 0, 0, 0, 0, 0;
%                        0, 0, 0, 0, 0, 0, 0, 0, 0];
%         resultado.pontoInicial{1} = inv(-resultado.Q) * resultado.C;
%         resultado.pontoInicial{2} = [500; 150; 200; 400; 100; 100; 50; 100; 200];
%         resultado.w = [0; 0; 0; 0; 0; 0; 0];
%         resultado.Sinicial = [1; 1; 1; 1; 1; 1; 1; 1; 1];

        resultado.Erro{1} = 10^-3;
        resultado.Erro{2} = 10^-3;
        resultado.Erro{3} = 10^-2;
        resultado.constant = [550; 309; 307; 240; 240; 240; 240; 240; 240; 126; 126; 126; 126; 0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0]; %->c
        resultado.C = [8.100000000000000;8.100000000000000;8.100000000000000;7.740000000000000;7.740000000000000;7.740000000000000;7.740000000000000;7.740000000000000;7.740000000000000;8.600000000000000; 8.600000000000000; 8.600000000000000;8.600000000000000;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0]; %->b
        resultado.C_quadratic = [0.00028; 0.00056; 0.00056; 0.00324; 0.00324; 0.00324; 0.00324; 0.00324; 0.00324; 0.00284; 0.00284; 0.00284; 0.00284; 0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0]; %->a  
        resultado.A =  [1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0; %<-------------Problema de Despacho Economico com 13 geradores
                        -1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
                        1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
                        0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
                        0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
                        0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
                        0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
                        0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
                        0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
                        0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
                        0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
                        0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
                        0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
                        0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0;
                        0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0;
                        0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0;
                        0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0;
                        0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0;
                        0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0;
                        0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0;
                        0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0;
                        0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0;
                        0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0;
                        0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0;
                        0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0;
                        0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0;
                        0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1];                                     
        resultado.b = [2520;0;680;0;360;0;360;-60;180;-60;180;-60;180;-60;180;-60;180;-60;180;-40;120;-40;120;-55;120;-55;120];
        resultado.Q = [0.000560000000000000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
                       0,0.00112000000000000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
                       0,0,0.00112000000000000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
                       0,0,0,0.00648000000000000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
                       0,0,0,0,0.00648000000000000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
                       0,0,0,0,0,0.00648000000000000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
                       0,0,0,0,0,0,0.00648000000000000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
                       0,0,0,0,0,0,0,0.00648000000000000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
                       0,0,0,0,0,0,0,0,0.00648000000000000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
                       0,0,0,0,0,0,0,0,0,0.00568000000000000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
                       0,0,0,0,0,0,0,0,0,0,0.00568000000000000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
                       0,0,0,0,0,0,0,0,0,0,0,0.00568000000000000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
                       0,0,0,0,0,0,0,0,0,0,0,0,0.00568000000000000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
                       0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
                       0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
                       0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
                       0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
                       0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
                       0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
                       0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
                       0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
                       0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
                       0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
                       0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
                       0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
                       0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
                       0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
                       0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
                       0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
                       0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
                       0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
                       0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
                       0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
                       0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
                       0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
                       0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
                       0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
                       0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
                       0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];                                   
        resultado.pontoInicial{1} = inv(-resultado.Q) * resultado.C;
        resultado.pontoInicial{2} = [659;321;325;150;145;155;165;155;155;65;65;75;85;600;80;300;60;300;60;90;30;90;30;90;30;90;30;90;30;90;30;65;15;65;15;50;15;50;15];
        resultado.w = [0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0];
        resultado.Sinicial = [30;30;30;30;30;30;30;30;30;30;30;30;30;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1];

          % Com Rosenbrock
%         resultado.Erro{1} = 10^-2;
%         resultado.Erro{2} = 10^-2;
%         resultado.Erro{3} = 10^-2;  
%         resultado.C = [-2; 0; 0; 0; 0];
%         resultado.C_quadratic = [0; 0; 0; 0; 0];
%         resultado.constant = [0; 0; 0; 0; 0];
%         resultado.A = [1, -1, 1, 0, 0; 0, -1, 0, 1, 0; 0, 1, 0, 0, 1];  %<------------- PROVA %Mehotra Exercicio 2
%         resultado.b = [15; -2; 15];
%         resultado.pontoInicial{2} = [2; 4; 17; 2; 11];%[10; 4; 9; 2; 11];%[10; 3; 8; 1; 12];
%         resultado.w = [0; 0; 0];
%         resultado.Sinicial = [1; 1; 1; 1; 1];
%         resultado.Q = Hessian(resultado.pontoInicial{2});
%         resultado.pontoInicial{1} = inv(-resultado.Q) * resultado.C;
      
        resultado.Mi = 1;
        resultado.e = constroiE(resultado.Sinicial);
     result = resultado;
 end
%_________________________________________________________________________% 

 function result = solveDualAfim(A, C, pontoInicial, Erro, e, SInicial, b, mi_ini, w_ini, Q_ini, C_quadratic, constant)
   resposta = ChooseSolutionMethod();
   if(resposta.answer1 == 2)
       if(resposta.answer2 == 5)
        MinMaxPrec = MinMaxBuscaUni();
        minBuscaUni = MinMaxPrec.min;
        maxBuscaUni = MinMaxPrec.max;
        precisao = MinMaxPrec.Precisao;
       else
        minBuscaUni = 0;
        maxBuscaUni = 0;
        precisao = 0;
       end
   else
    minBuscaUni = 0;
    maxBuscaUni = 0; 
    precisao = 0;
   end
   k = 1;
   flag = false;
   x{k} = pontoInicial{1};
   w{k} = w_ini;    
   BetaP = {};
   BetaD = {};
   Dw = {};
   Dx = {};
   Ds = {};
   PrevisorCorretor = {};
   Q{k} = Q_ini;
   mi{k} = mi_ini;
   Fo{k} = calculoFO(x{k}, C, C_quadratic, constant, resposta.answer1);
   grad{k} = CalculaGrad(resposta.answer1, Q{k}, x{k}, C);
   if((VerificaEquacao(A, x{k}, b) == true) && (verificanegatividade(x{k}) == false))
       resultado.Nome = VerificaNome(resposta);
       resultado.Fo = Fo;
       resultado.x = x;
       resultado.Iteracoe = k-1;
       flag = true;
   else
       x{k} = pontoInicial{2};
   end
%    if(verificanegatividade(Q{k}*x{k} + C - A'*w{k}) == false)
%       s{k} = Q{k}*x{k} + C - A'*w{k};
%    else
%        s{k} = SInicial;
%    end
   s{k} = SInicial;
   while (flag == false)           
       t{k} = constroiT(b, A, x{k});
       Q{k} = CalculaQ(Q_ini, x{k}, resposta.answer1);
       grad{k} = CalculaGrad(resposta.answer1, Q{k}, x{k}, C);
       hessian{k} = Hessian(x{k});
       u{k} = constroiU(resposta.answer1,C, A, w{k}, s{k}, Q{k}, x{k}, grad{k}); 
       X{k} = diag(x{k});
       S{k} = diag(s{k});
       v{k} = constroiV(mi{k}, e, X{k}, S{k});
       p{k} = constroiP(X{k}, v{k});
       teta{k} = constroiTETA(resposta.answer1 ,S{k}, X{k}, Q{k}, hessian{k});
       FO{k} =  calculoFO(x{k}, C, C_quadratic, constant, resposta.answer1);
       erro2{k} = Factbilidade(t{k}, b, Erro{2});
       erro3{k} = Factbilidade(u{k}, grad{k}, Erro{3});
       if((mi{k} < Erro{1})  && (erro2{k}.flag == true) && (erro3{k}.flag == true) || (isfinite(sum(x{k})) == false) || (k==200)) 
               resultado.Nome = VerificaNome(resposta);
               resultado.FO = FO;
               resultado.x = x;
               resultado.s = s;
               resultado.w = w;
               resultado.Dx = Dx;
               resultado.Ds = Ds;
               resultado.Dw = Dw; 
               resultado.mi = mi;
               resultado.BetaP = BetaP;
               resultado.BetaD = BetaD;
               resultado.Diagonal_S = S;
               resultado.Diagonal_X = X;
               resultado.Q = Q;
               resultado.t = t;
               resultado.u = u; 
               resultado.v = v;
               resultado.p = p;
               resultado.teta = teta;
               resultado.erro2 = erro2;
               resultado.erro3 = erro3;
               resultado.CtransposeDx = CtransposeDx;
               resultado.BtransposeDw = BtransposeDw;
               for(i = 1 : (length(x{k})/3))
                   soma(i) = x{k}(i);
               end
               resultado.soma = sum(soma);
               if(resposta.answer1 == 2)
                resultado.PrevisorCorretor =  PrevisorCorretor;
               end
               resultado.Iteracoe = k-1;
               resultado.TempoComputacional_segundos = cputime;
               flag = true;
       else
            Escolha{k} = EscolheMetodoResolucao(A, e, x{k}, w{k}, s{k}, t{k}, u{k}, v{k}, p{k}, teta{k}, X{k}, S{k}, resposta, k, Q{k}, C, mi{k}, C_quadratic, constant, minBuscaUni, maxBuscaUni, precisao);
            Dw{k} = Escolha{k}.Dw;
            Dx{k} = Escolha{k}.Dx;
            Ds{k} = Escolha{k}.Ds;
            CtransposeDx{k} = C'*Dx{k};
            BtransposeDw{k} = b'*Dw{k};
            if(resposta.answer1 == 2)
                PrevisorCorretor{k} = Escolha{k}.PrevCorre;
            end
            if((FactbilidadePrimal(t{k}, Dx{k}, C, 10^-12) == true) || (FactbilidadeDual(u{k}, Ds{k}, b, Dw{k}, 10^-12) == true))
                resultado.Nome = VerificaNome(resposta);
                resultado.x = 'Problema Ilimitado';
                resultado.Iteracoe = k-1;
                flag = true;
            else
                if((verificaErro(t{k}, Erro{2}) == true) && (verificaErro(Dx{k}, Erro{2}) == true) && (verificaErro(u{k}, Erro{3}) == true)&& (verificaErro(Ds{k}, Erro{3}) == true))
                   resultado.Nome = VerificaNome(resposta);
                   resultado.FO = FO;
                   resultado.x = x;
                   resultado.s = s;
                   resultado.w = w;
                   resultado.Dx = Dx;
                   resultado.Ds = Ds;
                   resultado.Dw = Dw; 
                   resultado.mi = mi;
                   resultado.BetaP = BetaP;
                   resultado.BetaD = BetaD;
                   resultado.Diagonal_S = S;
                   resultado.Diagonal_X = X;
                   resultado.Q = Q;
                   resultado.t = t;
                   resultado.u = u; 
                   resultado.v = v;
                   resultado.p = p;
                   resultado.teta = teta;
                   resultado.erro2 = erro2;
                   resultado.erro3 = erro3;
                   resultado.CtransposeDx = CtransposeDx;
                   resultado.BtransposeDw = BtransposeDw;
                   for(i = 1 : (length(x{k})/3))
                    soma(i) = x{k}(i);
                   end
                   resultado.soma = sum(soma);
                   if(resposta.answer1 == 2)
                    resultado.PrevisorCorretor = PrevisorCorretor;
                   end
                   resultado.Iteracoe = k-1;
                   resultado.TempoComputacional_segundos = cputime;
                   flag = true;
                else
                BetaP{k} = passoP(x{k}, Dx{k}, Q{k}, C, resposta.answer1);
                BetaD{k} = passoD(s{k}, Ds{k});
                x{k+1} = NextValue(x{k}, BetaP{k}, Dx{k});
                w{k+1} = NextValue(w{k}, BetaD{k}, Dw{k});
                s{k+1} = NextValue(s{k}, BetaD{k}, Ds{k});
                mi{k+1} = ProximaBarreira(mi{k});
                k = k + 1;
                end
            end    
       end      
   end
   result = resultado;
 end
 
function result = VerificaNome(resposta)
    if(resposta.answer1 == 1)
        resultado = 'Algoritmo Primal-Dual Afim Barreira Logaritma para PPQ';
    else
        if (resposta.answer1 == 2)
            if(resposta.answer2 == 1)
                resultado = 'Algoritmo Primal-Dual Afim com Previsor-Corretor de Mehotra para PPQ';
            else
                if(resposta.answer2 == 2)
                     resultado = 'Algoritmo Primal-Dual Afim com Previsor-Corretor de Tanabi para PPQ';
                else
                    if(resposta.answer2 == 3)
                        resultado = 'Algoritmo Primal-Dual Afim com Previsor-Corretor de Gondzio para PPQ';
                    else
                        if(resposta.answer2 == 4)
                            resultado = 'Algoritmo Primal-Dual Afim com Previsor-Corretor Misto para PPQ';
                        else
                            if(resposta.answer2 == 5)
                                resultado = 'Algoritmo Primal-Dual Afim com Previsor-Corretor Misto com Busca Unidimensional de Fibonacci para o Passo Previsor para PPQ';
                            end
                        end
                    end
                end
            end 
        else
            if (resposta.answer1 == 3)
                resultado = 'Algoritmo Primal-Dual Afim Barreira Logaritma para PNL (Rosenbrock)';
            end
        end
    end
    result = resultado;
end
% 
function result = CalculaGrad(answer, Q, x, c)
    if((answer == 1) || (answer == 2))
        resultado = Q*x+c;
    else
        if(answer == 3)
           resultado =  Grad(x);
        end
    end
    result = resultado;
end

function result = CalculaQ(Q, x, answer)
    if((answer == 1) || (answer == 2))
        resposta = Q;
    else if (answer == 3)
            resposta = Hessian(x);
        end
    end
    result = resposta;
end

function result = Hessian(x)
    result = [202, -400*x(2), 0, 0, 0; -400*x(2), 1200*x(2)^2-400*x(1), 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0];
end

function result = Grad(x)
    result =[202*x(1)-200*(x(2)^2)-2; 400*(x(2)^3)-400*x(1)*x(2);0; 0; 0];
end

function result = verificaErro(array, erro)
    flag = true;
    for(i = 1 : length(array))
       if(array(i) > erro)
           flag = false;
       end
    end
    result = flag;
end

function result = VerificaEquacao(A, x, b)
     aux1 = A*x;
     flag = true;
    for i = 1 : length(aux1)
       if(aux1(i) ~= b(i))
           flag = false;
       end
    end
    result = flag;
end

 function result = FactbilidadePrimal(t, Dx, c, epslon)
    flag = false;
    if((max(abs(t)) <= epslon) && (verificanegatividade(Dx) == false))
        flag = true;
    end
    result = flag;
 end

  function result = FactbilidadeDual(u, Ds, b, Dw, epslon)
    flag = false;
    if((max(abs(u)) <= epslon) && (verificanegatividade(Ds) == false))
        flag = true;
    end
    result = flag;
 end
 
 function result = Factbilidade(parametro1, parametro2, Erro)
    resultado.flag = false;
    resultado.valor = (max(abs(parametro1))/(max(abs(parametro2))+1));
    if((max(abs(parametro1))/(max(abs(parametro2))+1)) < Erro)
        resultado.flag = true;       
    end
    result = resultado;
 end
 
 function result = constroiT(b, A, x)
    result = b - A*x;
 end
 
 function result = constroiU(answer, c, A, w, s, Q, x, grad)
    if((answer == 1) || (answer == 2))
        resultado = Q*x + c - A'*w - s;
    else
        resultado = grad - A'*w - s;
    end
    result = resultado;
 end
 
 function result = constroiV(mi, e, X, S)
    result = mi*e - X*S*e;
 end
 
 function result = constroiP(X, v)
    result = inv(X)*v;
 end
 
 function result = constroiTETA(answer, S, X, Q, hessiana)
     if((answer == 1) || (answer == 2))
        resultado = inv(Q+inv(X)*S);
     else
         if(answer == 3)
             resultado = inv(hessiana + inv(X)*S);
         end
     end
     result = resultado;
 end
 
function result = calculoFO(X, C, C_quadratic, constant, answer)
    if((answer == 1) || (answer == 2))
        for i = 1 : length(X)
            fo(i) = X(i)*C(i) + (X(i)^2)*C_quadratic(i) + constant(i);
        end
        resposta = sum(fo);
    else if (answer == 3)
            resposta = 100*(X(1)- X(2)^2)^2 + (1- X(1))^2;
        end
    end

    result = resposta;
end
 
function result = constroiE(x)
    for i = 1 : length(x)
       aux(i) = 1; 
    end
    result = aux';
end

 function result = passoP(x, Dx, Q, c, answer)   
    betaP1_final = 0;
     for i = 1 : length(x)
        if Dx(i) < 0 
            betaP1(i) = -0.9995*x(i)/Dx(i);
        else
            betaP1(i) = 0;
        end
     end
     if(sum(abs(betaP1)) == 0)
         betaP1 = 1;
     else
        betaP1 = nonzeros(betaP1);
     end
     betaP1_final = min(betaP1);
     if((answer == 1) || (answer == 2))
        betaQ = - ((Dx') * (Q*x + c))/(Dx'*Q*Dx);
     else
         if(answer == 3)
             betaQ = - ((Dx') * (Grad(x)))/(Dx'*Hessian(x)*Dx);
         end
     end 
     if(betaQ > 0)
         resultado = min(betaP1_final, betaQ);
     else
         resultado = min(betaP1_final);
     end
    result = resultado;
 end
 
 function result = passoD(s, Ds)
     for i = 1 : length(s)
        if Ds(i) < 0 
            betaD(i) = -0.9995*s(i)/Ds(i);
        else
            betaD(i) = 0;
        end
     end
     if(sum(abs(betaD)) == 0)
         betaD = 1;
     else
        betaD = nonzeros(betaD);
     end
    result = min(betaD);
 end

 function result = Direcaow(A, teta, t, u, p)
    result = inv(A*teta*A')*(t + A*teta*(u - p));
 end
 
 function result = DirecaoX(teta, A, Dw, u, p) 
    result = teta*(A'*Dw - u + p);
 end

 function result = DirecaoS(X,p, S, Dx) 
    result = p - inv(X)*S*Dx;
 end
 
 function resultado = NextValue(valoratual, passo, direcao)
   resultado = valoratual + passo * direcao;
 end
 
 function result = verificanegatividade(S)
    flag = false;
    for i = 1 : length(S)
       if (S(i) <= 0)
           flag = true;
       end
    end
    result = flag;
 end
 
function result = ProximaBarreira(mi)
    result = mi/2;
end

function result = ProximaBarreiraPrevCorr(x, s)
    result = ((x'*s)/length(x));
end
 
function result = Previsor(A, e, x, w, s, t, u, v, teta, X, S, resposta, Q, C, mi, C_quadratic, constant, minBuscaUni, maxBuscaUni, precisao)
    if(resposta.answer2 == 1) %Mehotra
        resultado.Dw = inv(A*teta*A')*(t + A*teta*u);
        resultado.Dx = teta*(A'*resultado.Dw - u);
        resultado.Ds = -inv(X)*S*resultado.Dx;
        resultado.BetaP = passoP(x, resultado.Dx, Q, C, resposta.answer1);
        resultado.BetaD = passoD(s, resultado.Ds);
        resultado.x = NextValue(x, resultado.BetaP, resultado.Dx);
        resultado.w = NextValue(w, resultado.BetaD, resultado.Dw);
        resultado.s = NextValue(s, resultado.BetaD, resultado.Ds);
        resultado.mi = ProximaBarreiraPrevCorr(resultado.x, resultado.s);
        resultado.v = constroiV(resultado.mi, e, X, S) - diag(resultado.Dx) * diag(resultado.Ds) * e;
    end
    if(resposta.answer2 == 2) %Tanabi
        resultado.Dw = inv(A*teta*A')*(t + A*teta*u);
        resultado.Dx = teta*(A'*resultado.Dw - u);
        resultado.Ds = -inv(X)*S*resultado.Dx;
        resultado.BetaP = passoP(x, resultado.Dx, Q, C, resposta.answer1);
        resultado.BetaD = passoD(s, resultado.Ds);
        resultado.x = NextValue(x, resultado.BetaP, resultado.Dx);
        resultado.w = NextValue(w, resultado.BetaD, resultado.Dw);
        resultado.s = NextValue(s, resultado.BetaD, resultado.Ds);
        resultado.mi = ProximaBarreiraPrevCorr(resultado.x, resultado.s);
        resultado.v = constroiV(resultado.mi, e, X, S);
    end
    if(resposta.answer2 == 3) %Gondzio
        resultado.Dw = inv(A*teta*A')*(t + A*teta*u);
        resultado.Dx = teta*(A'*resultado.Dw - u);
        resultado.Ds = -inv(X)*S*resultado.Dx;
        resultado.BetaP = passoP(x, resultado.Dx, Q, C, resposta.answer1);
        resultado.BetaD = passoD(s, resultado.Ds);
        resultado.x = NextValue(x, resultado.BetaP, resultado.Dx);
        resultado.w = NextValue(w, resultado.BetaD, resultado.Dw);
        resultado.s = NextValue(s, resultado.BetaD, resultado.Ds);
        resultado.mi = resultado.x'*resultado.s;
        resultado.v = constroiV(resultado.mi, e, X, S) - diag(resultado.Dx) * diag(resultado.Ds) * e;
    end
    if(resposta.answer2 == 4) %Misto
        resultado.Dw = inv(A*teta*A')*(t + A*teta*(u - inv(X)*v));
        resultado.Dx = teta*(A'*resultado.Dw - u + inv(X)*v);
        resultado.Ds = inv(X)*(v - S*resultado.Dx);
        resultado.BetaP = passoP(x, resultado.Dx, Q, C, resposta.answer1);
        resultado.BetaD = passoD(s, resultado.Ds);
        resultado.x = NextValue(x, resultado.BetaP, resultado.Dx);
        resultado.w = NextValue(w, resultado.BetaD, resultado.Dw);
        resultado.s = NextValue(s, resultado.BetaD, resultado.Ds);
        resultado.mi = ProximaBarreiraPrevCorr(resultado.x, resultado.s);
        resultado.v = constroiV(resultado.mi, e, X, S) - diag(resultado.Dx) * diag(resultado.Ds) * e;
    end
    if(resposta.answer2 == 5) %Seminario    
        resultado.Dw = inv(A*teta*A')*(t + A*teta*(u - inv(X)*v));
        resultado.Dx = teta*(A'*resultado.Dw - u + inv(X)*v);
        resultado.Ds = inv(X)*(v - S*resultado.Dx);      
        resultado.BetaP = fibonacci(resultado.Dx,x, C, C_quadratic, constant, minBuscaUni, maxBuscaUni, precisao, resposta.answer1);
        resultado.BetaD = fibonacci(resultado.Ds,s, C, C_quadratic, constant, minBuscaUni, maxBuscaUni, precisao, resposta.answer1);
        resultado.x = NextValue(x, resultado.BetaP, resultado.Dx);
        resultado.w = NextValue(w, resultado.BetaD, resultado.Dw);
        resultado.s = NextValue(s, resultado.BetaD, resultado.Ds);
        resultado.mi = ProximaBarreiraPrevCorr(resultado.x, resultado.s);
        resultado.v = constroiV(resultado.mi, e, X, S) - diag(resultado.Dx) * diag(resultado.Ds) * e;        
    end
    result = resultado;
end

function result = Corretor(A, e, x, w, s, t, u, v, teta, X, S, resposta, k, p, Q, C, mi, C_quadratic, constant, minBuscaUni, maxBuscaUni, precisao)
    if(resposta.answer2 == 1) %Mehotra
        previsor = Previsor(A, e, x, w, s, t, u, v, teta, X, S, resposta, Q, C, mi, C_quadratic, constant, minBuscaUni, maxBuscaUni, precisao);
        resultado.Dw = inv(A*teta*A')*(t + A*teta*(u - inv(X)*previsor.v));
        resultado.Dx = teta*(A'*resultado.Dw-u+inv(X)*previsor.v);
        resultado.Ds = inv(X)*(previsor.v-S*resultado.Dx);
        resultado.Prev = previsor;
    end
    if(resposta.answer2 == 2)%Tanabi 
        previsor = Previsor(A, e, x, w, s, t, u, v, teta, X, S, resposta, Q, C, mi, C_quadratic, constant, minBuscaUni, maxBuscaUni, precisao);
        if(rem(k,2) ~= 0)       
            resultado.Dw = inv(A*teta*A')*(t + A*teta*(u - inv(X)*previsor.v));
            resultado.Dx = teta*(A'*resultado.Dw-u+inv(X)*previsor.v);
            resultado.Ds = inv(X)*(previsor.v-S*resultado.Dx);
            resultado.Prev = previsor;
        else
            resultado.Dw = previsor.Dw;
            resultado.Dx = previsor.Dx;
            resultado.Ds = previsor.Ds;
            resultado.Prev = previsor;
        end
        
    end
    if(resposta.answer2 == 3) %Gondzio
        previsor = Previsor(A, e, x, w, s, t, u, v, teta, X, S, resposta, Q, C, mi, C_quadratic, constant, minBuscaUni, maxBuscaUni, precisao);
        correction.Dw = inv(A*teta*A')*(t + A*teta*(u - inv(X)*previsor.v));
        correction.Dx = teta*(A'*correction.Dw-u+inv(X)*previsor.v);
        correction.Ds = inv(X)*(previsor.v-S*correction.Dx);      
        correction.BetaP = passoP(x, correction.Dx, Q, C, resposta.answer1);
        correction.BetaD = passoD(s, correction.Ds);
        correction.x = NextValue(x, correction.BetaP, correction.Dx);
        correction.w = NextValue(w, correction.BetaD, correction.Dw);
        correction.s = NextValue(s, correction.BetaD, correction.Ds);
        correction.mi = correction.x'*correction.s;
        if(previsor.mi < correction.mi)
            resultado.Dw = previsor.Dw;
            resultado.Dx = previsor.Dx;
            resultado.Ds = previsor.Ds;
            resultado.Prev = previsor;
        else
            resultado.Dw = correction.Dw;
            resultado.Dx = correction.Dx;
            resultado.Ds = correction.Ds;
            resultado.Prev = previsor;
        end
    end
    if(resposta.answer2 == 4) %Misto
        previsor = Previsor(A, e, x, w, s, t, u, v, teta, X, S, resposta, Q, C, mi, C_quadratic, constant, minBuscaUni, maxBuscaUni, precisao);
        resultado.Dw = inv(A*teta*A')*(t + A*teta*(u - inv(X)*previsor.v));
        resultado.Dx = teta*(A'*resultado.Dw-u+inv(X)*previsor.v);
        resultado.Ds = inv(X)*(previsor.v-S*resultado.Dx);
        resultado.Prev = previsor;
    end
    if(resposta.answer2 == 5) %Seminario
        previsor = Previsor(A, e, x, w, s, t, u, v, teta, X, S, resposta, Q, C, mi, C_quadratic, constant, minBuscaUni, maxBuscaUni, precisao);
        resultado.Dw = inv(A*teta*A')*(t + A*teta*(u - inv(X)*previsor.v));
        resultado.Dx = teta*(A'*resultado.Dw-u+inv(X)*previsor.v);
        resultado.Ds = inv(X)*(previsor.v-S*resultado.Dx);
        resultado.Prev = previsor;
    end
    result = resultado;
end

function result = MinMaxBuscaUni()
    resultado.min = input('Digite o valor de Minimo da Busca Unidimensional para o Passo Previsor:\n');
    resultado.max = input('Digite o valor de Maximo da Busca Unidimensional para o Passo Previsor:\n');
    resultado.Precisao = input('Digite o valor da Precisao da Busca Unidimensional para o Passo Previsor:\n');
    result = resultado;
end
 
function result = ChooseSolutionMethod()
    resultado.answer1 = input('Escolha o Metodo de Resolucao: \n 1-Primal-Dual PPQ \n 2-Primal Dual com Previsor-Corretor \n 3-Primal-Dual PNL (Rosenbrock) \n');
    resultado.answer1 = VerificaOpcao(resultado.answer1);
    if(resultado.answer1 == 2)
        resultado.answer2 = input('Escolha o Metodo Previsor-Corretor: \n 1-Mehotra \n 2-Tanabi \n 3-Gondzio \n 4-Misto \n 5-Misto com Busca Unidimensional de Fibonacci para o Passo Previsor \n ');
        resultado.answer2 = VerificaOpcaoMetodo(resultado.answer2); 
    end
    result = resultado;
end

function result = VerificaOpcao(answer)
    while((answer ~= 1) && (answer ~= 2) && (answer ~= 3))
       answer = input('Favor Escolher corretamente o metodo pelo qual deseja Executar o Metodo de Resolucao: \n 1-Primal-Dual PPQ \n 2-Primal Dual com Previsor-Corretor \n 3-Primal-Dual PNL (Rosenbrock) \n');
    end
    result = answer;
end

function result = VerificaOpcaoMetodo(answer)
    while((answer ~= 1) && (answer ~= 2) && (answer ~= 3) && (answer ~= 4) && (answer ~= 5))
       answer = input('Favor Escolher corretamente o metodo pelo qual deseja Executar o Metodo Previsor-Corretor: \n 1-Mehotra \n 2-Tanabi \n 3-Gondzio \n 4-Misto \n 5-Misto com Busca Unidimensional de Fibonacci para o Passo Previsor \n');
    end
    result = answer;
end 
 
function result = EscolheMetodoResolucao(A, e, x, w, s, t, u, v, p, teta, X, S, resposta, k, Q, C, mi, C_quadratic, constant, minBuscaUni, maxBuscaUni, precisao)
    if((resposta.answer1 == 1) || (resposta.answer1 == 3))    
        resultado.Dw = Direcaow(A, teta, t, u, p); %v, X
        resultado.Dx = DirecaoX(teta, A, resultado.Dw, u, p); %v, X
        resultado.Ds = DirecaoS(X, p, S, resultado.Dx);       
    else
        if(resposta.answer1 == 2)
           correct = Corretor(A, e, x, w, s, t, u, v, teta, X, S, resposta, k, p, Q, C, mi, C_quadratic, constant, minBuscaUni, maxBuscaUni, precisao);
           resultado.Dw = correct.Dw;
           resultado.Dx = correct.Dx;
           resultado.Ds = correct.Ds;
           resultado.PrevCorre = correct;
        end
    end
    result = resultado;
end

function result = fibonacci(d,x, C, C_quadratic, constant, min, max, precisao, answer)
% o intuito dessa parte e descobrir o alpha, a e b sao pontos que entre
% eles exista um alpha, o objetivo dessa parte da busca unidimensional e
% encontrar o minimo da funcao f(xk + alpha*dxk), alimentando o dxk como
% ogradiente da funcao principal.
%   a=0;
%   b=0.9995;
%   a=-1;
%   b=1;
%  l=0.2;
  a=min;
  b=max;
  l=precisao;  
  k = 1;
  metodo.nome = 'Fibonacci';
  metodo.a(k) = a;
  metodo.b(k) = b;
  metodo.fibonacci = determinar_N(metodo.a(k),metodo.b(k),l);
  metodo.n = length(metodo.fibonacci);
  metodo.lambda(k) = lambd(metodo.a(k),metodo.b(k),metodo.fibonacci,metodo.n,k);
  metodo.mi(k) = calculo_mi(metodo.a(k),metodo.b(k),metodo.fibonacci,metodo.n,k);  
  metodo.f_lambda(k) = objetivo(metodo.lambda(k),d,x, C, C_quadratic, constant, answer);  
  metodo.f_mi(k) = objetivo(metodo.mi(k),d,x, C, C_quadratic, constant, answer); 
  while (k<metodo.n-1)      
      if (objetivo(metodo.lambda(k),d,x, C, C_quadratic, constant, answer) > objetivo(metodo.mi(k),d,x, C, C_quadratic, constant, answer))
        metodo.a(k+1) = metodo.lambda(k);
        metodo.b(k+1) = metodo.b(k);
        metodo.lambda(k+1) = metodo.mi(k);
        metodo.mi(k+1) = calculo_mi(metodo.a(k+1),metodo.b(k+1),metodo.fibonacci,metodo.n,k);        
        metodo.f_mi(k+1) = objetivo(metodo.mi(k+1),d,x, C, C_quadratic, constant, answer);        
        metodo.f_lambda(k+1) =  metodo.f_mi(k+1);
      else
          metodo.a(k+1) = metodo.a(k);
          metodo.b(k+1) = metodo.mi(k);
          metodo.mi(k+1) = metodo.lambda(k);
          metodo.lambda(k+1) = lambd(metodo.a(k+1),metodo.b(k+1),metodo.fibonacci,metodo.n,k);      
          metodo.f_lambda(k+1) = objetivo(metodo.lambda(k+1),d,x, C, C_quadratic, constant, answer);          
          metodo.f_mi(k+1) =  metodo.f_lambda(k+1);
      end
    k = k+1;
  end 
  metodo.lambda(k) = metodo.lambda(k-1);
  metodo.mi(k) = metodo.mi(k-1)+l;
  if (objetivo(metodo.lambda(k),d,x, C, C_quadratic, constant, answer) < objetivo(metodo.mi(k),d,x, C, C_quadratic, constant, answer))
      metodo.a(k+1) = metodo.a(k);
      metodo.b(k+1) = metodo.mi(k);
  else
      metodo.a(k+1) = metodo.lambda(k);
      metodo.b(k+1) = metodo.b(k);
  end
  metodo.numero_iteracoes = k;
  if objetivo(metodo.a(k),d,x, C, C_quadratic, constant, answer) < objetivo(metodo.b(k),d,x, C, C_quadratic, constant, answer)
      saida = metodo.a(k);
  else
      saida = metodo.b(k);
  end
  result = saida;
end

function result = objetivo(alpha, direcao, z, C, C_quadratic, constant, answer)
    x = z + alpha*direcao;
    result = calculoFO(x, C, C_quadratic, constant, answer);  
end

function result = determinar_N(a,b,l)
    fn = (b-a)/l;
    n=1;
    fibonacci = sequencia_fibonacci(10);
    while (fibonacci(n)<=fn)
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

function result = calculo_mi(a,b,fibonacci,n,k)
    mi = a+(fibonacci(n-k)/fibonacci(n-k+1))*(b-a);
    result = mi;
end