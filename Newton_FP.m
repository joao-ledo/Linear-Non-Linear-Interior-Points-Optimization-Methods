%_________________________________________________________________________%
%              NEWTON METHOD TO FIND THE OPTIMAL POWER FLOW               %
%__*Developed by Joao Augusto Silva Ledo*_________________________________%

function result = Newton_FP()
  clear all;
  clc;
   k=1;
   E = 1.0*10^(-10);
   v = [1.05, 1.0, 1.0];
 % v = [1.1, 1.0, 1.0];
   teta = [0.0, 0.0, 0.0];
   pg = [0.0, 0.0, 0.3];
 % pg = [0.0, 0.0, 0.5];
   qg = [0.0, 0.0, 0.0];
   pc = [0.0, 0.5, 0.0];
 % pc = [0.0, 0.5, 0.5];
   qc = [0.0, 0.2, 0.0];
 % qc = [0.0, 0.2, 0.5];
   r = [0.2, 0.2, 0.1];
   x = [0.4, 0.4, 0.2];
   bShunt = [0.01, 0.01, 0.01];
   for i = 1 : length(r)
       b(i) = calculaB(r(i),x(i));
       g(i) = calculaG(r(i),x(i));  
   end
   newton(k).entrada = [v(2);teta(2);teta(3)];
 %  jacobiana = jacobiana_objetivo(pg,pc,qg,qc,teta,v,b,bShunt,g, newton(k).entrada);

   while(abs(max(calculoF(v, teta, g, b, bShunt, newton(k).entrada, pg, qg, pc, qc)))> E)         
        newton(k).jacobiana =  matriz_jacobiana(v, teta, g, b, bShunt, newton(k).entrada);
        newton(k+1).entrada = proximo_candidato(newton(k).jacobiana, v, teta, g, b, bShunt, newton(k).entrada, pg, qg, pc, qc);   
        k = k + 1;
   end 
   
   objective =  abs(max(calculoF(v, teta, g, b, bShunt, newton(k).entrada, pg, qg, pc, qc)));
   v(2)=newton(k).entrada(1);
   teta(2)=newton(k).entrada(2);
   teta(3)=newton(k).entrada(3);
   fluxos = CalculoFluxos(v, teta, g, b, bShunt, newton(k).entrada);
   perd =  perdas(v, teta, g, b, bShunt, newton(k).entrada);
   balanco = calculoF(v, teta, g, b, bShunt, newton(k).entrada, pg, qg, pc, qc);

   fp.nome = 'Power Flow';
   fp.modulo_tensao = v;
   fp.angulos_tensao = teta;
   
   fp.fluxo_P12 = fluxos(7);
   fp.fluxo_Q12 = fluxos(9);
   fp.fluxo_P21 = fluxos(1);
   fp.fluxo_Q21 = fluxos(3);
   
   fp.fluxo_P23 = fluxos(2);
   fp.fluxo_Q23 = fluxos(4);
   fp.fluxo_P32 = fluxos(6);
   fp.fluxo_Q32 = fluxos(10);
   
   fp.fluxo_P13 = fluxos(8);
   fp.fluxo_Q13 = fluxos(11);
   fp.fluxo_P31 = fluxos(5);
   fp.fluxo_Q31 = fluxos(12);

   fp.perda_ativa_1_2 = perd(1);
   fp.perda_ativa_2_3 = perd(2);
   fp.perda_ativa_3_1 = perd(3);
   fp.perda_reativa_1_2 = perd(4);
   fp.perda_reativa_2_3 = perd(5);
   fp.perda_reativa_3_1 = perd(6);
   
   fp.balanco_Equacao1 = balanco(1);  
   fp.balanco_Equacao2 = balanco(2); 
   fp.balanco_Equacao3 = balanco(3); 
   
   fp.objetivo = objective;
   fp.numero_de_iteracoes = k;   
   
   result = fp;
end

function result = matriz_jacobiana(v, teta, g, b, bShunt, entrada)
    a11 = g(1)*v(1)*cos(teta(1) - entrada(2)) - 2*g(3)*entrada(1) - 2*g(1)*entrada(1) + g(3)*v(3)*cos(entrada(2) - entrada(3)) - b(1)*v(1)*sin(teta(1) - entrada(2)) + b(3)*v(3)*sin(entrada(2) - entrada(3));            
    a12 = b(1)*v(1)*entrada(1)*cos(teta(1) - entrada(2)) + b(3)*entrada(1)*v(3)*cos(entrada(2) - entrada(3)) + g(1)*v(1)*entrada(1)*sin(teta(1) - entrada(2)) - g(3)*entrada(1)*v(3)*sin(entrada(2) - entrada(3));             
    a13 = g(3)*entrada(1)*v(3)*sin(entrada(2) - entrada(3)) - b(3)*entrada(1)*v(3)*cos(entrada(2) - entrada(3));             

    a21 = 2*entrada(1)*(b(1)+bShunt(1)) + 2*entrada(1)*(b(3) + bShunt(3)) - b(1)*v(1)*cos(teta(1) - entrada(2)) - b(3)*v(3)*cos(entrada(2) - entrada(3)) - g(1)*v(1)*sin(teta(1) - entrada(2)) + g(3)*v(3)*sin(entrada(2) - entrada(3));            
    a22 = g(1)*v(1)*entrada(1)*cos(teta(1) - entrada(2)) + g(3)*entrada(1)*v(3)*cos(entrada(2) - entrada(3)) - b(1)*v(1)*entrada(1)*sin(teta(1) - entrada(2)) + b(3)*entrada(1)*v(3)*sin(entrada(2) - entrada(3));            
    a23 = - g(3)*entrada(1)*v(3)*cos(entrada(2) - entrada(3)) - b(3)*entrada(1)*v(3)*sin(entrada(2) - entrada(3));           

    a31 = g(3)*v(3)*cos(entrada(2) - entrada(3)) - b(3)*v(3)*sin(entrada(2) - entrada(3));          
    a32 = - b(3)*entrada(1)*v(3)*cos(entrada(2) - entrada(3)) - g(3)*entrada(1)*v(3)*sin(entrada(2) - entrada(3));       
    a33 = b(2)*v(1)*v(3)*cos(teta(1) - entrada(3)) + b(3)*entrada(1)*v(3)*cos(entrada(2) - entrada(3)) + g(2)*v(1)*v(3)*sin(entrada(1) - entrada(3)) + g(3)*entrada(1)*v(3)*sin(entrada(2) - entrada(3));  
    
    result = [a11,a12,a13; a21,a22,a23; a31,a32,a33];
end

function result = calculoF(v, teta, g, b, bShunt, entrada, pg, qg, pc, qc)
    fluxos = CalculoFluxos(v, teta, g, b, bShunt, entrada);
    F1 = pg(2)-pc(2)-(fluxos(1)+fluxos(2)); 
    F2 = qg(2)-qc(2)-(fluxos(3)+fluxos(4)); 
    F3 = pg(3)-pc(3)-(fluxos(5)+fluxos(6)); 
    result = [F1; F2; F3];
end

function result = CalculoFluxos(v, teta, g, b, bShunt, entrada)
    p21=entrada(1)^2*g(1)-entrada(1)*v(1)*g(1)*cos(entrada(2)-teta(1))-entrada(1)*v(1)*b(1)*sin(entrada(2)-teta(1)); 
    p23=entrada(1)^2*g(3)-entrada(1)*v(3)*g(3)*cos(entrada(2)-entrada(3))-entrada(1)*v(3)*b(3)*sin(entrada(2)-entrada(3)); 
    
    q21=-entrada(1)^2*(b(1)+bShunt(1))+entrada(1)*v(1)*b(1)*cos(entrada(2)-teta(1))-entrada(1)*v(1)*g(1)*sin(entrada(2)-teta(1));
    q12=-v(1)^2*(b(1)+bShunt(1))+entrada(1)*v(1)*b(1)*cos(teta(1)-entrada(2))-entrada(1)*v(1)*g(1)*sin(teta(1)-entrada(2));
    
    q23=-entrada(1)^2*(b(3)+bShunt(3))+entrada(1)*v(3)*b(3)*cos(entrada(2)-entrada(3))-entrada(1)*v(3)*g(3)*sin(entrada(2)-entrada(3)); 
    q32=-v(3)^2*(b(3)+bShunt(3))+entrada(1)*v(3)*b(3)*cos(entrada(3)-entrada(2))-entrada(1)*v(3)*g(3)*sin(entrada(3)-entrada(2));
    
    q13=-v(1)^2*(b(2)+bShunt(2))+v(1)*v(3)*b(2)*cos(teta(1)-entrada(3))-v(1)*v(3)*g(2)*sin(teta(1)-entrada(3));
    q31=-v(3)^2*(b(2)+bShunt(2))+v(1)*v(3)*b(2)*cos(entrada(3)-teta(1))-v(1)*v(3)*g(2)*sin(entrada(3)-teta(1));
    
    p31=v(3)^2*g(2)-v(3)*v(1)*g(2)*cos(entrada(3)-teta(1))-v(3)*v(1)*b(2)*sin(entrada(3)-teta(1)); 
    p32=v(3)^2*g(3)-v(3)*entrada(1)*g(3)*cos(entrada(3)-entrada(2))-v(3)*entrada(1)*b(3)*sin(entrada(3)-entrada(2)); 
    
    p12=v(1)^2*g(1)-entrada(1)*v(1)*g(1)*cos(teta(1)-entrada(2))-entrada(1)*v(1)*b(1)*sin(teta(1)-entrada(2));
    p13=v(1)^2*g(2)-v(3)*v(1)*g(2)*cos(teta(1)-entrada(3))-v(3)*v(1)*b(2)*sin(teta(1)-entrada(3));
    
    result =[p21,p23,q21,q23,p31,p32, p12, p13, q12, q32, q13, q31];
end

function result = perdas(v, teta, g, b, bShunt, entrada)
    fluxos = CalculoFluxos(v, teta, g, b, bShunt, entrada);
    
    perda_ativa_12=fluxos(1)+fluxos(7);
    perda_ativa_23=fluxos(2)+fluxos(6);
    perda_ativa_31=fluxos(5)+fluxos(8);
    
    perda_reativa_12=fluxos(3)+fluxos(9);
    perda_reativa_23=fluxos(4)+fluxos(10);
    perda_reativa_31=fluxos(11)+fluxos(12);
    
    result = abs([perda_ativa_12,perda_ativa_23,perda_ativa_31,perda_reativa_12,perda_reativa_23,perda_reativa_31]);
end

function result = proximo_candidato(jacobiana, v, teta, g, b, bShunt, entrada, pg, qg, pc, qc)
    result = entrada - ((jacobiana)^(-1) * calculoF(v, teta, g, b, bShunt, entrada, pg, qg, pc, qc));
end

function result = calculaG(r,x)
    result = r/(r^2 + x^2);
end

function result = calculaB(r,x)
    result = (-x)/(r^2 + x^2);
end

% Used function, during the code development, to find the Jacobian matrix
function result = jacobiana_objetivo(pg,pc,qg,qc,teta,v,b,bShunt,g, entrada)
    syms pg1 pg2 pg3 pc1 pc2 pc3 qg1 qg2 qg3 qc1 qc2 qc3 teta1 teta2 teta3 v1 v2 v3 b1 b2 b3 bshunt1 bshunt2 bshunt3 g1 g2 g3;
    f1 = pg2 - pc2 - (v2^2*g1-v2*v1*g1*cos(teta2-teta1)-v2*v1*b1*sin(teta2-teta1)+v2^2*g3-v2*v3*g3*cos(teta2-teta3)-v2*v3*b3*sin(teta2-teta3));
    f2 = qg2 - qc2 - (-v2^2*(b1+bshunt1)+v2*v1*b1*cos(teta2-teta1)-v2*v1*g1*sin(teta2-teta1) -  v2^2*(b3+bshunt3)+v2*v3*b3*cos(teta2-teta3)-v2*v3*g3*sin(teta2-teta3));
    f3 = pg3 - pc3 - (v3^2*g2-v3*v1*g2*cos(teta3-teta1)-v3*v1*b2*sin(teta3-teta1) + v3^2*g3-v3*v2*g3*cos(teta3-teta2)-v3*v2*b3*sin(teta3-teta2));
    f = [f1;f2;f3];
    jacobiana = jacobian(f,[v2,teta2,teta3]);
  % pg1 = pg(1);
  % pg2 = pg(2);
  % pg3 = pg(3);
  % pc1 = pc(1);
  % pc2 = pc(2);
  % pc3 = pc(3);
  % qg1 = qg(1);
  % qg2 = qg(2);
  % qg3 = qg(3);
  % qc1 = qc(1);
  % qc2 = qc(2);
  % qc3 = qc (3);
  % teta1 = teta(1);
  % teta2 = entrada(2);
  % teta3 = entrada(3);
  % v1 = v(1);
  % %v2 = entrada(1);
  % v3 = v(3),
  % b1 = b(1);
  % b2 = b(2);
  % b3 = b(3);
  % bshunt1 = bShunt(1);
  % bshunt2 = bShunt(2);
  % bshunt3 = bShunt(3);
  % g1 = g(1);
  % g2 = g(2);
  % g3 = g(3);
  % jac = subs(jacobiana);
  % syms y;
  % f = y^2+2*y;
  % f_linha = diff(f,y);
  % f_2linhas = diff(f_linha,y);
  % y = x;
  % result = [subs(f_linha), subs(f_2linhas)];
    result = jacobiana;
end
