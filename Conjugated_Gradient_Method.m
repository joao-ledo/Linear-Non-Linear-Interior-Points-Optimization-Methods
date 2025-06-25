%_________________________________________________________________________%
%                        CONJUGATED GRADIENT METHOD                       %
%__*Developed by Joao Augusto Silva Ledo*_________________________________%

function result = GradientesConjugados()
    clear all;
    clc; 
    epslon = 0;
     k = 1;
    x{k} = [0;0]; 
    %x{k} = [8;9];
    %x{k} = [-1/2;1];
    g{k} = gradiente(x{k});
    d{k}= -g{k};
    
    while (max(abs(gradiente(x{k}))) >= epslon)
        alpha{k}=-((g{k}'*d{k})/(d{k}'*hessiana()*d{k}));
        x{k+1}=x{k}+alpha{k}*d{k};
        g{k+1}=gradiente(x{k+1});
        beta{k}=((g{k+1}'*hessiana()*d{k})/(d{k}'*hessiana()*d{k}));
        d{k+1}=-g{k+1}+beta{k}*d{k};
        k = k+1;    
    end

gc.nome = 'Conjugated Gradient Method';
gc.MelhorPonto = x{k-1}';


    result = gc;
end

 function result = objetivo(variave) 
    x=variavel(1);
    y=variavel(2);
    result = x^2 + y^2 + x*y - 3*x;
    %result = 4*(x-5)^2+(y-6)^2;
    %result = -12*y+4*x^2+4*y^2-4*x*y;
 end

function result = gradiente(variavel)
    x=variavel(1);
    y=variavel(2);
    result = [2*x+y-3;2*y+x] ;
    %result = [8*x-40;2*y-12] ;
    %result = [8*x-4*y;8*y-4*x-12];
end

function result = hessiana()
   result = [2, 1 ;1, 2];
   %result = [8,0;0,2];
   %result = [8,-4;-4,8];
end
