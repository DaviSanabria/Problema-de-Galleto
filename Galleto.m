clear all; clc;
global lisGenes;
global lisCalidades;
global lisDiversidades;
global pertdiversidades;
global genes;
global matrizC;
matrizC = zeros(9); % crea una matrizC de 9x9 inicialmente llena de ceros
%Crear la matrizC
for i = 1:9
    for j = 1:9
        if(i<6 && j<6)
            piramide(i,j)= i+j-1;
        elseif(i>4 && j>4)
            piramide(i,j)= 9-i+9-j+1;
        elseif(i>5 && j<5)
            piramide(i,j)= j+9-i;
        elseif(i<5 && j>5)
            piramide(i,j)= 9-j+i;
        end
    end
end

%Convertirla en foso
foso = piramide;
for i = 1:9
    for j = 1:9
        if(i>1 && j>1 && i<9 && j<9 && (i>6 || i<4 || j>6 || j<4))
            foso(i,j)= 0;
        end
    end
end

clear i, clear j;

matrizC = piramide;
genes=zeros(4,2); %matrizC de genes
genes(1,1)=1;
genes(1,2)=1;   
primero = [genes(1,1), genes(1,2)]; %se define la madre para la mutación inicial
calidades= zeros(1,4); %calidades de los genes
diversidades = zeros(1, 4); %diversidad de los genes
pertdiversidades = zeros(1,4); %genes seleccionados para evaluar las diversidades
pertdiversidades(1) = true;
calidades(1) = calidad(primero); %genes de calidades



numero = 100;%cantidad de iteraciones


iteraciones = zeros(1,numero);
% ciclo de mutuación
for i= 1:numero
    iteraciones(i) = mutacion();
end

mediaMutacion= mean(iteraciones)
modaMutacion= mode(iteraciones)

% ciclo de mutuación y cruce para la piramide
for i= 1:numero
    iteraciones(i) = mutacionCruce();
end

mediaMutacionCrucePiramide= mean(iteraciones)
modaMutacionCrucePiramide= mode(iteraciones)

matrizC= foso;
% ciclo de mutuación y cruce para el foso
for i= 1:numero
    iteraciones(i) = mutacionCruce();
end

mediaMutacionCruceFoso= mean(iteraciones)
modaMutacionCruceFoso= mode(iteraciones)

%función mutación 
function [result] = mutacion()
global genes;
genes(1,1)=1;
genes(1,2)=1; 
global matrizC;
% ciclo de mutuación
i=true;
generaciones = 1;
while i
    %# region de mutacion
        x = randi([1,2]); %aleatorio para la mutación
        if x>1
            c1=genes(1,1);
            c2=genes(1,2)+1;
            if c2 == 10
               c2 = 1;
            end
            hijo = [c1, c2];
        else
            c1=genes(1,1)+1;
            c2=genes(1,2);
            if c1 == 10
               c1 = 1;
            end
            hijo = [c1, c2];
        end    
    %# end
    
    %# region selección de genes
        seleccion(hijo);
    %# end
    
    %# región de evaluación de calidad buscada
    primero = [genes(1,1), genes(1,2)];
    calPrimero = calidad(primero);
    if calPrimero == 9
        i = false;
    else
        generaciones = generaciones + 1; 
    end
    %# end
    end
    result = generaciones;
end

%función mutación y cruce
function [result] = mutacionCruce()
global genes;
genes(1,1)=1;
genes(1,2)=1;
genes(2,1)=1;
genes(2,2)=1;
global matrizC;
% ciclo de mutuación
i=true;
generaciones = 1;
while i
    %# region de mutacion del primero
        x = randi([1,2]); %aleatorio para la mutación
        if x>1
            c1=genes(1,1);
            c2=genes(1,2)+1;
            if c2 == 10
               c2 = 1;
            end
            hijo1 = [c1, c2];
        else
            c1=genes(1,1)+1;
            c2=genes(1,2);
            if c1 == 10
               c1 = 1;
            end
            hijo1 = [c1, c2];
        end    
    %# end
    %# region de mutacion del segundo
       if genes(2,1) ~= 0 && genes(2,2) ~= 0
        x = randi([1,2]); %aleatorio para la mutación
        if x>1
            c1=genes(2,1);
            c2=genes(2,2)+1;
            if c2 == 10
               c2 = 1;
            end
            hijo2 = [c1, c2];
        else
            c1=genes(2,1)+1;
            c2=genes(2,2);
            if c1 == 10
               c1 = 1;
            end
            hijo2 = [c1, c2];
        end
       end
    %# end

    %# region selección de genes
    if genes(2,1) ~= 0 && genes(2,2) ~= 0
        %cruce
        temp= zeros(1,2);
        temp=hijo1;
        hijo1(2) = hijo2(2);
        hijo2(2) = temp(2);
        seleccionCruce(hijo1,hijo2);
    else
        seleccion(hijo1);
    end
    %# end
    
    %# región de evaluación de calidad buscada
    primero = [genes(1,1), genes(1,2)];
    segundo = [genes(2,1), genes(2,2)];
    calPrimero = calidad(primero);
    if segundo(1) ~=0 && segundo(2) ~= 0
        calSegundo = calidad(segundo);
        if calSegundo == 9
            i = false;
        end
    end
    if calPrimero == 9
        i = false;
    else
        generaciones = generaciones + 1; 
    end
    %# end
    end
    result = generaciones;
end




%función de calidad
function [result] = calidad(gen) 
    global matrizC;
    result = matrizC(gen(1,1),gen(1,2));
end

%función distancia cuadrática
function [result] = distanciaCua(gen1, gen2)
    result = ((gen1(1)-gen2(1))^2)+((gen1(2)-gen2(2))^2);
end

%función de diversidad
function [result] = diversidad(genes, perdiversidades)
    genesBase = zeros(1,4);
    gen1 =zeros (1,2);
    gen2 =zeros (1,2);
    lisDiversidades = zeros(1,length(perdiversidades));
    pos = 1;
    for i=1: length(perdiversidades)
        if perdiversidades(i) == true
           genesBase(pos)=i;
           pos =+ pos;
        end
    end
    genesBase = genesBase(1:pos);
    for i=1 : length(genes)
        if perdiversidades(i) ~= true && genes(i,1) ~= 0
            div=0;
            for v=1 : length(pos)
                geB = genesBase(v);
                gen1(1) = genes(geB,1);
                gen1(2) = genes(geB,2);
                gen2(1) = genes(i,1);
                gen2(2) = genes(i,2);
                dist = distanciaCua(gen1, gen2);
                div = div+(1/dist);     
            end
            lisDiversidades(i) = div;
        end
    end
    result = lisDiversidades;
end

%función rango calidades
function [result] = rangCal(lisCalidades)
    rangoCal = zeros(1,length(lisCalidades));
    for i = 1 : length(lisCalidades)
        mayor = 0;
        mayorAnt = 0;
        temp = 1;
        for j = 1 : length(lisCalidades)
            if i == 1
                if (lisCalidades(temp) <= lisCalidades(j)) && (lisCalidades(j) ~= 0)
                    mayor = j;
                    temp = j;
                end
            else
                if mayorAnt ~= 0
                    if lisCalidades(mayorAnt)-1 <= lisCalidades(j) && lisCalidades(j) < lisCalidades(mayorAnt) && lisCalidades(j) ~= 0
                        mayor = j;
                    end
                end
            end
        end
        if mayor ~= 0
            rangoCal(mayor) = i;
            mayorAnt = mayor;
        end
    end
    result = rangoCal;
end

%función calidades
function [result]= Calidades(lisGenes)
    global matrizC;
    lisCalid = zeros(1,length(lisGenes));
    for i = 1 :length(lisGenes)
        genX = zeros(1,2);
        genX(1) = lisGenes(i,1);
        genX(2) = lisGenes(i,2);
        if(genX(1)~=0)
            lisCalid(i) = calidad(genX);
        end
    end
    result = lisCalid;
end

%función rango diversidades
function [result] = rangDiv(lisDiversidades)
    rangoDiv = zeros(1,length(lisDiversidades));
    for i = 1 : length(lisDiversidades)
        menor = 0;
        for j = 1 : length(lisDiversidades)
            if i == 1
                maximo = max(lisDiversidades);
                if (maximo >= lisDiversidades(j)) && (lisDiversidades(j) ~= 0)
                    menor = j;
                end
            else
                if lisDiversidades(menorAnt)+1 >= lisDiversidades(j) && lisDiversidades(j) > lisDiversidades(menorAnt) && lisDiversidades(j) ~= 0
                    menor = j;
                end
            end
        end
        if menor ~= 0
            rangoDiv(menor) = i;
            menorAnt = menor;
        end
    end
    result = rangoDiv;
end

%función evaluacion
function [result] = evaluacion()
    global lisGenes
    global lisCalidades;
    global lisDiversidades;
    rangoCalidades = zeros(1,length(lisCalidades));
    rangoDiversidades = zeros(1,length(lisDiversidades));
    %rango calidades
    rangoCalidades = rangCal(lisCalidades);
    %rango diversidad
    rangoDiversidades = rangDiv(lisDiversidades);

    %combinar
    rangComb = zeros(1,length(lisCalidades));
    for i=1 : length(lisCalidades)
        if rangoCalidades(i) ~= 0 && rangoDiversidades(i) ~=0
            rangComb(i) =  (rangoCalidades(i)+rangoDiversidades(i))/2;
        elseif rangoCalidades(i) ~= 0
            rangComb(i) = rangoCalidades(i);
        else 
            rangComb(i) = rangoDiversidades(i);
        end
    end
    result = rangComb;
end

%función ordenar
function ordenar(rangoEvaluacion)
    global pertdiversidades;
    global lisGenes;
    tempperdiversidades = zeros(1, length(rangoEvaluacion));
    ListaFinal = zeros(length(rangoEvaluacion),2);
    for i = 1 : length(rangoEvaluacion)
        menor = 0;
        temp = 1;
        for j = 1 : length(rangoEvaluacion)
            if i == 1
                maximo = max(rangoEvaluacion);
                if (maximo >= rangoEvaluacion(j)) && (rangoEvaluacion(j) ~= 0)
                    menor = j;
                    temp = j;
                end
            else
                if rangoEvaluacion(menorAnt)+1 >= rangoEvaluacion(j) && rangoEvaluacion(j) > rangoEvaluacion(menorAnt) && rangoEvaluacion(j) ~= 0
                    menor = j;
                    temp = j;
                end
            end
        end
        if menor ~= 0
            ListaFinal(i,1) = lisGenes(menor,1);
            ListaFinal(i,2) = lisGenes(menor,2);
            if menor == length(rangoEvaluacion) || menor == length(rangoEvaluacion)-1;
                tempperdiversidades(i) = true;
            else
                tempperdiversidades(i) = pertdiversidades(menor);
            end
            menorAnt = menor;
        end
    end
    lisGenes = ListaFinal;
    pertdiversidades = tempperdiversidades(1:length(rangoEvaluacion)-1);
end

function [result] = seleccion(hijo)
    global pertdiversidades;
    global lisGenes;
    global lisCalidades;
    global lisDiversidades;
    global genes;
    global matrizC;
    calHijo=calidad(hijo);
    lisGenes = genes;
    cero = [0,0];
    lisGenes = cat(1,lisGenes, cero);
    lisGenes(5,1)=hijo(1);
    lisGenes(5,2)=hijo(2);
    lisCalidades = Calidades(lisGenes);
    lisCalidades(5) = calidad(hijo);
    lisPerdiversidades = pertdiversidades;
    lisPerdiversidades = cat(2,lisPerdiversidades, 0);
    lisDiversidades = zeros(1,5);
    lisDiversidades = diversidad(lisGenes, lisPerdiversidades);
    rangoEvaluacion = evaluacion();
    ordenar(rangoEvaluacion);
    lisGenes(end, :) = [];
    genes = lisGenes;
end

function [result] = seleccionCruce(hijo1, hijo2)
    global pertdiversidades;
    global lisGenes;
    global lisCalidades;
    global lisDiversidades;
    global genes;
    global matrizC;
    lisGenes = genes;
    cero = [0,0];
    lisGenes = cat(1,lisGenes, cero);
    lisGenes = cat(1,lisGenes, cero);
    lisGenes(5,1)=hijo1(1);
    lisGenes(5,2)=hijo1(2);
    lisGenes(6,1)=hijo2(1);
    lisGenes(6,2)=hijo2(2);
    lisCalidades = Calidades(lisGenes);
    lisCalidades(5) = calidad(hijo1);
    if hijo2(1) ~=0 && hijo2(2) ~= 0
        lisCalidades(6) = calidad(hijo2);
    end
    lisPerdiversidades = pertdiversidades;
    lisPerdiversidades = cat(2,lisPerdiversidades, 0);
    lisPerdiversidades = cat(2,lisPerdiversidades, 0);
    lisDiversidades = zeros(1,6);
    lisDiversidades = diversidad(lisGenes, lisPerdiversidades);
    rangoEvaluacion = evaluacion();
    ordenar(rangoEvaluacion);
    lisGenes(end, :) = [];
    lisGenes(end, :) = [];
    genes = lisGenes;
end