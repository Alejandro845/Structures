clear
clc
clf

%Ingresar el tipo de estructura tipoestru = [1 si la estructura es un
%portico o viga, 0 si es una cercha]
tipoestru = 1;

%SI SU ESTRUCTURA ES UNA CERCHA, ejecute la combinación crtl+f,y digite la 
%palabra CRCH, luego presione ENTER para ir a la sección correspondiente.

  if tipoestru==1

   %************Ingreso de datos (PORTICO O VIGA)***************
   
        %Las unidades de los valores de entrada deben ser KN y m obligatoriamente.
   
        %Ingrese las unidades de los resultados = [1 si son KN y m, 2 si son KN y mm
        %3 si son Kg y m, 4 si son Kg y mm, 5 si son Ton y m, 6 si son Ton y mm
        %7 si son N y m, 8 si son N y mm]
        UnidadesSalida = 8;
        
        
        %ingreso de coordenadas de los nudos
        % nudos = [#nudo, coordX, coodY;
        %     ...se repite para todos los nudos
        % ]
       

        nudos = [
        1 0 0;
        2 0 2;
        3 2 0;
        4 2 3;
        5 4 0;
        6 4 3;
        7 5 3
        ]

        
        %ingreso de datos de conexiones entre nudos
        % elementos = [#elemento, deNudo, aNudo;
        %           #elemento, deNudo, aNudo;
        %          ...se repite para todas las elementos
        % ]
        
        elementos = [
        1 1 2;
        2 2 4;
        3 3 4;
        4 4 6;
        5 5 6;
        6 6 7
        ]

        
        %Ingreso de tipo de estructura [0 para cercha, 1 para viga o
        %portico]
        
        %Ingreso de propiedades de sección de cada elemento
        % propSeccion = [#elemento, Area, Inercia, E;
        %                ...Se repite para todas las elementos
        %                ]
        A100 = ((100*100)-((100-2*3)*(100-2*3)))*1/1000^2;
        I100 = ((((100)*(100)^3)/12)-(((100-2*3)*(100-2*3)^3)/12))*1/1000^4;
        
        A90= ((90*90)-((90-2*2.5)*(90-2*2.5)))*1/1000^2;
        I90 = ((((90)*(90)^3)/12)-(((90-2*2.5)*(90-2*2.5)^3)/12))*1/1000^4;
        
        E = 2e8;
        
        propSeccion = [
        1 A100 I100 E;
        2 A90 I90 E;
        3 A100 I100 E;
        4 A90 I90 E;
        5 A100 I100 E;
        6 A90 I90 E
        ]
                      
        % Ingresar cargas puntuales
        % el nombre del vector se llamará puntuales
        % puntuales = [#nudo Fx Fy M...
        %               ....
        %               #nudo Fx Fy M]              
        
        puntuales = [
        7 0 -10 0;
        2 0 0 5
        ]

        
        % Ingresar cargas distribuidas
        
        % distribuidas = [ #noelemento fpara q1 q2
        %                 ...
        %                  #noelemento fpara q1 q2]
        % ** Nota q1 y q2 son positivas en la dirección local hacia abajo
        
        distribuidas = [
        2 -8 16 16;
        4 0 20 20;
        5 0 -20 -5
        ]
        
        % Ingresar la configuración de apoyos
        % apoyos = [ #nudoRestringido Dx Dy Giro
        %            ...
        %            #nudoRestringido Dx Dy Giro]
        % Dx, Dy, giro serán cero si el GDL está libre
        % y será 1 si el grado de libertad está restringido
        
        apoyos = [
        1 1 1 1;
        3 1 1 0;
        5 1 1 1
        ]
        
        % Escala de ploteado de deformaciones (1 para mostrar 0 para
        % ocultar)
        grafDef = 0;
        escDef = 15;
        
        grafM = 1;
        escM = 0.05;
        
        grafV = 0;
        escV = 0.05;
        
        grafN = 0;
        escN = 0.008;               
                       
        %**************************************************
        %  **********   INICIO DEL PROGRAMA    ********* 
        %**************************************************
        
        %Calculo de longitudes y ángulos de cada elemento
        temp = size(nudos); %variable temporal que guarda el tamaño de los nudos
        noNudos = temp(1) %extrae la casilla 1 de la variable temporal, es decir la cantidad de nudos
        temp = size(elementos); 
        noelementos = temp(1) %extrae la casilla 1 de la variable temporal reescrita para los elementos, es decir la cantidad de elementos
        
        for i=1:noelementos
           nudoInicio = elementos(i,2); %extrae la coordenada de inicio del nudo (fila2)
           nudoFin = elementos(i,3); %extrae la coordenada final del nudo (fila3)
           xi = nudos(nudoInicio,2); % asigna al nudo la coordenada  x de inicio respectiva
           yi = nudos(nudoInicio,3); % asigna al nudo la coordenada en y inicial
           xf = nudos(nudoFin,2); % asigna al nudo la coordenada de x final respectiva
           yf = nudos(nudoFin,3); % asigna al nudo la coordenada de  y final respectiva
           longelementos(i) = sqrt((xf-xi)^2+(yf-yi)^2); %calcula la distancia con la raiz de la diferencia de los cuadrados de las coordenadas y las almacena en un vector.
           anguloselementos(i) = acosd((xf-xi)/longelementos(i)); %calcula el angulo con arcocoseno, siempre medido desde la horizontal positiva
           if (yf-yi)<0
             anguloselementos(i) = -anguloselementos(i);
           end
        end
        longelementos
        anguloselementos   %como los valores arrojados por el arcoseno no pueden superar los 180 grados (rango de 0 a pi), puede arrojar un ángulo medido en el sentido negativo, pero sin el signo, por esa razón cuando la resta entre la coordenadas sea negativa, multiplico el valor del angulo por -1

        %Cálculo de matrices de rigidez locales de cada elemento
        
        for i=1:noelementos 
          elas = propSeccion(i,4); %extrae el módulo de elasticidad
          area = propSeccion(i,2);  %extrae el area
          inercia = propSeccion(i,3); %extrae la inercia
          long = longelementos(i); %extrae la longitud de la barra
          EAL = elas*area/long; %Calcula los valores que iran dentro de la matriz local
          EI12 = 12*elas*inercia/long^3;
          EI6 = 6*elas*inercia/long^2;
          EI4 = 4*elas*inercia/long;
          EI2 = 2*elas*inercia/long;
              
          kLocal(:,:,i) = [EAL 0 0 -EAL 0 0;
                    0 EI12 EI6 0 -EI12 EI6;
                    0 EI6 EI4 0 -EI6 EI2;
                    -EAL 0 0 EAL 0 0;
                    0 -EI12 -EI6 0 EI12 -EI6;
                    0 EI6 EI2 0 -EI6 EI4];    
        end   %ensambla la matriz de rigidez local para cada elemento, 
        % guarda todas las filas y todas las columnas en la “hoja i” para 
        % almacenar las matrices por separado y que la variable no se reescriba
        
        for i=1:noelementos 
          ang = anguloselementos(i); %extrae los angulos de los elementos
          cs = cosd(ang); %calculo de coseno y seno de cada angulo, con el d para que sea en grados
          sn = sind(ang);
          T = [cs sn 0 0 0 0; %montaje de la matriz de rotación
               -sn cs 0 0 0 0;
               0 0 1 0 0 0;
               0 0 0 cs sn 0;
               0 0 0 -sn cs 0;
               0 0 0 0 0 1];
          kGlobal(:,:,i) = T'*kLocal(:,:,i)*T;  %calculo de la matriz global para cada uno de los elementos, ya rotados
        
        end
        
        %ensamblaje de la matriz de rigidez de toda la estructura

        for i=1:noelementos
        GDLelementos(i,1) = elementos(i,1);
        GDLelementos(i,2) = elementos(i,2)*3-2;  
        GDLelementos(i,3) = elementos(i,2)*3-1;  
        GDLelementos(i,4) = elementos(i,2)*3;  
        GDLelementos(i,5) = elementos(i,3)*3-2;  
        GDLelementos(i,6) = elementos(i,3)*3-1;  
        GDLelementos(i,7) = elementos(i,3)*3;  
        end  %calcula los grados de libertad de cada elemento en sus extremos, partiendo de que noDelElemento*3=grado de libertad de giro, (noDelElemento*3)-2=grado de libertad horizontal, (noDelElemento*3)-1=grado de libertad vertical, ej. elemento 4. 4*3=12,  (4*3)-2=10, (4*3)-1=11, y así para todos los elementos.
        
        kEstructura = zeros(noNudos*3,noNudos*3); %Dimension de la matriz de rigidez de toda la estructura
        for i=1:noelementos
           gdl = GDLelementos(i,2:7); %Extrae los grados de libertad del elemento i
           kEstructura(gdl,gdl) = kGlobal(:,:,i)+kEstructura(gdl,gdl);
        end
        
        % Creacion de vector de cargas puntuales Q
        Q = zeros(noNudos*3,1); %Dimensiones del vector de cargas puntuales
        
        temp = size(puntuales); %Variable temportal que guarda el tamaño del vector
        noPuntuales = temp(1); %Extrae la cantidad de filas que hay en el vector
        for i=1:noPuntuales
          nudo = puntuales(i,1); %extrae el nudo en donde esta la carga
          Fx = puntuales(i,2); %extrae la fuerza en x que haya en ese nudo
          Fy = puntuales(i,3); %extrae la fuerza en y que haya en ese nudo
          M = puntuales(i,4); %extrae el momento que haya en ese nudo
          Q(nudo*3-2,1) = Fx + Q(nudo*3-2,1); %agrega la fuerza en x al vector de cargas puntuales
          Q(nudo*3-1,1) = Fy + Q(nudo*3-1,1); %agrega la fuerza en y al vector de cargas puntuales
          Q(nudo*3,1) = M + Q(nudo*3,1); %agrega el momento al vector de cargas puntuales
        end
        disp("Vector de cargas puntuales:");
        Q
        
        % Calculo  de vector de fuerzas distribuidas equivalentes
        f = zeros(noNudos*3,1); %Dimension del vector de fuerzas distribuidas
        temp = size(distribuidas); %extrae el tamaño del vector de distribuidas de entrada
        noDistribuidas = temp(1); %extrae la cantidad de filas de ese vector
        for i=1:noelementos
        fLocal(:,:,i) = [0;0;0;0;0;0]; %crea un vector de fuerzas locales para cada elemento
        
        end
        
        
        
        for i=1:noDistribuidas
          noelemento = distribuidas(i,1); %extrae el elemento
          fpara = distribuidas(i,2); %extrae la fuerza paralela
          q1 = distribuidas(i,3); %extrae la carga q1
          q2 = distribuidas(i,4); %extrae la carga q2
          nudoInicio = elementos(noelemento,2); %extrae el nudo incial del elemento
          nudoFin = elementos(noelemento,3); %extrae el nudo final del elemento
          L = longelementos(noelemento); %extrae la longitud del elemento
          
          f1 = fpara*L/2; %calcula la fuerza axial de los extremos
          f4 = fpara*L/2;
          
          if (q1<q2) %calcula las fuerzas horizontales y momentos de empotramiento perfecto en cada extremo
            f2 = -3/20*(q2-q1)*L-q1*L/2; 
            f3 = -(q2-q1)/30*L^2-q1*L^2/12;
            f5 = -7/20*(q2-q1)*L-q1*L/2;
            f6 = (q2-q1)/20*L^2+q1*L^2/12;
          elseif (q1==q2)
            f2 = -q1*L/2;
            f3 = -q1*L^2/12;
            f5 = f2;
            f6 = -f3;
          elseif (q1>q2)
            f2 = -7/20*(q1-q2)*L-q2*L/2;
            f3 = -(q1-q2)/20*L^2-q2*L^2/12;
            f5 = -3/20*(q1-q2)*L-q2*L/2;
            f6 = (q1-q2)/30*L^2+q2*L^2/12;
          end
          %armamos vector de fuerzas equivalentes locales
          fLocal(:,:,noelemento) = [f1 f2 f3 f4 f5 f6]';
        
          ang = anguloselementos(noelemento); %extrae el angulo del elemento
          cs = cosd(ang);
          sn = sind(ang);
          T = [cs sn 0 0 0 0;
               -sn cs 0 0 0 0;
               0 0 1 0 0 0;
               0 0 0 cs sn 0;
               0 0 0 -sn cs 0;
               0 0 0 0 0 1];
           
           fGlobal = T'*fLocal(:,:,noelemento); %calculo de las fuerzas globales para cada elemento
        
           f(nudoInicio*3-2,1) = fGlobal(1)+f(nudoInicio*3-2,1); %añade la fuerza global al vector de cargas distribuidas
           f(nudoInicio*3-1,1) = fGlobal(2)+f(nudoInicio*3-1,1);
           f(nudoInicio*3,1) = fGlobal(3)+f(nudoInicio*3,1);
        
           f(nudoFin*3-2,1) = fGlobal(4)+f(nudoFin*3-2,1);
           f(nudoFin*3-1,1) = fGlobal(5)+f(nudoFin*3-1,1);
           f(nudoFin*3,1) = fGlobal(6)+f(nudoFin*3,1);
           
        end
        
        disp("Vector de cargas equivalentes distribuidas f:")
        f
        
        %%Calculo de grados de libertad restringidos
        
        temp = size(apoyos); %extrae el tamaño de la matriz de apoyos
        noApoyos = temp(1); %extrae la cantidad de filas que hay en la matriz de apoyos
        contador = 0; %variable que contendra la posicion de una restriccion en un grado de libertad
        
        for i=1:noApoyos
          nudo = apoyos(i,1); %extrae el nodo donde hay alguna restriccion
          Dx = apoyos(i,2);
          Dy = apoyos(i,3);
          giro = apoyos(i,4);
        
          gdlx = nudo*3-2; %extrae los grados de libertad correspondientes a ese nodo
          gdly = nudo*3-1;
          gdlGiro = nudo*3;
          
          if (Dx ==1)
            contador = contador +1;
            gdlRestringidos(contador) = gdlx; %agrega el grado de libertad al vector de restricciones
          end
          if (Dy ==1)
            contador = contador +1;
            gdlRestringidos(contador) = gdly; %agrega el grado de libertad al vector de restricciones
          end
          if (giro ==1)
            contador = contador +1;
            gdlRestringidos(contador) = gdlGiro; %agrega el grado de libertad al vector de restricciones  
          end
        end
        
        %SOLUCION DEL SISTEMA DE ECUACIONES
        %creacion de matriz K reducida:
        kRed = kEstructura;
        kRed(:,gdlRestringidos) = []; %elimina las columnas en donde esten los grados de libertad restringidos
        kRed(gdlRestringidos,:) = []; %elimina las filas en donde esten los grados de libertad restringidos
        
        %creacion de vector de fuerzas reducidas
        fRed = f+Q;
        fRed(gdlRestringidos,:) = []; %elimina todas las filas donde esten los grados de libertad restringidos
        
        %resolviendo el sistema:
        uRed = linsolve(kRed,fRed); %linsolve utiliza gauss jordan para resolver sistemas de ecuaciones
        
        %%Creación de vector de grados de libertad libres
        gdlLibres = 1:1:noNudos*3;
        gdlLibres(gdlRestringidos) = []; %vector con los grados de libertad sin restricciones
          
        % Creacion de vector desplazamientos U con todos los ceros
        U = zeros(noNudos*3,1);
        
        temp = size(gdlLibres); %extrae el tamaño del vector de grados de libertad libres
        noGDLlibres = temp(2); %extrae la cantidad de columnas del vector de grados de libertad libres
        for i=1:noGDLlibres
          U(gdlLibres(i)) = uRed(i); %agrega el valor de desplazamientos en el grado de libertad correspondiente
        end
        
        %Resolución de reacciones del sistema:
        R = kEstructura*U-f-Q

        U
        
        %Conversiones de Unidades
        if (UnidadesSalida==1) %KN m - KN m
            
            disp("Desplazamientos de la estructura")
            U
            
            disp("Reacciones de la estructura")
            R
            
        elseif (UnidadesSalida==2) %KN m - KN mm
            
            for i=1:noNudos*3
                if (mod(i, 3) == 0) %si el grado de libertad es multiplo de 3, se trata de un giro
                    U1(i) = U(i); %guarda en una nueva variable para que la graficacion no se vea alterada
                else
                    U1(i) = U(i)*1000;
                end
            end  
            
            for i=1:noNudos*3
                if (mod(i, 3) == 0) %si el grado de libertad es multiplo de 3, se trata de un momento
                    R1(i) = R(i)*1000; %guarda en una nueva variable para que la graficacion no se vea alterada
                else
                    R1(i) = R(i);
                end
            end  
            
            
            disp("Desplazamientos de la estructura")
            U1 = U1' %imprime la transpuesta para que se vea como un vector columna
            
            disp("Reacciones de la estructura")
            R1 = R1' %imprime la transpuesta para que se vea como un vector columna
            
        elseif (UnidadesSalida==3) %KN m - Kg m
            
            for i=1:noNudos*3
                
               U1(i) = U(i); %guarda en una nueva variable para que la graficacion no se vea alterada
                
            end  
            
            for i=1:noNudos*3
                
               R1(i) = R(i)*101.97; %guarda en una nueva variable para que la graficacion no se vea alterada
   
            end  
            
            
            disp("Desplazamientos de la estructura")
            U1 = U1' %imprime la transpuesta para que se vea como un vector columna
            
            disp("Reacciones de la estructura")
            R1 = R1' %imprime la transpuesta para que se vea como un vector columna
            
        elseif (UnidadesSalida==4) %KN m - Kg mm
            
            for i=1:noNudos*3
                if (mod(i, 3) == 0) %si el grado de libertad es multiplo de 3, se trata de un giro
                    U1(i) = U(i); %guarda en una nueva variable para que la graficacion no se vea alterada
                else
                    U1(i) = U(i)*1000;
                end
            end  
            
            for i=1:noNudos*3
                if (mod(i, 3) == 0) %si el grado de libertad es multiplo de 3, se trata de un momento
                    R1(i) = R(i)*101.97*1000; %guarda en una nueva variable para que la graficacion no se vea alterada
                else
                    R1(i) = R(i)*101.97;
                end
            end  
            
            
            disp("Desplazamientos de la estructura")
            U1 = U1' %imprime la transpuesta para que se vea como un vector columna
            
            disp("Reacciones de la estructura")
            R1 = R1' %imprime la transpuesta para que se vea como un vector columna
            
        elseif (UnidadesSalida==5) %KN m - Ton m
            
            for i=1:noNudos*3
               
               U1(i) = U(i); %guarda en una nueva variable para que la graficacion no se vea alterada
               
            end  
            
            for i=1:noNudos*3
            
               R1(i) = R(i)*0.1; %guarda en una nueva variable para que la graficacion no se vea alterada
               
            end  
            
            
            disp("Desplazamientos de la estructura")
            U1 = U1' %imprime la transpuesta para que se vea como un vector columna
            
            disp("Reacciones de la estructura")
            R1 = R1' %imprime la transpuesta para que se vea como un vector columna
            
        elseif (UnidadesSalida==6) %KN m - Ton mm
            
            for i=1:noNudos*3
                if (mod(i, 3) == 0) %si el grado de libertad es multiplo de 3, se trata de un giro
                    U1(i) = U(i); %guarda en una nueva variable para que la graficacion no se vea alterada
                else
                    U1(i) = U(i)*1000;
                end
            end  
            
            for i=1:noNudos*3
                if (mod(i, 3) == 0) %si el grado de libertad es multiplo de 3, se trata de un momento
                    R1(i) = R(i)*0.1*1000; %guarda en una nueva variable para que la graficacion no se vea alterada
                else
                    R1(i) = R(i)*0.1;
                end
            end  
            
            
            disp("Desplazamientos de la estructura")
            U1 = U1' %imprime la transpuesta para que se vea como un vector columna
            
            disp("Reacciones de la estructura")
            R1 = R1' %imprime la transpuesta para que se vea como un vector columna
            
        elseif (UnidadesSalida==7) %KN m - N m
            
            for i=1:noNudos*3
            
               U1(i) = U(i); %guarda en una nueva variable para que la graficacion no se vea alterada
              
            end  
            
            for i=1:noNudos*3
            
               R1(i) = R(i)*1000; %guarda en una nueva variable para que la graficacion no se vea alterada
              
            end  
            
            
            disp("Desplazamientos de la estructura")
            U1 = U1' %imprime la transpuesta para que se vea como un vector columna
            
            disp("Reacciones de la estructura")
            R1 = R1' %imprime la transpuesta para que se vea como un vector columna
            
        elseif (UnidadesSalida==8) %KN m - N mm
            
            for i=1:noNudos*3
                if (mod(i, 3) == 0) %si el grado de libertad es multiplo de 3, se trata de un giro
                    U1(i) = U(i); %guarda en una nueva variable para que la graficacion no se vea alterada
                else
                    U1(i) = U(i)*1000;
                end
            end  
            
            for i=1:noNudos*3
                if (mod(i, 3) == 0) %si el grado de libertad es multiplo de 3, se trata de un momento
                    R1(i) = R(i)*1000*1000; %guarda en una nueva variable para que la graficacion no se vea alterada
                else
                    R1(i) = R(i)*1000;
                end
            end  
            
            
            disp("Desplazamientos de la estructura")
            U1 = U1' %imprime la transpuesta para que se vea como un vector columna
            
            disp("Reacciones de la estructura")
            R1 = R1' %imprime la transpuesta para que se vea como un vector columna
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%Postproceso. Cálculo de Solicitaciones%%%%
        for i=1:noelementos
          qq1(i) = 0;
          qq2(i) = 0;
        end
        
        for i=1:noDistribuidas
          noelemento = distribuidas(i,1); %extrae el elemento en donde hay alguna fuerza distribuida
          qq1(noelemento) = distribuidas(i,3); %agrega la fuerza distribuida del elemento a la variable
          qq2(noelemento) = distribuidas(i,4); %agrega la fuerza distribuida del elemento a la variable
        end
        
        
        for i=1:noelementos
          %calcular Desplazamientos de cada elemento
          gdlelemento = GDLelementos(i,2:7); %extrae los grados de libertad del elemento
          uelementoGlobal = U(gdlelemento,:); %extrae el valor de los desplazamientos en el grado de libertad correspondiente
          
          elas = propSeccion(i,4); %extrae el modulo de elasticidad del elemento
          inercia = propSeccion(i,3); %extrae la inercia del elemento
          EI = elas*inercia;
          ang = anguloselementos(i); %extrae el angulo del elemento
          cs = cosd(ang);
          sn = sind(ang);
          T = [cs sn 0 0 0 0; %matriz de transformacion
               -sn cs 0 0 0 0;
               0 0 1 0 0 0;
               0 0 0 cs sn 0;
               0 0 0 -sn cs 0;
               0 0 0 0 0 1];
          
          uLocal(:,:,i) = T*uelementoGlobal; %calcula los desplazamientos del elemento
          
          %calcular solicitaciones en extremos de cada elemento 
          fuerzaselemento = kLocal(:,:,i)*uLocal(:,:,i)-fLocal(:,:,i); %calcula las AIF del elemento
          N1(i) = fuerzaselemento(1,1)*-1; %se multiplica por menos 1 para aplicar el verdadero signo 
          V2(i) = fuerzaselemento(2,1);
          M3(i) = fuerzaselemento(3,1)*-1; %se multiplica por menos 1 para aplicar el verdadero signo 
          N4(i) = fuerzaselemento(4,1);
          V5(i) = fuerzaselemento(5,1)*-1; %se multiplica por menos 1 para aplicar el verdadero signo 
          M6(i) = fuerzaselemento(6,1);
          % Calcular constantes C1 y C2
          C2(i) = M3(i); %mediante las condiciones de borde, cuando x=0 M=M3, x=L M=M6
          L = longelementos(i);
          C1(i) = M6(i)/L+(qq2(i)-qq1(i))*L/6+qq1(i)*L/2-M3(i)/L;
          
          % Calculo de C3 y C4
          u3 = uLocal(3,1,i);
          u2 = uLocal(2,1,i);
          C3(i) = EI*u3;
          C4(i) = EI*u2;
          
          
          n1 = -(qq2(i)-qq1(i))/(6*L);
          n2 = -qq1(i)/2;
          n3 = C1(i);
          n4 = C2(i);
          
          cadena = [num2str(n1),"x^3+(",num2str(n2),")x^2+(",num2str(n3),")x+(",num2str(n4),")"]; %ecucacion de momento para cada elemento
          disp(["M(x) elemento ",num2str(i)," = ", cadena ])
        end
        
        
        
        
        
        
        
        
        
        
        %%+++++++++++++++++++++++++++++++++++++%%
        %%++++++++++++++GRAFICACIÓN++++++++++++%%
        %%+++++++++++++++++++++++++++++++++++++%%
        
        %%%%%%%%graficacion de portico sin deformar%%%%%%
        hold on %mantiene todas las graficas y no las reinicia
        for i=1:noelementos
          nudoInicio = elementos(i,2); %extrae el nodo inicial del elemento
          nudoFin = elementos(i,3); %extrae el nodo final del elemento
          
          
          xi = nudos(nudoInicio,2); %extrae las coordenadas del nodo inicial
          xf = nudos(nudoFin,2); %extrae las coordenadas del nodo final
          yi = nudos(nudoInicio,3); %extrae las coordenadas del nodo inicial
          yf = nudos(nudoFin,3); %extrae las coordenadas del nodo final
        
          plot([xi xf],[yi yf],"color","k","linewidth",3) %grafica todos los elementos con un colo negro y espesor de linea de 3
        end
        
        for i=1:noNudos
          xx = nudos(i,2);
          yy = nudos(i,3);
          scatter(xx,yy,"marker","o",...
              "markeredgecolor","k",...
              "markerfacecolor","r"); %grafica todos los nodos con un borde negro y contorno rojo
        end
        
        
        %%Graficacion de pórtico deformado%%%%
        if grafDef ==1
              for i=1:noNudos
                  
                  gdlX = i*3-2;
                  gdlY = i*3-1;
                  
                  nudosDesp(i,1) = nudos(i,1);
                  nudosDesp(i,2) = nudos(i,2)+U(gdlX)*escDef; %a la coordenada del nudo le suma el desplazamiento calculado
                  nudosDesp(i,3) = nudos(i,3)+U(gdlY)*escDef; %a la coordenada del nudo le suma el desplazamiento calculado
              end
        
              for i=1:noelementos
              %Plot de las elementos deformadas rectas
                  nudoInicio = elementos(i,2); %extrae el nodo inicial del elemento
                  nudoFin = elementos(i,3); %extrae el nodo final del elemento
                  
                  
                  xi = nudosDesp(nudoInicio,2);
                  xf = nudosDesp(nudoFin,2);
                  yi = nudosDesp(nudoInicio,3);
                  yf = nudosDesp(nudoFin,3);
        
          
            %plot de la gráfica deformada con curvas
            
                xo = nudosDesp(nudoInicio,2); %la coordenada inicial en x sera la coordenada desplazada
                yo = nudosDesp(nudoInicio,3); %la coordenada inicial en y sera la coordenada desplazada
                ang = anguloselementos(i); %extrae el angulo del elemento
                L = longelementos(i); %extrae la longitud del elemento
                elas = propSeccion(i,4); %extrae el modulo de elasticidad del elemento
                inercia = propSeccion(i,3); %extrae la inercia del elemento
                EI = elas*inercia;
        
                x = linspace(0,L,11); %parte en 11 la longitud del elemento
                u1 = uLocal(1,1,i); %desplazamiento local del extremo
                u4 = uLocal(4,1,i); %desplazamiento local del extremo
                deltaU = (u4-u1)*escDef;
                
                for j=1:11
                  xEsc(j) = x(j)+deltaU/10*(j-1);
                end
                
                y =1/EI*(-(qq2(i)-qq1(i))/L*x.^5/120-qq1(i)*x.^4/24+C1(i)*x.^3/6+C2(i)*x.^2/2+C3(i)*x+C4(i)); %ecuacion hallada con las condiciones de borde
                yEsc = (y-y(1))*escDef;
                
                xRot = xEsc*cosd(ang)-yEsc*sind(ang); %rota la coordenada x a la inclinacion del elemento
                yRot = xEsc*sind(ang)+yEsc*cosd(ang); %rota la coordenada y a la inclinacion del elemento        
                xTras = xRot+xo; %traslada la coordenada x hasta donde esta el elemento
                yTras = yRot+yo; %traslada la coordenada y hasta donde esta el elemento    
                plot(xTras,yTras,"linewidth",1,"color","b") %grafica con color azul y grosor de 1
        
            end
        end
        %Plotear diagrama de momento flector
        if grafM ==1
              for i=1:noelementos
                nudoInicio = elementos(i,2);
                nudoFin = elementos(i,3);
                
                
                xo = nudos(nudoInicio,2);
                yo = nudos(nudoInicio,3);
                ang = anguloselementos(i);
                L = longelementos(i);
                
                x = linspace(0,L,11); %parte la longitud del elemento en 11 partes
                y = -(qq2(i)-qq1(i))/L*x.^3/6-qq1(i)*x.^2/2+C1(i)*x+C2(i); %ecuacion hallada con condiciones de borde
                yEsc = -y*escM; %por convencion se cambia el signo
                
                xViga = x;
                yViga = linspace(0,0,11);
        
                %Rotacion de funcion de X
                xRot = x*cosd(ang)-yEsc*sind(ang); %rota la coordenada a la inclinacion del elemento
                yRot = x*sind(ang)+yEsc*cosd(ang); %rota la coordenada a la inclinacion del elemento
                %Rotacion de viga
                xVigaRot = xViga*cosd(ang)-yViga*sind(ang); %rota la coordenada a la inclinacion del elemento
                yVigaRot = xViga*sind(ang)+yViga*cosd(ang); %rota la coordenada a la inclinacion del elemento
        
                %Trasladar la funcion de X
                xTras = xRot+xo; %traslada la coordenada hasta donde esta el elemento
                yTras = yRot+yo; %traslada la coordenada hasta donde esta el elemento
                xVigaTras = xVigaRot+xo; %traslada la coordenada hasta donde esta el elemento
                yVigaTras = yVigaRot+yo; %traslada la coordenada hasta donde esta el elemento
                plot(xTras,yTras,"linewidth",1,"color",[0.5 0.5 1]) %grafica con color azul claro y grosor de 1
        
                for i=1:11
                  plot([xTras(i) xVigaTras(i)], [yTras(i) yVigaTras(i)],"linewidth",1,"color",[0.6 0.6 0.6]) %grafica las lineas perpendiculares al elemento para generar relleno a la grafica
                end
        
                text(xTras(1),yTras(1),num2str(y(1)),"rotation",ang) %pone valores numericos al inicio de la grafica
                text(xTras(11),yTras(11),num2str(y(11)),"rotation",ang) %pone valores numericos en la mitad de la grafica
                text(xTras(6),yTras(6),num2str(y(6)),"rotation",ang) %pone valores numericos al final de la grafica
                
              end %termina bucle for de Momentos
        end %termina el condicionante If de momento
        
        
        %Plotear diagrama de cortantes
        if grafV ==1
              for i=1:noelementos
                nudoInicio = elementos(i,2);
                nudoFin = elementos(i,3);
                
                
                xo = nudos(nudoInicio,2);
                yo = nudos(nudoInicio,3);
                ang = anguloselementos(i);
                L = longelementos(i);
                
                x = linspace(0,L,11); %parte la longitud del elemento en 11 partes
                y = -(qq2(i)-qq1(i))/L*x.^2/2-qq1(i)*x+C1(i); %ecuacion hallada con condiciones de borde
                yEsc = y*escV; %por convencion se mantiene el signo
                
                xViga = x;
                yViga = linspace(0,0,11);
        
                %Rotacion de funcion de X
                xRot = x*cosd(ang)-yEsc*sind(ang); %rota la coordenada a la inclinacion del elemento
                yRot = x*sind(ang)+yEsc*cosd(ang); %rota la coordenada a la inclinacion del elemento
                %Rotacion de viga
                xVigaRot = xViga*cosd(ang)-yViga*sind(ang); %rota la coordenada a la inclinacion del elemento
                yVigaRot = xViga*sind(ang)+yViga*cosd(ang); %rota la coordenada a la inclinacion del elemento
        
                %Trasladar la funcion de X
                xTras = xRot+xo; %traslada la coordenada hasta donde esta el elemento
                yTras = yRot+yo; %traslada la coordenada hasta donde esta el elemento
                xVigaTras = xVigaRot+xo; %traslada la coordenada hasta donde esta el elemento
                yVigaTras = yVigaRot+yo; %traslada la coordenada hasta donde esta el elemento
                plot(xTras,yTras,"linewidth",1,"color",[0.3 0.7 0.3]) %grafica con color claro y grosor de 1
        
                for i=1:11
                  plot([xTras(i) xVigaTras(i)], [yTras(i) yVigaTras(i)],"linewidth",1,"color",[0.6 0.6 0.6]) %grafica las lineas perpendiculares al elemento para generar relleno a la grafica
                end
        
                text(xTras(1),yTras(1),num2str(y(1)),"rotation",ang) %pone valores numericos al inicio de la grafica
                text(xTras(11),yTras(11),num2str(y(11)),"rotation",ang) %pone valores numericos al final de la grafica
                
              end %termina bucle for de Cortantes
        end %termina el condicionante If de Cortantes
        
        %Plotear diagrama de Axiales
        if grafN ==1
              for i=1:noelementos
                nudoInicio = elementos(i,2);
                nudoFin = elementos(i,3);
                
                
                xo = nudos(nudoInicio,2);
                yo = nudos(nudoInicio,3);
                ang = anguloselementos(i);
                L = longelementos(i);
                
                x = linspace(0,L,11); %parte la longitud del elemento en 11 partes
                y = (N4(i)-N1(i))/L*x+N1(i); %ecuacion hallada con condiciones de borde
                yEsc = y*escN; %por convencion se mantiene el signo
                
                xViga = x;
                yViga = linspace(0,0,11);
        
                %Rotacion de funcion de X
                xRot = x*cosd(ang)-yEsc*sind(ang); %rota la coordenada a la inclinacion del elemento
                yRot = x*sind(ang)+yEsc*cosd(ang); %rota la coordenada a la inclinacion del elemento
                %Rotacion de viga
                xVigaRot = xViga*cosd(ang)-yViga*sind(ang); %rota la coordenada a la inclinacion del elemento
                yVigaRot = xViga*sind(ang)+yViga*cosd(ang); %rota la coordenada a la inclinacion del elemento
        
                %Trasladar la funcion de X
                xTras = xRot+xo; %traslada la coordenada hasta donde esta el elemento
                yTras = yRot+yo; %traslada la coordenada hasta donde esta el elemento
                xVigaTras = xVigaRot+xo; %traslada la coordenada hasta donde esta el elemento
                yVigaTras = yVigaRot+yo; %traslada la coordenada hasta donde esta el elemento
                plot(xTras,yTras,"linewidth",1,"color",[0.7 0.3 0.3]) %grafica con color claro y grosor de 1
        
                for i=1:11
                  plot([xTras(i) xVigaTras(i)], [yTras(i) yVigaTras(i)],"linewidth",1,"color",[0.6 0.6 0.6]) %grafica las lineas perpendiculares al elemento para generar relleno a la grafica
                end
        
                text(xTras(1),yTras(1),num2str(y(1)),"rotation",ang) %pone valores numericos al inicio de la grafica
                text(xTras(11),yTras(11),num2str(y(11)),"rotation",ang) %pone valores numericos al final de la grafica
                
              end %termina bucle for de Normal
        end %termina el condicionante If de Normal
%%

%%
 


 else 
   
%%En caso de que la estructura a calcular no sea un pórtico ni una viga sino una cercha el código correspondiente empezará 
% desde acá, es decir, todos los datos deben ingresarse desde acá abajo en caso de tener datos ingresados en la parte 
% del código de pórticos y vigas estos no tendrán efecto en la ejecución del posterior código.
   %************Ingreso de datos (CRCH)***************
        
        %Desde este punto se ingresarán todos los datos necesarios para generar la estructura y sus valores desconocidos. 
        % Por favor seguir las indicaciones para cada lectura.
        
        %Ingrese las unidades de los resultados = [1 si son KN y m, 2 si son KN y mm
        %3 si son Kg y m, 4 si son Kg y mm, 5 si son Ton y m, 6 si son Ton y mm
        %7 si son N y m, 8 si son N y mm]
        UnidadesSalida = 8;
        
        %Ingreso de nodos de la forma
        % nodos=[No.Nodo, coordenada x, coordenada y]
        
        nodos =[1 0 0; 
                2 0.25 0;
                3 0.25 0.15;
                4 0.5 0.15;
                5 0.5 0.3;
                6 0.75 0.3;
                7 0.75 0.45;
                8 1 0.45;
                ];
        
        %Ingreso de barras de la forma
        % barras = [No.Barra, nodo de inicio, nodo de fin]
        
        barras = [1 1 2;
                  2 1 3;
                  3 2 3;
                  4 2 4;
                  5 3 4;
                  6 3 5;
                  7 4 5;
                  8 4 6;
                  9 5 6;
                  10 5 7;
                  11 6 7;
                  12 6 8;
                  13 7 8];
        % Defina las propiedades estas deben ser introducidas en (metros, Kpa)
        
        A330 = 0.00626;
        IPE330 = 0.00011768;
        
        A360= 0.00727;
        IPE360 = 0.00016267;
        
        A400 = 0.00845;
        IPE400 = 0.00023131;
        
        
        E = 2e8;

        %Defina las propiedades y materiales de cada barra, estas se deben introducir de la forma
        % propiedadesSeccion = [No.Barra, Área, Módulo de Elasticidad]
        
        propiedadesSeccion = [1 A330 E;
                              2 A360 E;
                              3 A400 E;
                              4 A330 E;
                              5 A400 E;
                              6 A360 E;
                              7 A400 E;
                              8 A330 E;
                              9 A400 E;
                              10 A360 E;
                              11 A400 E;
                              12 A330 E;
                              13 A360 E
                              ];
                                
                                
        
        %Ingreso de cargas, estas se deben ingresar de la forma
        %cargas = [No.Nodo Fx Fy]. La convención de signos debe ser positivo hacia la derecha, negativo hacia
        % la izquierda, positivo hacia arriba, negativo hacia abajo. Deben ser introducidas en (kN)
        cargas = [1 0 -17;
                  2 -20 0;
                  3 0 -17;
                  5 0 -17;
                  6 -5 0;
                  7 0 -17;
                  8 0 -20;
                 ];
        
        %Ingreso de apoyos estos se deben ingresar de la forma
        %apoyos = [No.Nodo, restricción x, restricción y] Restricción se refiere a la capacidad de desplazarse 
        % en esa dirección, si no se puede mover en esa dirección el valor a colocar será 1, si se puede mover es 0.
        apoyos = [1 1 1;
                  2 0 1;
                  8 1 1];
        
        %Graficación En este punto se decidirá qué quiere ser mostrado en el gráfico. Recuerde que para no tener errores en el código sólo un elemento de graficación debe estár activo
        %Si quiere que se grafique algo quite el 0 y ponga 1
        %Gráfica de la deformada de la cercha. Al ser los valores de desplazamiento tan pequeños es necesario usar escalas grandes.
        GrafD = 0;
        escala = 2000;
           %Gráfica de las fuerzas internas de la cercha. Al ser los valores de axial en cada barra tan grandes es necesario usar escalas pequeñas.
        GrafA = 1;
        escalaA = 0.002;
       
        %En este punto acaba el ingreso de datos. Al correr el programa obtendrá el vector de desplazamientos y de reacciones en consola. 
        %**************************************************************************
        %****************************INICIO DEL PROGRAMA***************************
        %**************************************************************************
        %No tocar
        
        %calculo de longitudes y ángulos de cada barra
        %”temp” es una variable temporal que obtiene el tamaño de la matriz en donde están los nodos y “Nonodos” traerá solo la primera columna, correspondiente al número de nodos de la cercha. De igual forma pasa con “temp1” y “Nobarras”.
        temp = size(nodos);
        Nonodos = temp(1);
        temp1 = size(barras);
        Nobarras = temp1(1);
        
        for i=1:Nobarras %este ciclo se repetirá para cada una de las barras de la cercha
            nodoinicio = barras(i,2); %extrae el valor del nodo de inicio de la matriz de barras (columna 2)
            nodofin = barras(i,3);%extrae el valor del nodo de fin de la matriz de barras (columna 3)
            xinicio = nodos(nodoinicio,2);%extrae el valor de la coordenada de inicio en x del nodo de inicio obtenido en las anteriores líneas
            yinicio = nodos(nodoinicio,3);%extrae el valor de la coordenada de inicio en y del nodo de inicio obtenido en las anteriores líneas
            xfin = nodos(nodofin,2);%extrae el valor de la coordenada de fin en x del nodo de inicio obtenido en las anteriores líneas
            yfin = nodos(nodofin,3);%extrae el valor de la coordenada de fin en y del nodo de inicio obtenido en las anteriores líneas
            longitudBarra(i) = sqrt((xfin-xinicio)^2+(yfin-yinicio)^2);%Al hacer uso del teorema de pitágoras se obtiene la distancia entre dos puntos, la cual será la longitud de la barra.
            angulosBarras(i) = acosd((xfin-xinicio)/longitudBarra(i));%Al hacer uso de la definición de coseno se halla el ángulo con respecto a la horizontal. 
            if (yfin-yinicio)<0
                angulosBarras(i) = -angulosBarras(i); %dado que algunos ángulos pueden ser mayores que 180° se tiene que tomar la parte negativa, es decir, el complemento.
            end
        end
        
        %calculo de matrices de rigidez locales de cada barra
        for i=1:Nobarras%este ciclo se repetirá para cada una de las barras de la cercha
            AELP=(propiedadesSeccion(i,2)*propiedadesSeccion(i,3))/longitudBarra(i);%para disminuir la cantidad de valores introducidos en la matriz de rigidez local se efectua la operación antes de introducir a la matriz. Esto trayendo los valores de las matrices digitadas en la entrada de datos.
            klocal(:,:,i)=[AELP -AELP;
                         -AELP AELP];%se guarda (:,:,i) dado que k será un tensor y para llenar todos los espacios se debe usar esta notación
        end
        
        %calculo de matrices de transformación para cada barra para calcular K
        %global de cada elemento
        for i=1:Nobarras%este ciclo se repetirá para cada una de las barras de la cercha
            cs = cosd(angulosBarras(i));%guarda el valor del coseno director de la barra i
            sn = sind(angulosBarras(i));%guarda el valor del seno director de la barra i
            Tlg = [cs 0;%genera la matriz de transformación de desplazamientos para cada barra
                   sn 0;
                   0 cs;
                   0 sn];
            Tgl = Tlg';%genera la matriz de transformación de cargas para cada barra al trasponer la matriz de transformación de desplazamientos.
            kglobal(:,:,i) = Tlg*klocal(:,:,i)*Tgl;%genera la matriz de rigidez global de cada barra
        end
        
        %calculo de matriz global de toda la estructura
        %antes de generar la matriz es necesario definir los grados de libertad de la cercha
        for i=1:Nobarras
            GDLbarras(i,1) = barras(i,1);%trae todas las barras de la cercha y la guarda en la primera columna
            GDLbarras(i,2) = barras(i,2)*2-1;%genera el grado de libertad en x para el primer nodo de la barra, el cual será 1 para el nodo 1.
            GDLbarras(i,3) = barras(i,2)*2;%genera el grado de libertad en y para el primer nodo de la barra, el cual será 2 para el nodo 1.
            GDLbarras(i,4) = barras(i,3)*2-1;%genera el grado de libertad en x para el ultimo nodo de la barra, el cual será 3 para el nodo 2.
            GDLbarras(i,5) = barras(i,3)*2;%genera el grado de libertad en y para el ultimo nodo de la barra, el cual será 4 para el nodo 2.
        end
        %se crea la matriz de tamaño 2*nodos x 2*nodos con todos los grados de libertad llena de 0 que sera llenada por la suma de las matrices de rigidez globales de cada elemento
        K = zeros(Nonodos*2,Nonodos*2);
        
        for i=1:Nobarras
            gdl = GDLbarras(i,2:5);%se extraen los valores de grados de libertad para cada i barra
            K(gdl,gdl) = kglobal(:,:,i) + K(gdl,gdl);%se llenan los valores de la matriz de la suma de las filas y columnas de los grados de libertad de cada matriz global por barra. Para que se sumen completamente, se suman de nuevo los valores que ya habían en la iteración anteriormente.
        end
        %creacion de vector de cargas
        %se crea el vector de cargas lleno de 0 que sera llenado con las cargas de
        %cada nodo
        Q = zeros(Nonodos*2,1);
        temp2 = size(cargas);%temporal que guardará el tamaño de la matriz  del vector de cargas introducido
        Nocargas = temp2(1);%obtendrá el valor de número de cargas de la variable temporal
        
        for i=1:Nocargas%se repite solo por la cantidad de cargas que se tengan en la estructura
            Q(cargas(i,1)*2-1,1) = cargas(i,2);%introduce el valor de la carga introducida en el vector de cargas buscando cada grado de libertad
            Q(cargas(i,1)*2,1) = cargas(i,3);
        end 
        
        %se tienen que crear las restricciones tanto en la matriz de rigidez como
        %en el vector de cargas, esto debido a los apoyos en donde el
        %desplazamiento es conocido. así se taparán las filas y columnas conocidas obteniendo solamente los valores de las desconocidas.
        
        temp3 = size(apoyos);%temporal que guardará el tamaño de la matriz  del vector de apoyos introducido
        Noapoyos = temp3(1);%obtendrá el valor de número de apoyos de la variable temporal

        contador = 0;%variable que guardará ciclo a ciclo el grado de libertad de los valores de desplazamiento conocidos
        for i=1:Noapoyos% se repetira por cada apoyo
            if (apoyos(i,2)==1)%observará si se tiene restricción en x
                contador = contador + 1;%el contador sumará 1
                gdl12(contador) = apoyos(i,1)*2-1;%buscará el grado de libertad donde el desplazamiento es 0 en x y lo guardará en gdl conocidos
            end
            if (apoyos(i,3)==1)%observará si se tiene restricción en y
                contador = contador + 1;
                gdl12(contador) = apoyos(i,1)*2;%buscará el grado de libertad donde el desplazamiento es 0 en y y lo guardará en gdl conocidos
            end
        end
        
        %Solucion sistema de ecuaciones
        %creacion de matriz K reducida y cargas reducida, sin filas ni columnas conocidas
        K11 = K;
        K11(:,gdl12) = [];%elimina las columnas de grados de libertad conocidas
        K11(gdl12,:) = [];%elimina las filas de grados de libertad conocidas
        Q11 = Q;
        Q11(gdl12,:) = [];%elimina las filas de grados de libertad conocidas
        
        %resolviendo el sistema K*delta=Q en donde solo se tienen valores desconocidos
        
        delta11 = linsolve(K11,Q11);%resuelve el sistema por GaussJordan y obtiene los desplazamientos para cada nodo
        
        %Vector de desplazamientos completo
        
        gdlLibres = 1:1:Nonodos*2;%crea un vector con los grados de libertad que no están restringidos
        gdlLibres(gdl12) = [];%los grados de libertad que están restringidos los llena en espacio vacío
        
        Delta = zeros(Nonodos*2,1);%genera un vector de 0 que será llenado con los valores de desplazamiento obtenidos
        temp4 = size(gdlLibres);%temporal que guarda el tamaño de la matriz de grados de liberad que no estan restringidos
        NogdlLibres = temp4(2);%guarda el tamaño de los grados de libertad que no están restringidos
        
        for i=1: NogdlLibres %se repetira para los grados de libertad que no están restringidos
            Delta(gdlLibres(i)) = delta11(i);%llena los espacios del vector con los valores de grados de libertad que no están restringidos
        end
        
        %Vector de reacciones
        Delta
        R = K*Delta-Q%dado que ya se tiene todos los datos generará un vector con todas las reacciones en todos los grados de libertad
        
        %Conversiones de Unidades
        if (UnidadesSalida==1) %KN m - KN m
            
            disp("Desplazamientos de la estructura")
            Delta
            
            disp("Reacciones de la estructura")
            R
            
        elseif (UnidadesSalida==2) %KN m - KN mm
            
            for i=1:Nonodos*2
                Delta1(i) = Delta(i)*1000;
                R1(i) = R(i)*1;
            end  
           
            disp("Desplazamientos de la estructura")
            Delta1 = Delta1'
            
            disp("Reacciones de la estructura")
            R1 = R1'
        
        elseif (UnidadesSalida==3) %KN m - Kg m
            
            for i=1:Nonodos*2
                Delta1(i) = Delta(i);
                R1(i) = R(i)*101.97;
            end  
           
            disp("Desplazamientos de la estructura")
            Delta1 = Delta1'
            
            disp("Reacciones de la estructura")
            R1 = R1'
            
        elseif (UnidadesSalida==4) %KN m - Kg mm
            
            for i=1:Nonodos*2
                Delta1(i) = Delta(i)*1000;
                R1(i) = R(i)*101.97;
            end  
           
            disp("Desplazamientos de la estructura")
            Delta1 = Delta1'
            
            disp("Reacciones de la estructura")
            R1 = R1'
            
        elseif (UnidadesSalida==5) %KN m - Ton m
            
            for i=1:Nonodos*2
                Delta1(i) = Delta(i);
                R1(i) = R(i)*0.1;
            end  
           
            disp("Desplazamientos de la estructura")
            Delta1 = Delta1'
            
            disp("Reacciones de la estructura")
            R1 = R1'
            
        elseif (UnidadesSalida==6) %KN m - Ton mm
            
            for i=1:Nonodos*2
                Delta1(i) = Delta(i)*1000;
                R1(i) = R(i)*0.1;
            end  
           
            disp("Desplazamientos de la estructura")
            Delta1 = Delta1'
            
            disp("Reacciones de la estructura")
            R1 = R1'
            
        elseif (UnidadesSalida==7) %KN m - N m
            
            for i=1:Nonodos*2
                Delta1(i) = Delta(i);
                R1(i) = R(i)*1000;
            end  
           
            disp("Desplazamientos de la estructura")
            Delta1 = Delta1'
            
            disp("Reacciones de la estructura")
            R1 = R1'
            
        elseif (UnidadesSalida==8) %KN m - N mm
            
            for i=1:Nonodos*2
                Delta1(i) = Delta(i)*1000;
                R1(i) = R(i)*1000;
            end  
           
            disp("Desplazamientos de la estructura")
            Delta1 = Delta1'
            
            disp("Reacciones de la estructura")
            R1 = R1'
            
        end
        
        %con los valores de cargas globales obtenidos, se pueden obtener las cargas locales
        %Calculo de cargas internas
        for i=1:Nobarras
            gdl = GDLbarras(i,2:5);%se extraen los valores de grados de libertad para cada i barra
            cs = cosd(angulosBarras(i));%guarda el valor del coseno director de la barra i 
            sn = sind(angulosBarras(i));%guarda el valor del seno director de la barra i
            Tlg = [cs 0;%genera la matriz de transformacion de desplazamientos
                   sn 0;
                   0 cs;
                   0 sn];
            Tgl = Tlg';%genera la matriz de transformacion de cargas
            cargaLocal = klocal(:,:,i)*Tgl*Delta(gdl,:);%obtiene los valores de carga local de inicio y final de la barra i
            deltaLocal = linsolve(cargaLocal,klocal(:,:,i));%obtiene los valores de desplazamiento local de inicio y final de la barra i
            N1(i) = cargaLocal(1,1)*-1;%guarda los valores de carga local de inicio de la barra
            N2(i) = cargaLocal(2,1);%guarda los valores de carga local de fin de la barra
	        d1(i) = deltaLocal(1,1);%guarda los valores de desplazamiento local de inicio de la barra
            d2(i) = deltaLocal(1,2);%guarda los valores de desplazamiento local de fin de la barra
        end
        N1;
        N2;
        d1;
        d2;

        %En este punto se han obtenido todos los valores necesarios al hacer uso de la matriz de rigidez. Desde este punto se graficaron estos valores.
        
        %***********************************************
        %**********GRAFICACION CARGAS INTERNAS**********
        %***********************************************
        %No tocar
        
        %graficación sin deformar
        
        hold on
        
        for i=1:Nobarras
            nodoinicio = barras(i,2);
            nodofin = barras(i,3);
            xinicio = nodos(nodoinicio,2);
            yinicio = nodos(nodoinicio,3);
            xfin = nodos(nodofin,2);
            yfin = nodos(nodofin,3);
            plot([xinicio xfin],[yinicio yfin],"Color","K","LineWidth",2.5)
        end
        
        for i=1:Nonodos
            xx = nodos(i,2);
            xy = nodos(i,3);
            scatter(xx,xy,"marker","o","MarkerEdgeColor","k","MarkerFaceColor","r")
        end
        
        %graficación deformada
        
        for i=1:Nonodos
            nodosdesplazados(i,1) = nodos(i,1);
            nodosdesplazados(i,2) = nodos(i,2) + Delta(i*2-1)*escala;
            nodosdesplazados(i,3) = nodos(i,3) + Delta(i*2)*escala;
        end
        
        if GrafD == 1
            for i=1:Nobarras
                nodoinicio = barras(i,2);
                nodofin = barras(i,3);
                xinicio = nodosdesplazados(nodoinicio,2);
                yinicio = nodosdesplazados(nodoinicio,3);
                xfin = nodosdesplazados(nodofin,2);
                yfin = nodosdesplazados(nodofin,3);
                plot([xinicio xfin],[yinicio yfin],"Color","b","LineWidth",1)
            end
        end
        
        %graficación axial
        
        if GrafA == 1
            for i=1:Nobarras
                nodoinicio = barras(i,2);
                nodofin = barras(i,3);
                xinicio = nodos(nodoinicio,2);
                yinicio = nodos(nodoinicio,3);
                xfin = nodos(nodofin,2);
                yfin = nodos(nodofin,3);
                longitudBarra(i) = sqrt((xfin-xinicio)^2+(yfin-yinicio)^2);
                angulosBarras(i) = acosd((xfin-xinicio)/longitudBarra(i));
                 if (yfin-yinicio)<0
                    angulosBarras(i) = -angulosBarras(i);
                end
        
                x = linspace(0,longitudBarra(i),11);
                y = (N2(i)-N1(i))/longitudBarra(i)*x+N1(i);
                yesc = y*escalaA;
                
                %rotacion funcion x
                xRot = x*cosd(angulosBarras(i))-yesc*sind(angulosBarras(i));
                yRot = x*sind(angulosBarras(i))+yesc*cosd(angulosBarras(i));
                %rotacion barra x
                xbarraRot = x*cosd(angulosBarras(i))-linspace(0,0,11)*sind(angulosBarras(i));
                ybarraRot = x*sind(angulosBarras(i))+linspace(0,0,11)*cosd(angulosBarras(i));
                %trasladar funcion x
                xTras = xRot+xinicio;
                yTras = yRot+yinicio;
                xbarraTras = xbarraRot+xinicio;
                ybarraTras = ybarraRot+yinicio;
                plot(xTras,yTras,"LineWidth",1,"Color",[0.3 0.7 0.3]);
        
                for i=1:11
                    plot([xTras(i) xbarraTras(i)], [yTras(i) ybarraTras(i)],"LineWidth",1,"color",[0.6 0.6 0.6])
                end
        
                text(xTras(1),yTras(1),num2str(y(1)),"Rotation",angulosBarras(i));
                text(xTras(11),yTras(11),num2str(y(11)),"Rotation",angulosBarras(i));
        
            end
        end
  end


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%% ANEXO %%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%

  %%%%%%%%%%%VIGA%%%%%%%%%%%%%

  %nudos = [
                %1 0 0;
                %2 4 0;
                %3 5 0;
                %4 6 0;
                %5 10 0;
                %6 12 0
                %]

 %elementos = [
                   % 1 1 2;
                   % 2 2 3;
                   % 3 3 4;
                   % 4 4 5;
                   % 5 5 6
                   % ]

 %A200=3054.65325411/1000^2;
  %A270=4495.21325411/1000^2;
   %IPE200=20956967.38621783/1000^4;
    %IPE270=56454954.51267295/1000^4;
     %E = 200000e6;

  %propSeccion = [
       % 1 A200 IPE200 E;
       % 2 A270 IPE270 E;
       % 3 A200 IPE200 E;
       % 4 A270 IPE270 E;
       % 5 A270 IPE270 E
        %]

% puntuales = [
                   % 1 0 -20 0;
                   % 3 0 0 50
                  %  ]

%distribuidas = [
                       % 4 0 0 40;
                       % 5 0 -30 -30
                       % ]

%apoyos = [
              %  2 1 1 0;
              %  4 0 1 0;
              %  5 0 1 0;
              %  6 1 1 1
               % ]

%grafDef = 1;
       % escDef = 700;
        
       % grafM = 0;
       % escM = 0.01;
        
       % grafV = 0;
       % escV = 0.01;
        
       % grafN = 0;
       % escN = 0.005;