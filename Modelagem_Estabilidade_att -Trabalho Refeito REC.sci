//Daniel Marins Silva D´Oliveira Martins : 12554571
//David Lopes F. de Carvalho :             12610305
//João Pedro Kawachi Chaves :              12745188
//Leonardo Fortunato de Carvalho :         9837910


clear; // Limpeza de variáveis
close(winsid()) // Fecha as janelas gráficas abertas anteriormente

//vetor de estados do estado inicial
theta_p_1_ini = 0   //rad/s
alfa_ini = 0    //rad
theta_p_2_ini = 0       //rad/s
v_ini =  25     //m/s
theta_2_ini = 0   //rad
//Estes valores também dependerão de cada situação inicial diferente
//Portanto, estes valores podem ter pequenas variações entre simulações
/////////////////////////////////////

//vetor de espaços do estado inicial
X0 = [theta_p_1_ini; alfa_ini; theta_p_2_ini; v_ini ; theta_2_ini] 
////////////////////////////////////

//parâmetros utilizados
//comprimentos (m):
a = 1.138
b = 2.412
c = -0.8
d = 4.67
e = 3.73

dist =  e - 1.25    //distância da carga ao centro de massa do reboque
//altera em cada situação de distribuição de massa, sendo a variável principal do fenômeno
//Deve levar em conta o comprimento total útil do reboque, e as dimensões da carga
//Cenário 1: dist = -d + 1.25 (carga na frente do reboque)
//Cenário 2: dist = e - 1.25 (carga no fundo do reboque)
////////////////////

//massas (kg):
mF = 5203 //massa carregada no eixo frontal (fonte: catálogo)
mT = 2454 //massa carregada no eixo traseiro (fonte: catálogo)
mc = mF + mT //massa do caminhão
mcarga = 6000
mreboque = 6000
mr = mcarga + mreboque
//////////////

//Momentos de inércia (kg m²):
Ic = 21015
Icarga = 1/12*mcarga*(2*2.5^2)
Ireboque = 20907
Ir = Icarga + mcarga*dist^2 + Ireboque
/////////////////////////////

g = 9.81 //m/s² , gravidade

//Forças verticais em cada eixo:
//Pesos em cada eixo:
PF = mF * g
PT = mT * g
PReb = mreboque * g
PE = PReb * d / (d+e)
PA = PReb * e / (d+e)
PC = mcarga * g

Fz_A = PA + PC * (e-dist)/(d+e) //Articulação
Fz_F = PF + Fz_A*(-c)/(a+b) //Frontal caminhão
Fz_T = PT + Fz_A*(a+b+c)/(a+b) //Traseiro caminhão
Fz_E = PE + PC*(d+dist)/(d+e) //Traseiro reboque (a parte frontal do reboque é apoiada ao caminhão)
////////////////////////////////


// intervalos de tempo
T = 20 //s - tempo total
Np = 10000 // numero de passos
t = linspace(0,T,Np); //vetor de tempo
//////////////////////

//função de cálculo de força lateral nos pneus
function F_lat = Fy (alfa, Fz)
    //Parâmetros fixados do pneu:
    a0 = 1.028 // -
    a1 = 2.014 //(*1000) 1/kN
    a2 = 710.501 // (*1000)-
    a3 = 5226.341 *180/%pi //(N/rad) 5226.341 // N/grau (gama = 0)
    a4 = 78.877 // kN
    a5 = 0.012 *180/%pi // (1/rad) = 0.012 1/grau
    a6 = -0.005 // 1/kN
    a7 = 0.67 // -
    
    gama = 0 // rad - ângulo de cambagem
    mi = 0.3 //Coef. atrito de operação
    //////////////////////////////
    
    
    Fz_kN = Fz/1000 // kN (Entra como kN, não N)
    mi_n = (a1*Fz_kN + a2)/1000 //Coef. atr
    
    // Equacionamento do pneu - Força lateral
    C = a0
    D = (a1*Fz_kN + a2) * Fz_kN
    BCD = a3 * sin(2* atan(Fz_kN/a4))*(1 - a5*abs(gama))
    B = BCD /( C * D )
    E = a6 * Fz_kN + a7
    
    alfa = asin(sin(alfa)) //correção de ângulos fora do intervalo de -90° a 90°
    alfa_eq = (mi_n / mi) * alfa  //correção do ângulo equivalente
    F_lat = -D*sin(C*atan(B*alfa_eq - E *(B*alfa_eq - atan(B*alfa_eq)))) * (mi / mi_n) //em N
    
endfunction
/////////////////////////////////////////////////////////////

///////////////////////////////// MODELO NÃO LINEAR /////////////////////////////////
//vetores para armazenar os valores das variáveis do vetor de estados NÃO LINEAR
theta_p_1= zeros(Np) //theta_p_1
alfa= zeros(Np) //alfa
theta_p_2= zeros(Np) //theta_p_2_ini
v= zeros(Np) //v
theta_2= zeros(Np) //theta_2
////////////////////////////////////////////////////////////////////////////////

function F_impulso = impulso(t)
    t_impulso = 1; // tempo de aplicação do impulso
    largura = 0.2; // Largura do "pulso" para fins de simulação.

    // Se o tempo atual estiver próximo do tempo de impulso, aplicar força
    if abs(t - t_impulso) < largura then
        F_impulso = 0.5*Fz_E;
    else
        F_impulso = 0;
    end
endfunction

//função não linear
function dy = penduloduplo(t,y)  //integração das equações diferenciais
    
    M11 = -d*mr*sin(y(5))
    M12 = -y(4)*(mr + mc)*sin(y(2))
    M13 = d*mr*sin(y(5))
    M14 = (mr + mc)*cos(y(2))
    M15 = 0
    
    M21 = -mr*(b+c+d*cos(y(5)))
    M22 = y(4)*(mr + mc)*cos (y(2))
    M23 = d*mr*cos(y(5))
    M24 = (mr + mc)*sin(y(2))
    M25 = 0
    
    M31 = mr*(b + c)*(b + c + d*cos(y(5))) + Ic
    M32 = -y(4)*mr*(b+c)*cos(y(2))
    M33 = -d*mr*(b+c)*cos(y(5))
    M34 = -mr*(b+c)*sin(y(2))
    M35 = 0
    
    M41 = Ir + d*mr*(d + (b+c)*cos(y(5)))
    M42 = -y(4)*d*mr*cos(y(2)+ y(5))
    M43 = -mr*d^2 - Ir
    M44 = -d*mr*sin(y(2) + y(5))
    M45 = 0
    
    M51 = 0
    M52 = 0
    M53 = 0
    M54 = 0
    M55 = 1
    
    M = [M11,M12,M13,M14,M15;M21,M22,M23,M24,M25;M31,M32,M33,M34,M35;M41,M42,M43,M44,M45;M51,M52,M53,M54,M55] //Matriz de massas
    
    //ângulos de deriva:
    angF = atan(  (y(4)*sin(y(2))+a*y(1))  /  (y(4)*cos(y(2)))  )
    angT = atan(  (y(4)*sin(y(2))-b*y(1))  /  (y(4)*cos(y(2)))  )
    angE = atan( (y(4)*sin(y(2) + y(5)) - (b+c)*y(1)*cos(y(5)) - (d+e)*(y(1) - y(3))  ) / ( y(4)*cos( y(2) + y(5) ) + (b+c)*y(1)*sin(y(5)) )  )
    
    Fy_F = Fy (angF,Fz_F) //Frontal caminhão
    Fy_T = Fy (angT,Fz_T) //Traseiro caminhão
    Fy_E = Fy (angE,Fz_E) //Eixos reboque
    
    f1 = (Fy_E+impulso(t)) * sin(y(5)) - (b+c)*y(1)^2*mr + (mc+mr)*y(4)*y(1)*sin(y(2)) - (y(1)-y(3))^2*d*mr*cos(y(5))
    f2 = Fy_T + Fy_F + (Fy_E+impulso(t))*cos(y(5)) - (mc+mr)*y(4)*y(1)*cos(y(2)) + ( y(1) - y(3) )^2*d*mr*sin(y(5))
    f3 = a*Fy_F - Fy_T*b - (b+c)*(  (Fy_E+impulso(t))*cos(y(5)) - y(4)*mr*cos(y(2))*y(1) + ( y(1)-y(3) )^2*d*mr*sin(y(5))  )
    f4 = d*( b*y(1)^2*mr*sin(y(5)) - (Fy_E-impulso(t)) + c*y(1)^2*mr*sin(y(5)) + y(4)*y(1)*mr*cos(y(2)+y(5))) - (Fy_E+impulso(t))*e
    f5 = y(3)
    
    f = [f1;f2;f3;f4;f5]
    Minv = inv(M)
    MAT = Minv*f
    
    
    dy(1)= MAT(1)
    dy(2)= MAT(2)
    dy(3)= MAT(3)
    dy(4)= MAT(4)
    dy(5)= MAT(5)
endfunction
///////////////////////////////////////////////////////////////////////
//X0 = [theta_p_1_ini; alfa_ini; theta_p_2_ini; v_ini ; theta_2_ini] 
//resolvendo a EDO
y = ode(X0,0,t,penduloduplo);

theta_p_1= y(1,:); //theta_p_1 (rad/s)
alfa= y(2,:); //alfa (rad)
theta_p_2= y(3,:); //theta_p_2 (rad/s)
v= y(4,:); //v (m/s)
theta_2= y(5,:); //theta_2 (rad)
///////////////////

//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ MODELO NÃO LINEAR ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

//Plotando gráficos da não linear:
//Ângulos

limalfa = 90*ones(t);
alfa_deg = alfa * 180 / %pi; // Convertendo radianos para graus
t_limite = length(t);

function t_limite=determinatlimite(t,vetor,limite) //Função que acha o t em que o valor do vetor extrapola o limite
    t_limite = t;
    i = 1;
    while i<=t_limite
        if vetor(i)>limite
            t_limite = i;
        end
        i = i+1;
    end
endfunction

t_limite = determinatlimite(t_limite,alfa_deg,limalfa);

function grafico=corrigeangulo(limite,angulo) //Função para analisar se o determinada variável extrapola o limite determinado. Se extrapolar, o valor é retirado do vetor. Se não, o vetor se mantém como antes.
    grafico = angulo;
    
    for i=1:length(angulo) //
    if abs(angulo(i)) > limite then
        grafico(i) = %nan;
    end
end
endfunction

alfa_deg = corrigeangulo(limalfa,alfa_deg);

scf(1); 
plot2d(t,limalfa,5);
plot2d(t,-limalfa,5);
plot2d(t,alfa_deg,2);
title("Ângulo alfa em função do tempo - Comparação com limite do modelo" ,'fontsize', 4);
xlabel("t(s)", 'fontsize', 3);
ylabel("Pos. angular(º)", 'fontsize', 3);

scf(2); 
plot2d(t,alfa_deg,2);
title("Ângulo alfa em função do tempo " ,'fontsize', 4);
xlabel("t(s)", 'fontsize', 3);
ylabel("Pos. angular(º)", 'fontsize', 3);

limtheta2 = 90*ones(t);
theta_2_deg = theta_2*180/%pi;
t_limite = determinatlimite(t_limite,theta_2_deg,limtheta2);
theta_2_deg = corrigeangulo(limtheta2,theta_2_deg);

scf(3); 
plot2d(t,limtheta2,5);
plot2d(t,-limtheta2,5);
plot2d(t,theta_2_deg,2);
title("Ângulo theta 2 em função do tempo - Comparação com limite do modelo" ,'fontsize', 4);
xlabel("t(s)", 'fontsize', 3);
ylabel("Pos. angular(º)", 'fontsize', 3);

scf(4); 
plot2d(t,theta_2_deg,2);;
title("Ângulo theta 2 em função do tempo " ,'fontsize', 4);
xlabel("t(s)", 'fontsize', 3);
ylabel("Pos. angular(º)", 'fontsize', 3);


function velocidade = corrigevelocidade(t,vetor) //Função para fazer com que as velocidades parem de ser exibidas quando atingirem o tempo em que pelo menos um dos ângulos extrapolou o limite determinado
    velocidade = vetor;
    for i=1:length(vetor)
        if i>t then
            velocidade(i) = %nan;
        end
    end
    
endfunction

theta_p_1_corrigido = corrigevelocidade(t_limite,theta_p_1);
theta_p_2_corrigido = corrigevelocidade(t_limite,theta_p_2);
v_corrigido = corrigevelocidade(t_limite,v);

//Velocidades angulares
scf(5);
plot2d(t,[theta_p_1_corrigido',theta_p_2_corrigido'],[2,3]);
title("Velocidades angulares em função do tempo" ,'fontsize', 4);
xlabel("t(s)", 'fontsize', 3);
ylabel("Vel. angular(rad/s)", 'fontsize', 3);
legend(['d Theta 1 /dt','d Theta 2 /dt'],4);

//Velocidade linear
scf(6);
plot2d(t,v_corrigido,2);
title("Velocidade linear em função do tempo" ,'fontsize', 4);
xlabel("t(s)", 'fontsize', 3);
ylabel("Velocidade (m/s)", 'fontsize', 3);

///////////////////////////////// MODELO LINEARIZADO /////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////

//////////////////////////CÁLCULO DAS MATRIZES\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
function x_dot=calculasemt(y) //Função criada sem necessidade de colocar o tempo, como na parte não linear
    theta_1_dot = y(1);
    alpha = y(2);
    theta_2_dot = y(3);
    v = y(4);
    theta_2 = y(5);
    
     //vetor de estados na condição de equilíbrio
    
    M11 = -d*mr*sin(y(5))
    M12 = -y(4)*(mr + mc)*sin(y(2))
    M13 = d*mr*sin(y(5))
    M14 = (mr + mc)*cos(y(2))
    M15 = 0
    
    M21 = -mr*(b+c+d*cos(y(5)))
    M22 = y(4)*(mr + mc)*cos (y(2))
    M23 = d*mr*cos(y(5))
    M24 = (mr + mc)*sin(y(2))
    M25 = 0
    
    M31 = mr*(b + c)*(b + c + d*cos(y(5))) + Ic
    M32 = -y(4)*mr*(b+c)*cos(y(2))
    M33 = -d*mr*(b+c)*cos(y(5))
    M34 = -mr*(b+c)*sin(y(2))
    M35 = 0
    
    M41 = Ir + d*mr*(d + (b+c)*cos(y(5)))
    M42 = -y(4)*d*mr*cos(y(2)+ y(5))
    M43 = -mr*d^2 - Ir
    M44 = -d*mr*sin(y(2) + y(5))
    M45 = 0
    
    M51 = 0
    M52 = 0
    M53 = 0
    M54 = 0
    M55 = 1
    
    M = [M11,M12,M13,M14,M15;M21,M22,M23,M24,M25;M31,M32,M33,M34,M35;M41,M42,M43,M44,M45;M51,M52,M53,M54,M55] //Matriz de massas
    
    //ângulos de deriva:
    angF = atan(  (y(4)*sin(y(2))+a*y(1))  /  (y(4)*cos(y(2)))  )
    angT = atan(  (y(4)*sin(y(2))-b*y(1))  /  (y(4)*cos(y(2)))  )
    angE = atan( (y(4)*sin(y(2) + y(5)) - (b+c)*y(1)*cos(y(5)) - (d+e)*(y(1) - y(3))  ) / ( y(4)*cos( y(2) + y(5) ) + (b+c)*y(1)*sin(y(5)) )  )
    
    Fy_F = Fy (angF,Fz_F) //Frontal caminhão
    Fy_T = Fy (angT,Fz_T) //Traseiro caminhão
    Fy_E = Fy (angE,Fz_E) //Eixos reboque
    
    f1 = Fy_E * sin(y(5)) - (b+c)*y(1)^2*mr + (mc+mr)*y(4)*y(1)*sin(y(2)) - (y(1)-y(3))^2*d*mr*cos(y(5))
    f2 = Fy_T + Fy_F + Fy_E*cos(y(5)) - (mc+mr)*y(4)*y(1)*cos(y(2)) + ( y(1) - y(3) )^2*d*mr*sin(y(5))
    f3 = a*Fy_F - Fy_T*b - (b+c)*(  Fy_E*cos(y(5)) - y(4)*mr*cos(y(2))*y(1) + ( y(1)-y(3) )^2*d*mr*sin(y(5))  )
    f4 = d*( b*y(1)^2*mr*sin(y(5)) - Fy_E + c*y(1)^2*mr*sin(y(5)) + y(4)*y(1)*mr*cos(y(2)+y(5))) - Fy_E*e
    f5 = y(3)
    
    f = [f1;f2;f3;f4;f5]
    Minv = inv(M)
    x_dot = Minv*f //Matriz de massas linearizada
endfunction

//Calculando o ponto de equilíbrio

// Palpite inicial para as variáveis
y0 = [0; 0; 0; v_ini; 0];

// Resolver o sistema de equações
[yequilibrio, info] = fsolve(y0,calculasemt);

//Função para calcular a Matriz A
function J = calcularJacobianNumerico(yequilibrio, delta)
    J = zeros(5, 5);
    for i = 1:5
        ytemp = yequilibrio;
        ytemp(i) = ytemp(i) + delta; // Perturbar a i-ésima variável
        xdot_perturbado = calculasemt(ytemp);
        J(:, i) = (xdot_perturbado - calculasemt(yequilibrio)) / delta;
    end
endfunction

//Funções para calcular a matriz B

function x_dot=calculacomimpulso(imp) //Função igual a calculasemt, mas para calcular o jacobiano na matriz B
    y = yequilibrio;
    
     //vetor de estados na condição de equilíbrio
    
    M11 = -d*mr*sin(y(5))
    M12 = -y(4)*(mr + mc)*sin(y(2))
    M13 = d*mr*sin(y(5))
    M14 = (mr + mc)*cos(y(2))
    M15 = 0
    
    M21 = -mr*(b+c+d*cos(y(5)))
    M22 = y(4)*(mr + mc)*cos (y(2))
    M23 = d*mr*cos(y(5))
    M24 = (mr + mc)*sin(y(2))
    M25 = 0
    
    M31 = mr*(b + c)*(b + c + d*cos(y(5))) + Ic
    M32 = -y(4)*mr*(b+c)*cos(y(2))
    M33 = -d*mr*(b+c)*cos(y(5))
    M34 = -mr*(b+c)*sin(y(2))
    M35 = 0
    
    M41 = Ir + d*mr*(d + (b+c)*cos(y(5)))
    M42 = -y(4)*d*mr*cos(y(2)+ y(5))
    M43 = -mr*d^2 - Ir
    M44 = -d*mr*sin(y(2) + y(5))
    M45 = 0
    
    M51 = 0
    M52 = 0
    M53 = 0
    M54 = 0
    M55 = 1
    
    M = [M11,M12,M13,M14,M15;M21,M22,M23,M24,M25;M31,M32,M33,M34,M35;M41,M42,M43,M44,M45;M51,M52,M53,M54,M55] //Matriz de massas
    
    //ângulos de deriva:
    angF = atan(  (y(4)*sin(y(2))+a*y(1))  /  (y(4)*cos(y(2)))  )
    angT = atan(  (y(4)*sin(y(2))-b*y(1))  /  (y(4)*cos(y(2)))  )
    angE = atan( (y(4)*sin(y(2) + y(5)) - (b+c)*y(1)*cos(y(5)) - (d+e)*(y(1) - y(3))  ) / ( y(4)*cos( y(2) + y(5) ) + (b+c)*y(1)*sin(y(5)) )  )
    
    Fy_F = Fy (angF,Fz_F) //Frontal caminhão
    Fy_T = Fy (angT,Fz_T) //Traseiro caminhão
    Fy_E = Fy (angE,Fz_E) //Eixos reboque
    
    f1 = (Fy_E+imp) * sin(y(5)) - (b+c)*y(1)^2*mr + (mc+mr)*y(4)*y(1)*sin(y(2)) - (y(1)-y(3))^2*d*mr*cos(y(5))
    f2 = Fy_T + Fy_F + (Fy_E+imp)*cos(y(5)) - (mc+mr)*y(4)*y(1)*cos(y(2)) + ( y(1) - y(3) )^2*d*mr*sin(y(5))
    f3 = a*Fy_F - Fy_T*b - (b+c)*(  (Fy_E+imp)*cos(y(5)) - y(4)*mr*cos(y(2))*y(1) + ( y(1)-y(3) )^2*d*mr*sin(y(5))  )
    f4 = d*( b*y(1)^2*mr*sin(y(5)) - (Fy_E-imp) + c*y(1)^2*mr*sin(y(5)) + y(4)*y(1)*mr*cos(y(2)+y(5))) - (Fy_E+imp)*e
    f5 = y(3)
    
    f = [f1;f2;f3;f4;f5]
    Minv = inv(M)
    x_dot = Minv*f //Matriz de massas linearizada
endfunction

//Calculando a entrada de equilíbrio
u0 = 0; 

testeequilibrio = calculacomimpulso(u0)

if testeequilibrio==0 then
    uequilibrio = u0;
else
    abort;
end

function J = JacobianoMatrizB(delta)

    imp_equilibrio = uequilibrio;
    imp_perturbado = imp_equilibrio + delta;
    
    xdot_equilibrio = calculacomimpulso(imp_equilibrio);
    xdot_perturbado = calculacomimpulso(imp_perturbado);
    
    J = (xdot_perturbado - xdot_equilibrio) *(1/delta);
endfunction


////////////////Matrizes para modelo linearizado em torno do equilíbrio
delta = 1e-5
MatA = calcularJacobianNumerico(yequilibrio, delta)
MatB = JacobianoMatrizB(delta)
//MatC = eye(5,5)
MatC = [1,0,0,0,0;0,1,0,0,0;0,0,1,0,0;0,0,0,0,1];
//MatD = zeros (5,1)
MatD = zeros(4,1)

//////////////////////Resolver sistema linear///////////////////////////////////////

sys = syslin( 'c' ,MatA,MatB,MatC,MatD)

//Função de transferência (colunas são as entradas, linhas são saídas)

G = ss2tf( sys )
disp ("Matriz das FT")
disp (G)

//Equação característica
p = poly(MatA,"s")
disp("Equação característica")
disp(p)

//Polos
polos = spec(MatA)
disp("Polos do sistema")
disp(polos)

//Routh-Hurwitz
[r,num] = routh_t(p)
disp("Tabela de Routh-Hurwitz")
disp(r)

if num == 0 then
    disp("Sistema é estável")
else
    disp("Há mudanças de sinal nas entradas da primeira coluna, portanto, o sistema é instável", num)
end
//// polos /////
f0 = scf(7)
for i = 1: length(polos);
    s = polos(i)
plot(real(s),imag(s),'marker','o','markerFaceColor','red','markerEdgeColor','red')
end
set(gca(),'grid',[1 1]);
a=gca();
xlabel('Real(Si)','fontsize' , 2 , 'fontstyle' , 6);
ylabel('Imaginario(Si)', 'fontsize', 2 , 'fontsizer', 6);
title("Diagrama de polos do sistema",'fontsize', 3.5);
h0= legend(['Polos'],1)
xsave(TMPDIR+"/polos.png",gcf())


/////bode///////
freq_range = logspace(-1,1,1000)
scf(8)
bode(G(1,1),freq_range)
title("Diag. Bode da saída theta ponto 1 em relação a entrada Fy_Impulso",'fontsize', 3.5)

scf(9)
bode(G(2,1),freq_range)
title("Diag. Bode da saída alpha em relação a entrada Fy_Impulso",'fontsize', 3.5)

scf(10)
bode(G(4,1),freq_range)
title("Diag. Bode da saída theta2 em relação a entrada Fy_Impulso",'fontsize', 3.5)

/////////////////////////// RODAR A LINEARIZADA //////////////////////////////////////

//vetor de espaços do estado inicial
XL0 = X0 

//tempo:
Npl = 100 // numero de passos
tl = linspace(0,T,Npl); //vetor de tempo

//Entrada impulso
u = zeros (1,Npl)

//Definindo módulo do sinal de entrada
for i = 1 : Npl
    u(1,i) = impulso(tl(i))
end

//Simulando linearmente

[Y_lin,x] = csim (u , tl , sys , XL0 )

//Y_lin: linha 1: theta ponto 1 ; linha 2: alfa; linha 3 : theta 2


//plotando gráficos de saída linear

scf(11) 
plot2d(tl,Y_lin(1,:),2);
title("Velocidade angular do caminhão em função do tempo - Linear" ,'fontsize', 3.5);
xlabel("t(s)", 'fontsize', 3);
ylabel("Vel. angular(rad/s)", 'fontsize', 3);

//plotando gráficos de saída linear

scf(12) 
plot2d(tl,Y_lin(2,:)*180/%pi,2);
title("Ângulo do caminhão com a direção de sua velocidade em função do tempo - Linear" ,'fontsize', 3.5);
xlabel("t(s)", 'fontsize', 3);
ylabel("Pos. angular(º)", 'fontsize', 3);

//plotando gráficos de saída linear

scf(13) 
plot2d(tl,Y_lin(4,:)*180/%pi,2);
title("Ângulo do reboque em relação ao caminhão em função do tempo - Linear" ,'fontsize', 3.5);
xlabel("t(s)", 'fontsize', 3);
ylabel("Pos. angular(º)", 'fontsize', 3);


