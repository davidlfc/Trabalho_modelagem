//Daniel Marins Silva D´Oliveira Martins : 12554571
//David Lopes F. de Carvalho :             12610305
//João Pedro Kawachi Chaves :              12745188
//Leonardo Fortunato de Carvalho :         9837910


clear; // Limpeza de variáveis
close(winsid()) // Fecha as janelas gráficas abertas anteriormente

//vetor de estados do estado inicial
theta_p_1_ini = 0.2   //rad/s
alfa_ini = 0    //rad
theta_p_2_ini = 0.2       //rad/s
v_ini =  30     //m/s
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

dist =  -d+1.25    //distância da carga ao centro de massa do reboque
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
    mi_n = (a1*Fz_kN + a2)/1000 //Coef. atrito nominal
    
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
    
    f1 = Fy_E * sin(y(5)) - (b+c)*y(1)^2*mr + (mc+mr)*y(4)*y(1)*sin(y(2)) - (y(1)-y(3))^2*d*mr*cos(y(5))
    f2 = Fy_T + Fy_F + Fy_E*cos(y(5)) - (mc+mr)*y(4)*y(1)*cos(y(2)) + ( y(1) - y(3) )^2*d*mr*sin(y(5))
    f3 = a*Fy_F - Fy_T*b - (b+c)*(  Fy_E*cos(y(5)) - y(4)*mr*cos(y(2))*y(1) + ( y(1)-y(3) )^2*d*mr*sin(y(5))  )
    f4 = d*( b*y(1)^2*mr*sin(y(5)) - Fy_E + c*y(1)^2*mr*sin(y(5)) + y(4)*y(1)*mr*cos(y(2)+y(5))) - Fy_E*e
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
////////////////////////////////////////////////////////////////////////

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
scf(1); 
plot2d(t,[alfa',theta_2'],[2,3]);
title("Ângulos em função do tempo" ,'fontsize', 5);
xlabel("t(s)", 'fontsize', 3);
ylabel("Pos. angular(rad)", 'fontsize', 3);
legend(['Alfa','Theta 2'],4);

//Velocidades angulares
scf(2);
plot2d(t,[theta_p_1',theta_p_2'],[2,3]);
title("Velocidades angulares em função do tempo" ,'fontsize', 5);
xlabel("t(s)", 'fontsize', 3);
ylabel("Vel. angular(rad)", 'fontsize', 3);
legend(['d Theta 1 /dt','d Theta 2 /dt'],4);

//Velocidade linear
scf(3);
plot2d(t,v,2);
title("Velocidade linear em função do tempo" ,'fontsize', 5);
xlabel("t(s)", 'fontsize', 3);
ylabel("Velocidade (m/s)", 'fontsize', 3);

///////////////////////////////// MODELO LINEARIZADO /////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////

//////////////////////////CÁLCULO DAS MATRIZES\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
function Matriz_A = M_A()
//Matriz A calculada em função para evitar armazenar as variáveis do código da mesma, que não possuem utilidade
    
    //Variáveis no equilíbrio estável
    theta_p_1_e = 0
    alfa_e = 0
    theta_p_2_e = 0
    v_e = v_ini //v equilíbrio não é zero
    theta_2_e = 0
    
    
    if dist < 0 then //dist < 0 -> carga na frente ; dist >= 0 -> carga no fundo
        FE = -168000
    else
        FE = -250000
    end
    
    
    Xe= [ theta_p_1_e ; alfa_e ; theta_p_2_e ; v_e ; theta_2_e ] //vetor de estados na condição de equilíbrio
    
    J11 = -2*(Xe(1)-Xe(3))*(b+c)*mr + (mc+mr)*Xe(4)*sin(Xe(2)) - 2*Xe(1)*d*mr*cos(Xe(5)) 
    J12 = (mc+mr)*Xe(4)*Xe(1)*cos(Xe(2))
    J13 = 2*(Xe(1)-Xe(3))*d*mr*cos(Xe(5))
    J14 = (mr+mc)*Xe(1)*sin(Xe(2))
    J15 = ( Xe(1) - Xe(3) )^2*d*mr*sin(Xe(5)) + FE*cos(Xe(5))
    
    J21 = -(mc+mr)* Xe(4)*cos(Xe(2)) - 2*(Xe(1)-Xe(3))*d*mr*cos(Xe(5))
    J22 = (mc+mr)*Xe(4)*Xe(1)*sin(Xe(2))
    J23 = -2*(Xe(1)-Xe(3))*d*mr*sin(Xe(5))
    J24 = -(mc+mr)*Xe(1)*cos(Xe(2))
    J25 = ( Xe(1) - Xe(3) )^2 * d*mr*cos(Xe(5)) + FE*(-sin(Xe(5)))
    
    J31 = -(b+c)*(  -Xe(4)*mr*cos(Xe(2)) + 2*(Xe(1)-Xe(3))*d*mr*sin(Xe(5))  )
    J32 = -(b+c)*( Xe(4)*mr*sin(Xe(2))*Xe(1) )
    J33 = -(b+c)*( -2*(Xe(1)-Xe(3))*d*mr*sin(Xe(5)) )
    J34 = -(b+c)*( -mr*cos(Xe(2))*Xe(1) )
    J35 = -(b+c)*( ( Xe(1) - Xe(3) )^2*mr*d*cos(Xe(5))  + FE*(-sin(Xe(5))) )
    
    J41 = d*( 2*(b+c)*Xe(1)*mr*sin(Xe(5)) + Xe(4)*mr*cos( Xe(2) + Xe(5) ) )
    J42 = d*Xe(4)*Xe(1)*mr*(-sin( Xe(2) + Xe(5) ))
    J43 = 0
    J44 = d*Xe(1)*mr*cos( Xe(2) + Xe(5) )
    J45 = d*( (b+c)*Xe(1)^2*mr*cos(Xe(5)) + Xe(4)*Xe(1)*mr*(-sin(Xe(2) + Xe(5))) )
    
    J51 = 0
    J52 = 0
    J53 = 1
    J54 = 0
    J55 = 0
    
    J = [J11,J12,J13,J14,J15;J21,J22,J23,J24,J25;J31,J32,J33,J34,J35;J41,J42,J43,J44,J45;J51,J52,J53,J54,J55] //Matriz de massas linearizada
    
    Matriz_A = J
    
endfunction

////////////////Matrizes para modelo linearizado em torno do equilíbrio
MatA = M_A()
MatB = [0,0,0;1,1,0;a,-b,0;0,0,-(d+e);0,0,0]
MatC = [1,0,0,0,0;0,1,0,0,0;0,0,0,0,1]
MatD = zeros (3,3)


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
f0 = scf(4)
for i = 1: length(polos);
    s = polos(i)
plot(real(s),imag(s),'marker','o','markerFaceColor','red','markerEdgeColor','red')
end
set(gca(),'grid',[1 1]);
a=gca();
xlabel('Real(Si)','fontsize' , 2 , 'fontstyle' , 6);
ylabel('Imaginario(Si)', 'fontsize', 2 , 'fontsizer', 6);
title("Diagrama de polos do sistema");
h0= legend(['Polos'],2)
xsave(TMPDIR+"/polos.png",gcf())


///////////////////////////////////////////////////////////////


/////bode///////
scf(5)
bode(G(1,1))
title("Diag. Bode da saída Theta ponto 1 em relação a entrada Fy_F")

scf(6)
bode(G(1,2))
title("Diag. Bode da saída Theta ponto 1 em relação a entrada Fy_T")

//Bode de G(1,3) não sairá pela natureza da entrada Fy_E

scf(7)
bode(G(2,1))
title("Diag. Bode da saída alfa em relação a entrada Fy_F")


scf(8)
bode(G(2,2))
title("Diag. Bode da saída alfa em relação a entrada Fy_T")

//Bode de G(1,3) não sairá pela natureza da entrada Fy_E


scf(9)
bode(G(3,1))
title("Diag. Bode da saída Theta 2 em relação a entrada Fy_F")


scf(10)
bode(G(3,2))
title("Diag. Bode da saída Theta 2 em relação a entrada Fy_T")

//Bode de G(1,3) não sairá pela natureza da entrada Fy_E

////////////////////////////////////////////////////////////////////


/////////////////////////// RODAR A LINEARIZADA //////////////////////////////////////

//vetor de estados do estado inicial
theta_p_1_ini_l = 0  //rad/s
alfa_ini_l = 0    //rad
theta_p_2_ini_l = 0       //rad/s
v_ini_l =  30     //m/s
theta_2_ini_l = 0   //rad
//Estes valores também dependerão de cada situação inicial diferente
//Portanto, estes valores podem ter pequenas variações entre simulações
/////////////////////////////////////

//vetor de espaços do estado inicial
XL0 = [theta_p_1_ini_l; alfa_ini_l; theta_p_2_ini_l; v_ini_l ; theta_2_ini_l] 

//tempo:
Npl = 1000 // numero de passos
tl = linspace(0,T,Npl); //vetor de tempo

//Entrada senoidal(natureza da Fy_E dentro da matriz A)
u = zeros (3,Npl)

//Definindo módulo dos sinais de entrada
//(valores calculados pela aproximação linear de pequenos ângulos nos gráficos de força dos pneus)
if dist < 0 then //dist < 0 -> carga na frente ; dist >= 0 -> carga no fundo
    FF = -232000
    FT = -185000
    FE = -168000
else
    FF = -220000
    FT = -196000
    FE = -250000
end
f = 50 //Hz, frequência
for i = 2 : Npl
    u(1,i) = FF * sin(2*%pi*f/tl(i))
    u(2,i) = FT * sin(2*%pi*f/tl(i))
    u(3,i) = FE * sin(2*%pi*f/tl(i))

end

//Simulando linearmente

[Y_lin,x] = csim (u , tl , sys , XL0 )

//Y_lin: linha 1: theta ponto 1 ; linha 2: alfa; linha 3 : theta 2

thetap1lin = Y_lin(1,:)

alfalin = Y_lin(2,:)

theta2lin = Y_lin(3,:)
//plotando gráficos de saída linear

scf(11) 
plot2d(tl,thetap1lin,2);
title("Velocidade angular do caminhão em função do tempo" ,'fontsize', 5);
xlabel("t(s)", 'fontsize', 3);
ylabel("Vel. angular(rad/s)", 'fontsize', 3);

//plotando gráficos de saída linear

scf(12) 
plot2d(tl,alfalin,2);
title("Ângulo do caminhão em função do tempo" ,'fontsize', 5);
xlabel("t(s)", 'fontsize', 3);
ylabel("Ângulo(rad)", 'fontsize', 3);

//plotando gráficos de saída linear

scf(13) 
plot2d(tl,theta2lin,2);
title("Ângulo do caminhão com o reboque em função do tempo" ,'fontsize', 5);
xlabel("t(s)", 'fontsize', 3);
ylabel("Ângulo(rad)", 'fontsize', 3);


