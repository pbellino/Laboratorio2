function [n, c, t] = cinetica_dir(rho,Q,param_num,param_cons)

% Calcula la cin�tica directa dando como dato la reactividad (rho) y el
% valor de la fuente de neutrones (Q).
%
% ENTRADA:
%         rho       Reactividad expresada en d�lares en funci�n del tiempo
%         Q         Fuente de neutrones en funci�n del tiempo. Si es
%                   constante se puede especificar como 'Q=cte'.
%       param_num   Un cell con los par�metros para los datos num�ricos:
%                   param_num={dt,tf,n_inicial}
%                   donde    dt: intervalo de discretizaci�n
%                            tf: tiempo final de la simulaci�n
%                     n_inicial: valor inicial de la densidad neutr�nica
%       param_cons  Un cell con los par�metros de las constantes cin�ticas:
%                   param_num={reactor,retardados,fotoneutrones,#fotoneutr}
%                   donde      reactor:    'RA3', 'RA1', 'CNA2', etc
%                           retardados: constantes de retardados ('Tuttle')
%                        fotoneutrones: constantes de fotoneutrones ('Keepin')
%                           #fotoneutr: cantidad de grupos a utilizar (los
%                                       primeros que m�s r�pido decaen)
% SALIDA:
%         n         Densidad neutr�nica
%         c         Concentraci�n de los precursores
%         t         Vector temporal
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Lectura de las constantes

reactor = param_cons{1};
nom_ret = param_cons{2};
nom_fot = param_cons{3};
n_fot   = param_cons{4};

[b, lam, beta_eff, Lstar, efot] = lee_constantes(reactor,nom_ret,nom_fot,n_fot);

dt    = param_num{1};
tf    = param_num{2};
ninic = param_num{3};

%% Definici�n num�rica

%dt   = 1;               % Intervalo temporal de integraci�n
%tf   = 1*24*3600;        % Tiempo total de la integraci�n
Kt   = floor(tf/dt);     % Cantidad de intervalos de la integraci�n
nd   = length(lam)+1;   % Cantidad de variables dle sistem (n_grupos + 1)

%ninic         = 4.1e13;    %Flujo del RA-3 operando a 8.5MW [1/cm^2 s]

t   = (0:Kt-1).*dt;

%% Defino la reactividad utilizada
% t_cor = 100;        % Tiempo en donde se produce el SCRAM
% rho = -15*ones(1,Kt);rho(t<=t_cor)=0;
% rho = 0.001.*sin(t.*2*pi*0.001);

% [rho_xe, t_xe] = func_xenon(60,0,tf,10,beta_eff,ninic);
% plot(t_xe/3600,rho_xe)

%% Construyo el vector de fuente
% Su tama�o ser� la cantidad de variables nd (a pesar de que s�lo importa
% la primera dimensi�n (associada a n(t)), y luego la dependencia temporal a
% trav�s de Kt
Qm  = zeros(nd,Kt);
% Si quiero definier la dependencia tempora, conviene hacerlo de forma
% vectorial:
if length(Q)==1
    Qvec = ones(1,Kt).*Q;
else
    Qvec = Q(:);
end
% Y lo meto en la matriz, s�lo en la parte de n(t)
Qm(1,:) = Qvec;

%% Resoluci�n de la cin�tica directa
% Busco resolver un sistema dy(t)/dt = A y(t) + Q(t)
% con y = [n(t) c_1(t) ... c_15(t)]
%
% Construyo la matr�z A con dimensi�n 16x16
u    = [lam' ; diag(-lam)];
frow = [-1 ; b]./Lstar;
A    = horzcat(frow,u);
% El valor A(1,1) es (rho-1)/Lstar. La condici�n inicial de reactor cr�tico
% impone que A(1,1) = -1/Lstar tal como se escribi�
clear u frow
%---------------------------------------------------------

% Condiciones iniciales

y         = zeros(nd,Kt);
% Condici�n inicial n(0) (definido al comienzo)
y(1,1)    = ninic;
% Condici�n inicial para los precursores (asumiendo reactor estacionario)
y(2:nd,1) = (y(1,1).* b./(lam.*Lstar))';

id = eye(size(A));              % Matriz identidad
for j=2:Kt
    A(1,1) = (rho(j)-1)/Lstar;     % El �nico elemento que cambia con el tiempo
    y(:,j) = (id-dt.*A)\(y(:,j-1)+dt.*Qm(:,j));  % M�todo de Euler impl�cito
end
clear id;

n = y(1,:)';
c = y(2:end,:)';

% % % % Gr�ficos
% % % figure

% % % 
% % % plot(t/3600,n,'k');
% % % semilogy(t/3600,n,'k');
% % % xlabel('Tiempo [h]');
% % % ylabel('Flujo [1/cm^2.s]');
% % % grid on

end