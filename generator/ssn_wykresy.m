function ssn_wykresy()
    % ustawienia srodowiska
    clear all;
    clc;
    clf;
    
    % masa podstawki, pierwszego wahadla i drugiego wahadla
    m = [5, 2, 3];
    % dlugosc pierwszego i drugiego wahadla
    l = [4, 4];
    % przesuniecie wahadla wymuszone w celu poprawienia jego ustawienia
    u = 0;
    % stala grawitacji
    g = 9.81;
    % krok czasowy symulacji
    dt = 0.02;
    
    % ustawienia macierzy na potrzeby rozwiazywania ukladu rownan
    % nieliniowych oraz metody LQR
    [T0, G, Q, R, M_, N_] = settings(m, l, g);
   
    % poczatkowe wartosci: polozenia podstawki, kata pierwszego wahadla
    % i kata drugiego wahadla
    y1_0 = [0, 0.1, 0];
    % poczatkowe wartosci: predkosci podstawki, predkosci katowej 
    % pierwszego wahadla, predkosci katowej drugiego wahadla
    y2_0 = [0, 0, 0];
    % warunki poczatkowe dla metody rozwiazujacej rownania nieliniowe
    y0 = [y1_0; y2_0];
       
    % wektor zawierajacy kolejne punkty czasowe
    t = [];
    % wektor danych dotyczacych polozenia i predkosci podstawki
    cart_data = [];
    % wektor danych dotyczacych katow pierwszego i drugiego wahadla
    phi_data = [];
    
    % petla wykonujaca symulacje
    i = 0;
    while 1
        % funkcja okreslajaca uklad rownan nieliniowych do rozwiazania
        lambda = @(t, y) odefun(t, y, M_, N_, u, m, l);
        % rozwiazywanie ukladu rownan nieliniowych
        [t_, y] = ode45(lambda, [(i-1)*dt, i*dt], y0);   
        % pobieranie najbardziej aktualnego zestawu parametrow
        % (dla i-tej chwili czasowej)
        y = y(end,:);
        % zamiana kolejnosci parametrow na potrzeby funkcji lambda
        y0 = [y(1) y(3) y(5); y(2) y(4) y(6)];

        % aktualizacja wektora punktow czasowych
        t = [t i*dt];
        % aktualizacja wartosci polozenia i predkosci podstawki
        cart_data = [cart_data [y(1); y(2)]];
        % aktualizacja wartosci katow obu wahadel
        phi_data = [phi_data [y(3); y(5)]];
        
        % wspolrzedne punktow dla podstawki
        cart_x = [y(1)-1, y(1)+1];
        cart_y = [0, 0];
        % wspolrzedne punktow dla pierwszego wahadla
        pendulum1_x = [y(1), y(1)+l(1)*sin(y(3))];
        pendulum1_y = [0, l(1)*cos(y(3))];
        % wspolrzedne punktow dla drugiego wahadla
        pendulum2_x = [pendulum1_x(2), pendulum1_x(2)+l(2)*sin(y(5))];       
        pendulum2_y = [pendulum1_y(2), pendulum1_y(2)+l(2)*cos(y(5))];
        
        % wykres przedstawiajacy podstawke oraz dwa wahadla
        subplot(3,1,1);
        plot(cart_x, cart_y, 'b', pendulum1_x, pendulum1_y, 'r', pendulum2_x, pendulum2_y, 'g');
        xlim([-8 8]);
        ylim([-4 12]);
        % wykres wartoœci polozenia i predkosci podstawki
        subplot(3,1,2);
        plot(t, cart_data);
        legend('x', 'v');
        % wykres wartoœci katow pierwszego i drugiego wahadla
        subplot(3,1,3);
        plot(t, phi_data);
        legend('\phi_1', '\phi_2')
        
        % obliczanie korekcji polozenia podstawki za pomoca metody LQR
        K = lqr(A(T0, M_, N_), B(T0, M_, G), Q, R);
        X = [y(1); y(3); y(5); y(2); y(4); y(6)-y(4)];
        u = K*X;
        
        % pauza w celu rysowania wykresu
        pause(0.001);
        % zwiekszanie licznika iteracji
        i = i+1;
    end
end

% funkcja okreslajaca rozwiazywany uklad rownan nieliniowych
function dydt = odefun(t, y, M_, N_, u, m, l)
    y1 = [y(1), y(3), y(5)];
    y2 = [y(2), y(4), y(6)];    
    dydt = [y2'; F(y1, y2, M_, N_, u, m, l)];
end

% funkcja F = F(y, y', u) okreslajaca jedno z rozwiazywanych rownan
% nieliniowych
function F = F(y, yp, M_, N_, u, m, l)
    F = (M_.*D(y)) \ f(y, yp, N_, u, m, l);
end

% macierz M
function M = M(m, l)
    M = [m(1)+m(2)+m(3), l(1)*(m(2)/2+m(3)), m(3)/2*l(2); ...
         l(1)*(m(2)/2+m(3)), (l(1)^2)*(m(2)/3+m(3)), l(1)*l(2)*m(3)/2; ...
         l(2)*m(3)/2, l(1)*l(2)*m(3)/2, (l(2)^2)*m(3)/3];
end

% macierz D
function D = D(y)
    D = [1, cos(y(2)), cos(y(3)); ...
         cos(y(2)), 1, cos(y(2)-y(3)); ...
         cos(y(3)), cos(y(2)-y(3)), 1];
end

% funkcja f = f(y, y', u)
function f = f(y, yp, N_, u, m, l)
    v1 = [l(1)*(m(2)/2+m(3))*(yp(2)^2)*sin(y(2)) + m(3)/2*l(2)*(yp(3)^2)*sin(y(3)); ...
          -l(1)*l(2)*m(3)/2*(yp(3)^2)*sin(y(2)-y(3)); ...
          l(1)*l(2)*m(3)/2*(yp(2)^2)*sin(y(2)-y(3))];
    v2 = N_ * [0; sin(y(2)); sin(y(3))];
    v3 = [u; 0; 0];
    f = v1 + v2 + v3;
end

% macierz N 
function N = N(m, l, g)
    N = [0, 0, 0; ...
         0, g*(m(2)/2+m(3))*l(1), 0; ...
         0, 0, g*l(2)*m(3)/2];
end

% macierz A w metodzie LQR
function A = A(T0, M, N)
    A = [zeros(3), eye(3); ((T0 / M) * N) / T0, zeros(3)];
end

% macierz B w metodzie LQR
function B = B(T0, M, G)
    B = [zeros(3,1); (T0 / M) * G];
end

% definicje macierzy na potrzeby rozwiazywania ukladu rownan
% nieliniowych oraz metody LQR
function [T0, G, Q, R, M_, N_] = settings(m, l, g)
    T0 = [1, 0, 0; 0, 1, 0; 0, -1, 1];
    G = [1; 0; 0];
    Q = [1000, 0, 0; 0, 100, 0; 0, 0, 300];
    Q = [Q, zeros(3); zeros(3,6)];
    R = 1;
    M_ = M(m, l);
    N_ = N(m, l, g);
end