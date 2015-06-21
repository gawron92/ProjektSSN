function ssn_projekt()
    % ustawienia srodowiska
    clear all;
    clc;
    
    % masa podstawki, pierwszego wahadla i drugiego wahadla
    m = [5, 2, 3];
    % dlugosc pierwszego i drugiego wahadla
    l = [4, 4];
    % stala grawitacji
    g = 9.81;
    % krok czasowy symulacji
    dt = 0.02;
    
    % poczatkowe wartosci: polozenia podstawki, kata pierwszego wahadla
    % i kata drugiego wahadla
    y1_0 = [0, 0.1, 0];
    % poczatkowe wartosci: predkosci podstawki, predkosci katowej 
    % pierwszego wahadla, predkosci katowej drugiego wahadla
    y2_0 = [0, 0, 0];
    
    % wybierz dostepne akcje
    c1 = 'Symulacja z użyciem LQR';
    c2 = 'Wygeneruj dane uczące';
    c3 = 'Naucz sieć neuronową';
    c4 = 'Symulacja z użuciem wytrenowanej sieci';
    
    choice = listdlg('Name', 'Wybierz akcję', 'PromptString', 'Dostępne akcję:', 'SelectionMode', 'single', 'ListString', {c1, c2, c3, c4}, 'ListSize', [300 100]);
    
    % domyslny warunek zakonczenia symulacji
    hasNext = @(i) i < 1000;
    
    % wykonanie procedury w zaleznosci od akcji wybranej przez uzytkownika
    if ~isempty(choice)
        switch choice
            case 1
                simulate(m, l, g, dt, y1_0, y2_0, hasNext, @forEachLqrAndPlot);
            case 2
                generateData(m, l, g, dt);
            case 3
                trainNetwork();
            case 4
                simulate(m, l, g, dt, y1_0, y2_0, hasNext, @forEachNeural);
        end

        % zamkniecie plikow jesli jakies zostaly otwarte
        fclose('all');
    end
end

% wytrenuj siec
function trainNetwork()
    % minimalne i maksymalne wartosci
    PR = [-8 8; -1 1; -1 1; -1 1; -1 1; -1 1];
    
    % rozmiar warsts
    S = [6 13 1];
    
    % funkcje przenoszenia
    TF = {'purelin', 'tansig','purelin'};
    
    % algorytm uczenia sie
    BTF = 'trainlm';
    
    % utworz siec
    net = newff(PR, S, TF, BTF);
    net.trainParam.epochs = 2000;

    % wczytaj dane uczace
    data = importdata('results.txt')';
    X = data(2:7, :);
    Y = data(8, :);
    
    % usun wartosci dla krotych algorytm jest rozbiezny
    cut = abs(Y)<500;
    X = X(:, cut);
    Y = Y(:, cut);
    
    % wytrenuj siec
    net = train(net, X, Y);
    
    % zapisz siec
    save net;
end

function generateData(m, l, g, dt)
    % definicja ile kroków czasowych użyć dla danych uczących
    hasNext = @(i) i < 200;
    
    % ile sumulacji przeprowadzić do wygenerowania danych uczących
    N = 20;
    
    % wyświetlenie paska postępu
    h = waitbar(0, 'Generowanie danych uczących');
    
    for i = 1:N
        % sprawdzanie czy losowe parametry początkowe nie są za duże
        test = Inf;
        while test > 0.2
            y1_0 = randn(1, 3)/10;
            y2_0 = randn(1, 3)/10;
            test = max(abs([y1_0 y2_0]));
        end
        
        % wykonanie symulacji
        simulate(m, l, g, dt, y1_0, y2_0, hasNext, @forEachGenerator);
        
        % jezeli uzytkownik nie kliknal "x" to zaktualizuj pasek posteku
        if ishandle(h)
            waitbar(i/N, h);
        % w przeciwnym razie zakoncz symulacje
        else
            break;
        end
    end
    
    % zamknij pasek postepu
    if ishandle(h)
        msgbox('Generowanie danych zakończone pomyślnie', 'Koniec');
        close(h);
    else
        msgbox('Generowanie przerwane na życzenie użytkownika', 'Koniec');
    end
end

% obliczanie korekcji polozenia podstawki za pomoca metody LQR
function u = forEachNeural(m, l, g, y, i, dt, T0, G, Q, R, M_, N_)
    persistent start;
    persistent net;
    
    % przy pierwszym uruchomieniu
    if isempty(start)
        start = [1];
        
        % wczytaj siec z pliku
        load net;
    end
    
    % wylicz wymuszenie u
    u = sim(net, [y(1); y(3); y(5); y(2); y(4); y(6)]);
    
    % narysuj wykresy
    forEachPlot(m, l, g, y, i, dt, T0, G, Q, R, M_, N_);
end

% obliczanie korekcji polozenia podstawki za pomoca metody LQR
function u = forEachLqr(m, l, g, y, i, dt, T0, G, Q, R, M_, N_)
    K = lqr(A(T0, M_, N_), B(T0, M_, G), Q, R);
    X = [y(1); y(3); y(5); y(2); y(4); y(6)-y(4)];
    u = K*X;
end

% zapis danych do pliku
function u = forEachGenerator(m, l, g, y, i, dt, T0, G, Q, R, M_, N_)
    persistent results
    
    if isempty(results) 
        results = fopen('results.txt', 'w+');
    end
    
    % wylicz wymuszenie u
    u = forEachLqr(m, l, g, y, i, dt, T0, G, Q, R, M_, N_);
    
    % zapisywanie wynikow
    saved_results = [i*dt y(1) y(3) y(5) y(2) y(4) y(6) u];
    fprintf(results, '%8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f\r\n', saved_results);
end

% rysowanie wykresow
function u = forEachPlot(m, l, g, y, i, dt, T0, G, Q, R, M_, N_)
    % wektor zawierajacy kolejne punkty czasowe
    persistent t;
    % wektor danych dotyczacych polozenia i predkosci podstawki
    persistent cart_data
    % wektor danych dotyczacych katow pierwszego i drugiego wahadla
    persistent phi_data;

    % inicjalizacja zmiennych
    if isempty(cart_data)
        t = [];
        cart_data = [];
        phi_data = [];
    end

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
    subplot('Position', [0.05, 0.62, 0.9, 0.3]);
    plot(cart_x, cart_y, 'b', pendulum1_x, pendulum1_y, 'r', pendulum2_x, pendulum2_y, 'g');
    xlim([-8 8]);
    ylim([-4 12]);
    title('Podgląd wachadła');
    
    % wykres warto?ci polozenia i predkosci podstawki
    subplot('Position', [0.05, 0.1, 0.4, 0.4]);
    plot(t, cart_data);
    legend('x [m]', 'v [m/s]', 'Location', 'southoutside');
    title('Podstawka wachadła');
    xlabel('Czas [s]');
    
    % wykres warto?ci katow pierwszego i drugiego wahadla
    subplot('Position', [0.55, 0.1, 0.4, 0.4]);
    plot(t, phi_data);
    legend('\phi1 [rad]', '\phi2 [rad]', 'Location', 'southoutside');
    title('Wychylenie członów');
    xlabel('Czas [s]');
    
    % przycisk zatrzymania symulacji
    uicontrol('Position',[5 5 70 25], ...
               'String', 'STOP', ...
               'Callback', @stopCb);
    
    % akcja po wcisnieciu przycisku zatrzymania symulacji
    function stopCb(popup, callbackdata)
        global stopSimulation
        stopSimulation = true;
        msgbox('Zatrzymano symulacje');
    end

    % rysuj wykres
    drawnow;
    
    u = Inf;
end

% rysowanie wykresow + liczenie wymuszenia z uzyciem lqr
function u = forEachLqrAndPlot(m, l, g, y, i, dt, T0, G, Q, R, M_, N_)
    % narysuj wykresy
    forEachPlot(m, l, g, y, i, dt, T0, G, Q, R, M_, N_);
    
    % wylicz nowa wartosc wymuszenia u
    u = forEachLqr(m, l, g, y, i, dt, T0, G, Q, R, M_, N_);
end

% wykonaj symulacje dla zadanych parametrow
function simulate(m, l, g, dt, y1_0, y2_0, hasNext, forEach)
    % zmienna globalna do zatrzymania symulacji
    global stopSimulation
    stopSimulation = false;
    
    % przesuniecie wahadla wymuszone w celu poprawienia jego ustawienia
    u = 0;
    
    % ustawienia macierzy na potrzeby rozwiazywania ukladu rownan
    % nieliniowych oraz metody LQR
    [T0, G, Q, R, M_, N_] = settings(m, l, g);
   
    % warunki poczatkowe dla metody rozwiazujacej rownania nieliniowe
    y0 = [y1_0; y2_0];
    
    % petla wykonujaca symulacje
    i = 0;
    while hasNext(i)
        % funkcja okreslajaca uklad rownan nieliniowych do rozwiazania
        lambda = @(t, y) odefun(t, y, M_, N_, u, m, l);
        
        % rozwiazywanie ukladu rownan nieliniowych
        [t_, y] = ode45(lambda, [(i-1)*dt, i*dt], y0);  
        
        % pobieranie najbardziej aktualnego zestawu parametrow
        % (dla i-tej chwili czasowej)
        y = y(end,:);
        
        % zamiana kolejnosci parametrow na potrzeby funkcji lambda
        % oraz kolejnej iteracji
        y0 = [y(1) y(3) y(5); y(2) y(4) y(6)];

        % nowa wartosc wymuszenia u
        u = forEach(m, l, g, y, i, dt, T0, G, Q, R, M_, N_);
        
        % zatrzymaj sumacje
        if stopSimulation
            break;
        end
        
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
