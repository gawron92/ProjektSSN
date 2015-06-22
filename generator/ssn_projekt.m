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
    c1 = 'Symulacja z uzyciem LQR (z wykresami)';
    c2 = 'Symulacja z uzyciem LQR (bez wykresow)';
    c3 = 'Wygeneruj dane uczace';
    c4 = 'Naucz siec neuronowa';
    c5 = 'Symulacja z uzyciem wytrenowanej sieci (z wykresami)';
    c6 = 'Symulacja z uzyciem wytrenowanej sieci (bez wykres�w)';
    c7 = 'Sprawdz jakosc wytrenowanej sieci (por�wnaj z LQR)';
    
    choice = listdlg('Name', 'Wybierz akcje', 'PromptString', 'Dostepne akcje:', 'SelectionMode', 'single', 'ListString', {c1, c2, c3, c4, c5, c6, c7}, 'ListSize', [400 100]);
    
    % warunek zakonczenia symulacji
    if ~isempty(choice) && choice ~= 3 && choice ~= 4
        % okno dialogowe do wprowadzania liczby krokow
        stepsInput = {'Maksymalna liczba krokow symulacji:'};
        steps = inputdlg(stepsInput, 'Liczba krokow symulacji', 1, {'1000'});
        hasNext = @(i) i < str2num(char(steps));
    else
        hasNext = @(i) i < 1000;
    end
    
    % usuwanie starych wynikow symulacji metoda LQR
    if ~isempty(choice) && choice ~= 1 && choice ~= 2
        if exist('lqr_data.txt', 'file')==2
            delete('lqr_data.txt');
        end
    end
    
    % usuwanie starych wynikow symulacji metoda SSN
    if ~isempty(choice) && choice ~= 5 && choice ~= 6
        if exist('neural_data.txt', 'file')==2
            delete('neural_data.txt');
        end
    end
    
    % wykonanie procedury w zaleznosci od akcji wybranej przez uzytkownika
    if ~isempty(choice)
        switch choice
            case 1
                [m l g y1_0 y2_0] = getParams();
                simulate(m, l, g, dt, y1_0, y2_0, hasNext, false, @forEachLqrAndPlot);
            case 2
                [m l g y1_0 y2_0] = getParams();
                simulate(m, l, g, dt, y1_0, y2_0, hasNext, true, @forEachLqrAndPlot);
            case 3
                [m l g y1_0 y2_0] = getParams();
                generateData(m, l, g, dt);
            case 4
                trainNetwork();
            case 5
                [m l g y1_0 y2_0] = getParams();
                simulate(m, l, g, dt, y1_0, y2_0, hasNext, false, @forEachNeural);
            case 6
                [m l g y1_0 y2_0] = getParams();
                simulate(m, l, g, dt, y1_0, y2_0, hasNext, true, @forEachNeural);
            case 7
                [m l g y1_0 y2_0] = getParams();
                simulate(m, l, g, dt, y1_0, y2_0, hasNext, true, @forEachTest);
        end

        % zamkniecie plikow jesli jakies zostaly otwarte
        fclose('all');
    end
end

function [m, l, g, y1_0, y2_0] = getParams()
    prompt = {'Masa podstawki', 'Masa czlonu #1', 'Masa czlonu #2', 'Dlugosc czlonu #1', 'Dlugosc czlonu #2', 'Przyspieszenie ziemskie'};
    answer = inputdlg(prompt, 'Parametry symulacji', 1, {'5', '2', '3', '4', '4', '9.81'});
    m = [str2num(char(answer(1))) str2num(char(answer(2))) str2num(char(answer(3)))];
    l = [str2num(char(answer(4))) str2num(char(answer(5)))];
    g = str2num(char(answer(6)));
    prompt = {'Polozenie podstawki', 'Wychylenie czlonu #1', 'Wychylenie czlonu #2', 'Predkosc podstawki', 'Predkosc katowa czlonu #1', 'Predkosc katowa czlonu #2'};
    answer = inputdlg(prompt, 'Warunki poczatkowe', 1, {'0', '0.1', '0', '0', '0', '0'});
    y1_0 = [str2num(char(answer(1))) str2num(char(answer(2))) str2num(char(answer(3)))];
    y2_0 = [str2num(char(answer(4))) str2num(char(answer(5))) str2num(char(answer(6)))];
end

% wytrenuj siec
function trainNetwork()
    % okno dialogowe do wprowadzania parametrow sieci neuronowej
    paramsInput = {'Liczba neuronow w warstwie ukrytej:','Liczba epok uczenia:','Akceptowalna wartosc bledu RMSE:'};
    params = inputdlg(paramsInput, 'Parametry sieci neuronowej oraz uczenia', 1, {'13','2000','1e-05'});
    
    % dostepme metody uczenia sieci
    m1 = 'trainlm';
    m2 = 'trainbfg'; 
    m3 = 'trainrp';
    m4 = 'traingd';
    
    % wybor metody uczenia sieci
    choice = listdlg('Name', 'Wybierz metode uczenia sieci', 'PromptString', 'Dostepne metody:', 'SelectionMode', 'single', 'ListString', {m1, m2, m3, m4}, 'ListSize', [300 100]);
    
    % ustawianie odpowiedniej metody
    learningMethod = m1;
    if ~isempty(choice)
        switch choice
            case 1
                learningMethod = m1;
            case 2
                learningMethod = m2;
            case 3
                learningMethod = m3;
            case 4
                learningMethod = m4;
        end
    end

    % minimalne i maksymalne wartosci
    PR = [-8 8; -1 1; -1 1; -1 1; -1 1; -1 1];
    
    % rozmiar warstw
    S = [6 str2num(char(params(1))) 1];
        
    % funkcje przenoszenia
    TF = {'purelin', 'tansig','purelin'};
    
    % algorytm uczenia sie
    BTF = learningMethod;
    
    % utworz siec
    net = newff(PR, S, TF, BTF);
    %net.trainParam.epochs = params(2);
    net.trainParam.epochs = str2num(char(params(2)));
    %net.trainParam.goal = params(3);
    net.trainParam.goal = str2num(char(params(3)));

    % wczytaj dane uczace
    data = importdata('generator_data.txt')';
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
    % definicja ile krokAlw czasowych uLLyć dla danych uczących
    hasNext = @(i) i < 200;
    
    % ile sumulacji przeprowadzić do wygenerowania danych uczących
    N = 20;
    
    % wyL�wietlenie paska postępu
    h = waitbar(0, 'Generowanie danych uczacych');
    
    for i = 1:N
        % sprawdzanie czy losowe parametry początkowe nie są za duLLe
        test = Inf;
        while test > 0.2
            y1_0 = randn(1, 3)/10;
            y2_0 = randn(1, 3)/10;
            test = max(abs([y1_0 y2_0]));
        end
        
        % wykonanie symulacji
        simulate(m, l, g, dt, y1_0, y2_0, hasNext, true, @forEachGenerator);
        
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
        msgbox('Generowanie danych zakonczone pomyslnie', 'Koniec');
        close(h);
    else
        msgbox('Generowanie przerwane na zyczenie uzytkownika', 'Koniec');
    end
end

function u = forEachTest(m, l, g, y, i, dt, T0, G, Q, R, M_, N_, isFast)
    persistent RMS;
    persistent U;
    persistent U2;
    persistent max;
    persistent max2;
    if isempty(RMS) 
        RMS = [];
        U = [];
        U2 = [];
        max = 10;
        max2 = 10;
    end
    
    u = forEachNeural(m, l, g, y, i, dt, T0, G, Q, R, M_, N_, true);
    u2 = forEachLqr(m, l, g, y, i, dt, T0, G, Q, R, M_, N_);
    
    
    
    RMS = [RMS rms(u-u2)];
    U = [U u];
    U2 = [U2 u2];
    
    maxCut = 50;
    if abs(u) > max & abs(u) < maxCut
        max = abs(u);
    end
    if abs(u2) > max & abs(u2) < maxCut
       max = abs(u2);
    end
    
    if abs(u) > max2
        max2 = abs(u);
    end
    if abs(u2) > max2
        max2 = abs(u2);
    end
     
    subplot(2, 2, 1);
    plot(U, U2, '+');
    axis square;
    xlim([-max max]);
    ylim([-max max]);
    xlabel('Wymuszenie u (SSN)');
    ylabel('Wymuszenie u (LQR)');
    title('Korelacja pomiedzy danymi [...]');
    hold on;
    plot([-max max], [-max max], 'g');
    hold off;
    
    subplot(2, 2, 2);
    plot(U, U2, '+');
    axis square;
    xlim([-max2 max2]);
    ylim([-max2 max2]);
    xlabel('Wymuszenie u (SSN)');
    ylabel('Wymuszenie u (LQR)');
    title('[...] generowanymi przez LQR i SSN');
    hold on;
    plot([-max2 max2], [-max2 max2], 'g');
    hold off;
    
    subplot(2, 1, 2);
    cut = 200;
    if i < cut
        x = 0:i;
        x = x*dt;
        plot(x, RMS);
    else
        x = (i-cut):(i);
        x = x*dt;
        plot(x, RMS(end-cut:end));
        xlim([x(1) x(end)]);
    end
    xlabel('Czat [s]');
    ylabel('RMS');
    title('Zmiana bledu sredniokwadratowego w czasie');
    
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
    
    drawnow;
end

% obliczanie korekcji polozenia podstawki za pomoca metody LQR
function u = forEachNeural(m, l, g, y, i, dt, T0, G, Q, R, M_, N_, isFast)
    persistent start;
    persistent net;
    persistent neural_data
    
    if isempty(neural_data) 
        neural_data = fopen('neural_data.txt', 'w+');
    end
    
    % przy pierwszym uruchomieniu
    if isempty(start)
        start = [1];
        
        % wczytaj siec z pliku
        load net;
    end
    
    % wylicz wymuszenie u
    u = sim(net, [y(1); y(3); y(5); y(2); y(4); y(6)]);
    
    % narysuj wykresy
    if ~isFast
        forEachPlot(m, l, g, y, i, dt, T0, G, Q, R, M_, N_);
    end
    
    % zapisywanie wynikow
    saved_results = [i*dt y(1) y(3) y(5) y(2) y(4) y(6) u];
    fprintf(neural_data, '%8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f\r\n', saved_results);
end

% obliczanie korekcji polozenia podstawki za pomoca metody LQR
function u = forEachLqr(m, l, g, y, i, dt, T0, G, Q, R, M_, N_)
    K = lqr(A(T0, M_, N_), B(T0, M_, G), Q, R);
    X = [y(1); y(3); y(5); y(2); y(4); y(6)-y(4)];
    u = K*X;
end

% zapis danych do pliku
function u = forEachGenerator(m, l, g, y, i, dt, T0, G, Q, R, M_, N_, isFast)
    persistent generator_data
    
    if isempty(generator_data) 
        generator_data = fopen('generator_data.txt', 'w+');
    end
    
    % wylicz wymuszenie u
    u = forEachLqr(m, l, g, y, i, dt, T0, G, Q, R, M_, N_);
    
    % zapisywanie wynikow
    saved_results = [i*dt y(1) y(3) y(5) y(2) y(4) y(6) u];
    fprintf(generator_data, '%8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f\r\n', saved_results);
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
    title('Podglad wachadla');
    
    % wykres warto?ci polozenia i predkosci podstawki
    subplot('Position', [0.05, 0.1, 0.4, 0.4]);
    plot(t, cart_data);
    legend('x [m]', 'v [m/s]', 'Location', 'southoutside');
    title('Podstawka wachadla');
    xlabel('Czas [s]');
    
    % wykres warto?ci katow pierwszego i drugiego wahadla
    subplot('Position', [0.55, 0.1, 0.4, 0.4]);
    plot(t, phi_data);
    legend('\phi1 [rad]', '\phi2 [rad]', 'Location', 'southoutside');
    title('Wychylenie czlonow');
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
function u = forEachLqrAndPlot(m, l, g, y, i, dt, T0, G, Q, R, M_, N_, isFast)
    persistent lqr_data
    
    if isempty(lqr_data) 
        lqr_data = fopen('lqr_data.txt', 'w+');
    end

    % narysuj wykresy
    if ~isFast
        forEachPlot(m, l, g, y, i, dt, T0, G, Q, R, M_, N_);
    end
    
    % wylicz nowa wartosc wymuszenia u
    u = forEachLqr(m, l, g, y, i, dt, T0, G, Q, R, M_, N_);
    
    % zapisywanie wynikow
    saved_results = [i*dt y(1) y(3) y(5) y(2) y(4) y(6) u];
    fprintf(lqr_data, '%8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f\r\n', saved_results);
end

% wykonaj symulacje dla zadanych parametrow
function simulate(m, l, g, dt, y1_0, y2_0, hasNext, isFast, forEach)
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
        u = forEach(m, l, g, y, i, dt, T0, G, Q, R, M_, N_, isFast);
        
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