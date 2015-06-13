# ProjektSSN
Sterowanie podwójnym odwróconym wahadłem z użyciem sztucznych sieci neuronowych

Materiały teoretyczne
============
Algorytm stabilizacji pozycji wahadła został opracowany w oparciu o poniższe materiały:
- http://www3.math.tu-berlin.de/Vorlesungen/SS12/Kontrolltheorie/matlab/inverted_pendulum.pdf
Głównie strona 3 - postać macierzy M(y), oraz f(y, y', u). Z drugiej macierzy usunięto drugi i trzeci
człon oznaczające odpowiednio hamowanie wahadła (np. w wyniku działania siły tarcia) oraz zewnętrzne wymuszenia.
Dwa ostatnie równania wykorzystano do stworzenia postaci układu równań nieliniowych, rozwiązywanego przez odpowiednią funkcję Matlaba
- http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.78.6105&rep=rep1&type=pdf
Wykorzystano postacie macierzy określone wzorami 3, 4 i 5 (na stronie 3) oraz niżej wymienione wartości współczynników w celu dostosowania
macierzy z pierwszego źródła w taki sposób, aby środek punktu materialnego był w środku każdej z dwóch części wahadła (a nie na końcu).
Dodatkowo wyżej wymienione macierze posłużyły do rozkładu macierzy M(y) oraz f(y, y', u) na iloczyny macierzy
- https://books.google.pl/books?id=rnUiAgAAQBAJ&pg=PA22&lpg=PA22&dq=state+space+equation+double+inverted+pendulum&source=bl&ots=EKB8xdgzjy&sig=RUL0ZShfNqiyh7tNfTzbXFGyR9M&hl=pl&sa=X&ved=0CDQQ6AEwAjgKahUKEwjynuqOrovGAhULPhQKHbQMALg#v=onepage&q=state%20space%20equation%20double%20inverted%20pendulum&f=true
Postaciemacierzy A i B wykorzystywane w metodzie LQR (strona 24) oraz stosowanie implementacji metody LQR w Matlabie (strona 25)
- http://www.ijarcsse.com/docs/papers/february2012/volume_2_issue_2/V2I2066A.pdf
Inne, mniej istotne informacje (porównanie z materiałami na powyższych stronach)
- http://eprints.lancs.ac.uk/20405/1/download.pdf
Inne, mniej istotne informacje (porównanie z materiałami na powyższych stronach)

Generator danych dla sieci
============
Generator (plik ssn_generator.m) umożliwia wygenerowanie zbioru danych uczących i testowych dla sieci neuronowej.
Wygenerowany plik zawiera następujące dane (w kolejnych kolumnach):
- wartość kroku czasowego
- pozycję podstawki na której znajduje się wahadło
- wartość kąta dla pierwszego wahadła
- wartość kąta dla drugiego wahadła
- wartość prędkości dla podstawki
- prędkość kątową pierwszego wahadła
- prędkość kątową drugiego wahadła
- wartość przesunięcia, którego użyto do zmiany pozycji podstawki w celu ustabilizowania wahadła

Wykres stanu wahadła i jego parametrów
============
Plik ssn_wykresy.m zawiera program umożliwiający obserwację zmian stanu wahadła (wizualizację jego stanu) oraz dwa
kolejne wykresy, z których pierwszy pokazuje zmiany pozycji i prędkości podstawki, natomiast drugi wyświetla
informacje o zmienie wartości kątów dla obu części wahadła
