% 1. Wczytanie pliku .csv z koordynatami x, y, z
filename = 'gt_seq_09_inertial_coord.csv'; % Zmień na nazwę Twojego pliku
data = readtable(filename);

% Wyciągnięcie kolumn x i y
x = data.x;
y = data.y;

% Dodanie stałej wartości z = 30 dla wszystkich punktów
z = 30 * ones(height(data), 1);

% 2. Określenie granicznych wymiarów płaszczyzny x, y
minX = min(x);
maxX = max(x);
minY = min(y);
maxY = max(y);

% 3. Rozmieszczenie 12 równomiernych punktów wewnątrz płaszczyzny
numPoints = 12;
xPoints = linspace(minX, maxX, sqrt(numPoints));
yPoints = linspace(minY, maxY, sqrt(numPoints));

[XGrid, YGrid] = meshgrid(xPoints, yPoints);
xResult = reshape(XGrid, [], 1);
yResult = reshape(YGrid, [], 1);

% Dodanie współrzędnej z = 30 dla nowych punktów
zResult = 30 * ones(size(xResult));

% 4. Zapisanie koordynatów x, y, z tych punktów do nowego pliku .csv
outputData = table(xResult, yResult, zResult, ...
                   'VariableNames', {'x', 'y', 'z'});
writetable(outputData, 'bs_coord.csv');
disp('Koordynaty (x, y, z) zostały zapisane w pliku output.csv');
