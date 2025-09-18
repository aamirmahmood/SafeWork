% 1. Pobranie pliku z koordynatami punktów x, y, z
pointsFile = 'gt_seq_09_inertial_coord.csv'; % Podaj nazwę pliku z punktami
pointsData = readtable(pointsFile);

% Wydobycie koordynatów punktów
pointsX = pointsData.x;
pointsY = pointsData.y;
pointsZ = pointsData.z;

% 2. Pobranie pliku z koordynatami stacji bazowych x, y, z
stationsFile = 'bs_coord.csv'; % Podaj nazwę pliku ze stacjami
stationsData = readtable(stationsFile);

% Wydobycie koordynatów stacji
stationsX = stationsData.x;
stationsY = stationsData.y;
stationsZ = stationsData.z;

% Liczba punktów i stacji
numPoints = height(pointsData);
numStations = height(stationsData);

% 3. Wyznaczenie 3 najbliższych stacji dla każdego punktu
nearestStations = zeros(numPoints, 3); % Macierz na numery najbliższych stacji

for i = 1:numPoints
    % Obliczanie odległości euklidesowych od punktu do każdej stacji
    distances = sqrt((stationsX - pointsX(i)).^2 + ...
                     (stationsY - pointsY(i)).^2 + ...
                     (stationsZ - pointsZ(i)).^2);
    
    % Sortowanie odległości i pobranie indeksów 3 najbliższych stacji
    [~, sortedIndices] = sort(distances);
    nearestStations(i, :) = sortedIndices(1:3);
end

% 4. Zapisanie wyników do nowego pliku .csv
outputTable = table(pointsX, pointsY, pointsZ, ...
                    nearestStations(:, 1), nearestStations(:, 2), nearestStations(:, 3), ...
                    'VariableNames', {'x', 'y', 'z', 'Station1', 'Station2', 'Station3'});

writetable(outputTable, 'new_input.csv');
disp('Wynik zapisano w pliku output.csv');
