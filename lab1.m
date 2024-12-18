clear all
%Задание 1
disp('Задание 1');
% Создание матрицы x1 размером 11x11
x1 = rand(11, 11);
disp('Матрица x1:');
disp(x1);

% Выделение подматрицы x2 размером 5x5
x2 = x1(end-4:end, end-4:end);
disp('Матрица x2:');
disp(x2);
%Задание 2

disp('Задание 2');
% Создание матрицы x1 размером 5x5 с произвольными числами от 0 до 1
x1 = rand(10, 10);

% Создание матрицы x2 с четными строками матрицы x1
% Используем индексацию с шагом 2, начиная со второй строки
x2 = x1(2:2:end, :);

% Создание матрицы x3 с нечетными столбцами матрицы x1
% Используем индексацию с шагом 2, начиная с первого столбца
x3 = x1(:, 1:2:end);

% Выводим матрицы для проверки
disp('Матрица x1:');
disp(x1);
disp('Матрица x2 (четные строки x1):');
disp(x2);
disp('Матрица x3 (нечетные столбцы x1):');
disp(x3);

%Задание 3
disp('Задание 3');
% Создание вектора-строки x1 с диапазоном значений от -1 до 1
% и разбиением на 100 элементов
x1 = linspace(-1, 1, 100);

% Зануление элементов вектора x1, которые меньше 0.5
% Используем логическую индексацию для нахождения и замены таких элементов
x1(x1 < 0.5) = 0;

% Выводим измененный вектор x1
disp('Измененный вектор x1:');
disp(x1);

%Задание 4
disp('Задание 4');
% Создание квадратных матриц x1 и x2 размером 5x5, заполненных
% произвольными числами от 0 до 1
x1 = rand(5,5);
x2 = rand(5,5);

% Перемножение матриц по правилам матричного умножения
resultMatrixMultiplication = x1 * x2;

% Перемножение матриц поэлементно
resultElementWiseMultiplication = x1 .* x2;

% Вывод результатов
disp('Результат матричного умножения матриц x1 и x2:');
disp(resultMatrixMultiplication);
disp('Результат поэлементного умножения матриц x1 и x2:');
disp(resultElementWiseMultiplication);

%Задание 5
disp('Задание 5');
% Создание вектор-строки x1 и вектор-столбца x2
% с одинаковым количеством элементов (например, 10) случайных чисел от 0 до 1
x1 = rand(1,10); % Вектор-строка с 10 случайными числами
x2 = rand(10,1); % Вектор-столбец с 10 случайными числами

% Перемножение векторов по правилам матричного умножения
% Результатом будет скаляр (число), так как умножаем вектор-строку на вектор-столбец
resultMatrixMultiplication = x1 * x2;

% Перемножение векторов поэлементно не применимо в строгом смысле,
% так как они имеют разные ориентации (строка и столбец).
% Однако, мы можем транспонировать вектор-столбец в вектор-строку
% или наоборот для поэлементного умножения. В этом примере транспонируем x2.
resultElementWiseMultiplication = x1 .* x2;

% Вывод результатов
disp('Результат матричного умножения векторов x1 и x2:');
disp(resultMatrixMultiplication);
disp('Результат поэлементного умножения векторов x1 и x2:');
disp(resultElementWiseMultiplication);

%Задание 6
disp('Задание 6');
% Задаем размеры матриц
n = 10; % Произвольное целое число для размера n
k = 10; %Произвольное целое число для размера k
m = 10; % Произвольное целое число для размера m

% Создаем матрицы x1 и x2 с произвольными числами от 0 до 1
x1 = rand(n, k);
x2 = rand(k, m);

% Перемножение матриц по правилам матричного умножения
tic; % Начало измерения времени
x3 = x1 * x2;
toc; % Конец измерения времени и вывод результата
disp('Результат матричного умножения матриц x1 и x2:');

% Реализация алгоритма поэлементного нахождения элементов матрицы x3
x3_elementwise = zeros(n, m); % Инициализация матрицы нулями
tic; % Начало измерения времени
for i = 1:n
    for j = 1:m
        for r = 1:k
            x3_elementwise(i, j) = x3_elementwise(i, j) + x1(i, r) * x2(r, j);
        end
    end
end
toc; % Конец измерения времени и вывод результата
disp('Результат поэлементного нахождения элементов матрицы x3:');

%Задание 7
disp('Задание 7');
% Создание вектора-строки x1 с случайными числами от 0 до 1
x1 = rand(1, 10); % Допустим, размер вектора - 10 элементов

% Поэлементное нахождение элементов вектора x2 (циклический подход)
tic; % Начало измерения времени
x2_loop = zeros(1, length(x1)-1); % Инициализация x2 с размером на 1 меньше x1
for i = 1:length(x1)-1
    x2_loop(i) = x1(i+1) - x1(i);
end
loopTime = toc; % Конец измерения времени

% Векторизованный подход для нахождения элементов x2
tic; % Начало измерения времени
x2_vectorized = x1(2:end) - x1(1:end-1);
vectorizedTime = toc; % Конец измерения времени

% Вывод результатов
disp('Вектор x1:');
disp(x1);
disp('Разности элементов (цикл):');
disp(x2_loop);
disp('Разности элементов (векторизованный):');
disp(x2_vectorized);

% Сравнение времени выполнения с использованием fprintf
fprintf('Время выполнения (цикл): %f сек.\n', loopTime);
fprintf('Время выполнения (векторизованный): %f сек.\n', vectorizedTime);

%Задание 8
disp('Задание 8');
% Создание матрицы x1 с размером 5x5, заполненной случайными числами от 0 до 1
x1 = rand(5, 5);

% Создание матрицы x2 для хранения разности строк матрицы x1
% Разность берется между последовательными строками,
% поэтому размер x2 будет на одну строку меньше, чем x1
x2 = diff(x1, 1, 1); % Второй аргумент '1' указывает на вычитание вдоль строк

% Создание матрицы x3 для хранения разности столбцов матрицы x1
% Разность берется между последовательными столбцами,
% поэтому размер x3 будет на один столбец меньше, чем x1
x3 = diff(x1, 1, 2); % Третий аргумент '2' указывает на вычитание вдоль столбцов

% Вывод результатов
disp('Матрица x1:');
disp(x1);
disp('Матрица x2 (разность строк x1):');
disp(x2);
disp('Матрица x3 (разность столбцов x1):');
disp(x3);

%Задание 9
disp('Задание 9');
% Создание вектора-строки x1 с 10 случайными числами от 0 до 1
x1 = rand(1, 10);

% Инициализация вектора x2
% Длина x2 будет на 2 меньше, чем длина x1,
% так как для крайних элементов не существует трех соседних элементов
x2 = zeros(1, length(x1) - 2);

% Вычисление среднего значения от трех соседних элементов вектора x1
for i = 1:length(x2)
    x2(i) = sum(x1(i:i+2)) / 3;
end

% Вывод векторов x1 и x2
disp('Вектор x1:');
disp(x1);
disp('Вектор x2 (среднее из трех соседних элементов x1):');
disp(x2);

%Задание 10
disp('Задание 10');
% Создание матрицы x1 размером 5x5, заполненной случайными числами от 0 до 1
x1 = rand(5, 5);

% Создание матрицы x2, содержащей только значения из x1 меньшие 0.5
% Используем логическую индексацию для выбора таких элементов
x2 = x1(x1 < 0.5);

% Преобразование x2 в матрицу того же размера, что и x1, где элементы, не удовлетворяющие условию, заменены на NaN (не число)
x2_fullSize = x1; % Создаем копию x1
x2_fullSize(x1 >= 0.5) = NaN; % Заменяем элементы, большие или равные 0.5, на NaN

% Вывод результатов
disp('Матрица x1:');
disp(x1);
disp('Элементы x1 меньшие 0.5 (вектор x2):');
disp(x2);
disp('Матрица x2 с элементами меньше 0.5 и NaN для остальных:');
disp(x2_fullSize);

%Задание 11
disp('Задание 11');
% Создание матрицы x1 размером 5x5, заполненной случайными числами от 0 до 1
x1 = rand(5, 5);

% Инициализация матрицы x2 с теми же размерами, что и x1
x2 = zeros(size(x1));

% Заполнение матрицы x2 единицами на позициях, где значения x1 меньше 0.5
x2(x1 < 0.5) = 1;

% Вывод матриц x1 и x2
disp('Матрица x1:');
disp(x1);
disp('Матрица x2 (единицы на позициях значений x1 меньше 0.5):');
disp(x2);

%Задание 12
disp('Задание 12');
% Задаем размер матрицы
N = 5; % Пример размера матрицы 5x5

% Создаем векторы координат для осей x и y
% Линейное пространство от -1 до 1, количество точек соответствует размеру матрицы N
x = linspace(-1, 1, N);
y = linspace(-1, 1, N);

% Генерируем две матрицы координатных сеток
[X, Y] = meshgrid(x, y);

% Вычисляем расстояние от центра до каждой точки матрицы
% Центр в координатах x и y предполагается равным 0, поэтому просто используем X и Y
Distance = sqrt(X.^2 + Y.^2);

% Результат в матрице x1
x1 = Distance;

% Выводим результат
disp('Матрица x1 с расстояниями от центра до края:');
disp(x1);

%Задание 13
disp('Задание 13');
% Создание вектора-строки x1 с диапазоном значений от -1 до 1, разбитым на 1000 элементов
x1 = linspace(-1, 1, 1000);

% Создание вектора-строки x2, содержащего значения синуса от соответствующих элементов x1
x2 = sin(2*pi*x1);

% Создание вектора-строки x3, содержащего значения косинуса от соответствующих элементов x1
x3 = cos(2*pi*x1);

figure;

% Построение графика для синуса
plot(x1, x2, 'r'); % 'r' означает красный цвет линии
hold on; % Удержание текущего графика для добавления второй функции

% Добавление на график значения косинуса
plot(x1, x3, 'b'); % 'b' означает синий цвет линии
hold off;
% Добавление легенды для различия графиков
legend('sin(2\pi x)', 'cos(2\pi x)');

% Название графика
title('Графики функций синуса и косинуса');

% Метки осей
xlabel('x');
ylabel('y');

% Включение сетки для удобства восприятия
grid on;

%Задание 14
disp('Задание 14');
% Создание вектора-строки x1
x1 = linspace(-1, 1, 1000);

% Векторизованный подход
tic; % Начало измерения времени для векторизованного подхода
x2_vectorized = sin(x1); % Вычисление синуса для каждого элемента x1
vectorizedTime = toc; % Замер времени выполнения

% Поэлементный подход
x2_loop = zeros(1, length(x1)); % Инициализация вектора x2 для поэлементного подхода
tic; % Начало измерения времени для поэлементного подхода
for i = 1:length(x1)
    x2_loop(i) = sin(x1(i)); % Вычисление синуса для каждого элемента x1 поэлементно
end
loopTime = toc; % Замер времени выполнения

% Вывод результатов времени выполнения
fprintf('Время выполнения векторизованного подхода: %f сек.\n', vectorizedTime);
fprintf('Время выполнения поэлементного подхода: %f сек.\n', loopTime);

%Задание 15
disp('Задание 15');
% Создание вектора-строки x1 с 10000 случайными числами от 0 до 1
x1 = rand(1, 10000);

% Использование встроенной функции для нахождения суммы
tic; % Начало измерения времени
sumUsingBuiltIn = sum(x1);
timeUsingBuiltIn = toc; % Замер времени выполнения

% Реализация поэлементного алгоритма для нахождения суммы
tic; % Начало измерения времени
sumUsingLoop = 0; % Инициализация переменной для хранения суммы
for i = 1:length(x1)
    sumUsingLoop = sumUsingLoop + x1(i);
end
timeUsingLoop = toc; % Замер времени выполнения

% Вывод сумм и времени выполнения
fprintf('Сумма с использованием встроенной функции: %f\n', sumUsingBuiltIn);
fprintf('Время использования встроенной функции: %f сек.\n', timeUsingBuiltIn);

fprintf('Сумма с использованием поэлементного алгоритма: %f\n', sumUsingLoop);
fprintf('Время использования поэлементного алгоритма: %f сек.\n', timeUsingLoop);


%Задание 16
disp('Задание 16');
% Создание вектора-строки x1 с диапазоном значений от -1 до 1
x1 = linspace(-1, 1, 100);

% Создание вектора-строки x2 с диапазоном значений от -10 до 10
x2 = linspace(-10, 10, 100);

% Вычисление произведения значений синуса от x1 на значения косинуса от x2
x3 = sin(x1) .* cos(x2);
figure;
% Построение графика x3
plot(x3);
title('График произведения синуса x1 на косинус x2');
xlabel('Индекс элемента');
ylabel('Значение x3');
grid on; % Включение сетки для лучшей визуализации


%Задание 17
disp('Задание 17');
% Создание матрицы x1 размером 5x5 с комплексными числами
% Реальная и мнимая части чисел выбираются случайно от 0 до 1
x1 = rand(5) + 1i * rand(5)

% Вычисление модуля комплексных чисел в x1
x2 = abs(x1);

% Вычисление фазы комплексных чисел в x1
x3 = angle(x1);

% Вывод результатов
disp('Матрица x1 (комплексные числа):');
disp(x1);

disp('Матрица x2 (модуль чисел из x1):');
disp(x2);

disp('Матрица x3 (фаза чисел из x1):');
disp(x3);


%Задание 18
disp('Задание 18');
% Задание размера матрицы
N = 100; % Пример размера для квадратной матрицы

% Создание координатных сеток для осей x и y
[x, y] = meshgrid(linspace(-1, 1, N), linspace(-1, 1, N));

% Вычисление расстояния от центра до каждой точки в матрице
% Центр координат предполагается в (0,0)
x1 = sqrt(x.^2 + y.^2);

% Создание второй матрицы на основе значений exp(-x1.*x1/0.125)
x2 = exp(-x1 .* x1 / 0.125);

% Визуализация второй матрицы с помощью imagesc()
figure;
imagesc(x2);
colorbar; % Добавление шкалы цвета для интерпретации значений
title('Визуализация матрицы с exp(-x1.*x1/0.125)');
xlabel('X');
ylabel('Y');


%Задание 19
disp('Задание 19');
% Создание квадратной матрицы 10x10, заполненной нулями
x1 = zeros(10, 10);

% Генерация диапазона значений от 0 до 1, разбитого на 10 элементов
diagonalValues = linspace(0, 1, 10);

% Вставка значений на главную диагональ матрицы
for i = 1:10
    x1(i, i) = diagonalValues(i);
end

% Альтернативный и более эффективный способ без цикла:
% x1 = diag(linspace(0, 1, 10));

% Вывод матрицы
disp('Матрица с элементами от 0 до 1 на главной диагонали:');
disp(x1);




















































































