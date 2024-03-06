close all;
clear all;
% Параметры задачи
lambda = 0.532; % Длина волны в метрах
n = 1; % Показатель преломления среды
k = 2*pi/lambda; % Волновое число


for NA_values = [0.65, 0.8, 0.95]

  alpha = asin(NA_values)

  T = @(theta) cos(theta).^0.5;

  r = 0;
  z = linspace(-2 * lambda, 2 * lambda, 1000);
  psi = 0;

  calculate_Ex_1 = @(theta, r, z) 1 .* T(theta) * sin(theta) .* (1 - cos(theta)) .* ...
  exp(1i * k * z * cos(theta)) .* besselj(2, k * r * sin(theta));

  calculate_Ex_2 = @(theta, r, z) 1 .* T(theta) .* sin(theta) .* (1 + cos(theta)) * ...
  exp(1i * k * z * cos(theta)) .* besselj(0, k * r * sin(theta));

  calculate_Ey = @(theta, r, z) 1 .* T(theta) .* sin(theta) .* (1 - cos(theta)) * ...
  exp(1i * k * z * cos(theta)) .* besselj(2, k * r * sin(theta));

  calculate_Ez = @(theta, r, z) 1 .* T(theta) * sin(theta)^2 .* ...
  exp(1i * k * z * cos(theta)) .* besselj(1, k * r * sin(theta));

  Ex = -1i * cos(2 * psi) * integral(@(theta) calculate_Ex_1(theta, r, z), 0, alpha, 'ArrayValued', true) - ...
  1i * integral(@(theta) calculate_Ex_2(theta, r, z), 0, alpha, 'ArrayValued', true)

  Ey = -1i * sin(2*psi) * integral(@(theta) calculate_Ey(theta, r, z), 0, alpha, 'ArrayValued', true)

  Ez = -2 * cos(psi) * integral(@(theta) calculate_Ez(theta, r, z), 0, alpha, 'ArrayValued', true)

  I = abs(Ex).^2 + abs(Ey).^2 + abs(Ez).^2;

  % Визуализация результатов
  figure;
  plot(z, I);
  xlabel('z (м)');
  ylabel('Интенсивность');
  title(['Распределение интенсивности для NA = ', num2str(NA_values)]);

end
