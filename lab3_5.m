close all;
clear all;

% Параметры
clarity = 15;
lambda = 0.532 * clarity;
k = (2 * pi) / lambda;
screen_size = 4;
R = 2.5 * lambda;
z_values = [0.5*lambda, lambda, 5*lambda, 10*lambda, 50*lambda];
C_values = -(1i*k*z_values) / (2*pi);
U = cell(1, length(z_values));

% Вычисление результирующих полей
tic
for z_idx = 1:length(z_values)
    U{z_idx} = zeros(2*screen_size*clarity, 2*screen_size*clarity);
    for i = 1:(2*screen_size*clarity + 1)
        for j = 1:(2*screen_size*clarity + 1)
            X = i - (screen_size*clarity + 1);
            Y = j - (screen_size*clarity + 1);
            func = @(x,y) exp(1i*2*pi*(z_values(z_idx)/lambda)) * ...
                (exp(1i*k*sqrt((X-x).^2+(Y-y).^2+z_values(z_idx).^2)) ./ ...
                ((X-x).^2+(Y-y).^2+z_values(z_idx).^2)) .* ...
                (1+(1i./(k.*sqrt((X-x).^2+(Y-y).^2+z_values(z_idx).^2))));
            U{z_idx}(i,j) = integral2(func,-R,R, @(x)-sqrt(R^2 - x.^2), @(x)sqrt(R^2 - x.^2));
        end
    end
    U{z_idx} = C_values(z_idx) .* U{z_idx};
end
toc

function plotIntensity(U, screen_size, titleStr, figureNumber)
    mid_row = round(size(U, 1) / 2);
    intensity = U(mid_row, :) .* conj(U(mid_row, :));
    x_coords = linspace(-screen_size, screen_size, length(intensity));
    figure;
    plot(x_coords, real(intensity));
    xlabel('x');
    ylabel('Intensity');
    title(titleStr);

    saveas(gcf, sprintf('figure_%d.png', figureNumber));
end


for i = 1:length(U)
    plotIntensity(U{i}, screen_size, sprintf('Intensity profile along x-axis (Figure %d)', i), i);
end


lambda = 0.532;
k = (2 * pi) / lambda;
D = 5 * lambda;

z10 = 10 * lambda;
z50 = 50 * lambda;

r_prime = linspace(-screen_size, screen_size, 100);

I10 = (besselj(1, k * D * r_prime / (2 * z10)) ./ (k * D * r_prime / (2 * z10))).^2;
I50 = (besselj(1, k * D * r_prime / (2 * z50)) ./ (k * D * r_prime / (2 * z50))).^2;

I10(r_prime == 0) = 1;
I50(r_prime == 0) = 1;

function plotNormalizedIntensity(U, I, screen_size, z, figure_number)
    mid_row = round(size(U, 1) / 2);

    intensity = U(mid_row, :) .* conj(U(mid_row, :));

    r_prime = linspace(0, screen_size, length(I));

    original_x = linspace(0, screen_size, length(intensity));

    interpolated_intensity = interp1(original_x, intensity, r_prime, 'linear', 'extrap');

    normalized_I = I / max(I);
    normalized_intensity = interpolated_intensity / max(interpolated_intensity);

    figure;
    plot(r_prime, normalized_intensity, 'b', 'DisplayName', 'RS-I');
    hold on;
    plot(r_prime, normalized_I, 'g', 'DisplayName', 'Airy');

    xlabel('r, \mum');
    ylabel('Intensity');
    title(sprintf('Intensity Distribution', z));
    legend show;
    grid on;

    axis([0, screen_size, 0, 1]);

    saveas(gcf, sprintf('figure_%d.png', figure_number));
end


plotNormalizedIntensity(U{4}, I10, 3, z10, 6);

plotNormalizedIntensity(U{5}, I50, 15, z50, 7);

