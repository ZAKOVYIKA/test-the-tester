classdef AerodynamicParadox < handle                                
                                                                           
% Данный класс относится к моделированию т.н. "аэродинамического парадокса 
% спутника" — попадая в верхние слои атмосферы, космический аппарат, 
% испытывая торможение в разреженном газе, увеличивает при этом скорость 
% своего движения.
% https://kvant.mccme.ru/pdf/2011/04/Travin.pdf
% 
    
    properties 

        % Мировые константы                                                
        gravitiConstMultMassEarth = 4.0020e+14; % Гравитационная постоянная умноженная на массу Земли в СИ
        rEarth                    = 6400000; % Радиус Земли в метрах

        % Параметры спутника                                               
        squareOfMidel             = 1; % Площадь Миделя в отн.ед.
        heightOfOrbit             = 500000; % Высота орбиты спутника в метрах
        phiPolarAngleRadian       = 0; % Полярный угол положения спутника в рад
        phiPolarAngleSpeed        % Скорость изменения(производная) полярного угла рад/сек

        % Параметры окружающей среды                                       
        coeffOfDrag               = 2; % Коэффициент лобового сопротивления
        densityOfAtmosphere       = 5.*10.^(-8); % плотность атмосферы на орбите кг/м3
        heightOfOrbitRederence    = 40000; % Констатнта для высоты орбиты спутника в метрах

        % Константы задачи
        tauStepTime               = 0.1; % Шаг по сетке времени

    end
    
    methods
        function obj = AerodynamicParadox(name, value)              
                                                                           
            % Конструктор
            arguments (Repeating)                                          
                name
                value
            end

            % Вычисление производной полярного угла
            obj.phiPolarAngleSpeed = sqrt(obj.gravitiConstMultMassEarth / (obj.heightOfOrbit + obj.rEarth) ^ 3);         
            
            for valueIdx = 1:numel(value)                                         
                obj.(name{valueIdx}) = value{valueIdx};
            end
        end
        
        function [figureHandle] = plotVelocityTrajectory(obj)        
            % Рисунок

            [radius, angle, radiisTimeDerivative, anglesTimeDerivative] = obj.solve_diff_equation(10); 
            
            figureHandle = figure;
            tileHandle = tiledlayout(1, 2);

            nexttile;                                                                                                                                             
            plot(sqrt(anglesTimeDerivative.^2.*radius.^2 ...     
                                + radiisTimeDerivative.^2));      
            title('Скорость КА');
            xlabel('Время за 10 периодов КА');
            ylabel('Скорость КА в метрах в секунду');

            nexttile;
            polarplot( 0:0.01:2*pi, obj.rEarth.*ones(size(0:0.01:2*pi)), 'r');
            title('Положение КА');
           
            hold on;
            polarplot(angle, radius, 'b'); 
            legend('В первый период','В последний период');
        end

        function [radius, phi, radiusDerivative, phiDerivative] = solve_diff_equation(obj, periodsAmount)
            % Решение разностной схемы для движения спутника по орбите с трением. 

            timesteps = [0:obj.tauStepTime:periodsAmount.*2.*pi./obj.phiPolarAngleSpeed]';     

            radius           = zeros(numel(timesteps), 1);
            phi              = zeros(numel(timesteps), 1);
            radiusDerivative = zeros(numel(timesteps), 1);
            phiDerivative    = zeros(numel(timesteps), 1);
            
            % 0-ой узел
            radius(1)        = obj.heightOfOrbit + obj.rEarth;
            phi(1)           = obj.phiPolarAngleRadian;

            % 1-ый узел
            radius(2)        = obj.heightOfOrbit + obj.rEarth;
            phi(2)           = obj.phiPolarAngleRadian + obj.phiPolarAngleSpeed.*obj.tauStepTime;
            
            % Вычисление произведения постоянных в силе трения
            prefactorForce = obj.coeffOfDrag.*obj.squareOfMidel.*obj.densityOfAtmosphere./2;

            for timeIdx = 2:numel(timesteps)-1                                   
                % 1-ые производные и угол
                findDeriviate = @(val2, val1) (val2 - val1)./obj.tauStepTime;              
                findAngleBeta = @(radiusDerivative, phiDerivative, radius) atan(radiusDerivative./(phiDerivative.*radius));        
                
                % Сила трения 
                calcFrictionForce = @(radiusDerivative, phiDerivative, radius) prefactorForce.*(radiusDerivative.^2 + phiDerivative.^2.*radius.^2)...
                    .* exp(- (radius - obj.rEarth)./obj.heightOfOrbitRederence); 
                
                % Производные
                radiusDerivative(timeIdx) = findDeriviate(radius(timeIdx), radius(timeIdx-1));
                phiDerivative(timeIdx) = findDeriviate(phi(timeIdx), phi(timeIdx-1));
                    
                frictionForce = calcFrictionForce(radiusDerivative(timeIdx), phiDerivative(timeIdx), radius(timeIdx));
            
                angleBeta = findAngleBeta(radiusDerivative(timeIdx), phiDerivative(timeIdx), radius(timeIdx));
                
                % Координаты
                radius(timeIdx+1) = 2.*radius(timeIdx) - radius(timeIdx-1) + obj.tauStepTime.^2.*(radius(timeIdx).*phiDerivative(timeIdx).^2 ...
                    - obj.gravitiConstMultMassEarth./radius(timeIdx).^2 - frictionForce.*sin(angleBeta));
                phi(timeIdx+1) = 2.*phi(timeIdx) - phi(timeIdx-1) + obj.tauStepTime.^2./radius(timeIdx) .* (-frictionForce.*cos(angleBeta)...
                    - 2.*phiDerivative(timeIdx).*radiusDerivative(timeIdx));
            end
            
            % Производные последнего узла
            radiusDerivative(timeIdx+1) = findDeriviate(radius(timeIdx+1), radius(timeIdx));
            phiDerivative(timeIdx+1)    = findDeriviate(phi(timeIdx+1), phi(timeIdx));
            
            % Удалить 0-ой узел
            radius(1) = [];
            phi(1) = [];
            radiusDerivative(1) = [];
            phiDerivative(1) = [];

        end
    end

end












































