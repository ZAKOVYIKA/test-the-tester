classdef aerodynamicParadox < handle                                
                                                                           
% Данный класс относится к моделированию т.н. "аэродинамического парадокса 
% спутника" — попадая в верхние слои атмосферы, космический аппарат, 
% испытывая торможение в разреженном газе, увеличивает при этом скорость 
% своего движения.
% https://kvant.mccme.ru/pdf/2011/04/Travin.pdf
% 
    
    properties 
        
        % Мировые константы                                                
        G = 6.67e-11; % м3/(кг*с2) 
        earthMass = 5.97e24; % кг 
        earthRadius = 6371e3; % м

        % Параметры окружающей среды                                       
        Cx = 2; 
        p0 = 5e-8; % кг/м3
        H = 40000; % м?

        % Параметры спутника                                               
        satelliteArea = 1; % отн. м2
        satelliteMass = 100; % отн. кг

        startHeight = 500e3; % м

        startPhi = 0; % рад
        startRadius % м
        
        startPhiDer  % рад/сек 
        startOrbitDer = 0; % м/сек 
        
        % Константы задачи
        dT = 0.1; % с

    end
    
    methods
        function obj = aerodynamicParadox(name, value) 
            % Конструктор           
                                                                           
            arguments (Repeating)                                          
                name
                value
            end     
            
            for i = 1:numel(value)                                         
                obj.(name{i}) = value{i};
            end
            
            obj.startRadius = obj.earthRadius + obj.startHeight;
            obj.startPhiDer = sqrt(obj.G * obj.earthMass / (obj.startRadius)^3);

        end
        
        function [figureHandle] = plotVelocityAndTrajectory(obj, periodAmount) 
            % Рисунок     
            
            [radiusArray, phiArray, radiusDerArray, phiDerArray, linearSpeed, timesteps] = ...
                obj.solveDiffEquation(periodAmount); 
            
            figureHandle = figure;
            tileHandle = tiledlayout(1, 2);

            ax1 = nexttile;                                                
            ax2 = nexttile;                                                
            
                                                                           
            h1 = plot(ax1, timesteps, linearSpeed);      
            
            axes(ax2)
            h2 = polarplot( [0:0.01:2.*pi], ...                             
                           obj.earthRadius.*ones(1, numel([0:0.01:2.*pi])), 'r');  

            ax2 = gca;
            hold(ax2, 'on')
            h3 = polarplot(ax2, phiArray, radiusArray, 'b');
            
        end
    end

        %%
    methods (Access = private)
        function [radiusArray, phiArray, radiusDerArray, phiDerArray, linearSpeed, timesteps] = ...
                solveDiffEquation(obj, periodAmount)
            % Решение разностной схемы для движения спутника по орбите с трением. 

            dummy = periodAmount * 2 * pi / obj.startPhiDer;
            timesteps = (0:obj.dT:dummy)';     

            radiusArray = zeros(numel(timesteps), 1);
            phiArray = zeros(numel(timesteps), 1);
            radiusDerArray = zeros(numel(timesteps), 1);
            phiDerArray = zeros(numel(timesteps), 1);
            linearSpeed = zeros(numel(timesteps), 1);
            
            % 0-ой узел
            radiusArray(1) = obj.startRadius;
            phiArray(1) = obj.startPhi;
            radiusDerArray(1) = obj.startOrbitDer;
            phiDerArray(1) = obj.startPhiDer;
            linearSpeed(1) = sqrt(radiusDerArray(1)^2 +...
                radiusArray(1)^2 * phiDerArray(1)^2);

%             % 1-ый узел
%             radiusArray(2) = obj.Hic + obj.Rearth;
%             phiArray(2) = obj.Fic + obj.dFdtic.*obj.tau;

            for idx = 2:numel(timesteps) 
                % расчет в idx-1 точке
                % расчет плотности воздуха
                p = obj.p0 * exp(-(radiusArray(idx-1) - obj.earthRadius) / obj.H );
                % сила трения
                Ffr = obj.Cx * obj.satelliteArea * (p * linearSpeed(idx-1)^2)/2;
                % угол между векторами скоростей, рад
                beta = atan(radiusDerArray(idx-1) / (phiDerArray(idx-1) * radiusArray(idx-1)));
                % расчет второй производной по радиусу
                radiusSecDer = radiusArray(idx-1) * phiDerArray(idx-1)^2 ...
                    - obj.G * obj.earthMass / radiusArray(idx-1)^2 ...
                    - Ffr * sin(beta) / obj.satelliteMass;
                % вторая производная угла положения
                phiSecDer = -((Ffr * cos(beta) / obj.satelliteMass) ...
                    - (2 * phiDerArray(idx-1) * radiusDerArray(idx-1))) ...
                    / radiusArray(idx-1);

                radiusArray(idx) = radiusArray(idx-1) + obj.dT * radiusDerArray(idx-1) ...
                    + radiusSecDer * obj.dT^2 / 2;

                if (radiusArray(idx) < obj.earthRadius)
                    radiusArray(idx:end) = [];
                    phiArray(idx:end) = [];
                    radiusDerArray(idx:end) = [];
                    phiDerArray(idx:end) = [];
                    linearSpeed(idx:end) = [];
                    timesteps(idx:end) = [];
                    break
                end
                phiArray(idx) = phiArray(idx-1) + obj.dT * phiDerArray(idx-1) ...
                    + phiSecDer * obj.dT^2 / 2;
                radiusDerArray(idx) = radiusDerArray(idx-1) + radiusSecDer * obj.dT;
                phiDerArray(idx) = phiDerArray(idx-1) + phiSecDer * obj.dT;
                linearSpeed(idx) = sqrt(radiusDerArray(idx)^2 +...
                    radiusArray(idx)^2 * phiDerArray(idx)^2);
 
            end

        end
    end

end
