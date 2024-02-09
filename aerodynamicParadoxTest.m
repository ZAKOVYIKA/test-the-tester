classdef aerodynamicParadoxTest < matlab.unittest.TestCase
    
    properties
        %
    end
    
    methods (TestClassTeardown)
        function teardown(~)
            clearvars;
            close all;
        end
    end

    methods (Test)
        function testClassCreation(testCase)
            % Проверка возможности создания класса
            instance = aerodynamicParadox;
            testCase.verifyClass(instance, 'aerodynamicParadox', 'Объект не соответствует классу');

        end
        
        function testFigureReturn(testCase)
            % Проверка что функция возращает объект figure с дочерними
            % элементами
            instance = aerodynamicParadox;
            figureInstance = instance.plotVelocityAndTrajectory(1);
            testCase.verifyClass(figureInstance, 'matlab.ui.Figure', 'Объект результата не соответствует классу');
            testCase.verifyClass(figureInstance.Children, 'matlab.graphics.layout.TiledChartLayout');
            testCase.verifyNotEmpty(figureInstance.Children.Children);

        end

        function testAttributesChangePosibility(testCase)
            % проверка возможности изменения значений атрибутов класса объекта
            instance = aerodynamicParadox;
            dummy = instance.satelliteMass;
            instance.satelliteMass = dummy + 10;
            testCase.verifyNotEqual(instance.satelliteMass, dummy)

        end


    end
end

