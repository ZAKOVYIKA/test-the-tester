classdef AerodynamicParadoxResultTest < matlab.unittest.TestCase
    methods (Test)
         function testRadiusEarth(testCase)
             actSolution = AerodynamicParadox;
             actRadius = actSolution.solve_diff_equation(10);
             actRadius = actRadius(end);
             startRadius = 6.9e+6;
             testCase.verifyLessThan(actRadius, startRadius)
         end
         function testSpeedSat(testCase)
             actSolution = AerodynamicParadox;
             [actRadius, actAngle, actRadiisTimeDerivative, actAnglesTimeDerivative] = actSolution.solve_diff_equation(10);
             speedStart = sqrt(actAnglesTimeDerivative(1).^2.*actRadius(1).^2 + actRadiisTimeDerivative(1).^2);
             speedFinish = sqrt(actAnglesTimeDerivative(end).^2.*actRadius(end).^2 + actRadiisTimeDerivative(end).^2);
             testCase.verifyLessThan(speedStart, speedFinish);
         end
    end
end