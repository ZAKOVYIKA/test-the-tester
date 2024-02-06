classdef AerodynamicParadoxNegativeTest < matlab.unittest.TestCase
    methods (Test)
         function testRadiusEarth(testCase)
             actSolution = AerodynamicParadox;
             actRadius = actSolution.solve_diff_equation(0);
             actRadius = actRadius(end);
             startRadius = 6.9e+6;
             testCase.verifyLessThan(actRadius, startRadius)
         end
    end
end