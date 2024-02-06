classdef AerodynamicParadoxpropertiesTest < matlab.unittest.TestCase

    properties
        TestFigure
    end

    methods (Test)
         function testRadiusEarth(testCase)
             actSolution = AerodynamicParadox;
             expSolution = 6400000;
             testCase.verifyEqual(actSolution.rEarth, expSolution)
         end

         function testGravitiConstMultMassEarth(testCase)
             actSolution = AerodynamicParadox;
             expSolution = 4.0020e+14;
             testCase.verifyEqual(actSolution.gravitiConstMultMassEarth, expSolution)
         end

         function createFigure(testCase)
            actSolution = AerodynamicParadox;
            testCase.TestFigure = actSolution.plotVelocityTrajectory;
            testCase.verifyNotEmpty(testCase.TestFigure)
        end
    end

    methods (TestMethodTeardown)
        function closeFigure(testCase)
            close(testCase.TestFigure)
        end
    end
end