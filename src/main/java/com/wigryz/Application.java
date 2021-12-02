package com.wigryz;

import com.wigryz.algorithms.Algorithms;
import com.wigryz.algorithms.IntegrationScheme;
import com.wigryz.structures.Element4x2D;
import com.wigryz.structures.Grid;
import com.wigryz.structures.Side;

import java.util.Arrays;

public class Application {

    public static void main(String[] args) {
//        Function function1 = x -> 5 * Math.pow(x[0], 2) + 3 * x[0] + 6;
//        Function function2 = x -> 5 * Math.pow(x[0], 2) * Math.pow(x[1], 2) + 3 * x[0] * x[1] + 6;
//
//        System.out.println(Algorithms.gauss1D(function1, IntegralScheme.INTEGRAL_SCHEME_1N));
//        System.out.println(Algorithms.gauss1D(function1, IntegralScheme.INTEGRAL_SCHEME_2N));
//        System.out.println(Algorithms.gauss2D(function2, IntegralScheme.INTEGRAL_SCHEME_1N));
//        System.out.println(Algorithms.gauss2D(function2, IntegralScheme.INTEGRAL_SCHEME_2N));

        Element4x2D element4x2D = new Element4x2D(IntegrationScheme.INTEGRATION_SCHEME_1N);
//
//        for (int i = 0; i < 4; i++) {
//            System.out.println("\nSide: " + i + "\n");
//            double[][] result = Algorithms.calculateHBC(0.0125, 25, (short) i, element4x2D);
//            System.out.println(Arrays.deepToString(result).replace("], ", "]\n"));
//        }
//
//        Side sideLeft = new Side(IntegrationScheme.INTEGRATION_SCHEME_1N, Side.LEFT);
//        Side sideBottom = new Side(IntegrationScheme.INTEGRATION_SCHEME_1N, Side.BOTTOM);
//        Side sideRight = new Side(IntegrationScheme.INTEGRATION_SCHEME_1N, Side.RIGHT);
//        Side sideTop = new Side(IntegrationScheme.INTEGRATION_SCHEME_1N, Side.TOP);

        Grid grid = new Grid(0.1, 0.1, 4,  4);
        grid.calculate(element4x2D);
        grid.agregate();
        grid.calculateT();
        grid.getHG();
    }
}
