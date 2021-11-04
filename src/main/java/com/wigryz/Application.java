package com.wigryz;

import com.wigryz.algorithms.IntegralScheme;
import com.wigryz.structures.Element4x2D;
import com.wigryz.structures.Grid;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;

public class Application {

    public static void main(String[] args) {
//        Function function1 = x -> 5 * Math.pow(x[0], 2) + 3 * x[0] + 6;
//        Function function2 = x -> 5 * Math.pow(x[0], 2) * Math.pow(x[1], 2) + 3 * x[0] * x[1] + 6;
//
//        System.out.println(Algorithms.gauss1D(function1, IntegralScheme.INTEGRAL_SCHEME_1N));
//        System.out.println(Algorithms.gauss1D(function1, IntegralScheme.INTEGRAL_SCHEME_2N));
//        System.out.println(Algorithms.gauss2D(function2, IntegralScheme.INTEGRAL_SCHEME_1N));
//        System.out.println(Algorithms.gauss2D(function2, IntegralScheme.INTEGRAL_SCHEME_2N));

        Element4x2D element4x2D = new Element4x2D(IntegralScheme.INTEGRAL_SCHEME_1N);
//        System.out.println(element4x2D);

        Grid grid = new Grid(0.025, 0.025, 2, 2);
        grid.calculate(element4x2D);
    }
}
