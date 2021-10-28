package com.wigryz;

import com.wigryz.algorithms.Function;
import com.wigryz.algorithms.Gauss;
import com.wigryz.algorithms.IntegralScheme;
import com.wigryz.structures.Element4x2D;

import java.util.Arrays;

public class Application {

    public static void main(String[] args) {
//        Function function1 = x -> 5 * Math.pow(x[0], 2) + 3 * x[0] + 6;
//        Function function2 = x -> 5 * Math.pow(x[0], 2) * Math.pow(x[1], 2) + 3 * x[0] * x[1] + 6;
//
//        System.out.println(Gauss.gauss1D(function1, IntegralScheme.INTEGRAL_SCHEME_1N));
//        System.out.println(Gauss.gauss1D(function1, IntegralScheme.INTEGRAL_SCHEME_2N));
//        System.out.println(Gauss.gauss2D(function2, IntegralScheme.INTEGRAL_SCHEME_1N));
//        System.out.println(Gauss.gauss2D(function2, IntegralScheme.INTEGRAL_SCHEME_2N));

        Element4x2D element4x2D = new Element4x2D(IntegralScheme.INTEGRAL_SCHEME_1N);
        System.out.println(element4x2D);
//        System.out.println(Arrays.deepToString(element4x2D.getEtaArray()).replace("], ", "]\n"));
//        System.out.println();
//        System.out.println(Arrays.deepToString(element4x2D.getKsiArray()).replace("], ", "]\n"));
    }
}
