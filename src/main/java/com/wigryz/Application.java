package com.wigryz;

import com.wigryz.algorithms.Function;
import com.wigryz.algorithms.Gauss;

public class Application {
    public static void main(String[] args) {
        Function function1 = x -> 5 * Math.pow(x[0], 2) + 3 * x[0] + 6;
        Function function2 = x -> 5 * Math.pow(x[0], 2) * Math.pow(x[1], 2) + 3 * x[0] * x[1] + 6;

        System.out.println(Gauss.gauss1D(function1, Gauss.INTEGRAL_SCHEME_1N));
        System.out.println(Gauss.gauss1D(function1, Gauss.INTEGRAL_SCHEME_2N));
        System.out.println(Gauss.gauss2D(function2, Gauss.INTEGRAL_SCHEME_1N));
        System.out.println(Gauss.gauss2D(function2, Gauss.INTEGRAL_SCHEME_2N));
    }
}
