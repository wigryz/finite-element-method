package com.wigryz.algorithms;

import java.util.List;

public class Gauss {

    public static IntegralScheme INTEGRAL_SCHEME_1N =
        new IntegralScheme(
            List.of(0, 1),
            List.of(-1d/Math.sqrt(3), 1d/Math.sqrt(3)),
            List.of(1d, 1d));

    public static IntegralScheme INTEGRAL_SCHEME_2N =
        new IntegralScheme(
            List.of(0, 1, 2),
            List.of(-Math.sqrt(3)/Math.sqrt(5), 0d, Math.sqrt(3)/Math.sqrt(5)),
            List.of(5d/9d, 8d/9d, 5d/9d));

    public static double gauss1D(Function fun, IntegralScheme scheme) {
        double result = 0d;
        for (int i=0 ; i < scheme.k.size() ; i++) {
            result += scheme.coefficients.get(i) * fun.apply(scheme.nodes.get(i));
        }
        return result;
    }

    public static double gauss2D(Function fun, IntegralScheme scheme) {
        double result = 0d;
        for(int i=0 ; i < scheme.k.size() ; i++) {
            for(int j=0 ; j < scheme.k.size() ; j++) {
                result += scheme.coefficients.get(i)
                    * scheme.coefficients.get(j)
                    * fun.apply(scheme.nodes.get(i), scheme.nodes.get(j));
            }
        }
        return result;
    }

}
