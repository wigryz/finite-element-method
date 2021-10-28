package com.wigryz.algorithms;

import java.util.List;

public class Gauss {

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
