package com.wigryz.algorithms;

import lombok.Getter;

import java.util.List;

@Getter
public enum IntegrationScheme {

    INTEGRATION_SCHEME_1N(
        List.of(0, 1),
        List.of(-1.0/Math.sqrt(3), 1.0/Math.sqrt(3)),
        List.of(1.0, 1.0)),

    INTEGRATION_SCHEME_2N(
        List.of(0, 1, 2),
        List.of(-Math.sqrt(3)/Math.sqrt(5), 0d, Math.sqrt(3)/Math.sqrt(5)),
        List.of(5d/9d, 8d/9d, 5d/9d));

    IntegrationScheme(List<Integer> k, List<Double> nodes, List<Double> coefficients) {
        this.k = k;
        this.nodes = nodes;
        this.coefficients = coefficients;
    }

    List<Integer> k;
    List<Double> nodes;
    List<Double> coefficients;
}