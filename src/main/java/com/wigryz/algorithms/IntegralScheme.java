package com.wigryz.algorithms;

import lombok.AllArgsConstructor;
import lombok.Data;
import lombok.NoArgsConstructor;

import java.util.List;

@Data
@NoArgsConstructor
@AllArgsConstructor
public class IntegralScheme {

    List<Integer> k;
    List<Double> nodes;
    List<Double> coefficients;
}
