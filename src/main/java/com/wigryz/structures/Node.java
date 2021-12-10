package com.wigryz.structures;


import lombok.AllArgsConstructor;
import lombok.Getter;
import lombok.Setter;
import lombok.ToString;

@Getter
@Setter
@AllArgsConstructor
@ToString
public class Node {

    private int id;
    private double x;
    private double y;
    private short boundaryCondition; // flaga warunku brzegowego
    private double temperature;
}
