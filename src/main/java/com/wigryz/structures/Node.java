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
    private short BC; // flaga warunku brzegowego
    private double temperature;

//    public Node(int id, double x, double y, short BC, double initialTemperature) {
//        this.id = id;
//        this.x = x;
//        this.y = y;
//        this.BC = BC;
//        this.temperature = initialTemperature;
//    }
}
