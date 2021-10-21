package com.wigryz.structures;

import lombok.AllArgsConstructor;
import lombok.Getter;
import lombok.Setter;

import java.util.ArrayList;
import java.util.List;

@Getter
@Setter
@AllArgsConstructor
public class Grid {
    private double h; //height of table
    private double b; //width of table
    private int nH;  //number of elements on height
    private int nB; //number of elements on width
    private int nN; //number of nodes
    private int nE; //number of elements
    private List<Node> nodes;
    private List<Element> elements;

    private double dx; //delta x
    private double dy; // delta y

    public Grid(double h, double b, int nH, int nB) {
        this.h = h;
        this.b = b;
        this.nH = nH;
        this.nB = nB;
        nodes = new ArrayList<>();
        elements = new ArrayList<>();
        nN = nH * nB;
        nE = (nH - 1) * (nB - 1);
        dx = b / (nB - 1);
        dy = h / (nH - 1);
        createElement();
        createNodes();
    }

    private void createNodes() {
        int id = 1;
        for (double x = 0; x <= b; x += dx) {
            for (double y = 0; y <= h; y += dy) {
                nodes.add(new Node(id++, x, y));
            }
        }
    }

    private void createElement() {
        int id = 1;
        for (int n = 1 ; n +nH + 1 <= nN ; n++) {
            int id2 = n + nH;
            int id3 = id2 + 1;
            int id4 = n + 1;
            if((n % nH) == 0)
                continue;
            elements.add(new Element(id++, n, id2, id3, id4));
        }
    }
}
