package com.wigryz.structures;

import com.wigryz.algorithms.Algorithms;
import lombok.AllArgsConstructor;
import lombok.Getter;
import lombok.Setter;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

@Getter
@Setter
@AllArgsConstructor
public class Grid {
    private double h; //height of table
    private double b; //width of table
    private int nH;  //number of nodes on height
    private int nB; //number of nodes on width
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
        for (int n = 1; n + nH + 1 <= nN; n++) {
            int id2 = n + nH;
            int id3 = id2 + 1;
            int id4 = n + 1;
            if ((n % nH) == 0)
                continue;
            elements.add(new Element(id++, n, id2, id3, id4));
        }
    }

    public void calculate(Element4x2D element) {
        for (int i = 0; i < getNE(); i++) {
            double[][] H = new double[4][4];
            for (int j = 0; j < element.getNumberOfPoints(); j++) {
                double[][] J = new double[2][2];
                double[][] JInv = new double[2][2];
                Algorithms.jacobian(i, j, J, JInv, element, this);
//                System.out.println(Arrays.deepToString(J).replace("], ", "]\n"));
//                System.out.println(Arrays.deepToString(JInv).replace("], ", "]\n"));
                double[][] HOfIntegralPoint = calculateH(JInv, j, element);
                System.out.println(Arrays.deepToString(HOfIntegralPoint).replace("], ", "]\n") + "\n");
                for(int g=0 ; g < H.length ; g++) {
                    for (int h=0 ; h < H[g].length ; h++) {
                        H[h][g] += HOfIntegralPoint[h][g];
                    }
                }
            }
            getElements().get(i).setH(H);
            System.out.println(Arrays.deepToString(H).replace("], ", "]\n"));
        }
    }

    public double[][] calculateH(double[][] JInv, int i, Element4x2D element) {
        double k_t = 30.0;
        double dV = 0.000156;
        double[] xArray = new double[4];
        double[] yArray = new double[4];

        double[][] etaArray = element.getEtaArray();
        double[][] ksiArray = element.getKsiArray();

        for (int j = 0; j < element.getNumberOfPoints(); j++) {
            xArray[j] += JInv[0][0] * ksiArray[i][j];
            xArray[j] += JInv[0][1] * etaArray[i][j];

            yArray[j] += JInv[1][0] * ksiArray[i][j];
            yArray[j] += JInv[1][1] * etaArray[i][j];
        }
        //cos nie dziala z obliczaniem macierzy H dla kolejnych punktów całkowania
        RealMatrix xMatrix = new Array2DRowRealMatrix(xArray);
        RealMatrix resultX = xMatrix.multiply(xMatrix.transpose());
        RealMatrix yMatrix = new Array2DRowRealMatrix(yArray);
        RealMatrix resultY = yMatrix.multiply(yMatrix.transpose());
        return resultX.add(resultY).scalarMultiply(k_t).scalarMultiply(dV).getData();
    }
}
