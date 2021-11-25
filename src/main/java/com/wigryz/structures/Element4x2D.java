package com.wigryz.structures;

import com.wigryz.algorithms.IntegrationScheme;
import lombok.Data;

import java.util.Arrays;

@Data
public class Element4x2D {

    IntegrationScheme integrationScheme;

    private double coord = 1d/Math.sqrt(3);
    private double[][] etaArray;
    private double[][] ksiArray;
    private int numberOfPoints;

    private Side[] sides;

    public Element4x2D(IntegrationScheme integrationScheme) {
        this.integrationScheme = integrationScheme;
//        coord = integrationScheme.getNodes().get(1); // TODO jak to powinno byc
        numberOfPoints = (int)Math.pow(integrationScheme.getK().size(), 2);
        this.etaArray = new double[4][numberOfPoints];
        this.ksiArray = new double[4][numberOfPoints];
        sides = new Side[4];
        fillSides();
        fillArrays();
    }

    private void fillArrays() {
        double ksi = 0.0;
        double eta = 0.0;
        for(int i = 0; i < 4; i++)
        {
            switch (i) {
                case 0 -> { // pc1
                    ksi = -coord;
                    eta = -coord;
                }
                case 1 -> { // pc2
                    ksi = coord;
                    eta = -coord;
                }
                case 2 -> { // pc3
                    ksi = coord;
                    eta = coord;
                }
                case 3 -> { // pc4
                    ksi = -coord;
                    eta = coord;
                }
                default -> {
                    System.out.println("SHOULD NEVER GET HERE.");
                }
            }

            ksiArray[i][0] = -0.25*(1-eta);
            ksiArray[i][1] =  0.25*(1-eta);
            ksiArray[i][2] =  0.25*(1+eta);
            ksiArray[i][3] = -0.25*(1+eta);

            etaArray[i][0] = -0.25*(1-ksi);
            etaArray[i][1] = -0.25*(1+ksi);
            etaArray[i][2] =  0.25*(1+ksi);
            etaArray[i][3] =  0.25*(1-ksi);
        }
    }

    private void fillSides() {
        for(short i=0 ; i < 4 ; i++) {
            this.sides[i] = new Side(integrationScheme, i);
        }
    }

    @Override
    public String toString() {
        return "Element4x2D{" +
            "\netaArray=\n" + Arrays.deepToString(etaArray).replace("], ", "]\n") +
            "\nksiArray=\n" + Arrays.deepToString(ksiArray).replace("], ", "]\n") +
            '}';
    }
}