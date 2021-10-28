package com.wigryz.structures;

import com.wigryz.algorithms.IntegralScheme;
import lombok.Data;
import lombok.ToString;

import java.util.Arrays;

@Data
public class Element4x2D {

    private final double coord = 1d/Math.sqrt(3);
    private double[][] etaArray;
    private double[][] ksiArray;
    private int numberOfNodes;

    public Element4x2D(IntegralScheme integralScheme) {
        numberOfNodes = (int)Math.pow(integralScheme.getK().size(), 2);
        this.etaArray = new double[4][numberOfNodes];
        this.ksiArray = new double[4][numberOfNodes];
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

    @Override
    public String toString() {
        return "Element4x2D{" +
            "\netaArray=\n" + Arrays.deepToString(etaArray).replace("], ", "]\n") +
            "\nksiArray=\n" + Arrays.deepToString(ksiArray).replace("], ", "]\n") +
            '}';
    }
}
