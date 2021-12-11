package com.wigryz.structures;

import com.wigryz.algorithms.IntegrationScheme;
import lombok.Data;

import java.util.Arrays;

@Data
public class Element4x2D {

    IntegrationScheme integrationScheme;

    private double coord = 1d / Math.sqrt(3);
    private double[][] etaArray;
    private double[][] ksiArray;
    private double[][] array;
    private int numberOfPoints;


    private Side[] sides;

    public Element4x2D(IntegrationScheme integrationScheme) {
        this.integrationScheme = integrationScheme;
//        coord = integrationScheme.getNodes().get(1); // TODO jak to powinno byc
        numberOfPoints = (int) Math.pow(integrationScheme.getK().size(), 2);
        this.etaArray = new double[numberOfPoints][4];
        this.ksiArray = new double[numberOfPoints][4];
        this.array = new double[numberOfPoints][4];
        sides = new Side[4];
        if (integrationScheme.equals(IntegrationScheme.INTEGRATION_SCHEME_1N)) {
            fillArrays1N();
        } else {
            fillArrays2N();
        }
        fillSides();
    }

    private void fillArrays1N() {
        double ksi = 0.0;
        double eta = 0.0;
        for (int i = 0; i < numberOfPoints; i++) {
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
                default -> System.out.println("SHOULD NEVER GET HERE.");
            }

            ksiArray[i][0] = -0.25 * (1 - eta);
            ksiArray[i][1] = 0.25 * (1 - eta);
            ksiArray[i][2] = 0.25 * (1 + eta);
            ksiArray[i][3] = -0.25 * (1 + eta);

            etaArray[i][0] = -0.25 * (1 - ksi);
            etaArray[i][1] = -0.25 * (1 + ksi);
            etaArray[i][2] = 0.25 * (1 + ksi);
            etaArray[i][3] = 0.25 * (1 - ksi);

            array[i][0] = 0.25 * (1 - ksi) * (1 - eta);
            array[i][1] = 0.25 * (1 + ksi) * (1 - eta);
            array[i][2] = 0.25 * (1 + ksi) * (1 + eta);
            array[i][3] = 0.25 * (1 - ksi) * (1 + eta);
        }
    }

    private void fillArrays2N() {
        double ksi = 0.0;
        double eta = 0.0;
        for (int i = 0; i < numberOfPoints; i++) {
            // KSI
            if(i == 0 || i == 3 || i == 7) { //lewa
                ksi = integrationScheme.getNodes().get(0);
            }

            if(i == 4 || i == 8 || i == 6) { //srodek
                ksi = integrationScheme.getNodes().get(1);
            }

            if(i == 1 || i == 5 || i == 2) { // prawa
                ksi = integrationScheme.getNodes().get(2);
            }
            // ETA
            if(i == 0 || i == 4 || i == 1) { //dol
                eta = integrationScheme.getNodes().get(0);
            }

            if(i == 7 || i == 8 || i == 5) { //srodek
                eta = integrationScheme.getNodes().get(1);
            }

            if(i == 3 || i == 6 || i == 2) { // gora
                eta = integrationScheme.getNodes().get(2);
            }

            ksiArray[i][0] = -0.25 * (1 - eta);
            ksiArray[i][1] = 0.25 * (1 - eta);
            ksiArray[i][2] = 0.25 * (1 + eta);
            ksiArray[i][3] = -0.25 * (1 + eta);

            etaArray[i][0] = -0.25 * (1 - ksi);
            etaArray[i][1] = -0.25 * (1 + ksi);
            etaArray[i][2] = 0.25 * (1 + ksi);
            etaArray[i][3] = 0.25 * (1 - ksi);

            array[i][0] = 0.25 * (1 - ksi) * (1 - eta);
            array[i][1] = 0.25 * (1 + ksi) * (1 - eta);
            array[i][2] = 0.25 * (1 + ksi) * (1 + eta);
            array[i][3] = 0.25 * (1 - ksi) * (1 + eta);
        }
    }

    private void fillSides() {
        for (short i = 0; i < 4; i++) {
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

/*
    N1 = 0.25*(1-ksi)*(1-eta);
    N2 = 0.25*(1+ksi)*(1-eta);
    N3 = 0.25*(1+ksi)*(1+eta);
    N4 = 0.25*(1-ksi)*(1+eta);
 */