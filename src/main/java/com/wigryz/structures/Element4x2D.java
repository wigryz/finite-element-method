package com.wigryz.structures;

import com.wigryz.algorithms.IntegrationScheme;
import lombok.Data;

import java.util.Arrays;

@Data
public class Element4x2D {

    IntegrationScheme integrationScheme;
    private double[][] dNdEta;
    private double[][] dNdKsi;
    private double[][] n;
    private int numberOfPoints;


    private Side[] sides;

    public Element4x2D(IntegrationScheme integrationScheme) {
        this.integrationScheme = integrationScheme;
        numberOfPoints = (int) Math.pow(integrationScheme.getK().size(), 2);
        this.dNdEta = new double[numberOfPoints][4];
        this.dNdKsi = new double[numberOfPoints][4];
        this.n = new double[numberOfPoints][4];
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
            // KSI
            if (i == 0 || i == 3) { //lewa
                ksi = integrationScheme.getNodes().get(0);
            }
            if (i == 1 || i == 2) { //prawa
                ksi = integrationScheme.getNodes().get(1);
            }
            // ETA
            if (i == 0 || i == 1) { //dol
                eta = integrationScheme.getNodes().get(0);
            }
            if (i == 2 || i == 3) { //gora
                eta = integrationScheme.getNodes().get(1);
            }
            fillRow(i, ksi, eta);
        }
    }

    private void fillArrays2N() {
        double ksi = 0.0;
        double eta = 0.0;
        for (int i = 0; i < numberOfPoints; i++) {
            // KSI
            if (i == 0 || i == 3 || i == 7) { //lewa
                ksi = integrationScheme.getNodes().get(0);
            }
            if (i == 4 || i == 8 || i == 6) { //srodek
                ksi = integrationScheme.getNodes().get(1);
            }
            if (i == 1 || i == 5 || i == 2) { // prawa
                ksi = integrationScheme.getNodes().get(2);
            }
            // ETA
            if (i == 0 || i == 4 || i == 1) { //dol
                eta = integrationScheme.getNodes().get(0);
            }
            if (i == 7 || i == 8 || i == 5) { //srodek
                eta = integrationScheme.getNodes().get(1);
            }
            if (i == 3 || i == 6 || i == 2) { // gora
                eta = integrationScheme.getNodes().get(2);
            }
            fillRow(i, ksi, eta);
        }
    }

    private void fillRow(int i, double ksi, double eta) {
        dNdKsi[i][0] = -0.25 * (1 - eta);
        dNdKsi[i][1] = 0.25 * (1 - eta);
        dNdKsi[i][2] = 0.25 * (1 + eta);
        dNdKsi[i][3] = -0.25 * (1 + eta);

        dNdEta[i][0] = -0.25 * (1 - ksi);
        dNdEta[i][1] = -0.25 * (1 + ksi);
        dNdEta[i][2] = 0.25 * (1 + ksi);
        dNdEta[i][3] = 0.25 * (1 - ksi);

        n[i][0] = 0.25 * (1 - ksi) * (1 - eta);
        n[i][1] = 0.25 * (1 + ksi) * (1 - eta);
        n[i][2] = 0.25 * (1 + ksi) * (1 + eta);
        n[i][3] = 0.25 * (1 - ksi) * (1 + eta);
    }

    private void fillSides() {
        for (short i = 0; i < 4; i++) {
            this.sides[i] = new Side(integrationScheme, i);
        }
    }

    @Override
    public String toString() {
        return "Element4x2D{" +
            "\netaArray=\n" + Arrays.deepToString(dNdEta).replace("], ", "]\n") +
            "\nksiArray=\n" + Arrays.deepToString(dNdKsi).replace("], ", "]\n") +
            '}';
    }
}