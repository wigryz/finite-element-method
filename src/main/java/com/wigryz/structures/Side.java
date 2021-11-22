package com.wigryz.structures;

import com.wigryz.algorithms.IntegrationScheme;
import lombok.Getter;
import lombok.Setter;
import lombok.ToString;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

@Getter
@Setter
@ToString
public class Side {

    public static short BOTTOM = 0;
    public static short RIGHT = 1;
    public static short TOP = 2;
    public static short LEFT = 3;

    double ksi[];
    double eta[];
    double N[][];

    public Side(IntegrationScheme integralScheme, short wallId) {
        int numberOfPoints = integralScheme.getK().size();
        this.ksi = new double[numberOfPoints];
        this.eta = new double[numberOfPoints];
        N = new double[numberOfPoints][4];

        List<Double> integrationPoints = new ArrayList<>(integralScheme.getNodes());

        if (wallId == BOTTOM || wallId == TOP) {
            for (int i = 0; i < numberOfPoints; i++) {
                ksi[i] = integrationPoints.get(i);
                eta[i] = (wallId == BOTTOM) ? -1.0 : 1.0;
            }
        } else if (wallId == LEFT || wallId == RIGHT) {
            Collections.reverse(integrationPoints);
            for (int i = 0; i < numberOfPoints; i++) {
                ksi[i] = (wallId == LEFT) ? -1.0 : 1.0;
                eta[i] = integrationPoints.get(i);
            }
        }

        for (int i = 0; i < integralScheme.getK().size(); i++) {
            for (int j = 0; j < 4; j++) { // co jesli mamy 3 punktowy schemat calkowania
                N[i][0] = (0.25) * (1 - ksi[i]) * (1 - eta[i]);
                N[i][1] = (0.25) * (1 + ksi[i]) * (1 - eta[i]);
                N[i][2] = (0.25) * (1 + ksi[i]) * (1 + eta[i]);
                N[i][3] = (0.25) * (1 - ksi[i]) * (1 + eta[i]);
            }
        }
    }
}
