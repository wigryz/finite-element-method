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

    public static final short BOTTOM = 0;
    public static final short RIGHT = 1;
    public static final short TOP = 2;
    public static final short LEFT = 3;

    double[] ksi;
    double[] eta;
    double[][] n;

    public Side(IntegrationScheme integralScheme, short wallId) {
        int numberOfPoints = integralScheme.getK().size();
        this.ksi = new double[numberOfPoints];
        this.eta = new double[numberOfPoints];
        n = new double[numberOfPoints][4];

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
            n[i][0] = (0.25) * (1 - ksi[i]) * (1 - eta[i]);
            n[i][1] = (0.25) * (1 + ksi[i]) * (1 - eta[i]);
            n[i][2] = (0.25) * (1 + ksi[i]) * (1 + eta[i]);
            n[i][3] = (0.25) * (1 - ksi[i]) * (1 + eta[i]);
        }
    }
}
