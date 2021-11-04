package com.wigryz.algorithms;

import com.wigryz.structures.Element4x2D;
import com.wigryz.structures.Grid;
import java.util.List;

public class Algorithms {

    public static double gauss1D(Function fun, IntegralScheme scheme) {
        double result = 0d;
        for (int i=0 ; i < scheme.k.size() ; i++) {
            result += scheme.coefficients.get(i) * fun.apply(scheme.nodes.get(i));
        }
        return result;
    }

    public static double gauss2D(Function fun, IntegralScheme scheme) {
        double result = 0d;
        for(int i=0 ; i < scheme.k.size() ; i++) {
            for(int j=0 ; j < scheme.k.size() ; j++) {
                result += scheme.coefficients.get(i)
                    * scheme.coefficients.get(j)
                    * fun.apply(scheme.nodes.get(i), scheme.nodes.get(j));
            }
        }
        return result;
    }


    /*
    4 - 3
    |   |
    1 - 2
    Liczenie jakobianu:
     - pobranie listy id-ków z danego elementu,
     - iterowanie po liscie idków,
      * pobranie wspolrzednych nodea o pierwszym id,
      * pobranie eta i ksi z tablicy o indeksie [numer punktu calkowania][numer id]
      * sumowanie iloczynu eta i wspolrzednej (x/y) w danym punkcie,
     - liczenie wyznacznika macierzy,
     - liczenie jakobianu odwroconego,

    Wynikiem jest zawartosc tablicy J i JInv.
     */
    public static void jacobian(int i, int j, double[][] J, double[][] JInv, Element4x2D element,
                                Grid grid) {
        List<Integer> nodeIdList = grid.getElements().get(i).getIdList();

        for(int k=0 ; k < nodeIdList.size() ; k++) {
            double x = grid.getNodes().get(nodeIdList.get(k) - 1).getX();
            double y = grid.getNodes().get(nodeIdList.get(k) - 1).getY();

            double eta = element.getEtaArray()[j][k];
            double ksi = element.getKsiArray()[j][k];
            // dYdEta
            J[0][0] += eta * y;
            // -dYdKsi
            J[0][1] += -ksi * y;
            // -dXdEta
            J[1][0] += -eta * x;
            // dXdKsi
            J[1][1] += ksi * x;
        }
        double det = J[0][0] * J[1][1] - J[1][0] * J[0][1];

        for(int k=0 ; k < J.length; k++)
            for(int h=0 ; h < J[0].length ; h++)
                JInv[k][h] = 1d/det * J[k][h];
    }
}
