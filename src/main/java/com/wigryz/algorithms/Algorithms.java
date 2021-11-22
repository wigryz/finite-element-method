package com.wigryz.algorithms;

import com.wigryz.structures.Element;
import com.wigryz.structures.Element4x2D;
import com.wigryz.structures.Grid;
import com.wigryz.structures.Node;
import com.wigryz.structures.Side;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;

import java.util.ArrayList;
import java.util.List;

public class Algorithms {

    private Algorithms() {
    }

    public static double gauss1D(Function fun, IntegrationScheme scheme) {
        double result = 0d;
        for (int i = 0; i < scheme.k.size(); i++) {
            result += scheme.coefficients.get(i) * fun.apply(scheme.nodes.get(i));
        }
        return result;
    }

    public static double gauss2D(Function fun, IntegrationScheme scheme) {
        double result = 0d;
        for (int i = 0; i < scheme.k.size(); i++) {
            for (int j = 0; j < scheme.k.size(); j++) {
                result += scheme.coefficients.get(i)
                    * scheme.coefficients.get(j)
                    * fun.apply(scheme.nodes.get(i), scheme.nodes.get(j));
            }
        }
        return result;
    }

    public static void jacobian(int i, int j, double[][] J, double[][] JInv, Element4x2D element,
                                Grid grid) {
        List<Integer> nodeIdList = grid.getElements().get(i).getIdList();

        for (int k = 0; k < nodeIdList.size(); k++) {
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
        //obliczanie jakobianu (wyznacznika macierzy jakobiego)
        double det = J[0][0] * J[1][1] - J[1][0] * J[0][1];

        // transponowanie macierzy jakobiego
        JInv[0][0] = J[1][1];
        JInv[0][1] = -J[0][1];
        JInv[1][0] = -J[1][0];
        JInv[1][1] = J[0][0];

        for (int k = 0; k < J.length; k++)
            for (int h = 0; h < J[0].length; h++)
                JInv[k][h] = 1d / det * J[k][h];
    }

    public static double[][] calculateHOfIntPoint(double[][] JInv, int integrationPoint,
                                                  Element4x2D element) {
        double k_t = 30.0; // wspolczynnik przewodzenia ciepla
        double dV = 0.000156; // zmiana objetosci //to jest chyba wyznacznik jakobianu

        double[][] etaArray = element.getEtaArray();
        double[][] ksiArray = element.getKsiArray();

        double[] dNi_dx = new double[etaArray.length];
        double[] dNi_dy = new double[etaArray.length];


        for (int j = 0; j < etaArray.length; j++) {
            dNi_dx[j] += JInv[0][0] * ksiArray[integrationPoint][j];
            dNi_dx[j] += JInv[0][1] * etaArray[integrationPoint][j];

            dNi_dy[j] += JInv[1][0] * ksiArray[integrationPoint][j];
            dNi_dy[j] += JInv[1][1] * etaArray[integrationPoint][j];
        }

        RealMatrix xMatrix = new Array2DRowRealMatrix(dNi_dx);
        RealMatrix resultX = xMatrix.multiply(xMatrix.transpose());
        RealMatrix yMatrix = new Array2DRowRealMatrix(dNi_dy);
        RealMatrix resultY = yMatrix.multiply(yMatrix.transpose());
        return resultX.add(resultY).scalarMultiply(k_t).scalarMultiply(dV).getData();
    }

    public static double[][] calculateHBC(Grid grid, int elementId, double alpha,
                                          Element4x2D universalElement) {
        Element element = grid.getElements().get(elementId);
        List<Node> nodes = new ArrayList<>(element.getIdList().size());
        element.getIdList().forEach(id -> nodes.add(grid.getNodes().get(id - 1)));

        RealMatrix HBC = new Array2DRowRealMatrix(4, 4);

        double detJ = 0.0;

        if (nodes.get(0).getBC() == nodes.get(1).getBC() && nodes.get(0).getBC() != 0) { //dolna
            detJ = (nodes.get(1).getX() - nodes.get(0).getX()) / 2.0;
            HBC = HBC.add(
                new Array2DRowRealMatrix(
                    calculateSideHBC(detJ, alpha, Side.BOTTOM, universalElement)));
        }
        if (nodes.get(1).getBC() == nodes.get(2).getBC() && nodes.get(1).getBC() != 0) { //prawa
            detJ = (nodes.get(2).getY() - nodes.get(1).getY()) / 2.0;
            HBC = HBC.add(
                new Array2DRowRealMatrix(
                    calculateSideHBC(detJ, alpha, Side.RIGHT, universalElement)));
        }
        if (nodes.get(2).getBC() == nodes.get(3).getBC() && nodes.get(0).getBC() != 0) { //gorna
            detJ = (nodes.get(2).getX() - nodes.get(3).getX()) / 2.0;
            HBC = HBC.add(
                new Array2DRowRealMatrix(
                    calculateSideHBC(detJ, alpha, Side.TOP, universalElement)));
        }
        if (nodes.get(3).getBC() == nodes.get(0).getBC() && nodes.get(0).getBC() != 0) { //lewa
            detJ = (nodes.get(3).getY() - nodes.get(0).getY()) / 2.0;
            HBC = HBC.add(
                new Array2DRowRealMatrix(
                    calculateSideHBC(detJ, alpha, Side.LEFT, universalElement)));
        }
        return HBC.getData();
    }

    /*
        PYTANIE
        czy nie łatwiej byłoby obliczyć raz cztery macierze Hbc dla każdej ze ścian
        a później tylko mnożyć je przez alfę i detJ?
    */


    // alfa to współczynnik przewodzenia ciepła czy jakos tak
    // SKĄD WZIĄĆ DET J ?!?
    public static double[][] calculateSideHBC(double detJ, double alpha, short side,
                                              Element4x2D element) {
        IntegrationScheme scheme = element.getIntegrationScheme();
        double[][] NArray = element.getSides()[side].getN();

        RealMatrix result = new Array2DRowRealMatrix(element.getNumberOfPoints(),
            element.getNumberOfPoints());

        for (int i = 0; i < scheme.k.size(); i++) {
            RealMatrix nRow = new Array2DRowRealMatrix(NArray[i]);
            result = result.add(nRow.multiply(nRow.transpose()));
        }
        return result.scalarMultiply(detJ).scalarMultiply(alpha).getData();
    }
}
