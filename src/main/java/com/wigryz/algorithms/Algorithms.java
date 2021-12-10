package com.wigryz.algorithms;

import com.wigryz.structures.Element;
import com.wigryz.structures.Element4x2D;
import com.wigryz.structures.Grid;
import com.wigryz.structures.Node;
import com.wigryz.structures.Side;
import com.wigryz.utilities.Configuration;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

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

    public static double jacobian(int i, int j, double[][] jacobian, double[][] inverseJacobian,
                                  Element4x2D element, Grid grid) {
        List<Integer> nodeIdList = grid.getElements().get(i).getIdList();

        for (int k = 0; k < nodeIdList.size(); k++) {
            double x = grid.getNodes().get(nodeIdList.get(k) - 1).getX();
            double y = grid.getNodes().get(nodeIdList.get(k) - 1).getY();

            double eta = element.getEtaArray()[j][k];
            double ksi = element.getKsiArray()[j][k];
            // dXdKsi
            jacobian[0][0] += ksi * x;
            // dYdKsi
            jacobian[0][1] += ksi * y;
            // dXdEta
            jacobian[1][0] += eta * x;
            // dYdEta
            jacobian[1][1] += eta * y;
        }
        //obliczanie jakobianu (wyznacznika macierzy jakobiego)
        double detJ = jacobian[0][0] * jacobian[1][1] - jacobian[1][0] * jacobian[0][1];

        // transponowanie macierzy jakobiego
        inverseJacobian[0][0] = jacobian[1][1];
        inverseJacobian[0][1] = -jacobian[0][1];
        inverseJacobian[1][0] = -jacobian[1][0];
        inverseJacobian[1][1] = jacobian[0][0];

        for (int k = 0; k < jacobian.length; k++)
            for (int h = 0; h < jacobian[0].length; h++)
                inverseJacobian[k][h] = 1d / detJ * jacobian[k][h];

        return detJ;
    }

    public static Map<String, double[][]> calculateHAndCOfIntPoint(double[][] inverseJacobian,
                                                                   int integrationPoint,
                                                                   double detJ,
                                                                   Element4x2D element) {
        double kT = Configuration.getInstance().conductivity(); // wspolczynnik przewodzenia ciepla
        double dV = detJ;

        double[][] etaArray = element.getEtaArray();
        double[][] ksiArray = element.getKsiArray();

        double[] dNidx = new double[etaArray.length];
        double[] dNidy = new double[etaArray.length];

        for (int j = 0; j < etaArray.length; j++) {
            dNidx[j] += inverseJacobian[0][0] * ksiArray[integrationPoint][j];
            dNidx[j] += inverseJacobian[0][1] * etaArray[integrationPoint][j];

            dNidy[j] += inverseJacobian[1][0] * ksiArray[integrationPoint][j];
            dNidy[j] += inverseJacobian[1][1] * etaArray[integrationPoint][j];
        }

        RealMatrix xMatrix = new Array2DRowRealMatrix(dNidx);
        RealMatrix resultX = xMatrix.multiply(xMatrix.transpose());
        RealMatrix yMatrix = new Array2DRowRealMatrix(dNidy);
        RealMatrix resultY = yMatrix.multiply(yMatrix.transpose());
        double[][] h = resultX.add(resultY).scalarMultiply(kT).scalarMultiply(dV).getData();

        // obliczanie macierzy C
        double specificHeat = Configuration.getInstance().specificHeat();
        double density = Configuration.getInstance().density();

        double[] array = element.getArray()[integrationPoint];
        RealVector vector = new ArrayRealVector(array);
        RealMatrix resultMatrix = vector.outerProduct(vector);

        double[][] c = resultMatrix
            .scalarMultiply(specificHeat)
            .scalarMultiply(density)
            .scalarMultiply(detJ)
            .getData();
        return Map.of("H", h, "C", c);
    }

    public static Map<String, Object> calculateHbcAndP(Grid grid, int elementId, double alpha,
                                                       double t, Element4x2D universalElement) {
        Element element = grid.getElements().get(elementId);
        List<Node> nodes = new ArrayList<>(element.getIdList().size());
        element.getIdList().forEach(id -> nodes.add(grid.getNodes().get(id - 1)));

        RealMatrix hbcMatrix = new Array2DRowRealMatrix(4, 4);
        RealVector pMatrix = new ArrayRealVector(4);

        double detJ;
        if (nodes.get(0).getBoundaryCondition() == nodes.get(1).getBoundaryCondition() &&
            nodes.get(0).getBoundaryCondition() != 0) { //dolna
            detJ = calculateDetJ(nodes.get(0), nodes.get(1));
            hbcMatrix = hbcMatrix.add(
                new Array2DRowRealMatrix(
                    calculateSideHBC(detJ, alpha, Side.BOTTOM, universalElement)));
            pMatrix = pMatrix.add(new ArrayRealVector(
                calculateSideP(detJ, alpha, t, Side.BOTTOM, universalElement)));
        }
        if (nodes.get(1).getBoundaryCondition() == nodes.get(2).getBoundaryCondition() &&
            nodes.get(1).getBoundaryCondition() != 0) { //prawa
            detJ = calculateDetJ(nodes.get(1), nodes.get(2));
            hbcMatrix = hbcMatrix.add(
                new Array2DRowRealMatrix(
                    calculateSideHBC(detJ, alpha, Side.RIGHT, universalElement)));
            pMatrix = pMatrix.add(new ArrayRealVector(
                calculateSideP(detJ, alpha, t, Side.RIGHT, universalElement)));
        }
        if (nodes.get(2).getBoundaryCondition() == nodes.get(3).getBoundaryCondition() &&
            nodes.get(2).getBoundaryCondition() != 0) { //gorna
            detJ = calculateDetJ(nodes.get(3), nodes.get(2));
            hbcMatrix = hbcMatrix.add(
                new Array2DRowRealMatrix(
                    calculateSideHBC(detJ, alpha, Side.TOP, universalElement)));
            pMatrix = pMatrix.add(new ArrayRealVector(
                calculateSideP(detJ, alpha, t, Side.TOP, universalElement)));
        }
        if (nodes.get(3).getBoundaryCondition() == nodes.get(0).getBoundaryCondition() &&
            nodes.get(3).getBoundaryCondition() != 0) { //lewa
            detJ = calculateDetJ(nodes.get(0), nodes.get(3));
            hbcMatrix = hbcMatrix.add(
                new Array2DRowRealMatrix(
                    calculateSideHBC(detJ, alpha, Side.LEFT, universalElement)));
            pMatrix = pMatrix.add(new ArrayRealVector(
                calculateSideP(detJ, alpha, t, Side.LEFT, universalElement)));
        }
        return Map.of("HBC", hbcMatrix.getData(), "P", pMatrix.toArray());
    }

    /*
        PYTANIE
        czy nie łatwiej byłoby obliczyć raz cztery macierze Hbc dla każdej ze ścian
        a później tylko mnożyć je przez alfę i detJ?
    */

    // alfa to współczynnik przewodzenia ciepła czy jakos tak
    public static double[][] calculateSideHBC(double detJ, double alpha, short side,
                                              Element4x2D element) {
        IntegrationScheme scheme = element.getIntegrationScheme();
        double[][] nArray = element.getSides()[side].getN();

        RealMatrix result = new Array2DRowRealMatrix(element.getNumberOfPoints(),
                                                     element.getNumberOfPoints());

        for (int i = 0; i < scheme.k.size(); i++) {
            RealMatrix nRow = new Array2DRowRealMatrix(nArray[i]);
            result = result.add(nRow.multiply(nRow.transpose()));
        }
        return result.scalarMultiply(detJ).scalarMultiply(alpha).getData();
    }

    private static double[] calculateSideP(double detJ, double alpha, double t, short side,
                                           Element4x2D element) {
        IntegrationScheme scheme = element.getIntegrationScheme();
        double[][] nArray = element.getSides()[side].getN();

        RealVector result = new ArrayRealVector(element.getNumberOfPoints());

        for (int i = 0; i < scheme.k.size(); i++) {
            RealVector nRow = new ArrayRealVector(nArray[i]);
            result = result.add(nRow.mapMultiply(t));
        }
        return result.mapMultiply(detJ).mapMultiply(alpha).toArray();
    }

    private static double calculateDetJ(Node x1, Node x2) {
        return Math.sqrt(Math.pow(x2.getX() - x1.getX(), 2) +
                             Math.pow((x1.getY() - x2.getY()), 2)) / 2.0;
    }
}
