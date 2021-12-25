package com.wigryz.algorithms;

import com.wigryz.structures.Element;
import com.wigryz.structures.Element4x2D;
import com.wigryz.structures.Grid;
import com.wigryz.structures.Node;
import com.wigryz.structures.Side;
import com.wigryz.utilities.Configuration;
import com.wigryz.utilities.MatrixUtils;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.function.ToDoubleBiFunction;
import java.util.function.ToDoubleFunction;

public class Algorithms {

    private static final int NUMBER_OF_SHAPE_FUNCTIONS = 4;
    private static final Logger log = LogManager.getLogger();

    private Algorithms() {
    }

    public static double gauss1D(ToDoubleFunction<Double> fun, IntegrationScheme scheme) {
        double result = 0d;
        for (int i = 0; i < scheme.k.size(); i++) {
            result += scheme.coefficients.get(i) * fun.applyAsDouble(scheme.nodes.get(i));
        }
        return result;
    }

    public static double gauss2D(ToDoubleBiFunction<Double, Double> fun, IntegrationScheme scheme) {
        double result = 0d;
        for (int i = 0; i < scheme.k.size(); i++) {
            for (int j = 0; j < scheme.k.size(); j++) {
                result += scheme.coefficients.get(i)
                          * scheme.coefficients.get(j)
                          * fun.applyAsDouble(scheme.nodes.get(i), scheme.nodes.get(j));
            }
        }
        return result;
    }

    public static double jacobian(int elementId, int integrationPoint, double[][] jacobian,
                                  double[][] inverseJacobian,
                                  Element4x2D element, Grid grid) {
        List<Integer> nodeIdList = grid.getElements().get(elementId).getIdList();

        for (int i = 0; i < NUMBER_OF_SHAPE_FUNCTIONS; i++) {
            double x = grid.getNodes().get(nodeIdList.get(i) - 1).getX();
            double y = grid.getNodes().get(nodeIdList.get(i) - 1).getY();

            double dNdEta = element.getDNdEta()[integrationPoint][i];
            double dNdKsi = element.getDNdKsi()[integrationPoint][i];
            // dXdKsi
            jacobian[0][0] += dNdKsi * x;
            // dYdKsi
            jacobian[0][1] += dNdKsi * y;
            // dXdEta
            jacobian[1][0] += dNdEta * x;
            // dYdEta
            jacobian[1][1] += dNdEta * y;
        }
        //obliczanie jakobianu (wyznacznika macierzy jakobiego)
        double detJ = (jacobian[0][0] * jacobian[1][1]) - (jacobian[1][0] * jacobian[0][1]);

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
        record CoefficientIdsOfIntegrationPoint(int IP, int ksi, int eta) {
        }
        List<CoefficientIdsOfIntegrationPoint> points;

        if (element.getIntegrationScheme() == IntegrationScheme.INTEGRATION_SCHEME_1N) {
            points = List.of(new CoefficientIdsOfIntegrationPoint(0, 0, 0),
                             new CoefficientIdsOfIntegrationPoint(1, 1, 0),
                             new CoefficientIdsOfIntegrationPoint(2, 1, 1),
                             new CoefficientIdsOfIntegrationPoint(3, 0, 1));
        } else {
            points = List.of(new CoefficientIdsOfIntegrationPoint(0, 0, 0),
                             new CoefficientIdsOfIntegrationPoint(1, 2, 0),
                             new CoefficientIdsOfIntegrationPoint(2, 2, 2),
                             new CoefficientIdsOfIntegrationPoint(3, 0, 2),
                             new CoefficientIdsOfIntegrationPoint(4, 1, 0),
                             new CoefficientIdsOfIntegrationPoint(5, 2, 1),
                             new CoefficientIdsOfIntegrationPoint(6, 1, 2),
                             new CoefficientIdsOfIntegrationPoint(7, 0, 1),
                             new CoefficientIdsOfIntegrationPoint(8, 1, 1));
        }

        double conductivity = Configuration.getInstance().conductivity();

        double[][] dNdEta = element.getDNdEta();
        double[][] dNdKsi = element.getDNdKsi();

        double[] dNidx = new double[dNdEta[0].length];
        double[] dNidy = new double[dNdEta[0].length];

        for (int j = 0; j < NUMBER_OF_SHAPE_FUNCTIONS; j++) {
            dNidx[j] += inverseJacobian[0][0] * dNdKsi[integrationPoint][j];
            dNidx[j] += inverseJacobian[0][1] * dNdEta[integrationPoint][j];

            dNidy[j] += inverseJacobian[1][0] * dNdKsi[integrationPoint][j];
            dNidy[j] += inverseJacobian[1][1] * dNdEta[integrationPoint][j];
        }

        RealMatrix xMatrix = new Array2DRowRealMatrix(dNidx);
        RealMatrix resultX = xMatrix.multiply(xMatrix.transpose());
        RealMatrix yMatrix = new Array2DRowRealMatrix(dNidy);
        RealMatrix resultY = yMatrix.multiply(yMatrix.transpose());
        double[][] h = resultX.add(resultY)
                              .scalarMultiply(conductivity)
                              .scalarMultiply(detJ)
                              .scalarMultiply(element.getIntegrationScheme().getCoefficients()
                                                     .get(points.get(integrationPoint).eta))
                              .scalarMultiply(element.getIntegrationScheme().getCoefficients()
                                                     .get(points.get(integrationPoint).ksi))
                              .getData();

        // obliczanie macierzy C
        double specificHeat = Configuration.getInstance().specificHeat();
        double density = Configuration.getInstance().density();

        double[] nArray = element.getN()[integrationPoint];
        RealVector vector = new ArrayRealVector(nArray);
        RealMatrix resultMatrix = vector.outerProduct(vector);

        double[][] c = resultMatrix
            .scalarMultiply(specificHeat)
            .scalarMultiply(density)
            .scalarMultiply(detJ)
            .scalarMultiply(element.getIntegrationScheme().getCoefficients()
                                   .get(points.get(integrationPoint).eta))
            .scalarMultiply(element.getIntegrationScheme().getCoefficients()
                                   .get(points.get(integrationPoint).ksi))
            .getData();
        return Map.of("H", h, "C", c);
    }

    public static Map<String, Object> calculateHbcAndP(Grid grid, int elementId, double alpha,
                                                       double t, Element4x2D universalElement) {
        Element element = grid.getElements().get(elementId);
        List<Node> nodes = new ArrayList<>(element.getIdList().size());
        element.getIdList().forEach(id -> nodes.add(grid.getNodes().get(id - 1)));

        RealMatrix hbcMatrix = new Array2DRowRealMatrix(4, 4);
        RealVector pVector = new ArrayRealVector(4);

        double detJ;
        if (nodes.get(0).getBoundaryCondition() == nodes.get(1).getBoundaryCondition() &&
            nodes.get(0).getBoundaryCondition() != 0) { //dolna
            detJ = calculateDetJ(nodes.get(0), nodes.get(1));
            hbcMatrix = hbcMatrix.add(
                new Array2DRowRealMatrix(
                    calculateSideHBC(detJ, alpha, Side.BOTTOM, universalElement)));
            pVector = pVector.add(new ArrayRealVector(
                calculateSideP(detJ, alpha, t, Side.BOTTOM, universalElement)));
            log.debug("Hbc matrix of bottom side:\n{}\n",
                      MatrixUtils.matrixToString(hbcMatrix.getData()));
            log.debug("pMatrix of bottom side:\n{}\n",
                      Arrays.toString(pVector.toArray()));
        }
        if (nodes.get(1).getBoundaryCondition() == nodes.get(2).getBoundaryCondition() &&
            nodes.get(1).getBoundaryCondition() != 0) { //prawa
            detJ = calculateDetJ(nodes.get(1), nodes.get(2));
            hbcMatrix = hbcMatrix.add(
                new Array2DRowRealMatrix(
                    calculateSideHBC(detJ, alpha, Side.RIGHT, universalElement)));
            pVector = pVector.add(new ArrayRealVector(
                calculateSideP(detJ, alpha, t, Side.RIGHT, universalElement)));
            log.debug("Hbc matrix of right side:\n{}\n",
                      MatrixUtils.matrixToString(hbcMatrix.getData()));
            log.debug("pMatrix of right side:\n{}\n",
                      Arrays.toString(pVector.toArray()));
        }
        if (nodes.get(2).getBoundaryCondition() == nodes.get(3).getBoundaryCondition() &&
            nodes.get(2).getBoundaryCondition() != 0) { //gorna
            detJ = calculateDetJ(nodes.get(3), nodes.get(2));
            hbcMatrix = hbcMatrix.add(
                new Array2DRowRealMatrix(
                    calculateSideHBC(detJ, alpha, Side.TOP, universalElement)));
            pVector = pVector.add(new ArrayRealVector(
                calculateSideP(detJ, alpha, t, Side.TOP, universalElement)));
            log.debug("Hbc matrix of top side:\n{}\n",
                      MatrixUtils.matrixToString(hbcMatrix.getData()));
            log.debug("pMatrix of top side:\n{}\n",
                      Arrays.toString(pVector.toArray()));
        }
        if (nodes.get(3).getBoundaryCondition() == nodes.get(0).getBoundaryCondition() &&
            nodes.get(3).getBoundaryCondition() != 0) { //lewa
            detJ = calculateDetJ(nodes.get(0), nodes.get(3));
            hbcMatrix = hbcMatrix.add(
                new Array2DRowRealMatrix(
                    calculateSideHBC(detJ, alpha, Side.LEFT, universalElement)));
            pVector = pVector.add(new ArrayRealVector(
                calculateSideP(detJ, alpha, t, Side.LEFT, universalElement)));
            log.debug("Hbc matrix of left side:\n{}\n",
                      MatrixUtils.matrixToString(hbcMatrix.getData()));
            log.debug("pMatrix of left side:\n{}\n",
                      Arrays.toString(pVector.toArray()));
        }
        return Map.of("HBC", hbcMatrix.getData(), "P", pVector.toArray());
    }

    public static double[][] calculateSideHBC(double detJ, double alpha, short side,
                                              Element4x2D element) {
        IntegrationScheme scheme = element.getIntegrationScheme();
        double[][] nArray = element.getSides()[side].getN();

        RealMatrix result = new Array2DRowRealMatrix(NUMBER_OF_SHAPE_FUNCTIONS,
                                                     NUMBER_OF_SHAPE_FUNCTIONS);

        for (int i = 0; i < scheme.k.size(); i++) {
            RealMatrix nRow = new Array2DRowRealMatrix(nArray[i]);
            result =
                result.add(nRow.multiply(nRow.transpose())
                               .scalarMultiply(element.getIntegrationScheme()
                                                      .getCoefficients()
                                                      .get(i)));
        }
        return result.scalarMultiply(detJ)
                     .scalarMultiply(alpha)
                     .getData();
    }

    private static double[] calculateSideP(double detJ, double alpha, double t, short side,
                                           Element4x2D element) {
        IntegrationScheme scheme = element.getIntegrationScheme();
        double[][] nArray = element.getSides()[side].getN();

        RealVector result = new ArrayRealVector(NUMBER_OF_SHAPE_FUNCTIONS);

        for (int i = 0; i < scheme.k.size(); i++) {
            RealVector nRow = new ArrayRealVector(nArray[i]);
            result = result.add(nRow.mapMultiply(t).mapMultiply(
                element.getIntegrationScheme().getCoefficients().get(i)));
        }
        return result.mapMultiply(detJ).mapMultiply(alpha).toArray();
    }

    private static double calculateDetJ(Node x1, Node x2) {
        return Math.sqrt(Math.pow(x2.getX() - x1.getX(), 2) +
                         Math.pow((x1.getY() - x2.getY()), 2)) / 2.0;
    }
}
