package com.wigryz.structures;

import com.wigryz.algorithms.Algorithms;
import com.wigryz.utilities.Configuration;
import com.wigryz.utilities.MatrixUtils;
import lombok.AllArgsConstructor;
import lombok.Getter;
import lombok.Setter;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.DecompositionSolver;
import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.DoubleSummaryStatistics;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

@Getter
@Setter
@AllArgsConstructor
public class Grid {

    private static final Logger log = LogManager.getLogger();

    private double height; //height of table
    private double width; //width of table
    private int nH;  //number of nodes on height
    private int nB; //number of nodes on width
    private int nN; //number of nodes
    private int nE; //number of elements
    private double dx; //delta x
    private double dy; // delta y
    private List<Node> nodes;
    private List<Element> elements;
    private double[][] globalH;
    private double[] globalP;
    private double[][] globalC;
    private int currentTime = 0;


    public Grid(double height, double width, int nH, int nB, Element4x2D element) {
        this.height = height;
        this.width = width;
        this.nH = nH;
        this.nB = nB;
        nodes = new ArrayList<>();
        elements = new ArrayList<>();
        nN = nH * nB;
        nE = (nH - 1) * (nB - 1);
        dx = width / (nB - 1);
        dy = height / (nH - 1);
        globalH = new double[nN][nN];
        globalC = new double[nN][nN];
        globalP = new double[nN];
        createElement();
        createNodes();

        calculate(element);
        aggregate();
    }

    @Deprecated
    public Grid(Configuration configuration, boolean fromFile) {
        log.info("Configuration:\n{}", configuration);
        this.height = configuration.heightOfGrid();
        this.width = configuration.widthOfGrid();
        this.nH = configuration.numberOfNodesOnHeight();
        this.nB = configuration.numberOfNodesOnWidth();
        nodes = new ArrayList<>();
        elements = new ArrayList<>();
        nN = nH * nB;
        nE = (nH - 1) * (nB - 1);
        dx = width / (nB - 1);
        dy = height / (nH - 1);
        globalH = new double[nN][nN];
        globalC = new double[nN][nN];
        globalP = new double[nN];
        createElement();
        createNodes();

        calculate(new Element4x2D(configuration.integrationScheme()));
        aggregate();
    }

    public Grid(Configuration configuration) {
        log.info("Configuration:\n{}", configuration);
        nodes = configuration.nodes();
        elements = configuration.elements();
        nN = configuration.numberOfNodes();
        nB = configuration.numberOfNodesOnWidth();
        nE = configuration.numberOfElements();

        globalH = new double[nN][nN];
        globalC = new double[nN][nN];
        globalP = new double[nN];

        Element4x2D element4x2D = new Element4x2D(configuration.integrationScheme());
        calculate(element4x2D);
        aggregate();
    }

    private void createNodes() {
        double initialTemperature = Configuration.getInstance().initialTemperature();
        int id = 1;
        double x = 0.0;
        double y;
        for (double i = 0; i < nB; i++) {
            y = 0.0;
            for (double j = 0; j < nH; j++) {
                short bc = 0;
                if (i == 0 || i == nB - 1 || j == 0 || j == nH - 1)
                    bc = 1;
                nodes.add(new Node(id++, x, y, bc, initialTemperature));
                y += dy;
            }
            x += dx;
        }
    }

    private void createElement() {
        int id = 1;
        for (int n = 1; n + nH + 1 <= nN; n++) {
            int id2 = n + nH;
            int id3 = id2 + 1;
            int id4 = n + 1;
            if ((n % nH) == 0)
                continue;
            elements.add(new Element(id++, n, id2, id3, id4));
        }
    }

    public void calculate(Element4x2D element) {
        log.info("Started calculating grid.");
        for (int i = 0; i < getNE(); i++) {
            log.info("\nElement with id {}:\n{}\n", i, elements.get(i));
            log.info("Nodes:\n{}",
                     elements.get(i).getIdList().stream()
                             .map(id -> getNodes().get(id - 1)).map(Node::toString)
                             .collect(Collectors.joining("\n")));
            double[][] hLocal = new double[4][4];
            double[][] cLocal = new double[4][4];
            for (int j = 0; j < element.getNumberOfPoints(); j++) {
                log.info("Integration point number {}.", j);
                double[][] jacobianMatrix = new double[2][2];
                double[][] inverseJacobianMatrix = new double[2][2];
                double detJ = Algorithms.jacobian(i, j, jacobianMatrix, inverseJacobianMatrix,
                                                  element, this);
                log.info("detJ: {}\n", detJ);
                log.info("Jacobian matrix:\n{}\n",
                         MatrixUtils.matrixToString(jacobianMatrix));
                log.info("Inverted jacobian matrix:\n{}\n",
                         MatrixUtils.matrixToString(inverseJacobianMatrix));
                Map<String, double[][]> hAndC =
                    Algorithms.calculateHAndCOfIntPoint(inverseJacobianMatrix, j, detJ, element);
                double[][] hOfIntegralPoint = hAndC.get("H");
                double[][] cOfIntegralPoint = hAndC.get("C");
                log.info("H matrix:\n{}\n", MatrixUtils.matrixToString(hOfIntegralPoint));
                log.info("C matrix:\n{}\n", MatrixUtils.matrixToString(cOfIntegralPoint));
                for (int k = 0; k < hLocal.length; k++) {
                    for (int l = 0; l < hLocal[k].length; l++) {
                        hLocal[k][l] += hOfIntegralPoint[k][l];
                        cLocal[k][l] += cOfIntegralPoint[k][l];
                    }
                }
            }
            Map<String, Object> hbcAndP =
                Algorithms.calculateHbcAndP(this, i, Configuration.getInstance().alfa(),
                                            Configuration.getInstance().ambientTemperature(),
                                            element);
            log.info("Hbc matrix:\n{}\n",
                     MatrixUtils.matrixToString((double[][]) hbcAndP.get("HBC")));
            log.info("P vector:\n{}\n", Arrays.toString((double[]) hbcAndP.get("P")));

            getElements().get(i).setH(hLocal);
            getElements().get(i).setC(cLocal);
            getElements().get(i).setHbc((double[][]) hbcAndP.get("HBC"));
            getElements().get(i).setP((double[]) hbcAndP.get("P"));
        }
    }

    public void aggregate() {
        for (int i = 0; i < getNE(); i++) {
            Element element = elements.get(i);
            int[] id = element.getIdList().stream().mapToInt(val -> val).toArray();
            for (int h = 0; h < 4; h++) {
                for (int g = 0; g < 4; g++) {
                    globalH[id[h] - 1][id[g] - 1] += (element.getH())[h][g] + (element.getHbc())[h][g];
                    globalC[id[h] - 1][id[g] - 1] += (element.getC())[h][g];
                }
                globalP[id[h] - 1] += (element.getP())[h];
            }
        }
    }

    public void iterate() {
        double dT = Configuration.getInstance().simulationStepTime();

        RealMatrix hGlobalMatrix = new Array2DRowRealMatrix(this.getGlobalH(), false);
        RealMatrix cGlobalMatrix = new Array2DRowRealMatrix(this.getGlobalC(), false);
//        RealMatrix cGlobalMatrix = new Array2DRowRealMatrix(this.nN, this.nN);
        RealVector pGlobalVector = new ArrayRealVector(this.getGlobalP());

        RealVector t0 = new ArrayRealVector(getNodes().stream()
                                                      .mapToDouble(Node::getTemperature)
                                                      .toArray());
        RealMatrix hDash = hGlobalMatrix.add(cGlobalMatrix.scalarMultiply(1 / dT));
        RealVector pDash = pGlobalVector.add(cGlobalMatrix.scalarMultiply(1 / dT)
                                                          .operate(t0));

        DecompositionSolver solver = new LUDecomposition(hDash).getSolver();
        RealVector result = solver.solve(pDash);

        IntStream.range(0, nN)
                 .forEach(n -> this.getNodes().get(n).setTemperature(result.getEntry(n)));
        currentTime += Configuration.getInstance().simulationStepTime();
        printTemperatures();
    }

    public void printTemperatures() {
        log.warn(String.format("Temperature after %d seconds:", currentTime));
        double[] temperatures = this.getNodes()
                                    .stream()
                                    .mapToDouble(Node::getTemperature)
                                    .toArray();
        int i = 1;

        StringBuilder temperaturesGrid = new StringBuilder();
        temperaturesGrid.append("\n");
        for (double temp : temperatures) {
            temperaturesGrid.append(String.format("%.5f ", temp));
            if (i % nB == 0)
                temperaturesGrid.append("\n");
            i++;
        }
        DoubleSummaryStatistics summary = this.getNodes()
                                              .stream()
                                              .mapToDouble(Node::getTemperature)
                                              .summaryStatistics();
        log.info(temperaturesGrid.toString());
        log.warn(String.format("Min temp: %.5f Max temp: %.5f%n",
                               summary.getMin(), summary.getMax()));
    }
}