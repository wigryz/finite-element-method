package com.wigryz.structures;

import com.wigryz.algorithms.Algorithms;
import com.wigryz.utilities.Configuration;
import lombok.AllArgsConstructor;
import lombok.Getter;
import lombok.Setter;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.DecompositionSolver;
import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

import java.util.ArrayList;
import java.util.DoubleSummaryStatistics;
import java.util.List;
import java.util.Map;
import java.util.stream.IntStream;

@Getter
@Setter
@AllArgsConstructor
public class Grid {
    private double h; //height of table
    private double b; //width of table
    private int nH;  //number of nodes on height
    private int nB; //number of nodes on width
    private int nN; //number of nodes
    private int nE; //number of elements
    private List<Node> nodes;
    private List<Element> elements;
    private double[][] HG;
    private double[] PG;
    private double[][] CG;

    private int currentTime = 0;

    private double dx; //delta x
    private double dy; // delta y

    public Grid(double h, double b, int nH, int nB) {
        this.h = h;
        this.b = b;
        this.nH = nH;
        this.nB = nB;
        nodes = new ArrayList<>();
        elements = new ArrayList<>();
        nN = nH * nB;
        nE = (nH - 1) * (nB - 1);
        dx = b / (nB - 1);
        dy = h / (nH - 1);
        HG = new double[nN][nN];
        CG = new double[nN][nN];
        PG = new double[nN];
        createElement();
        createNodes();
    }

    private void createNodes() {
        double initialTemperature = Configuration.getInstance().initialTemperature();
        int id = 1;
        double x = 0.0;
        double y = 0.0;
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
        for (int i = 0; i < getNE(); i++) {
            double[][] H = new double[4][4];
            double[][] C = new double[4][4];
            for (int j = 0; j < element.getNumberOfPoints(); j++) {
                double[][] J = new double[2][2];
                double[][] JInv = new double[2][2];
                double detJ = Algorithms.jacobian(i, j, J, JInv, element, this);
                Map<String, double[][]> HAndC = Algorithms.calculateHOfIntPoint(JInv, j, detJ,
                                                                                element);
                double[][] HOfIntegralPoint = HAndC.get("H");
                double[][] COfIntegralPoint = HAndC.get("C");
                for (int g = 0; g < H.length; g++) {
                    for (int h = 0; h < H[g].length; h++) {
                        H[h][g] += HOfIntegralPoint[h][g];
                        C[h][g] += COfIntegralPoint[h][g];
                    }
                }
            }
            Map<String, Object> HbcAndP =
                Algorithms.calculateHbcAndP(this, i, Configuration.getInstance().alfa(),
                                            Configuration.getInstance().ambientTemperature(),
                                            element);
            getElements().get(i).setH(H);
            getElements().get(i).setC(C);
            getElements().get(i).setHbc((double[][]) HbcAndP.get("HBC"));
            getElements().get(i).setP((double[]) HbcAndP.get("P"));
//            System.out.println(Arrays.deepToString(Hbc).replace("], ", "]\n"));
        }
    }

    public void agregate() {
        for (int i = 0; i < getNE(); i++) {
            Element element = elements.get(i);
            int[] id = element.getIdList().stream().mapToInt(val -> val).toArray();
            for (int h = 0; h < 4; h++) {
                for (int g = 0; g < 4; g++) {
                    HG[id[h] - 1][id[g] - 1] += (element.getH())[h][g] + (element.getHbc())[h][g];
                    CG[id[h] - 1][id[g] - 1] += (element.getC())[h][g];
                }
                PG[id[h] - 1] += (element.getP())[h];
            }
        }
    }

    public void calculateT() {
        RealMatrix HG = new Array2DRowRealMatrix(getHG(), false);
        RealVector PG = new ArrayRealVector(getPG(), false);
        DecompositionSolver solver = new LUDecomposition(HG).getSolver();
        RealVector result = solver.solve(PG);
        //przypisać do każdego node`a

        System.out.println(result);
        System.out.println("Calculated temperature!");
    }

    public void calculateHandP() {
        double dT = Configuration.getInstance().simulationStepTime();

        RealMatrix HGMatrix = new Array2DRowRealMatrix(getHG(), false);
        RealMatrix CGMatrix = new Array2DRowRealMatrix(getCG(), false);
        RealMatrix HNew = HGMatrix.add(CGMatrix.scalarMultiply(1 / dT));
        RealVector PGMatrix = new ArrayRealVector(getPG());

        RealMatrix t0 = new Array2DRowRealMatrix(getNodes().stream()
                                                           .mapToDouble(Node::getTemperature)
                                                           .toArray());
        RealVector PNew = PGMatrix.add(
            CGMatrix.scalarMultiply(1 / dT).multiply(t0).getColumnVector(0));

//        this.HG = HNew.getData();
//        this.PG = PNew.toArray();

        DecompositionSolver solver = new LUDecomposition(HNew).getSolver();
        RealVector result = solver.solve(PNew);
        IntStream.range(0, nN)
                 .forEach(n -> this.getNodes().get(n).setTemperature(result.getEntry(n)));
        currentTime += 50;
        printTemperatures();
    }

    public void printTemperatures() {
        System.out.println("\n\n\nTemperature after " + currentTime + " seconds:");
        double[] temperatures = this.getNodes()
                                    .stream()
                                    .mapToDouble(Node::getTemperature)
                                    .toArray();
        int i=1;
        for(double temp : temperatures) {
            System.out.printf("%.5f ", temp);
            if(i % nB == 0)
                System.out.println();
            i++;
        }
        DoubleSummaryStatistics summary = this.getNodes()
                                              .stream()
                                              .mapToDouble(Node::getTemperature)
                                              .summaryStatistics();
        System.out.printf("%nMin temp: %.5f Max temp: %.5f", summary.getMin(), summary.getMax());
    }
}