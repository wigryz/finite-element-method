package com.wigryz;

import com.wigryz.algorithms.IntegrationScheme;
import com.wigryz.structures.Element4x2D;
import com.wigryz.structures.Grid;
import com.wigryz.structures.Node;
import com.wigryz.utilities.Configuration;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.MethodOrderer;
import org.junit.jupiter.api.Order;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.TestMethodOrder;

import static org.junit.jupiter.api.Assertions.*;

@TestMethodOrder(MethodOrderer.OrderAnnotation.class)
class ApplicationTest {

    Configuration conf = Configuration.getInstance();

    double[] tempAfterOneStep =
        {365.81547, 249.01534, 249.01534, 365.81547,
            249.01534, 110.03798, 110.03798, 249.01534,
            249.01534, 110.03798, 110.03798, 249.01534,
            365.81547, 249.01534, 249.01534, 365.81547};
    double[] tempAfterTwoSteps =
        {502.59171, 353.09987, 353.09987, 502.59171,
            353.09987, 168.83702, 168.83702, 353.09987,
            353.09987, 168.83702, 168.83702, 353.09987,
            502.59171, 353.09987, 353.09987, 502.59171};

    @BeforeEach
    void setUp() {
        conf.loadConfigurationHardcoded();
    }

    @Test
    @Order(1)
    void testFor1NIntegralScheme() {
        Element4x2D element4x2D = new Element4x2D(IntegrationScheme.INTEGRATION_SCHEME_1N);

        Grid grid = new Grid(0.1, 0.1,
                             4, 4,
                             element4x2D);
        grid.printTemperatures();
        grid.iterate();
        double[] temperatures = grid.getNodes()
                                    .stream()
                                    .mapToDouble(Node::getTemperature)
                                    .toArray();
        assertTrue(compareTwoArraysWithGivenPrecision(temperatures, tempAfterOneStep, 0.0001));
        grid.iterate();
        temperatures = grid.getNodes()
                           .stream()
                           .mapToDouble(Node::getTemperature)
                           .toArray();
        assertTrue(compareTwoArraysWithGivenPrecision(temperatures, tempAfterTwoSteps, 0.0001));
    }

    @Test
    @Order(2)
    void testFor2NIntegralScheme() {

        Element4x2D element4x2D = new Element4x2D(IntegrationScheme.INTEGRATION_SCHEME_2N);

        Grid grid = new Grid(0.1, 0.1,
                             4, 4,
                             element4x2D);
        grid.printTemperatures();
        grid.iterate();
        double[] temperatures = grid.getNodes()
                                    .stream()
                                    .mapToDouble(Node::getTemperature)
                                    .toArray();
        assertTrue(compareTwoArraysWithGivenPrecision(temperatures, tempAfterOneStep, 0.0001));
        grid.iterate();
        temperatures = grid.getNodes()
                           .stream()
                           .mapToDouble(Node::getTemperature)
                           .toArray();
        assertTrue(compareTwoArraysWithGivenPrecision(temperatures, tempAfterTwoSteps, 0.0001));
    }

    @Test
    void compareN1ToN2() {
        System.out.println("1N");
        Element4x2D element4x2D1N = new Element4x2D(IntegrationScheme.INTEGRATION_SCHEME_1N);
        Grid grid1 = new Grid(0.025, 0.025,
                              2, 2,
                              element4x2D1N);

        System.out.println("2N");

        Element4x2D element4x2D2N = new Element4x2D(IntegrationScheme.INTEGRATION_SCHEME_2N);
        Grid grid2 = new Grid(0.025, 0.025,
                              2, 2,
                              element4x2D2N);
        assertTrue(compareTwoNestedArraysWithGivenPrecision(grid1.getGlobalH(), grid2.getGlobalH(),
                                                            0.001));
        assertTrue(compareTwoNestedArraysWithGivenPrecision(grid1.getGlobalC(), grid2.getGlobalC(),
                                                            0.001));
        assertTrue(compareTwoArraysWithGivenPrecision(grid1.getGlobalP(),
                                                      grid2.getGlobalP(),
                                                      0.001));
    }

    boolean compareTwoNestedArraysWithGivenPrecision(double[][] A, double[][] B, double precision) {
        for (int i = 0; i < A.length; i++) {
            for (int j = 0; j < A[0].length; j++) {
                if (Math.abs(A[i][j] - B[i][j]) > precision)
                    return false;
            }
        }
        return true;
    }


    boolean compareTwoArraysWithGivenPrecision(double[] A, double[] B, double precision) {
        for (int i = 0; i < A.length; i++) {
            if (Math.abs(A[i] - B[i]) > precision)
                return false;
        }
        return true;
    }
}