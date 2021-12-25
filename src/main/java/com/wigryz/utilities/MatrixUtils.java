package com.wigryz.utilities;

import java.util.Arrays;

public class MatrixUtils {

    public static String matrixToString(double[][] matrix) {
        return Arrays.deepToString(matrix).replace("], ", "]\n");
    }
}
