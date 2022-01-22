package com.wigryz;

import com.wigryz.structures.Grid;
import com.wigryz.utilities.Configuration;

public class Application {

    public static void main(String[] args) {

        Configuration conf = Configuration.getInstance();
//        conf.loadConfigurationFromFile();
//        conf.loadConfigurationHardcoded();
//        conf.loadInitDataFromFile("Test1_4_4.txt");
//        conf.loadInitDataFromFile("Test2_4_4_MixGrid.txt");
//        conf.loadInitDataFromFile("Test3_31_31_kwadrat.txt");
        conf.loadInitDataFromFile("Test4_31_31_trapez.txt");

//        Grid grid = new Grid(conf, true);
        Grid grid = new Grid(conf);
        grid.printTemperatures();
        while(grid.getCurrentTime() < conf.simulationTime()) {
            grid.iterate();
        }
    }
}
