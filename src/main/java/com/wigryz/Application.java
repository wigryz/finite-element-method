package com.wigryz;

import com.wigryz.structures.Element4x2D;
import com.wigryz.structures.Grid;
import com.wigryz.utilities.Configuration;

public class Application {

    public static void main(String[] args) {

        Configuration conf = Configuration.getInstance();
        conf.loadConfigurationFromFile();

        Element4x2D element4x2D = new Element4x2D(conf.integrationScheme());
        Grid grid = new Grid(conf.heightOfGrid(), conf.widthOfGrid(),
                             conf.numberOfNodesOnHeight(), conf.numberOfNodesOnWidth(),
                             element4x2D);
        grid.printTemperatures();
        for(int i=0 ; i < 10 ; i++) {
            grid.iterate();
        }
    }
}
