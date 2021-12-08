package com.wigryz.utilities;

import java.io.File;
import java.io.FileNotFoundException;
import java.net.URISyntaxException;
import java.net.URL;
import java.util.NoSuchElementException;
import java.util.Objects;
import java.util.Scanner;

public class Configuration {
    private static Configuration instance;

    private double initialTemperature;
    private int simulationTime;             // [s]
    private int simulationStepTime;         // [s]
    private double ambientTemperature;
    private double alfa;
    private double heightOfGrid;
    private double widthOfGrid;
    private int numberOfNodesOnHeight;
    private int numberOfNodesOnWidth;
    private double specificHeat;            //unused
    private double conductivity;
    private double density;                 //unused

    private Configuration() {
    }

    public static synchronized Configuration getInstance() {
        if (instance == null) {
            instance = new Configuration();
        }
        return instance;
    }

    public void loadConfigurationFromFile() {
        URL url = getClass().getClassLoader().getResource("init.txt");
        try (Scanner in = new Scanner(new File(Objects.requireNonNull(url).toURI()))) {
            in.next();
            initialTemperature = in.nextDouble();
            in.next();
            simulationTime = in.nextInt();
            in.next();
            simulationStepTime = in.nextInt();
            in.next();
            ambientTemperature = in.nextDouble();
            in.next();
            alfa = in.nextDouble();
            in.next();
            heightOfGrid = in.nextDouble();
            in.next();
            widthOfGrid = in.nextDouble();
            in.next();
            numberOfNodesOnHeight = in.nextInt();
            in.next();
            numberOfNodesOnWidth = in.nextInt();
            in.next();
            specificHeat = in.nextDouble();
            in.next();
            conductivity = in.nextDouble();
            in.next();
            density = in.nextDouble();
        } catch (FileNotFoundException e) {
            System.out.println("Failed during opening init file.");
            System.exit(1);
        } catch (NoSuchElementException e) {
            System.out.println("Failed during reading init file.");
            System.exit(1);
        } catch (URISyntaxException e) {
            e.printStackTrace();
        }
    }

    public void initializeConfiguration() {
        initialTemperature = 100.0;
        simulationTime = 500;
        simulationStepTime = 50;
        ambientTemperature = 1200.0;
        alfa = 300.0;
        heightOfGrid = 0.100;
        widthOfGrid = 0.100;
        numberOfNodesOnHeight = 4;
        numberOfNodesOnWidth = 4;
        specificHeat = 700.0;
        conductivity = 25.0;
        density = 7800.0;
    }

    public double initialTemperature() {
        return initialTemperature;
    }

    public int simulationTime() {
        return simulationTime;
    }

    public int simulationStepTime() {
        return simulationStepTime;
    }

    public double ambientTemperature() {
        return ambientTemperature;
    }

    public double alfa() {
        return alfa;
    }

    public double heightOfGrid() {
        return heightOfGrid;
    }

    public double widthOfGrid() {
        return widthOfGrid;
    }

    public int numberOfNodesOnHeight() {
        return numberOfNodesOnHeight;
    }

    public int numberOfNodesOnWidth() {
        return numberOfNodesOnWidth;
    }

    public double specificHeat() {
        return specificHeat;
    }

    public double conductivity() {
        return conductivity;
    }

    public double density() {
        return density;
    }
}
