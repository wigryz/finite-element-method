package com.wigryz.utilities;

import com.wigryz.algorithms.IntegrationScheme;
import com.wigryz.structures.Element;
import com.wigryz.structures.Node;

import java.io.File;
import java.io.FileNotFoundException;
import java.net.URISyntaxException;
import java.net.URL;
import java.util.ArrayList;
import java.util.List;
import java.util.NoSuchElementException;
import java.util.Objects;
import java.util.Scanner;

public class Configuration {
    private static Configuration instance;

    private double initialTemperature;      // temperatura poczatkowa w wezlach
    private int simulationTime;             // czas symulacji [s]
    private int simulationStepTime;         // skok czasu w symulacji [s]
    private double ambientTemperature;      // temperatura otoczenia
    private double alfa;                    // wspolczynnik konwekcyjnej wymiany ciepla
    private double heightOfGrid;
    private double widthOfGrid;
    private int numberOfNodesOnHeight;
    private int numberOfNodesOnWidth;
    private double specificHeat;            // cieplo wlasciwe
    private double conductivity;            // wspolczynnik przewodzenia ciepla|przewodnosc cieplna
    private double density;                 // gestosc
    private IntegrationScheme integrationScheme;

    //init data
    private int numberOfNodes;
    private int numberOfElements;
    private List<Node> nodes = new ArrayList<>();
    private List<Element> elements = new ArrayList<>();

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
            in.next();
            integrationScheme = in.nextInt() == 1 ? IntegrationScheme.INTEGRATION_SCHEME_1N :
                IntegrationScheme.INTEGRATION_SCHEME_2N;
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

    public void loadConfigurationHardcoded() {
        initialTemperature = 100.0;
        simulationTime = 500;
        simulationStepTime = 50;
        ambientTemperature = 1200.0;
        alfa = 300.0;
        heightOfGrid = 0.1;
        widthOfGrid = 0.1;
        numberOfNodesOnHeight = 4;
        numberOfNodesOnWidth = 4;
        specificHeat = 700.0;
        conductivity = 25.0;
        density = 7800.0;
        integrationScheme = IntegrationScheme.INTEGRATION_SCHEME_1N;
    }

    public void loadInitDataFromFile(String path) {
        URL url = getClass().getClassLoader().getResource(path);
        try (Scanner in = new Scanner(new File(Objects.requireNonNull(url).toURI()))) {
            in.next();
            simulationTime = in.nextInt();
            in.next();
            simulationStepTime = in.nextInt();
            in.next();
            conductivity = in.nextDouble();
            in.next();
            alfa = in.nextDouble();
            in.next();
            ambientTemperature = in.nextDouble();
            in.next();
            initialTemperature = in.nextDouble();
            in.next();
            density = in.nextDouble();
            in.next();
            specificHeat = in.nextDouble();
            in.next();
            in.next();
            numberOfNodes = in.nextInt();
            in.next();
            in.next();
            numberOfElements = in.nextInt();
            in.next();
            integrationScheme = in.nextInt() == 1 ? IntegrationScheme.INTEGRATION_SCHEME_1N :
                IntegrationScheme.INTEGRATION_SCHEME_2N;
            in.nextLine();
            in.nextLine();
            readNodes(numberOfNodes, in);
            in.nextLine();
            readElements(numberOfElements, in);
            in.nextLine();
            readBoundaryConditions(in);
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

    private void readNodes(int numberOfNodes, Scanner in) {
        short bc = 0;
        for(int i=0 ; i < numberOfNodes ; i++) {
            String[] line = in.nextLine().split(",");
            int index = Integer.parseInt(line[0].strip());
            double x = Double.parseDouble(line[1].strip());
            double y = Double.parseDouble(line[2].strip());
            nodes.add(new Node(index, x, y, bc, initialTemperature));
        }
    }

    private void readElements(int numberOfElements, Scanner in) {
        for(int i=0 ; i < numberOfElements ; i++) {
            String[] line = in.nextLine().split(",");
            int index = Integer.parseInt(line[0].strip());
            int id1 = Integer.parseInt(line[1].strip());
            int id2 = Integer.parseInt(line[2].strip());
            int id3 = Integer.parseInt(line[3].strip());
            int id4 = Integer.parseInt(line[4].strip());
            if(i == 0)
                numberOfNodesOnWidth = id4 - 1;
            elements.add(new Element(index, id1, id2, id3, id4));
        }
    }

    private void readBoundaryConditions(Scanner in) {
        short bc = 1;
        while(in.hasNext()) {
            String[] line = in.next().split(",");
            int nodeIndex = Integer.parseInt(line[0].strip());
            nodes.get(nodeIndex - 1).setBoundaryCondition(bc);
        }
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

    public IntegrationScheme integrationScheme() {
        return integrationScheme;
    }

    public int numberOfNodes() {
        return numberOfNodes;
    }

    public int numberOfElements() {
        return numberOfElements;
    }

    public List<Node> nodes() {
        return nodes;
    }

    public List<Element> elements() {
        return elements;
    }

}
