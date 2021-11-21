import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Random;
import java.util.Scanner;

public class Main {
    public static void readFileItems(Scanner sc, int NumOfPoints, ArrayList<Point> points) {
        sc.nextLine();
        for (int i = 0; i < NumOfPoints; i++) {
            String[] splitLine;
            splitLine = sc.nextLine().split(" ");
            points.add(new Point(Double.valueOf(splitLine[0]), Double.valueOf(splitLine[1])));
        }

    }


    public static void readFile(Scanner sc, Parameters parameters, ArrayList<Point> points) {
        int NumOfPoints = 0;
        String[] splitLine;
        if (sc.hasNextLine()) {
            splitLine = sc.nextLine().split(" ");
            NumOfPoints = Integer.valueOf(splitLine[0]);
            parameters.setPolynomialDegree(Integer.valueOf(splitLine[1]));
            readFileItems(sc, NumOfPoints, points);
        }
    }

    public static Gene nonUniformMutation(Gene gene, Parameters parameters, int currentGeneration) {
        Random rand = new Random();
        boolean lower = true;
        double y = 0.0;
        double delta = 0.0;
        if (rand.nextDouble() <= 0.5)
            y = gene.getCoefficient() - (parameters.getLB());
        else {
            y = (parameters.getUB()) - gene.getCoefficient();
            lower = false;
        }
        delta = y * (1 - Math.pow(rand.nextDouble(), Math.pow((1 - (currentGeneration / parameters.getNumOfGenerations())), parameters.getB())));
        if (lower)
            gene.setCoefficient(gene.getCoefficient() - delta);
        else
            gene.setCoefficient(gene.getCoefficient() + delta);

        return gene;

    }

    public static void mutation(Chromosome chromosome, Parameters parameters, ArrayList<Point> points, int currentGeneration) {
        Random rand = new Random();
        int index = 0;
        for (Gene i : chromosome.genes) {
            if (rand.nextDouble() <= parameters.getMutationRate()) {
                i = nonUniformMutation(chromosome.genes.get(index), parameters, currentGeneration);
            }
            index++;

        }

        chromosome.calcFitness(parameters, points);
    }

    public static Boolean chromosomeExists(Chromosome chromosome, ArrayList<Chromosome> population) {
        for (Chromosome i : population) {
            if (chromosome.toString().equals(i.toString())) {
                return true;
            }
        }
        return false;
    }

    public static void initializePopulation(Parameters parameters, ArrayList<Point> points, ArrayList<Chromosome> population) {

        Random rand = new Random();
        Chromosome chromosome;
        ArrayList<Double> coefficients = new ArrayList<Double>();
        int randWholeNum = 0;
        double randDoubleNum = 0.0;


        for (int j = 0; j < parameters.getPopulationSize(); j++) {

            for (int k = 0; k <= parameters.getPolynomialDegree(); k++) {
                randWholeNum = rand.nextInt(parameters.getUB());
                randDoubleNum = rand.nextFloat();
                randDoubleNum = randDoubleNum * Math.pow(10, 2);
                randDoubleNum = Math.floor(randDoubleNum);
                randDoubleNum = randDoubleNum / Math.pow(10, 2);

                coefficients.add(((double) randWholeNum + (randWholeNum != 10 ? randDoubleNum : 0)) * (rand.nextBoolean() ? 1 : -1));
            }

            chromosome = new Chromosome(coefficients);
            coefficients.clear();

            if (!chromosomeExists(chromosome, population)) {
                //chromosome.calcFitness(parameters, points);
                population.add(chromosome);
            } else {
                j--;
            }
        }

    }

    public static boolean crossover(Chromosome chromosome1, Chromosome chromosome2, int currentGen, Parameters parameters, ArrayList<Point> points, ArrayList<Chromosome> population) {
        Random rand = new Random();
        ArrayList<Integer> splitPoints = new ArrayList<>();
        int randNum = 0;
        boolean crossoverFlag = false;

        if (!(parameters.getNumOfSplitPoints() <= parameters.getPolynomialDegree()))
            parameters.setNumOfSplitPoints(1);


        for (int i = 0; i < parameters.getNumOfSplitPoints(); i++) {
            randNum = rand.nextInt(chromosome1.genes.size()) + 1;
            if (randNum == chromosome1.genes.size())
                randNum--;

            if (!splitPoints.contains(randNum))
                splitPoints.add(randNum);
            else {
                i--;
            }
        }

        Collections.sort(splitPoints);

        for (int i : splitPoints) {
            System.out.println(i);
        }

        Chromosome newChromo1 = new Chromosome();
        Chromosome newChromo2 = new Chromosome();
        int j = 0;

        Boolean crossoverSwitch = false;

        for (int i = 0; i < chromosome1.genes.size(); i++) {
            if (j < splitPoints.size() && i >= splitPoints.get(j)) {
                crossoverSwitch = !crossoverSwitch;
                j++;
            }
            if (!crossoverSwitch) {
                newChromo1.genes.add(chromosome1.genes.get(i));
                newChromo2.genes.add(chromosome2.genes.get(i));
            } else {
                newChromo1.genes.add(chromosome2.genes.get(i));
                newChromo2.genes.add(chromosome1.genes.get(i));
            }
        }

        mutation(newChromo1, parameters, points, currentGen);
        mutation(newChromo2, parameters, points, currentGen);
        if (!chromosomeExists(newChromo1, population)) {
            population.add(newChromo1);
            population.remove(chromosome1);
            population.remove(chromosome2);
            crossoverFlag = true;
        }
        if (!chromosomeExists(newChromo2, population)) {
            population.add(newChromo2);
            population.remove(chromosome1);
            population.remove(chromosome2);
            crossoverFlag = true;
        }

        return crossoverFlag;

    }

    //TODO: tournamentSelection
    /*public static ArrayList<Chromosome> tournamentSelection() {

        ArrayList<Chromosome> selected = new ArrayList<>();
    }*/

    public static void setSuitablePopulationSize(Parameters parameters) {
        if (parameters.isFixedPopulationSize() == false) {
            if (parameters.getChromosomeSize() < 9)
                parameters.setPopulationSize((int) Math.pow(2, parameters.getChromosomeSize() - 1));
            else
                parameters.setPopulationSize(750);
        }
    }

    public static void solve(Parameters parameters, Scanner sc, ArrayList<Point> points, ArrayList<Chromosome> population) throws IOException {
        // TODO: 11/17/2021   int currentGeneration=0;

        FileWriter fileWriter = new FileWriter("Output.txt");
        int generationRepetitionCounter = 0;
        String bestChromosome = " ";
        for (int q = 0; q < parameters.getTCsNum(); q++) {
            readFile(sc, parameters, points);

            parameters.setChromosomeSize(points.size());
            setSuitablePopulationSize(parameters);

            System.out.println("===============================================");
            System.out.println("TC number: " + (q + 1));
            System.out.println("Population size: " + parameters.getPopulationSize());
            System.out.println("------------------------");

            initializePopulation(parameters, points, population);

            ArrayList<Chromosome> selectedChromosomes = new ArrayList<>();

            for (int i = 0; i < parameters.getNumOfGenerations(); i++) {
                System.out.println("At generation:" + (i + 1));
                selectedChromosomes = selectChromosomes(parameters, population);

                if (selectedChromosomes.size() >= 2) {
                    System.out.println("Number of New Chromosomes: " + (selectedChromosomes.size()));
                    for (int j = 0; j < selectedChromosomes.size(); j += 2) {
                        crossover(selectedChromosomes.get(j), selectedChromosomes.get(j + 1), parameters, population);
                    }
                }
                // for (Chromosome k : population) {
                //  System.out.println(k.toString());
                //}
                if (bestChromosome.equals(getBestFitness(population))) {
                    generationRepetitionCounter++;
                } else {
                    generationRepetitionCounter = 0;
                    bestChromosome = getBestFitness(population);
                }

                if (generationRepetitionCounter > 2) {
                    i = parameters.getNumOfGenerations() + 1;
                }

            }

            fileWriter.write("Case: " + (q + 1) + " " + bestChromosome + "\n");
            System.out.println(bestChromosome);
            System.out.println("===========================================");

            points.clear();
            population.clear();

        }

        fileWriter.close();
    }


    public static void main(String[] args) throws IOException {

        Parameters parameters = new Parameters(2, false, 9, 2, 0.00, 0.5,
                -10, 10, 1, 2);
        String readFilePath = "";
        File file = new File(readFilePath);
        //Scanner sc = new Scanner(file);

        ArrayList<Chromosome> population = new ArrayList<>();
        ArrayList<Point> points = new ArrayList<>();

        initializePopulation(parameters, points, population);

        for (Chromosome i : population) {
            System.out.println(i.toString());
        }
        System.out.println("--------------------");
        crossover(population.get(0), population.get(1), 1, parameters, points, population);

        for (Chromosome i : population) {
            System.out.println(i.toString());
        }

        /*if (sc.hasNextLine()) {
            parameters.setNumberOfSets(sc.nextInt());
            //  solve(parameters, sc, points, population);
        }*/


    }
}
