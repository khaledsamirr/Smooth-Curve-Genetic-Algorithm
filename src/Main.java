import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Random;
import java.util.Scanner;
import java.util.concurrent.TimeUnit;

public class Main {
    public static void readFileItems(Scanner sc, int NumOfPoints, ArrayList<Point> points) {
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
            parameters.setNumberOfSets(NumOfPoints);
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

        setCoeffPrecision(chromosome, 4);
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

    public static <T> double setPrecisionNum(T num, int numAfterPoint) {
        Double dNum = (double) num;
        dNum = dNum * Math.pow(10, numAfterPoint >= 1 ? numAfterPoint : 0);
        dNum = Math.floor(dNum);
        dNum = dNum / Math.pow(10, numAfterPoint >= 1 ? numAfterPoint : 0);
        return dNum;
    }

    public static void setCoeffPrecision(Chromosome chromosome, int numAfterPoint) {
        for (int i = 0; i < chromosome.genes.size(); i++) {
            chromosome.genes.get(i).setCoefficient(setPrecisionNum(chromosome.genes.get(i).getCoefficient(), numAfterPoint));
        }
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
                randDoubleNum = setPrecisionNum(randDoubleNum, 4);

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
    public static ArrayList<Chromosome> tournamentSelection(ArrayList<Chromosome> population, Parameters parameters, ArrayList<Point> points) {
        ArrayList<Chromosome> matingPool = new ArrayList<>();
        double fit1 = 0.0, fit2 = 0.0;
        for (int i = 0; i < population.size(); i += 2) {
            fit1 = population.get(i).getFitness();
            fit2 = population.get(i + 1).getFitness();
            if (fit1 < fit2) {
                matingPool.add(population.get(i));
            } else {
                matingPool.add(population.get(i + 1));
            }
        }
        return matingPool;
    }

    public static Chromosome getBestChromosome(ArrayList<Chromosome> population) {
        double bestFitness = 0.0;
        Chromosome bestChromosome = new Chromosome();

        for (int i = 0; i < population.size(); i++) {
            if (i == 0) {
                bestFitness = population.get(i).getFitness();
                bestChromosome = population.get(i);
            }
            if (population.get(i).getFitness() < bestFitness) {
                bestFitness = population.get(i).getFitness();
                bestChromosome = population.get(i);
            }
        }

        return bestChromosome;
    }

    public static void setSuitablePopulationSize(Parameters parameters) {
        if (!parameters.isFixedPopulationSize())
            parameters.setPopulationSize((int) Math.pow(2, parameters.getChromosomeSize() - 1));
    }

    public static void solve(Parameters parameters, Scanner sc, ArrayList<Point> points, ArrayList<Chromosome> population) throws IOException {
        StopWatch sw = new StopWatch();
        StopWatch swWholeTime = new StopWatch();

        swWholeTime.start();

        FileWriter fileWriter = new FileWriter("Output.txt");
        int generationRepetitionCounter = 0;
        String bestChromosomeFitness = " ";
        Chromosome TCsChromosome = new Chromosome();
        ArrayList<Integer> bestChromosomes = new ArrayList<>();

        for (int q = 0; q < parameters.getTCsNum(); q++) {
            readFile(sc, parameters, points);
            System.out.println("Polynomial Degree: " + parameters.getPolynomialDegree());

            setSuitablePopulationSize(parameters);

            System.out.println("===============================================");
            System.out.println("TC number: " + (q + 1));
            System.out.println("------------------------");

            sw.reset();
            sw.start();
            initializePopulation(parameters, points, population);
            sw.stop();
            System.out.println(sw.getElapsedTime().toMillis() + "ms on Initialization");
            System.out.println(sw.getElapsedTime().toSeconds() + "s on Initialization");
            System.out.println(sw.getElapsedTime().toMinutes() + "m on Initialization");
            System.out.println(sw.getElapsedTime().toHours() + "hr on Initialization");
            System.out.println();

            sw.reset();
            sw.start();

            for (int i = 0; i < population.size(); i++) {
                population.get(i).calcFitness(parameters, points);
            }

            ArrayList<Chromosome> selectedChromosomes = new ArrayList<>();

            for (int i = 0; i < parameters.getNumOfGenerations(); i++) {

                selectedChromosomes = tournamentSelection(population, parameters, points);

                if (selectedChromosomes.size() >= 2) {

                    for (int j = 0; j < selectedChromosomes.size(); j += 2) {

                        crossover(selectedChromosomes.get(j), selectedChromosomes.get(j + 1), i + 1, parameters, points, population);

                    }
                }

                TCsChromosome = getBestChromosome(population);

                /*if (bestChromosomeFitness.equals(String.valueOf(getBestChromosome(population).getFitness()))) {
                    generationRepetitionCounter++;
                } else {
                    generationRepetitionCounter = 0;
                    TCsChromosome = getBestChromosome(population);
                    bestChromosomeFitness = String.valueOf(TCsChromosome.getFitness());
                }*/

                if (generationRepetitionCounter > 2) {
                    i = parameters.getNumOfGenerations() + 1;
                } else {
                    System.out.println("Generation:" + (i + 1) + " Finished of TC Num: " + (q + 1));
                }
            }
            String chromosomeSeq = "";

            setCoeffPrecision(TCsChromosome, 2);
            chromosomeSeq = TCsChromosome.toString();
            double err = setPrecisionNum(TCsChromosome.getFitness(), 4);
            fileWriter.write(chromosomeSeq.substring(1, chromosomeSeq.length() - 1) + ", Error = " + err + "\n");


            System.out.println(bestChromosomeFitness);
            //bestChromosomes.add(Integer.valueOf(bestChromosomeFitness));

            points.clear();
            population.clear();

            sw.stop();
            System.out.println(sw.getElapsedTime().toMillis() + "ms TC: " + (q + 1));
            System.out.println(sw.getElapsedTime().toSeconds() + "s TC: " + (q + 1));
            System.out.println(sw.getElapsedTime().toMinutes() + "m TC: " + (q + 1));
            System.out.println(sw.getElapsedTime().toHours() + "hr TC: " + (q + 1));

            System.out.println("===========================================");

        }

        fileWriter.close();
        swWholeTime.stop();
        System.out.println(swWholeTime.getElapsedTime().toMillis() + "ms");
        System.out.println(swWholeTime.getElapsedTime().toSeconds() + "s");
        System.out.println(swWholeTime.getElapsedTime().toMinutes() + "m");
        System.out.println(swWholeTime.getElapsedTime().toHours() + "hr");

    }


    public static void main(String[] args) throws IOException {

        Parameters parameters = new Parameters(1024, true, 1024, 2, 0.05, 0.6,
                -10, 10, 1, 2);
        String readFilePath = "D:\\College\\Soft Computing\\Assignment\\Assignment #2\\input-2.txt";
        File file = new File(readFilePath);
        Scanner sc = new Scanner(file);

        ArrayList<Chromosome> population = new ArrayList<>();
        ArrayList<Point> points = new ArrayList<>();

        if (sc.hasNextLine()) {
            parameters.setTCsNum(Integer.valueOf(sc.nextLine()));
            solve(parameters, sc, points, population);
        }

    }
}