import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
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
            splitLine=sc.nextLine().split(" ");
            NumOfPoints = Integer.valueOf(splitLine[0]);
            parameters.setPolynomialDegree(Integer.valueOf(splitLine[1]));
            readFileItems(sc, NumOfPoints, points);
        }
    }
    public static Gene nonUniformMutation(Gene gene, Parameters parameters,int currentGeneration) {
        Random rand = new Random();
        boolean lower=true;
        double y=0.0;
        double delta=0.0;
        if(rand.nextDouble()<=0.5)
            y=gene.getCoefficient()-(parameters.getLB());
        else {
            y = (parameters.getUB()) - gene.getCoefficient();
            lower=false;
        }
        delta=y*(1-Math.pow(rand.nextDouble(),Math.pow((1-(currentGeneration/parameters.getNumOfGenerations())),parameters.getB())));
        if(lower)
            gene.setCoefficient(gene.getCoefficient()-delta);
        else
            gene.setCoefficient(gene.getCoefficient()+delta);

        return gene;

    }
    public static void mutation(Chromosome chromosome, Parameters parameters,ArrayList<Point>points,int currentGeneration) {
        Random rand = new Random();
        int index=0;
        for (Gene i : chromosome.genes) {
            if (rand.nextDouble() <= parameters.getMutationRate()) {
                i=nonUniformMutation(chromosome.genes.get(index),parameters,currentGeneration);
            }
            index++;

        }

        chromosome.calcFitness(parameters,points);
    }

    /*
    public static void solve(Parameters parameters, Scanner sc, ArrayList<Item> items, ArrayList<Chromosome> population) throws IOException {
        // TODO: 11/17/2021   int currentGeneration=0;

        FileWriter fileWriter = new FileWriter("Output.txt");
        int generationRepetitionCounter = 0;
        String bestChromosome = " ";
        for (int q = 0; q < parameters.getTCsNum(); q++) {
            readFile(sc, parameters, items);

            parameters.setChromosomeSize(items.size());
            setSuitablePopulationSize(parameters);

            System.out.println("===============================================");
            System.out.println("TC number: " + (q + 1));
            System.out.println("Population size: " + parameters.getPopulationSize());
            System.out.println("------------------------");

            initializePopulation(parameters, items, population);

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

            items.clear();
            population.clear();

        }

        fileWriter.close();
    }
    */

    public static void main(String[] args) throws IOException {

        Parameters parameters = new Parameters(24, false, 9, 3, 0.05, 0.5,
                -10,10,1);
        String readFilePath = "";
        File file = new File(readFilePath);
        Scanner sc = new Scanner(file);

        ArrayList<Chromosome> population = new ArrayList<>();
        ArrayList<Point> points = new ArrayList<>();

        if (sc.hasNextLine()) {
            parameters.setNumberOfSets(sc.nextInt());
          //  solve(parameters, sc, points, population);
        }


    }
}
