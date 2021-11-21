import java.util.ArrayList;

public class Chromosome {

    double fitness = 0.0;
    double error = 0.0;
    ArrayList<Gene> genes = new ArrayList<>();

    Chromosome() {
    }

    public Chromosome(ArrayList<Double> coefficients) {
        for (Double i : coefficients) {
            genes.add(new Gene(i));
        }
    }

    public double getError() {
        return error;
    }

    public void setError(double error) {
        this.error = error;
    }

    /*public void addGene(Gene gene) {
        genes.add(gene);
    }*/

    public void print() {
        String sequence = "[";
        for (int i = 0; i < genes.size(); i++) {
            sequence += genes.get(i).getCoefficient();
            if (i + 1 != genes.size())
                sequence += ",";
        }
        sequence += "]";

        System.out.println(sequence);
    }

    public void calcFitness(Parameters parameters, ArrayList<Point> points) {
        double fitness = 0.0;
        double err = 0.0;
        double Ycalc = 0.0;
        for (int i = 0; i < parameters.getNumberOfSets(); i++) {
            for (int j = 0; j < genes.size(); j++) {
                Ycalc += genes.get(j).getCoefficient() * (Math.pow(points.get(i).getX(), j));
            }
            err = Math.pow(Ycalc - points.get(i).getY(), 2);
            fitness += err;
        }
        this.fitness = fitness / parameters.getNumberOfSets();
    }

    public String toString() {
        String sequenceActivated = "[";

        for (int i = 0; i < genes.size(); i++) {
            sequenceActivated += genes.get(i).toString();
            if (i < genes.size() - 1) {
                sequenceActivated += ",";
            }
        }
        sequenceActivated += "]";
        return sequenceActivated;
    }

    public double getFitness() {
        return fitness;
    }
}
