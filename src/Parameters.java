public class Parameters {
    private int TCsNum = 0;
    private int LB = 0;
    private int UB = 0;
    private double b = 0;
    private int NumberOfSets = 0;
    private int populationSize = 10;
    private boolean fixedPopulationSize = true;
    private double mutationRate = 0.05;
    private double crossoverRate = 0.9;
    private int numOfSplitPoints = 1;
    private int chromosomeSize = 2;
    private int numOfGenerations = 0;
    private int polynomialDegree = 1;

    public Parameters(int populationSize, boolean fixedPopulationSize, int numOfGenerations, int numOfSplitPoints,
                      double mutationRate, double crossoverRate, int lowerBound, int upperBound, double b, int polynomialDegree) {
        this.numOfSplitPoints = numOfSplitPoints;
        this.numOfGenerations = numOfGenerations;
        this.fixedPopulationSize = fixedPopulationSize;
        this.populationSize = populationSize;
        this.mutationRate = mutationRate;
        this.crossoverRate = crossoverRate;
        this.LB = lowerBound;
        this.UB = upperBound;
        this.b = b;
        this.polynomialDegree = polynomialDegree;
    }

    public int getTCsNum() {
        return TCsNum;
    }

    public void setTCsNum(int TCsNum) {
        this.TCsNum = TCsNum;
    }

    public double getB() {
        return b;
    }

    public void setB(double b) {
        this.b = b;
    }

    public int getNumberOfSets() {
        return NumberOfSets;
    }

    public void setNumberOfSets(int numberOfSets) {
        NumberOfSets = numberOfSets;
    }

    public int getLB() {
        return LB;
    }

    public void setLB(int LB) {
        this.LB = LB;
    }

    public int getUB() {
        return UB;
    }

    public void setUB(int UB) {
        this.UB = UB;
    }

    public int getPolynomialDegree() {
        return polynomialDegree;
    }

    public void setPolynomialDegree(int polynomialDegree) {
        this.polynomialDegree = polynomialDegree;
    }

    public int getPopulationSize() {
        return populationSize;
    }

    public void setPopulationSize(int populationSize) {
        this.populationSize = populationSize;
    }

    public double getMutationRate() {
        return mutationRate;
    }

    public void setMutationRate(double mutationRate) {
        this.mutationRate = mutationRate;
    }

    public double getCrossoverRate() {
        return crossoverRate;
    }

    public void setCrossoverRate(double crossoverRate) {
        this.crossoverRate = crossoverRate;
    }

    public int getNumOfSplitPoints() {
        return numOfSplitPoints;
    }

    public void setNumOfSplitPoints(int numOfSplitPoints) {
        if (numOfSplitPoints < chromosomeSize)
            this.numOfSplitPoints = numOfSplitPoints;
        else
            this.numOfSplitPoints = 1;
    }

    public int getChromosomeSize() {
        return chromosomeSize;
    }

    public void setChromosomeSize(int chromosomeSize) {
        this.chromosomeSize = chromosomeSize;
    }

    public int getNumOfGenerations() {
        return numOfGenerations;
    }

    public void setNumOfGenerations(int numOfGenerations) {
        this.numOfGenerations = numOfGenerations;
    }

    public boolean isFixedPopulationSize() {
        return fixedPopulationSize;
    }

    public void setFixedPopulationSize(boolean fixedPopulationSize) {
        this.fixedPopulationSize = fixedPopulationSize;
    }
}
