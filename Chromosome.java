package assignment2;

import java.util.ArrayList;


public class Chromosome  {
    ArrayList<Float> genes  = new ArrayList<>();
    float fitness;
    boolean discarded =false;
    boolean selcted =false;


    public Chromosome(){
        discarded=false;
        selcted =false;
    }
    public Chromosome(ArrayList<Float> genes, float fitness, boolean discarded, boolean selcted) {
        this.genes = genes;
        this.fitness = fitness;
        this.discarded=discarded;
        this.selcted=selcted;

    }

}

