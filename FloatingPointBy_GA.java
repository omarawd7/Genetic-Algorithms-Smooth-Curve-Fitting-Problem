package assignment2;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Random;

public class FloatingPointBy_GA {

    static final Double pc = 0.6; // probability of crossover
    static final Double pm = 0.5; // probability of Mutation
    static Random rand = new Random(); // instance of random class

    static float get_error_number(Chromosome[] c, point[] p, int number_of_points, int num_of_bits,
            int population_size, PrintWriter myWriter) throws IOException {
        float C_error = 0;
        float square_error;
        int ind_of_best_fitness = 0;

        for (int k = 0; k < population_size; k++) {
            if (c[k].fitness > c[ind_of_best_fitness].fitness) {
                ind_of_best_fitness = k;
            }
        }

        myWriter.println("coefficients: ");
        for (int i = 0; i < num_of_bits; i++) {
            myWriter.println(c[ind_of_best_fitness].genes.get(i));
        }
    //    myWriter.println( "best fitness: "+ c[ind_of_best_fitness].fitness);

        square_error = (1 / c[ind_of_best_fitness].fitness);
        return square_error;
    }
//elitist replacement.
    static Chromosome[] replace(Chromosome[] c, Chromosome[] c1, int populationSize) {
        Chromosome cMin = new Chromosome();
        Chromosome cTemp = new Chromosome();
        int index_min = 0;
       
        for (int i = 0; i < populationSize; i++) {
            cTemp = c[i];
            if ( i==0 || cMin.fitness > cTemp.fitness) {
                cMin = cTemp;
                index_min = i;
            }

        }
        if (c[index_min].fitness < c1[0].fitness) {
            c[index_min] = c1[0];
        }
        // replace the second Chromosome
        index_min = 0;
        for (int i = 0; i < populationSize; i++) {
            cTemp = c[i];
            if (i==0 || cMin.fitness > cTemp.fitness ) {
                cMin = cTemp;
                index_min = i;
            }
        }

        if (c[index_min].fitness < c1[1].fitness) {
            c[index_min] = c1[1];
        }
        return c;
    }
//non-uniform mutation
    static Chromosome mutation(Chromosome c, int bits_num, int t,
            int T) {
        for (int i = 0; i < bits_num; i++) {

            Random random = new Random();
            float rm = random.nextFloat(0, 1);
            float y;
            float U_delta = 0;
            float L_delta = 0;
            float TY_delta = 0;
            float degree_of_non_uniformity = 1;
            if (rm <= pm) {
                L_delta = c.genes.get(i) - (-10);
                U_delta = 10 - c.genes.get(i);
                Random random1 = new Random();
                float r1 = random1.nextFloat(0, 1);
                if (r1 <= 0.5) {
                    y = L_delta;
                    Random random2 = new Random();
                    float r = random2.nextFloat(0, 1);
                    TY_delta = y * (float) (Math.pow((1 - r), Math.pow((1 - (t + 1)) / T, degree_of_non_uniformity)));
                    c.genes.set(i, c.genes.get(i) - TY_delta);
                } else {
                    y = U_delta;
                    Random random3 = new Random();
                    float r = random3.nextFloat(0, 1);
                    TY_delta = y * (float) (Math.pow((1 - r), Math.pow((1 - t) / T, degree_of_non_uniformity)));
                    c.genes.set(i, c.genes.get(i) + TY_delta);
                }

            }

        }

        return c;
    }
//2-point crossover
    static Chromosome[] crossOver(Chromosome[] c, int first_candidates_index, int sec_candidates_index,
            int bitsNum) {

        Double r3 = rand.nextDouble(0, 1);
        Chromosome[] c1 = new Chromosome[2];
        c1[0] = new Chromosome(c[first_candidates_index].genes, c[first_candidates_index].fitness,
                c[first_candidates_index].discarded, c[first_candidates_index].selcted);
        c1[1] = new Chromosome(c[sec_candidates_index].genes, c[sec_candidates_index].fitness,
                c[sec_candidates_index].discarded, c[sec_candidates_index].selcted);

        if (r3 > pc) {
            return c1;
        }
        int r1 = rand.nextInt(bitsNum - 1);
        int r2 = rand.nextInt(bitsNum - 1);

        if (r1 > r2) {
            for (int i = r1; i < r2; i++) {
                float tmp = c1[0].genes.get(i);
                c1[0].genes.set(i, c1[1].genes.get(i));
                c1[1].genes.set(i, tmp);
            }
        } else {
            for (int i = r2; i < r1; i++) {
                float tmp = c1[0].genes.get(i);
                c1[0].genes.set(i, c1[1].genes.get(i));
                c1[1].genes.set(i, tmp);
            }
        }

        return c1;
    }
   // tournament selection
    static Chromosome[] tournament_selection(Chromosome[] c, int population_size, int k, int bits_num, int t, int T,point[] p, int number_of_points ) {
        ArrayList<Chromosome> candidates = new ArrayList<>();
        ArrayList<Integer> candidates_ind = new ArrayList<>();
        Chromosome[] ca = new Chromosome[2];

        int candidates_index = 0;

        for (int i = 0; i < k; i++) {// select K Candidates
            candidates_index = selectCandidates(c, population_size, k);
            c[candidates_index].selcted = true;
            candidates.add(c[candidates_index]);
            candidates_ind.add(candidates_index);

            // System.out.print(candidates_index +" c[candidates_index]: ");
            // System.out.println(c[candidates_index].fitness);
        }
        float best_fitness = -1;
        int best_fitness_ind = 0;

        for (int i = 0; i < k; i++) {
            if (best_fitness == -1 || candidates.get(i).fitness > best_fitness) {
                best_fitness = candidates.get(i).fitness;
                best_fitness_ind = i;
            }
            c[candidates_ind.get(i)].selcted = false;

        }

        c[candidates_ind.get(best_fitness_ind)].selcted = true;

        ca[0] = new Chromosome(c[candidates_ind.get(best_fitness_ind)].genes,
                c[candidates_ind.get(best_fitness_ind)].fitness, c[candidates_ind.get(best_fitness_ind)].discarded,
                c[candidates_ind.get(best_fitness_ind)].selcted);

        // ca[0]=c[candidates_ind.get(best_fitness_ind)];

        // System.out.print(candidates_ind.get(best_fitness_ind)+" first selected
        // chromosome");
        // System.out.println(c[candidates_ind.get(best_fitness_ind)].fitness);
        ArrayList<Chromosome> candidates2 = new ArrayList<>();
        ArrayList<Integer> candidates_ind2 = new ArrayList<>();
        int candidates_index2 = 0;
        for (int i = 0; i < k; i++) {
            candidates_index2 = selectCandidates(c, population_size, k);
            c[candidates_index2].selcted = true;
            candidates2.add(c[candidates_index2]);
            candidates_ind2.add(candidates_index2);
        }

        float best_fitness2 = -1;
        int best_fitness_ind2 = 0;

        for (int i = 0; i < k; i++) {
            if (best_fitness2 == -1 || candidates2.get(i).fitness >= best_fitness2) {

                best_fitness2 = candidates2.get(i).fitness;
                best_fitness_ind2 = i;

            }
            c[candidates_ind2.get(i)].selcted = false;

        }

        c[candidates_ind2.get(best_fitness_ind2)].selcted = true;

        ca[1] = new Chromosome(c[candidates_ind2.get(best_fitness_ind2)].genes,
                c[candidates_ind2.get(best_fitness_ind2)].fitness, c[candidates_ind2.get(best_fitness_ind2)].discarded,
                c[candidates_ind2.get(best_fitness_ind2)].selcted);

        // ca[1]=c[candidates_ind2.get(best_fitness_ind2)];
        // System.out.print(candidates_ind2.get(best_fitness_ind2)+" second selected
        // chromosome");
        // System.out.println(c[candidates_ind2.get(best_fitness_ind2)].fitness);

        // free the selected choromosom
        c[candidates_ind.get(best_fitness_ind)].selcted = false;
        c[candidates_ind2.get(best_fitness_ind2)].selcted = false;
        ca[0].selcted = false;
        ca[1].selcted = false;

        ca = crossOver(c, candidates_ind.get(best_fitness_ind), candidates_ind2.get(best_fitness_ind2), bits_num);
        ca[0] = mutation(ca[0], bits_num, t, T);
        ca[1] = mutation(ca[1], bits_num, t, T);
        ca=Calculate_fitness( ca,  p,  number_of_points, bits_num, 2);
        c = replace(c, ca, population_size);

        return c;

    }

    static int selectCandidates(Chromosome[] c, int population_size, int k) {
        Random random = new Random(); // instance of random class
        float totalFitness = 0;
        int candidates_index = 0;
        float cumulativeFitness = 0;

        for (int i = 0; i < population_size; i++) {

            if (!c[i].selcted) {
                totalFitness += c[i].fitness;

            }

        }
        // System.out.print("totalFitness: ");
        // System.out.println(totalFitness);
        float randomFitness = random.nextFloat(totalFitness);

        // System.out.print("randomFitness: ");
        // System.out.println(randomFitness);
        for (int j = 0; j < population_size; j++) {
            if (!c[j].selcted) {
                cumulativeFitness += c[j].fitness;

                if (cumulativeFitness >= randomFitness) {
                    c[j].selcted = true;
                    candidates_index = j;
                    break;
                }
            }
        }
        // System.out.print("cumulativeFitness: ");
        // System.out.println(cumulativeFitness);
        // System.out.print("candidates_index: ");
        // System.out.println(candidates_index);
        return candidates_index;
    }

    static Chromosome[] Calculate_fitness(Chromosome[] c, point[] p, int number_of_points, int num_of_bits,
            int population_size) {
        for (int k = 0; k < population_size; k++) {
            float C_error = 0;

            for (int i = 0; i < number_of_points; i++) {
                float error = 0;

                for (int j = 0; j < num_of_bits; j++) {
                    error += c[k].genes.get(j) * Math.pow(p[i].x, j);
                }
                error -= p[i].y;
                error = (float) Math.pow(error, 2);
                C_error += error;// cummulative error for all point errors
            }
            C_error = C_error / number_of_points;
            c[k].fitness = (1 / C_error);

            // System.out.println("c[k].fitness =");

            // System.out.println(c[k].fitness);

        }

        return c;
    }

    static Chromosome[] Population_initialization(Chromosome[] c, int polynomial_degree, int population_size) {

        for (int i = 0; i < population_size; i++) {
            c[i] = new Chromosome();
            for (int j = 0; j < polynomial_degree + 1; j++) {
                Random rand = new Random(); // instance of random class

                c[i].genes.add(rand.nextFloat(-10, 10));
            }
        }
        return c;

    }

    public static void main(String[] args) throws IOException {
        FileWriter fileWriter = new FileWriter(
                "C:/Users/User/Desktop/Fcai-4-1st term/Soft Computing/Codes/assignment2/out.txt");
        PrintWriter myWriter = new PrintWriter(fileWriter);
        File file = new File("C:/Users/User/Desktop/Fcai-4-1st term/Soft Computing/Codes/assignment2/input.txt");
        BufferedReader br = new BufferedReader(new FileReader(file));
        int Number_of_datasets = 1;
        int population_size = 10000;
        int number_of_points = 0;
        int polynomial_degree = 0;
        int numberOfgenerations = 1000;

        String st = br.readLine();
        Number_of_datasets = Integer.parseInt(st);
        for (int j = 0; j < Number_of_datasets; j++) {
            myWriter.println("   dataset Number: " + j);

            st = br.readLine();
            String[] arrOfStr = st.split(" ");
            number_of_points = Integer.parseInt(arrOfStr[0]);
            polynomial_degree = Integer.parseInt(arrOfStr[1]);
            // System.out.println(number_of_points);
            // System.out.println(polynomial_degree);

            point[] p = new point[number_of_points];
            Chromosome[] c = new Chromosome[population_size];
            for (int i = 0; i < number_of_points; i++) {
                st = br.readLine();
                String[] arrOfStr1 = st.split(" ");
                p[i] = new point();
                p[i].x = Float.parseFloat(arrOfStr1[0]);
                p[i].y = Float.parseFloat(arrOfStr1[1]);
                // System.out.println(p[i].x);
                // System.out.println(p[i].y);

            }

            float square_error = 0;
            c = Population_initialization(c, polynomial_degree, population_size);



            for (int jj = 0; jj < numberOfgenerations; jj++) {
                c = Calculate_fitness(c, p, number_of_points, polynomial_degree + 1, population_size);
                c = tournament_selection(c, population_size, 50, polynomial_degree + 1, jj, numberOfgenerations ,p,number_of_points);
            }


            c = Calculate_fitness(c, p, number_of_points, polynomial_degree + 1, population_size);
            square_error = get_error_number(c, p, number_of_points, polynomial_degree + 1, population_size, myWriter);
            myWriter.println("MSE:");
            myWriter.println(square_error);

    /*             myWriter.println("last gen fitness:");

        for(int jj=0;jj<population_size;jj++){

             myWriter.println(c[jj].fitness);

            }*/

        }
        myWriter.close();
    }

}



