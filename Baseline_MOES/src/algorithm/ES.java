package algorithm;

import java.util.ArrayList;
import java.util.LinkedList;
import java.util.Random;
import org.apache.commons.math3.distribution.MultivariateNormalDistribution;
import util.Graph;
import util.Path;
import util.Point;
import java.io.BufferedReader;
import java.io.InputStreamReader;

public class ES {
    public final int numOfGeneration = 50; // number of generation
    public final int elite = 30;
    public final int children = 80;
    public Path particles[] = new Path[children];
    public Path gBest, gBestDistance, gBestSmooth, gBestSafety;
    public double startPopulation[];
    public double candidate[];
    public int rank;
    public int rank0Count;
    public Path initialCandidate;
    public Graph graph;
    public double R; // radius
    public double maxPointy = 15;
    public double minPointy = -15;
    public double maxVariance = 15;
    public double minVariance = -15;
    public double mean[];
    public static Point startPoint;
    public static Point endPoint;
    public int numR; // number of R in map
    Random random = new Random();
    public String numberTeString;
    public double identityMatrix[][];
    public double[] standardDevi;
    public Path[] paretoFront;
    public double AB;
    public LinkedList<Point> resultDistance = new LinkedList<Point>();
    public LinkedList<Point> resultSafety = new LinkedList<Point>();
    public LinkedList<Point> resultSmooth = new LinkedList<Point>();
    public ArrayList<LinkedList<Point>> resultPareto = new ArrayList<LinkedList<Point>>();

    public ES(int numR, Point start, Point end, Graph graph, String numberTeString) {
        startPoint = start;
        endPoint = end;
        this.numR = numR;
        this.graph = graph;
        this.AB = Math.hypot(end.x - start.x, end.y - start.y);
        this.R = AB / (numR + 1);
        this.numberTeString = numberTeString;
    }

    public void initialize(int numR) {
        standardDevi = new double[numR];
        mean = new double[numR];

        double tmp = Path.convertPointToPointToBeginning(16.635425876805183,
        58.090775125257885, R,
        startPoint, endPoint);
        System.out.println(tmp);
        Point pointstest = Path.convertPointToPoint(1.1224148098978217, R,
        startPoint, endPoint);
        System.out.println(pointstest.x + " " + pointstest.y);

        do {
        double pointy[] = new double[numR];
        Point points[] = new Point[numR];
        for (int j = 0; j < numR; j++) {
        standardDevi[j] = random.nextDouble() * (maxVariance - minVariance) +
        minVariance;
        do {
        pointy[j] = random.nextDouble() * (maxPointy - minPointy) + minPointy;
        points[j] = Path.convertPointToPoint(pointy[j], (j + 1) * R, startPoint,
        endPoint);
        } while (!points[j].inCoordinate());
        }
        initialCandidate = new Path(numR, R, pointy, points);
        } while (pathCollision(initialCandidate) == true);

        identityMatrix = new double[numR][numR];
        for (int i = 0; i < numR; i++) {
            identityMatrix[i][i] = 1;
        }
        gBest = initialCandidate;
        gBestDistance = initialCandidate;
        gBestSafety = initialCandidate;
        gBestSmooth = initialCandidate;
    }

    public boolean pathCollision(Path path) {
        for (int i = 0; i < numR; i++) {
            if (i == 0) {
                if (graph.isIntersectLine(startPoint, path.points[i])) {
                    return true;
                }
                if (numR == 1) {
                    return graph.isIntersectLine(endPoint, path.points[i]);
                }
            } else if (i == numR - 1) {
                if (graph.isIntersectLine(endPoint, path.points[i])
                        || graph.isIntersectLine(path.points[i], path.points[i - 1])) {
                    return true;
                }
            } else if (i != 0) {
                if (graph.isIntersectLine(path.points[i], path.points[i - 1])) {
                    return true;
                }
            }
        }
        return false;
    }

    public static double[] add(double[] first, double[] second) {
        int length = first.length < second.length ? first.length : second.length;
        double[] result = new double[length];
        for (int i = 0; i < length; i++) {
            result[i] = first[i] + second[i];
        }
        return result;
    }

    public static double[] minusSquare(double[] first, double[] second) {
        int length = first.length < second.length ? first.length : second.length;
        double[] result = new double[length];
        for (int i = 0; i < length; i++) {
            result[i] = (first[i] - second[i]) * (first[i] - second[i]);
        }
        return result;
    }

    public static double[] multiple(double[] first, double[] second) {
        int length = first.length < second.length ? first.length : second.length;
        double[] result = new double[length];
        for (int i = 0; i < length; i++) {
            result[i] = first[i] * second[i];
        }
        return result;
    }

    public void bubbleSort(Path arr[], int n) {
        for (int i = 0; i < n - 1; i++)
            for (int j = 0; j < n - i - 1; j++)
                if (arr[j].distance > arr[j + 1].distance) {
                    Path temp = arr[j];
                    arr[j] = arr[j + 1];
                    arr[j + 1] = temp;
                }
    }

    public boolean checkDominate(Path particle1, Path particle2) {
        if ((particle1.distance <= particle2.distance
                && particle1.pathSafety(graph) <= particle2.pathSafety(graph)
                && particle1.pathSmooth() <= particle2.pathSmooth())
                && (particle1.distance < particle2.distance
                        && particle1.pathSafety(graph) < particle2.pathSafety(graph)
                        && particle1.pathSmooth() < particle2.pathSmooth())) {
            return true;
        } else
            return false;
    }

    public boolean checkRank(Path[] particle) {
        for (int i = 0; i < particle.length; i++) {
            if (particle[i].rank == -1)
                return false;
        }
        return true;
    }

    // Start crowding distance
    // Return rank of each particle
    public int[] particleRank(Path[] particles, int type) {
        int len = particles.length;
        int[] rank = new int[len];
        double[] obj = new double[len];
        int count;
        // Sap xep cac particle theo tieu chi
        if (type == 1) {
            for (int i = 0; i != len; i++) {
                if (particles[i].points[0] != null) {
                    obj[i] = particles[i].distance;
                } else {
                    obj[i] = Double.POSITIVE_INFINITY;
                }
            }
        } else if (type == 2) {
            for (int i = 0; i != len; i++) {
                if (particles[i].points[0] != null) {
                    obj[i] = particles[i].pathSafety(graph);
                } else {
                    obj[i] = Double.POSITIVE_INFINITY;
                }
            }
        } else if (type == 3) {
            for (int i = 0; i != len; i++) {
                if (particles[i].points[0] != null) {
                    obj[i] = particles[i].pathSmooth();
                } else {
                    obj[i] = Double.POSITIVE_INFINITY;
                }
            }
        }

        for (int i = 0; i != len; i++) {
            count = 0;
            for (int j = 0; j != len; j++) {
                if (j != i && obj[j] >= obj[i]) {
                    count++; // Count number of particle worse than obj[i]
                }
            }
            rank[i] = len - count - 1;
            for (int k = 0; k != i; k++) {
                if (rank[k] == rank[i]) {
                    rank[i] += 1;
                }
            }
        }
        return rank;
    }

    public int[] particleRankCD(double[] CDparticles) {
        int len = CDparticles.length;
        int[] rank = new int[len];
        int count;
        for (int i = 0; i != len; i++) {
            count = 0;
            for (int j = 0; j != len; j++) {
                if (j != i && CDparticles[j] <= CDparticles[i]) {
                    count++; // Count number of particle worse than obj[i]
                }
            }
            rank[i] = len - count - 1;
            for (int k = 0; k != i; k++) {
                if (rank[k] == rank[i]) {
                    rank[i] += 1;
                }
            }
        }
        return rank;
    }

    // Return index of each particle from highest rank to lowest rank
    public int[] indexRank(int[] rank) {
        int length = rank.length;
        int index[] = new int[length];
        for (int i = 0; i < length; i++) {
            index[rank[i]] = i;
        }
        return index;
    }

    public double[] crowdingDistance(Path[] particles) {
        int len = particles.length;
        double[] CD = new double[len];
        double[] dis = new double[len];
        double[] safety = new double[len];
        double[] smooth = new double[len];
        int[] rankDistance = new int[len];
        int[] rankSafety = new int[len];
        int[] rankSmooth = new int[len];
        int[] rerankDistance = new int[len];
        int[] rerankSafety = new int[len];
        int[] rerankSmooth = new int[len];
        rankDistance = particleRank(particles, 1);
        rankSafety = particleRank(particles, 2);
        rankSmooth = particleRank(particles, 3);
        rerankDistance = indexRank(rankDistance);
        rerankSafety = indexRank(rankSafety);
        rerankSmooth = indexRank(rankSmooth);

        for (int i = 0; i != len; i++) {
            CD[i] = 0;
            if (particles[i].points[0] != null) {
                dis[i] = particles[i].distance;
                safety[i] = particles[i].pathSafety(graph);
                smooth[i] = particles[i].pathSmooth();
            }
        }
        int index = 0; // Count number of particle that is null
        for (int i = 0; i != len; i++) {
            if (particles[i].points[0] == null) {
                index++;
            }
        }

        for (int i = 0; i != len; i++) {
            if (rankDistance[i] == 0 || rankDistance[i] == (len - 1 - index)) {
                CD[i] += Double.POSITIVE_INFINITY;
            }
            if (rankSafety[i] == 0 || rankSafety[i] == (len - 1 - index)) {
                CD[i] += Double.POSITIVE_INFINITY;
            }
            if (rankSmooth[i] == 0 || rankSmooth[i] == (len - 1 - index)) {
                CD[i] += Double.POSITIVE_INFINITY;
            }

            if (particles[i].points[0] == null) {
                CD[i] = 0;
            } else if (rankDistance[i] != 0 && rankDistance[i] != (len - 1 - index) && rankSmooth[i] != 0
                    && rankSmooth[i] != (len - 1 - index) && rankSafety[i] != 0 && rankSafety[i] != (len - 1 - index)) {
                // For example i = 0, particle 0 has rank distance = 4 -> need to find particle with
                // rank distance = 3 and 5 and their distance
                CD[i] = CD[i] + (dis[rerankDistance[rankDistance[i] + 1]] - dis[rerankDistance[rankDistance[i] - 1]])
                        / (dis[rerankDistance[len - 1 - index]] - dis[rerankDistance[0]]);
                CD[i] = CD[i] + (safety[rerankSafety[rankSafety[i] + 1]] - safety[rerankSafety[rankSafety[i] - 1]])
                        / (safety[rerankSafety[len - 1 - index]] - safety[rerankSafety[0]]);
                CD[i] = CD[i] + (smooth[rerankSmooth[rankSmooth[i] + 1]] - smooth[rerankSmooth[rankSmooth[i] - 1]])
                        / (smooth[rerankSmooth[len - 1 - index]] - smooth[rerankSmooth[0]]);
            }
        }
        return CD;
    }
    // End crowding distance

    public void run() {
        initialize(numR);

        MultivariateNormalDistribution mnd = new MultivariateNormalDistribution(mean, identityMatrix);
        startPopulation = initialCandidate.pointy;

        // Run numOfGeneration generation
        for (int iter = 0; iter < numOfGeneration; iter++) {
            // System.out.println();
            System.out.println("Iteration: " + iter);
            // Generate n children in 1 generation, only generate child that doesnt collide
            // ArrayList<Path> particlesArrayList = new ArrayList<>();
            for (int i = 0; i < children; i++) {
                do {
                    double pointy[] = new double[numR];
                    Point points[] = new Point[numR];
                    pointy = add(startPopulation, multiple(standardDevi, mnd.sample()));
                    for (int j = 0; j < numR; j++) {
                        points[j] = Path.convertPointToPoint(pointy[j], (j + 1) * R, startPoint, endPoint);
                    }
                    particles[i] = new Path(numR, R, pointy, points);
                    particles[i].distance();
                } while (pathCollision(particles[i]) == true);
                // particlesArrayList.add(particles[i]);
            }

            // Begin multi-objective
            ArrayList<Integer> chosenParticleIndex = new ArrayList<Integer>();
            rank = 0;
            int indexOfLastRank = 0;
            // Select only elite child, if chosenParticleIndex has more element than elite,
            // use crowding distance to sort best child

            Path particlestmp[] = particles.clone();
            while (chosenParticleIndex.size() < elite) {
                int i1 = 0;
                indexOfLastRank = chosenParticleIndex.size();
                while (i1 < particlestmp.length) {
                    int i2 = 0;
                    boolean cont = true;
                    while (i2 < particlestmp.length && cont == true) {
                        if ((i2 == (particlestmp.length) - 1) && (!checkDominate(particles[i2], particles[i1]))) {
                            chosenParticleIndex.add(i1);
                        }
                        if (checkDominate(particles[i2], particles[i1])) {
                            cont = false;
                        }
                        i2++;
                    }
                    i1++;
                }

                // Choose unrated particle
                particlestmp = new Path[children - chosenParticleIndex.size()];
                int j = 0;
                for (int i = 0; i < particlestmp.length; i++) {
                    if (chosenParticleIndex.contains(j))
                        j++;
                    particlestmp[i] = particles[j];
                    j++;
                }
                if (rank == 0) {
                    rank0Count = chosenParticleIndex.size();
                }
                rank++;
            }
            // Choose all partcle of rank n for crowding distance sort
            Path[] crowdingDistanceSort = new Path[chosenParticleIndex.size() - indexOfLastRank];
            for (int i = 0; i < chosenParticleIndex.size() - indexOfLastRank; i++) {
                crowdingDistanceSort[i] = particles[chosenParticleIndex.get(i + indexOfLastRank)];
            }

            double[] selectedCD = crowdingDistance(crowdingDistanceSort);

            // Sort crowding distance
            int[] rankCD = new int[selectedCD.length];
            int[] rerankCD = new int[selectedCD.length];

            rankCD = particleRankCD(selectedCD);

            rerankCD = indexRank(rankCD);
            Path[] elitePaths = new Path[elite];
            for (int i = 0; i < indexOfLastRank; i++) {
                elitePaths[i] = particles[chosenParticleIndex.get(i)];
            }
            for (int i = indexOfLastRank; i < elite; i++) {
                elitePaths[i] = particles[rerankCD[i - indexOfLastRank]];
            }
            // End multiobjective

            // Calculate new standard deviation
            standardDevi = new double[numR];
            for (int i = 0; i < elite; i++) {
                standardDevi = add(standardDevi, minusSquare(elitePaths[i].pointy,
                        startPopulation));
            }
            for (int i = 0; i < numR; i++) {
                standardDevi[i] = standardDevi[i] / elite;
                standardDevi[i] = Math.sqrt(standardDevi[i]);
            }

            // Calculate new mean
            startPopulation = new double[numR];
            for (int i = 0; i < elite; i++) {
                startPopulation = add(startPopulation, elitePaths[i].pointy);
            }
            for (int i = 0; i < numR; i++) {
                startPopulation[i] = startPopulation[i] / elite;
            }

            // End

            if (iter == numOfGeneration - 1) {
                paretoFront = new Path[rank0Count];
                for (int i = 0; i < rank0Count; i++) {
                    paretoFront[i] = particles[chosenParticleIndex.get(i)];
                }

                int[] rankPareto = new int[paretoFront.length];
                int[] rerankParetoDistance = new int[paretoFront.length];
                int[] rerankParetoSafety = new int[paretoFront.length];
                int[] rerankParetoSmooth = new int[paretoFront.length];

                // Choose gBest by distance
                rankPareto = particleRank(paretoFront, 1);
                rerankParetoDistance = indexRank(rankPareto);
                rankPareto = particleRank(paretoFront, 2);
                rerankParetoSafety = indexRank(rankPareto);
                rankPareto = particleRank(paretoFront, 3);
                rerankParetoSmooth = indexRank(rankPareto);

                gBestDistance = particles[rerankParetoDistance[0]];
                gBestSafety = particles[rerankParetoSafety[0]];
                gBestSmooth = particles[rerankParetoSmooth[0]];
            }
        }

        resultDistance.add(startPoint);
        for (int i = 0; i < numR; i++) {
            resultDistance.add(gBestDistance.points[i]);
        }
        resultDistance.add(endPoint);
        resultDistance.removeLast();
        resultDistance.removeFirst();

        resultSafety.add(startPoint);
        for (int i = 0; i < numR; i++) {
            resultSafety.add(gBestSafety.points[i]);
        }
        resultSafety.add(endPoint);
        resultSafety.removeLast();
        resultSafety.removeFirst();

        resultSmooth.add(startPoint);
        for (int i = 0; i < numR; i++) {
            resultSmooth.add(gBestSmooth.points[i]);
        }
        resultSmooth.add(endPoint);
        resultSmooth.removeLast();
        resultSmooth.removeFirst();
    }
}
