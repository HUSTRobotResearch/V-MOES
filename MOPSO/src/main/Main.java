package main;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.Scanner;

import algorithm.PSO;
import gui.GUIRobotics;
import util.Graph;
import util.Path;
import util.Point;

public class Main {
	public static void main(String[] args) throws IOException {
		int count = 0;
		boolean flag = false;
		while (count < 20) {
			long time = System.currentTimeMillis();
			// Change the test number here
			String FILE_URL = "exp2/mopso_test30.txt";
			File file = new File(FILE_URL);
			String numberTeString = "30";
			// Create GUI
            GUIRobotics gui = new GUIRobotics(800, 100, 10);
            gui.generateEnvironment("input/obstacle_" + numberTeString + ".txt");

			// Read input data
			Graph graph = new Graph("input/obstacle_" + numberTeString + ".txt");
			LinkedList<Point> pointsToVisit = readPointData("input/input_" + numberTeString + ".txt");

			ArrayList<Point> result = new ArrayList<Point>();
			PSO pso = new PSO(10, pointsToVisit.get(0), pointsToVisit.get(1), graph);
			try {

				pso.run();
				result.add(pointsToVisit.get(0));
				for (int j = 0; j < pso.result.size(); j++) {
					result.add(pso.result.get(j));
				}
				result.add(pointsToVisit.get(1));

				// To write all point into a file
				try {
					String FILE_URL_TESTING = "testing/mopso_test_point_" + numberTeString + ".txt";
					File file_test = new File(FILE_URL_TESTING);
					if (!file_test.exists()) {
						file_test.createNewFile();
					}
					FileWriter fw = new FileWriter(file_test.getAbsoluteFile(), true);
					BufferedWriter bw = new BufferedWriter(fw);

					for (Point point : result) {
						bw.write(point.x + " " + point.y + "\n");
					}
					bw.write("-1\n");
					bw.close();
				} catch (Exception e) {
					System.out.println("Something went wrong!");
					e.printStackTrace();
				}
				// --------------------------------------------------------
			} catch (Exception e) {
				System.out.println("Something went wrong!");
				e.printStackTrace();
			}

			time = System.currentTimeMillis() - time;
			try {
				if (!file.exists()) {
					file.createNewFile();
				}
				FileWriter fw = new FileWriter(file.getAbsoluteFile(), true);
				BufferedWriter bw = new BufferedWriter(fw);
				for (Path path : pso.NaParticles) {
					if (path.points[0] != null) {
						flag = true;

						bw.write("Path :" + "  " + (double) Math.round(path.distance * 10000) / 10000 + " "
								+ (double) Math.round(path.pathSafety(graph) * 10000) / 10000 + "  "
								+ (double) Math.round(path.pathSmooth() * 10000) / 10000 + "\n");
					}
				}
				if (flag == true)
					count++;
				flag = false;
				bw.write("Total execution time: " + (time) + "\n");
				bw.write("----------------------\n");
				bw.close();
			} catch (Exception e) {
				System.out.println("Something went wrong!");
				e.printStackTrace();
			}
			System.out.println("Time:\t" + time + " ms");
			System.out.println("result" + result);
			System.out.println("Done!");
		}
	}

	public static LinkedList<Point> readPointData(String filename) throws FileNotFoundException {
		Scanner scan = new Scanner(new File(filename));
		LinkedList<Point> pointsToVisit = new LinkedList<Point>();
		double x = scan.nextDouble();
		while (x != -1) {
			double y = scan.nextDouble();
			pointsToVisit.addLast(new Point(x, y));
			x = scan.nextDouble();
		}
		scan.close();

		return pointsToVisit;
	}
}
