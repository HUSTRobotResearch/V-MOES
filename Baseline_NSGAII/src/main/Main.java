package main;

import java.awt.Color;
import java.io.*;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.Scanner;

import algorithm.NSGAII;
import graph.GraphDivision;
import gui.GUIRobotics;
import util.Graph;
import util.Path;
import util.Point;

public class Main {
	public static void main(String[] args) throws IOException {
		int count = 0;
		boolean flag = false;
		while (count < 1) {
			final long startTime = System.currentTimeMillis();
			// Change the test number here
			String FILE_URL = "output/nsgaii_test29.txt";
			File file = new File(FILE_URL);
			String numberTeString = "9";

			// Create GUI
			 GUIRobotics gui = new GUIRobotics(800, 100, 10);
			 gui.generateEnvironment("input/obstacle_" + numberTeString + ".txt");
			Graph graph = new Graph("input/obstacle_" + numberTeString + ".txt");
			// Read input data
			ArrayList<Point> result = new ArrayList<Point>();

			LinkedList<Point> startEndPoints = readPointData("input/input_" + numberTeString + ".txt");
			GraphDivision graphDivision = new GraphDivision(graph, startEndPoints);
			NSGAII nsgaii = new NSGAII(graph, startEndPoints.getFirst(), startEndPoints.getLast());
			try {
				LinkedList<Point> findPathAfterFixed = nsgaii.getPathAfterFixed();
				result.add(startEndPoints.get(0));
				for (int j = 0; j < findPathAfterFixed.size(); j++) {
					result.add(findPathAfterFixed.get(j));
				}
				result.add(startEndPoints.get(1));

				 for (int i = 0; i < (findPathAfterFixed.size() - 1); i++) {
				 gui.canvas.drawLines(result, startEndPoints);
				 }

				System.out.println("output" + result);

				System.out.println("Done");
				System.out.println("Path Distance: " + (double) Math.round(nsgaii.pathDistance * 10000) / 10000);
				System.out.println("Path Smooth: " + (double) Math.round(nsgaii.pathSmooth * 10000) / 10000);
				System.out.println("Path Safety: " + (double) Math.round(nsgaii.pathSafety * 10000) / 10000);
				for (int j = 0; j < findPathAfterFixed.size(); j++) {
					System.out.println(nsgaii.pointsToVisitAfterFixed.get(j).x + " " + nsgaii.pointsToVisitAfterFixed.get(j).y);
				}
				
			} catch (Exception e) {
				System.out.println("Something went wrong!");
				e.printStackTrace();
			}

			final long endTime = System.currentTimeMillis();
			try {
				if (!file.exists()) {
					file.createNewFile();
				}
				FileWriter fw = new FileWriter(file.getAbsoluteFile(), true);
				BufferedWriter bw = new BufferedWriter(fw);
				for (Path path : nsgaii.NEWPOP) {
					if (((double) Math.round(path.pathSmooth() * 10000) / 10000 != 0.0)
							&& ((double) (path.pathDistance()) < 200)) {
						try {
							String FILE_URL_PATH_POINT = "output/path_point/nsgaii_test_point_" + numberTeString + ".txt";
							File file_test = new File(FILE_URL_PATH_POINT);
							if (!file_test.exists()) {
								file_test.createNewFile();
							}
							FileWriter fw_path_point = new FileWriter(file_test.getAbsoluteFile(), true);
							BufferedWriter bw_path_point = new BufferedWriter(fw_path_point);

							for (Point point : path.points) {
								bw_path_point.write(point.x + " " + point.y + "\n");
							}
							bw_path_point.write("----------------------\n");
							bw_path_point.close();
						} catch (Exception e) {
							System.out.println("Something went wrong!");
							e.printStackTrace();
						}
						flag = true;

						bw.write("Path :" + "  " + (double) Math.round(path.pathDistance() * 10000) / 10000 + " "
								+ (double) Math.round(path.pathSafety(graph) * 10000) / 10000 + "  "
								+ (double) Math.round(path.pathSmooth() * 10000) / 10000 + "\n");
					}
				}
				if (flag == true)
					count++;
				flag = false;
				bw.write("Total execution time: " + (endTime - startTime) + "\n");
				bw.write("----------------------\n");
				bw.close();
			} catch (Exception e) {
				System.out.println("Something went wrong!");
				e.printStackTrace();
			}
			System.out.println("Total execution time: " + (endTime - startTime));
		
			
			System.out.println("End!");
			count++;
		}
	}

	public static LinkedList<Point> readPointData(String filename) throws IOException {
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
