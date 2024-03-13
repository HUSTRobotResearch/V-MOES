// package gui;

// import java.awt.*;
// import java.awt.event.*;
// import java.awt.geom.Ellipse2D;
// import java.awt.geom.Rectangle2D;
// import java.awt.geom.Line2D;
// import java.io.FileInputStream;
// import java.io.FileNotFoundException;
// import java.util.ArrayList;
// import java.util.LinkedList;
// import java.util.Scanner;

// import graph.ObstaclesGraph;
// import gui.GUIRobotics.MyCanvas;
// import util.Point;

// public class GUIRobotics {

// 	private Frame mainFrame;
// 	private Panel controlPanel;
// 	public MyCanvas canvas;
// 	private static int size = 600;
// 	private static double range = 110;
// 	private static int numOfRange = 11;

// 	public GUIRobotics(int size, double range, int numOfRange) {
// 		GUIRobotics.size = size;
// 		GUIRobotics.range = range;
// 		GUIRobotics.numOfRange = numOfRange;
// 		prepareGUI();
// 	}

// 	private void prepareGUI() {
// 		mainFrame = new Frame("GUI for Robot Path Planning");
// 		mainFrame.setSize(size + 60, size + 60);
// 		mainFrame.setResizable(false);
// 		mainFrame.setLayout(new GridLayout(1, 1));
// 		mainFrame.setLocationRelativeTo(null);

// 		mainFrame.addWindowListener(new WindowAdapter() {
// 			public void windowClosing(WindowEvent windowEvent) {
// 				System.exit(0);
// 			}
// 		});

// 		controlPanel = new Panel();
// 		controlPanel.setLayout(new FlowLayout());

// 		mainFrame.add(controlPanel);
// 		mainFrame.setVisible(true);
// 	}

// 	public void generateEnvironment(String obtacles_file) {
// 		canvas = new MyCanvas();

// 		controlPanel.add(canvas);

// 		// Draw obstacles
// 		Graphics2D g2 = (Graphics2D) canvas.getGraphics();
// 		ArrayList<ObstaclesGraph> obstacles = ObstaclesGraph.getObstacles(obtacles_file);

// 		for (ObstaclesGraph obstacle : obstacles) {
// 			Polygon polygon = new Polygon();
// 			for (int i = 0; i < obstacle.points.size(); i++) {
// 				polygon.addPoint((int) (MyCanvas.OX + obstacle.points.get(i).x * size / range),
// 						(int) (MyCanvas.OY - obstacle.points.get(i).y * size / range));
// 			}
// 			g2.setColor(Color.darkGray);
// 			g2.fill(polygon);
// 			for (int i = 0; i < obstacle.points.size() - 1; i++) {
// 				canvas.drawLine(obstacle.points.get(i), obstacle.points.get(i + 1), Color.darkGray);
// 			}
// 			canvas.drawLine(obstacle.points.get(0), obstacle.points.get(obstacle.points.size() - 1), Color.darkGray);
// 		}

// 		// draw Oxy

// 		g2.drawLine(MyCanvas.OX, MyCanvas.OY, MyCanvas.OY, MyCanvas.OY);
// 		g2.drawLine(MyCanvas.OX, MyCanvas.OY, MyCanvas.OX, MyCanvas.OX);
// 		g2.drawString("O", MyCanvas.OX - 10, MyCanvas.OY + 10);
// 		g2.drawString("x", MyCanvas.OY + 5, MyCanvas.OY);
// 		g2.drawString("y", MyCanvas.OX, MyCanvas.OX);

// 		// draw gird
// 		float[] dash1 = { 2f, 0f, 2f };
// 		BasicStroke bs1 = new BasicStroke(1, BasicStroke.CAP_BUTT, BasicStroke.JOIN_ROUND, 1.0f, dash1, 2f);
// 		g2.setStroke(bs1);
// 		for (int i = 0; i < numOfRange; i++) {
// 			g2.draw(new Line2D.Double(MyCanvas.OX, MyCanvas.OY - size / range * 10 * (i + 1), MyCanvas.OY,
// 					MyCanvas.OY - size / range * 10 * (i + 1)));

// 			g2.draw(new Line2D.Double(MyCanvas.OX + size / range * 10 * (i + 1), MyCanvas.OY,
// 					MyCanvas.OX + size / range * 10 * (i + 1), MyCanvas.OX));

// 			g2.drawString(String.valueOf(range / numOfRange * (i + 1)),
// 					(int) (MyCanvas.OX + size / range * 10 * (i + 1) - 10), MyCanvas.OY + 10);

// 			g2.drawString(String.valueOf(range / numOfRange * (i + 1)), MyCanvas.OX,
// 					(int) (MyCanvas.OY - size / range * 10 * (i + 1)));
// 		}

// 		mainFrame.setVisible(true);
// 	}

// 	public static class MyCanvas extends Canvas {

// 		static int OX = 10;
// 		static int OY = size - OX;

// 		private double alpha = size / range;

// 		public MyCanvas() {
// 			setSize(size, size);
// 		}

// 		public void drawPoint(Point p) {
// 			Graphics2D g2 = (Graphics2D) getGraphics();
// 			g2.setColor(Color.BLUE);
// 			g2.fill(new Ellipse2D.Double(OX + p.x * alpha - 3, OY - p.y * alpha - 3, 6, 6));
// 		}

// 		public void drawPoint(Point p, Color color) {
// 			Graphics2D g2 = (Graphics2D) getGraphics();
// 			g2.setColor(color);
// 			g2.fill(new Ellipse2D.Double(OX + p.x * alpha - 3, OY - p.y * alpha - 3, 10, 10));
// 		}

// 		public void drawLine(Point p1, Point p2) {
// 			Graphics2D g2 = (Graphics2D) getGraphics();
// 			g2.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
// 			g2.setStroke(new BasicStroke(2));
// 			g2.setColor(Color.RED);
// 			g2.draw(new Line2D.Double(OX + p1.x * alpha, OY - p1.y * alpha, OX + p2.x * alpha, OY - p2.y * alpha));
// 		}

// 		public void drawLine(Point p1, Point p2, Color color) {
// 			Graphics2D g2 = (Graphics2D) getGraphics();
// 			g2.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
// 			g2.setStroke(new BasicStroke(2));
// 			g2.setColor(color);
// 			g2.draw(new Line2D.Double(OX + p1.x * alpha, OY - p1.y * alpha, OX + p2.x * alpha, OY - p2.y * alpha));
// 		}

// 		public void drawLines(ArrayList<Point> points, LinkedList<Point> pointsToVisit, Color color) {
// 			Graphics2D g2 = (Graphics2D) getGraphics();
// 			g2.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
// 			g2.setColor(color);
// 			g2.setStroke(new BasicStroke(2.0f));
// 			for (int i = 0; i < points.size() - 1; i++) {
// 				g2.draw(new Line2D.Double(OX + points.get(i).x * alpha, OY - points.get(i).y * alpha,
// 						OX + points.get(i + 1).x * alpha, OY - points.get(i + 1).y * alpha));
// 			}

// 			g2.setColor(Color.RED);
// 			for (int i = 0; i < points.size(); i++) {
// 				Point pt = new Point(points.get(i).x, points.get(i).y);
// 				if (pt.indexInSet(pointsToVisit) == -1) {
// 					g2.setColor(Color.black);
// 					g2.fill(new Rectangle2D.Double(OX + points.get(i).x * alpha - 3, OY - points.get(i).y * alpha - 3,
// 							6,
// 							6));
// 					g2.setColor(Color.RED);
// 				} else
// 					g2.fill(new Ellipse2D.Double(OX + points.get(i).x * alpha - 4, OY - points.get(i).y * alpha - 4, 8,
// 							8));
// 			}
// 		}

// 		public void drawLineWithMiddle(Point pt1, Point pt2) {
// 			Graphics2D g2 = (Graphics2D) getGraphics();
// 			g2.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
// 			g2.setColor(Color.BLUE);
// 			g2.draw(new Line2D.Double(OX + pt1.x * alpha, OY - pt1.y * alpha, OX + pt2.x * alpha, OY - pt2.y * alpha));

// 			g2.setColor(Color.RED);
// 			Point pt = new Point((pt1.x + pt2.x) / 2, (pt1.y + pt2.y) / 2);
// 			g2.fill(new Ellipse2D.Double(OX + pt.x * alpha - 3, OY - pt.y * alpha - 3, 6, 6));
// 		}

// 		public void drawLineWithoutMiddle(Point pt1, Point pt2) {
// 			Graphics2D g2 = (Graphics2D) getGraphics();
// 			g2.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
// 			g2.setColor(Color.BLUE);
// 			g2.draw(new Line2D.Double(OX + pt1.x * alpha, OY - pt1.y * alpha, OX + pt2.x * alpha, OY - pt2.y * alpha));
// 		}
// 	}

// 	public static void main(String[] args) throws FileNotFoundException {

// 		GUIRobotics gui = new GUIRobotics(600, 110, 11);
// 		gui.generateEnvironment("obstacles.txt");

// 		System.setIn(new FileInputStream("maklink.txt"));
// 		Scanner sc = new Scanner(System.in);

// 		sc.nextLine();
// 		String string = sc.nextLine();
// 		while (!string.equals("-1")) {

// 			// Read path planning
// 			String numbers[] = string.replaceAll(",", "").replaceAll("\\(", "").replaceAll("\\)", "").split("\\s+");

// 			ArrayList<Point> points = new ArrayList<>();
// 			for (int i = 0; i < numbers.length / 2; i++) {
// 				points.add(new Point(Double.parseDouble(numbers[2 * i]), Double.parseDouble(numbers[2 * i + 1])));
// 			}

// 			for (int i = 0; i < points.size() / 2; i++) {
// 				gui.canvas.drawLineWithoutMiddle(points.get(2 * i), points.get(2 * i + 1));
// 			}

// 			string = sc.nextLine();
// 		}
// 		sc.close();
// 	}

// }
package gui;

import java.awt.*;
import java.awt.event.*;
import java.awt.geom.Ellipse2D;
import java.awt.geom.Rectangle2D;
import java.awt.geom.Line2D;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.Scanner;

import javax.swing.JLabel;

import graph.ObstaclesGraph;
import gui.GUIRobotics.MyCanvas;
import util.Point;

public class GUIRobotics {

	private Frame mainFrame;
	private Panel controlPanel;
	public MyCanvas canvas;
	private static int size = 600;
	private static double range = 100;
	private static int numOfRange = 10;	

	public GUIRobotics(int size, double range, int numOfRange) {
		GUIRobotics.size = size;
		GUIRobotics.range = range;
		GUIRobotics.numOfRange = numOfRange;
		prepareGUI();
	}

	private void prepareGUI() {
		mainFrame = new Frame("GUI for Robot Path Planning");
		mainFrame.setSize(size + 100, size + 100);
		mainFrame.setResizable(false);
		mainFrame.setLayout(new GridLayout(1, 1));
		mainFrame.setLocationRelativeTo(null);

		mainFrame.addWindowListener(new WindowAdapter() {
			public void windowClosing(WindowEvent windowEvent) {
				System.exit(0);
			}
		});

		controlPanel = new Panel();
		controlPanel.setLayout(new FlowLayout());

		mainFrame.add(controlPanel);
		mainFrame.setVisible(true);
		
	}

	// void generateEnvironment(Graph graph) {
	// canvas = new MyCanvas();
	// controlPanel.add(canvas);
	// Graphics2D g2 = (Graphics2D) canvas.getGraphics();
	// for (Polygons obstacle : graph.obstacles) {
	// Polygon polygon = new Polygon();
	// for (int i = 0; i < obstacle.cornum; i++) {
	// polygon.addPoint((int) (MyCanvas.OX + obstacle.point[i].x * size / range),
	// (int) (MyCanvas.OY - obstacle.point[i].y * size / range));
	// }
	// g2.fill(polygon);
	// for (int i = 0; i < obstacle.cornum - 1; i++) {
	// canvas.drawLine(obstacle.point[i], obstacle.point[i + 1], Color.BLACK);
	// }
	// canvas.drawLine(obstacle.point[0], obstacle.point[obstacle.cornum - 1],
	// Color.BLACK);
	// }

	// // draw Oxy
	// // TODO: Make text larger
	// g2.drawLine(MyCanvas.OX, MyCanvas.OY, MyCanvas.OY, MyCanvas.OY);
	// g2.drawLine(MyCanvas.OX, MyCanvas.OY, MyCanvas.OX, MyCanvas.OX);
	// g2.drawString("O", MyCanvas.OX - 10, MyCanvas.OY + 10);
	// g2.drawString("x", MyCanvas.OY, MyCanvas.OY);
	// g2.drawString("y", MyCanvas.OX, MyCanvas.OX);

	// // draw gird
	// float[] dash1 = { 2f, 0f, 2f };
	// BasicStroke bs1 = new BasicStroke(1,
	// BasicStroke.CAP_BUTT,
	// BasicStroke.JOIN_ROUND,
	// 1.0f,
	// dash1,
	// 2f);
	// g2.setStroke(bs1);
	// for (int i = 0; i < numOfRange; i++) {
	// g2.draw(new Line2D.Double(MyCanvas.OX, MyCanvas.OY - size / range * 10 * (i +
	// 1),
	// MyCanvas.OY, MyCanvas.OY - size / range * 10 * (i + 1)));

	// g2.draw(new Line2D.Double(MyCanvas.OX + size / range * 10 * (i + 1),
	// MyCanvas.OY,
	// MyCanvas.OX + size / range * 10 * (i + 1), MyCanvas.OX));

	// g2.drawString(String.valueOf(range / numOfRange * (i + 1)),
	// (int) (MyCanvas.OX + size / range * 10 * (i + 1) - 10),
	// MyCanvas.OY + 10);

	// g2.drawString(String.valueOf(range / numOfRange * (i + 1)),
	// MyCanvas.OX,
	// (int) (MyCanvas.OY - size / range * 10 * (i + 1)));
	// }
	// }

	public void generateEnvironment(String obtacles_file, double ratio, boolean draw_borders) {
		canvas = new MyCanvas();
//		canvas.setSize(size + 30 , size + 30);

		controlPanel.add(canvas);

		// Draw obstacles
		ArrayList<ObstaclesGraph> obstacles = ObstaclesGraph.getObstacles(obtacles_file, Math.sqrt(ratio));
		Graphics2D g2 = (Graphics2D) canvas.getGraphics();

		for (ObstaclesGraph obstacle : obstacles) {
			if (draw_borders) {
				Polygon polygon = new Polygon();
				for (int i = 0; i < obstacle.borders.size(); i++) {
					polygon.addPoint((int) (MyCanvas.OX + obstacle.borders.get(i).x * size / range),
							(int) (MyCanvas.OY - obstacle.borders.get(i).y * size / range));
				}
				for (int i = 0; i < obstacle.borders.size() - 1; i++) {
					canvas.drawLine(obstacle.borders.get(i), obstacle.borders.get(i + 1), Color.RED);
				}
				canvas.drawLine(obstacle.borders.get(0), obstacle.borders.get(obstacle.points.size() - 1), Color.RED);
			}
			Polygon polygon = new Polygon();
			for (int i = 0; i < obstacle.points.size(); i++) {
				polygon.addPoint((int) (MyCanvas.OX + obstacle.points.get(i).x * size / range + 20),
						(int) (MyCanvas.OY - obstacle.points.get(i).y * size / range));
			}
			g2.setColor(Color.darkGray);
			g2.fill(polygon);
			for (int i = 0; i < obstacle.points.size() - 1; i++) {
				canvas.drawLine(obstacle.points.get(i), obstacle.points.get(i + 1), Color.darkGray);
			}
			canvas.drawLine(obstacle.points.get(0), obstacle.points.get(obstacle.points.size() - 1), Color.darkGray);
		}

		// draw Oxy

		g2.drawLine(MyCanvas.OX + 20, MyCanvas.OY, MyCanvas.OY, MyCanvas.OY);
		g2.drawLine(MyCanvas.OX + 20, MyCanvas.OY, MyCanvas.OX + 20, MyCanvas.OX);
		g2.drawString("O", MyCanvas.OX + 10, MyCanvas.OY + 10);
		g2.drawString("x", MyCanvas.OY + 5, MyCanvas.OY);
		g2.drawString("y", MyCanvas.OX  + 20, MyCanvas.OX);

		// draw gird
		float[] dash1 = { 2f, 0f, 2f };
		BasicStroke bs1 = new BasicStroke(1, BasicStroke.CAP_BUTT, BasicStroke.JOIN_ROUND, 1.0f, dash1, 2f);
		g2.setStroke(bs1);
		for (int i = 0; i < numOfRange - 1; i++) {
			g2.setFont(new Font("VnArial Bold", Font.PLAIN, 20));
			g2.draw(new Line2D.Double(MyCanvas.OX + 20, MyCanvas.OY - size / range * 10 * (i + 1), MyCanvas.OY,
					MyCanvas.OY - size / range * 10 * (i + 1)));

			g2.draw(new Line2D.Double(MyCanvas.OX + size / range * 10 * (i + 1) + 20, MyCanvas.OY,
					MyCanvas.OX + size / range * 10 * (i + 1) + 20, MyCanvas.OX));

			g2.drawString(String.valueOf((int)(range / numOfRange * (i + 1))),
					(int) (MyCanvas.OX + size / range * 10 * (i + 1) - 10) + 20,(int) MyCanvas.OY + 20);
		
			g2.drawString(String.valueOf((int)(range / numOfRange * (i + 1))), MyCanvas.OX - 5,
					(int) (MyCanvas.OY - size / range * 10 * (i + 1)) + 5) ;
		
		}

		mainFrame.setVisible(true);
	}

	public static class MyCanvas extends Canvas {

		static int OX = 10;
		static int OY = size - OX;

		private double alpha = size / range;

		public MyCanvas() {
			setSize(size + 30, size + 30);
		}

		public void drawLine(ObstaclesGraph.Point p1, ObstaclesGraph.Point p2) {
			Graphics2D g2 = (Graphics2D) getGraphics();
			g2.setRenderingHint(RenderingHints.KEY_ANTIALIASING,
					RenderingHints.VALUE_ANTIALIAS_ON);
			g2.setStroke(new BasicStroke(2));
			g2.setColor(Color.RED);
			g2.draw(new Line2D.Double(OX + p1.x * alpha + 20, OY - p1.y * alpha, OX + p2.x * alpha + 20, OY - p2.y * alpha));
			g2.setColor(Color.red);
		}

		public void drawLine(ObstaclesGraph.Point p1, ObstaclesGraph.Point p2, Color color) {
			Graphics2D g2 = (Graphics2D) getGraphics();
			g2.setRenderingHint(RenderingHints.KEY_ANTIALIASING,
					RenderingHints.VALUE_ANTIALIAS_ON);
			g2.setStroke(new BasicStroke(2));
			g2.setColor(color);
			g2.draw(new Line2D.Double(OX + p1.x * alpha + 20, OY - p1.y * alpha, OX + p2.x * alpha + 20, OY - p2.y * alpha));
			g2.setColor(Color.red);
		}

		public void drawPoint(Point p) {
			Graphics2D g2 = (Graphics2D) getGraphics();
			g2.setColor(Color.BLUE);
			g2.fill(new Ellipse2D.Double(OX + p.x * alpha - 3 + 20, OY - p.y * alpha - 3, 6, 6));
		}

		public void drawPoint(Point p, Color color) {
			Graphics2D g2 = (Graphics2D) getGraphics();
			g2.setColor(color);
			g2.fill(new Ellipse2D.Double(OX + p.x * alpha - 3 + 20, OY - p.y * alpha - 3, 10, 10));
		}

		public void drawLine(Point p1, Point p2) {
			Graphics2D g2 = (Graphics2D) getGraphics();
			g2.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
			g2.setStroke(new BasicStroke(2));
			g2.setColor(Color.RED);
			g2.draw(new Line2D.Double(OX + p1.x * alpha + 20, OY - p1.y * alpha, OX + p2.x * alpha + 20, OY - p2.y * alpha));
		
		}

		public void drawLine(Point p1, Point p2, Color color) {
			Graphics2D g2 = (Graphics2D) getGraphics();
			g2.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
			g2.setStroke(new BasicStroke(2));
			g2.setColor(color);
			g2.draw(new Line2D.Double(OX + p1.x * alpha + 20, OY - p1.y * alpha, OX + p2.x * alpha + 20, OY - p2.y * alpha));
		}

		public void drawLines(ArrayList<Point> points, LinkedList<Point> pointsToVisit, Color lineColor,
				Color dotColor) {
			Graphics2D g2 = (Graphics2D) getGraphics();
			g2.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
			g2.setColor(lineColor);
			Stroke dashed = new BasicStroke(3, BasicStroke.CAP_BUTT, BasicStroke.JOIN_BEVEL, 0, new float[]{21, 9 ,3, 9}, 10);
//			g2.setStroke(dashed);
			g2.setStroke(new BasicStroke(2.0f));
			
			for (int i = 0; i < points.size() - 1; i++) {
				
				g2.draw(new Line2D.Double(OX + points.get(i).x * alpha + 20, OY - points.get(i).y * alpha,
						OX + points.get(i + 1).x * alpha + 20, OY - points.get(i + 1).y * alpha));
				
			}

			g2.setColor(Color.RED);
			g2.setFont(new Font("VnArial Bold", Font.PLAIN, 20));
			g2.drawString("S", (int) (MyCanvas.OX + points.get(0).x * alpha),(int) (MyCanvas.OY - points.get(0).y * alpha + 5));
				
			for (int i = 0; i < points.size(); i++) {
				Point pt = new Point(points.get(i).x, points.get(i).y);
				if (pt.indexInSet(pointsToVisit) == -1) {
					g2.setColor(dotColor);
					g2.fill(new Rectangle2D.Double(OX + points.get(i).x * alpha - 3 + 20, OY - points.get(i).y * alpha - 3, 5, 5));
					g2.setColor(Color.RED);
					g2.drawString("T", (int)(MyCanvas.OX + points.get(12).x * alpha - 3 + 30), (int)(MyCanvas.OY - points.get(12).y * alpha + 5));
				} else
					g2.fill(new Ellipse2D.Double(OX + points.get(i).x * alpha - 4 + 20, OY - points.get(i).y * alpha - 4, 8,
							8));
			}
		}
		public void drawLinesVMOES(ArrayList<Point> points, LinkedList<Point> pointsToVisit, Color lineColor,
				Color dotColor) {
			Graphics2D g2 = (Graphics2D) getGraphics();
			g2.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
			g2.setColor(lineColor);
			Stroke dashed = new BasicStroke(3, BasicStroke.CAP_SQUARE, BasicStroke.JOIN_MITER, 10.0f, new float[]{9}, 0);
//			g2.setStroke(dashed);
			g2.setStroke(new BasicStroke(2));
			for (int i = 0; i < points.size() - 1; i++) {
				g2.draw(new Line2D.Double(OX + points.get(i).x * alpha + 20, OY - points.get(i).y * alpha,
						OX + points.get(i + 1).x * alpha + 20, OY - points.get(i + 1).y * alpha));
			}

			g2.setColor(Color.RED);
			
			for (int i = 0; i < points.size(); i++) {
				Point pt = new Point(points.get(i).x, points.get(i).y);
				if (pt.indexInSet(pointsToVisit) == -1) {
					g2.setColor(dotColor);
					g2.fill(new Polygon( new int[] {(int)(OX + points.get(i).x * alpha + 20),(int)(OX + points.get(i).x * alpha-3 + 20),(int)(OX + points.get(i).x * alpha+3 + 20)},
							 			 new int[] {(int)(OY - points.get(i).y * alpha+3),(int)(OY - points.get(i).y * alpha-3),(int)(OY - points.get(i).y * alpha-3)},
							 			 3));
					g2.setColor(Color.RED);
				} else
					g2.fill(new Ellipse2D.Double(OX + points.get(i).x * alpha - 4 + 20, OY - points.get(i).y * alpha - 4, 8,
							8));
			}
		}
		public void drawLinesMOPSO(ArrayList<Point> points, LinkedList<Point> pointsToVisit, Color lineColor,
				Color dotColor) {
			Graphics2D g2 = (Graphics2D) getGraphics();
			g2.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
			g2.setColor(lineColor);
			Stroke dashed = new BasicStroke(3, BasicStroke.CAP_ROUND, BasicStroke.JOIN_ROUND, 0, new float[]{9}, 10);
//			g2.setStroke(dashed);
			g2.setStroke(new BasicStroke(2));
			for (int i = 0; i < points.size() - 1; i++) {
				g2.draw(new Line2D.Double(OX + points.get(i).x * alpha + 20, OY - points.get(i).y * alpha,
						OX + points.get(i + 1).x * alpha + 20, OY - points.get(i + 1).y * alpha));
			}

			g2.setColor(Color.RED);
			
			for (int i = 0; i < points.size(); i++) {
				Point pt = new Point(points.get(i).x, points.get(i).y);
				if (pt.indexInSet(pointsToVisit) == -1) {
					g2.setColor(dotColor);
					g2.fillOval((int)(OX + points.get(i).x * alpha + 20), (int)(OY - points.get(i).y * alpha - 2), 5, 5);
					g2.setColor(Color.RED);
				} else
					g2.fill(new Ellipse2D.Double(OX + points.get(i).x * alpha - 4 + 20, OY - points.get(i).y * alpha - 4, 8,
							8));
			}
		}
		public void drawLinesNSGAII(ArrayList<Point> points, LinkedList<Point> pointsToVisit, Color lineColor,
				Color dotColor) {
			Graphics2D g2 = (Graphics2D) getGraphics();
			g2.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
			g2.setColor(lineColor);
			g2.setStroke(new BasicStroke(2.0f));
			for (int i = 0; i < points.size() - 1; i++) {
				g2.draw(new Line2D.Double(OX + points.get(i).x * alpha + 20, OY - points.get(i).y * alpha,
						OX + points.get(i + 1).x * alpha + 20, OY - points.get(i + 1).y * alpha));
			}

			g2.setColor(Color.RED);
			
			for (int i = 0; i < points.size(); i++) {
				Point pt = new Point(points.get(i).x, points.get(i).y);
				if (pt.indexInSet(pointsToVisit) == -1) {
					g2.setColor(dotColor);
					g2.drawLine((int)(OX + points.get(i).x * alpha-3 + 20), (int)(OY - points.get(i).y * alpha-3), (int)(OX + points.get(i).x * alpha + 3 + 20), (int)(OY - points.get(i).y * alpha+3));
					g2.setColor(Color.RED);
				} else
					g2.fill(new Ellipse2D.Double(OX + points.get(i).x * alpha - 4 + 20, OY - points.get(i).y * alpha - 4, 8,
							8));
			}
		}
		

		public void drawLineWithMiddle(Point pt1, Point pt2) {
			Graphics2D g2 = (Graphics2D) getGraphics();
			g2.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
			g2.setColor(Color.BLUE);
			g2.draw(new Line2D.Double(OX + pt1.x * alpha, OY - pt1.y * alpha, OX + pt2.x * alpha, OY - pt2.y * alpha));

			g2.setColor(Color.RED);
			Point pt = new Point((pt1.x + pt2.x) / 2, (pt1.y + pt2.y) / 2);
			g2.fill(new Ellipse2D.Double(OX + pt.x * alpha - 3, OY - pt.y * alpha - 3, 6, 6));
		}

		public void drawLineWithoutMiddle(Point pt1, Point pt2) {
			Graphics2D g2 = (Graphics2D) getGraphics();
			g2.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
			g2.setColor(Color.BLUE);
			g2.draw(new Line2D.Double(OX + pt1.x * alpha, OY - pt1.y * alpha, OX + pt2.x * alpha, OY - pt2.y * alpha));
		}
	}

	public static void main(String[] args) throws FileNotFoundException {
		double ratio = 0.9;
		boolean draw_borders = true;
		GUIRobotics gui = new GUIRobotics(600, 110, 11);
		gui.generateEnvironment("obstacles.txt", ratio, draw_borders);

		System.setIn(new FileInputStream("maklink.txt"));
		Scanner sc = new Scanner(System.in);

		sc.nextLine();
		String string = sc.nextLine();
		while (!string.equals("-1")) {

			// Read path planning
			String numbers[] = string.replaceAll(",", "").replaceAll("\\(", "").replaceAll("\\)", "").split("\\s+");

			ArrayList<Point> points = new ArrayList<>();
			for (int i = 0; i < numbers.length / 2; i++) {
				points.add(new Point(Double.parseDouble(numbers[2 * i]), Double.parseDouble(numbers[2 * i + 1])));
			}

			for (int i = 0; i < points.size() / 2; i++) {
				gui.canvas.drawLineWithoutMiddle(points.get(2 * i), points.get(2 * i + 1));
			}

			string = sc.nextLine();
		}
		sc.close();
	}

}
