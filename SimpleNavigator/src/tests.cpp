#include <gtest/gtest.h>

#include "s21_graph.h"
#include "s21_graph_algorithms.h"

TEST(graph, loadGraphFromFile1) {
  Graph graph;
  EXPECT_ANY_THROW(graph.loadGraphFromFile("gen0.txt"));
}

TEST(graph, loadGraphFromFile2) {
  Graph graph;
  EXPECT_NO_THROW(graph.loadGraphFromFile("gen1.txt"));
}

TEST(graph, loadGraphFromFile3) {
  Graph graph;
  Graph graph1 = graph;
  EXPECT_NO_THROW(auto graph3 = std::move(graph1));
}

TEST(graph, exportGraphToDot1) {
  Graph graph;
  graph.loadGraphFromFile("gen5.txt");
  EXPECT_NO_THROW(graph.exportGraphToDot("gen5"));
}

TEST(graph, exportGraphToDot2) {
  Graph graph;
  graph.loadGraphFromFile("Agen3.txt");
  EXPECT_NO_THROW(graph.exportGraphToDot("Agen3"));
}

TEST(graph, exportGraphToDot3) {
  Graph graph;
  EXPECT_NO_THROW(graph.loadGraphFromFile("gen5.txt"));
  EXPECT_NO_THROW(graph.exportGraphToDot(""));
}

TEST(graph_algorithms, depthFirstSearch) {
  Graph graph;
  graph.loadGraphFromFile("gen5.txt");
  GraphAlgorithms ga;
  std::vector<int> path = ga.depthFirstSearch(graph, 1);
  ASSERT_FALSE(path.empty());
}

TEST(graph_algorithms, breadthFirstSearch) {
  Graph graph;
  graph.loadGraphFromFile("gen5.txt");
  GraphAlgorithms ga;
  std::vector<int> path = ga.breadthFirstSearch(graph, 1);
  ASSERT_FALSE(path.empty());
}

TEST(graph_algorithms, getShortestPathBetweenVertices) {
  Graph graph;
  graph.loadGraphFromFile("gen5.txt");
  GraphAlgorithms ga;
  double W = ga.getShortestPathBetweenVertices(graph, 1, 5);
  ASSERT_GE(W, 0);
}

TEST(graph_algorithms, getShortestPathsBetweenAllVertices) {
  Graph graph;
  graph.loadGraphFromFile("gen5.txt");
  GraphAlgorithms ga;
  S21Matrix<double> matrix = ga.getShortestPathsBetweenAllVertices(graph);
  ASSERT_EQ(matrix.get_cols(), 5);
  ASSERT_EQ(matrix.get_rows(), 5);
}

TEST(graph_algorithms, getLeastSpanningTree) {
  Graph graph;
  graph.loadGraphFromFile("gen5.txt");
  GraphAlgorithms ga;
  S21Matrix<double> matrix = ga.getLeastSpanningTree(graph);
  ASSERT_EQ(matrix.get_cols(), 5);
  ASSERT_EQ(matrix.get_rows(), 5);
}
TEST(graph_algorithms, solveTravelingSalesmanProblem) {
  Graph graph;
  graph.loadGraphFromFile("gen5.txt");
  GraphAlgorithms ga;
  TsmResult result;
  result = ga.solveTravelingSalesmanProblem(graph);
  ASSERT_GE(result.distance, 5);
  ASSERT_EQ(result.visit.size(), 6);
}
TEST(graph_algorithms, annealingAlgorithm) {
  Graph graph;
  graph.loadGraphFromFile("gen5.txt");
  GraphAlgorithms ga;
  TsmResult result;
  result = ga.annealingAlgorithm(graph);
  ASSERT_GE(result.distance, 0);
  ASSERT_EQ(result.visit.size(), 5);
}
TEST(graph_algorithms, SolveTspBranchAndBound) {
  Graph graph;
  graph.loadGraphFromFile("gen5.txt");
  GraphAlgorithms ga;
  TsmResult result;
  result = ga.SolveTspBranchAndBound(graph);
  ASSERT_GE(result.distance, 0);
  ASSERT_EQ(result.visit.size(), 6);
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
