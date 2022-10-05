#pragma once

#include <random>
#include <set>
#include <vector>

#include "s21_graph.h"

struct TsmResult {
  std::vector<int> visit;
  double distance;
};

class Ant {
 public:
  explicit Ant(int startVertex) { tabu.insert(startVertex); }
  std::set<int> tabu;

  double rand() {
    std::random_device rd;
    std::default_random_engine eng(rd());
    std::uniform_real_distribution<double> distr(0, 1);
    return distr(eng);
  }
};

class GraphAlgorithms {
 public:
  //  Functions
  std::vector<int> depthFirstSearch(Graph &graph, int startVertex);
  std::vector<int> breadthFirstSearch(Graph &graph, int startVertex);
  std::vector<double> algoD(Graph &graph, int vertex1, int vertex2);
  double getShortestPathBetweenVertices(Graph &graph, int vertex1, int vertex2);
  S21Matrix<double> getShortestPathsBetweenAllVertices(Graph &graph);
  S21Matrix<double> getLeastSpanningTree(Graph &graph);
  TsmResult solveTravelingSalesmanProblem(Graph &graph);
  TsmResult annealingAlgorithm(Graph &graph);
  TsmResult SolveTspBranchAndBound(Graph &graph);
  TsmResult helperSolveTravel(Graph &graph);

 private:
  S21Matrix<double> oneAntPath(Graph &graph,
                               const S21Matrix<double> &countFeromn,
                               int curVertex);
  std::vector<double> calculateProbability(Graph &graph, int curVertex,
                                           double *sumWish,
                                           const S21Matrix<double> &countFeromn,
                                           Ant ant);
  double checkDistance(Graph &graph, const double &curDistance,
                       const S21Matrix<double> &countFeromn,
                       std::vector<int> *path);

  int randomVal(int val);
  double distanceLen(const std::vector<int> &path, Graph &graph);
  double ReduceMatrix(S21Matrix<double> *matrix);
  double BranchAndBound(S21Matrix<double> matrix, double total_cost,
                        int where_am_i, const int kCities, int *been_there,
                        int *travel_path, int pos);
  const S21Matrix<double> matrixTransform(S21Matrix<double> new_matrix,
                                          S21Matrix<double> matrix,
                                          const int &where_am_i, const int &i,
                                          std::vector<double> *row,
                                          std::vector<double> *col,
                                          double *travel_cost, int step);
  //
};
