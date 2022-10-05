#pragma once

#include <iostream>

#include "s21_graph.h"
#include "s21_graph_algorithms.h"
#include "s21_matrix_oop.h"

class ConsoleInterface {
 public:
  void start();

 private:
  void read();
  void menu();
  void DFS();
  void BFS();
  void nextChange();
  template <typename T>
  void printVector(std::vector<T> path);
  template <typename T>
  void printMatrix(S21Matrix<T> res);
  void dijkstra();
  void floydWarshallAlgo();
  void primAlgorithmAlgo();
  void exportDot();
  void ant();
  void antVsOtherAlgo();

  //   void antVsOther();

  Graph graph_;
  GraphAlgorithms algo_;

  std::string end = "\u001b[0m";
  std::string end1 = "\u001b[0m\n";
  std::string style1 = "\u001b[1;35;5;117m";
  std::string style2 = "\u001b[1;33;5;117m";
  std::string style3 = "\u001b[1;32;5;117m";
  std::string style4 = "\u001b[1;31;5;117m";
};
