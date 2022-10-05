#pragma once

#include <fstream>
#include <iostream>

#include "s21_matrix_oop.h"

class Graph {
 public:
  Graph();
  Graph(const Graph& other);
  Graph(Graph&& other);
  ~Graph();

  Graph& operator=(const Graph& other);
  Graph& operator=(Graph&& other);
  double& operator()(int row, int column) {
    return matrix_->operator()(row, column);
  }

  void loadGraphFromFile(const std::string& filename);
  void exportGraphToDot(const std::string& filename);
  int getSize() { return size_; }
  S21Matrix<double>* getMatrix() { return matrix_; }

 private:
  S21Matrix<double>* matrix_;
  int size_;
  bool isWeightedGraph_;  //  Взвешенный граф
  bool isDirectedGraph_;  //  Ориентированный граф
};
