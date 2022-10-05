#include "s21_graph_algorithms.h"

#include <cmath>

#include "s21_queue.h"
#include "s21_stack.h"

std::vector<int> GraphAlgorithms::depthFirstSearch(Graph &graph,
                                                   int startVertex) {
  if (startVertex <= 0 || graph.getSize() <= 0 ||
      startVertex > graph.getSize()) {
    throw std::invalid_argument("Ivalid data on DFS function");
  }
  s21::stack<int> graphVertex;
  std::vector<bool> vertexVisits(graph.getSize());
  std::vector<int> rez;
  graphVertex.push(startVertex - 1);
  vertexVisits[startVertex - 1] = true;
  while (!graphVertex.empty()) {
    int curVert = graphVertex.pop();
    rez.push_back(curVert + 1);
    for (int i = 0; i < graph.getSize(); ++i) {
      if (graph(curVert, i) > 0.0 && vertexVisits[i] == false) {
        graphVertex.push(i);
        vertexVisits[i] = true;
      }
    }
  }
  return rez;
}

std::vector<int> GraphAlgorithms::breadthFirstSearch(Graph &graph,
                                                     int startVertex) {
  if (startVertex <= 0 || graph.getSize() <= 0 ||
      startVertex > graph.getSize()) {
    throw std::invalid_argument("Ivalid data on BFS function");
  }
  s21::queue<int> graphVertex;
  std::vector<bool> vertexVisits(graph.getSize());
  std::vector<int> rez;
  graphVertex.push(startVertex - 1);
  vertexVisits[startVertex - 1] = true;
  while (!graphVertex.empty()) {
    int curVert = graphVertex.pop();
    rez.push_back(curVert + 1);
    for (int i = 0; i < graph.getSize(); ++i) {
      if (graph(curVert, i) > 0.0 && vertexVisits[i] == false) {
        graphVertex.push(i);
        vertexVisits[i] = true;
      }
    }
  }
  return rez;
}

std::vector<double> GraphAlgorithms::algoD(Graph &graph, int vertex1,
                                           int vertex2) {
  if (vertex1 < 1 || vertex2 < 1 || graph.getSize() < 1 ||
      vertex1 > graph.getSize() || vertex2 > graph.getSize()) {
    throw std::invalid_argument("Ivalid data on algoD function");
  }
  std::vector<double> label, prevVert;
  std::vector<bool> constLabel;
  for (int i = 0; i < graph.getSize(); ++i) {
    label.push_back(INFINITY);
    constLabel.push_back(false);
    prevVert.push_back(NAN);
  }
  int curVert = 0, count = 1;
  label[0] = curVert;
  constLabel[0] = true;
  prevVert[0] = 0.0;
  while (count != graph.getSize()) {
    double minPos;
    double minVal = __DBL_MAX__;
    for (int i = 0; i < graph.getSize(); ++i) {
      if (graph(curVert, i) > 0.0 && constLabel[i] == false) {
        if (label[i] > (label[curVert] + graph(curVert, i))) {
          label[i] = label[curVert] + graph(curVert, i);
          prevVert[i] = curVert;
        }
        if (label[i] < minVal) {
          minVal = label[i];
          minPos = i;
        }
      }
    }
    for (int i = 0; i < graph.getSize(); ++i) {
      if (constLabel[i] == false) {
        if (label[i] < minVal) {
          minVal = label[i];
          minPos = i;
        }
      }
    }
    constLabel[minPos] = true;
    curVert = minPos;
    ++count;
  }
  if (vertex1 > vertex2) {
    std::swap(vertex1, vertex2);
  }
  std::vector<double> rezPath;
  rezPath.insert(rezPath.begin(), vertex2);
  double z = prevVert[vertex2 - 1];
  while (z != vertex1 - 1) {
    rezPath.insert(rezPath.begin(), z + 1);
    z = prevVert[z];
  }
  rezPath.insert(rezPath.begin(), z + 1);
  return rezPath;
}

double GraphAlgorithms::getShortestPathBetweenVertices(Graph &graph,
                                                       int vertex1,
                                                       int vertex2) {
  std::vector<double> rez = algoD(graph, vertex1, vertex2);
  double distance = 0.0;
  for (int i = 0; i < (int)rez.size() - 1; ++i) {
    distance += graph((int)(rez[i] - 1), (int)(rez[i + 1] - 1));
  }
  return distance;
}

S21Matrix<double> GraphAlgorithms::getShortestPathsBetweenAllVertices(
    Graph &graph) {
  S21Matrix<double> *myMatrix = graph.getMatrix();
  int size = graph.getSize();
  S21Matrix<double> newMatrix(*myMatrix);
  for (int i = 0; i < size; i++) {
    for (int j = 0; j < size; j++) {
      if (newMatrix(i, j) == 0 && i != j) newMatrix(i, j) = INFINITY;
    }
  }
  for (int k = 0; k < size; k++) {
    for (int i = 0; i < size; i++) {
      for (int j = 0; j < size; j++) {
        newMatrix(i, j) =
            std::min(newMatrix(i, j), (newMatrix(i, k) + newMatrix(k, j)));
      }
    }
  }
  return newMatrix;
}

S21Matrix<double> GraphAlgorithms::getLeastSpanningTree(Graph &graph) {
  if (graph.getSize() < 1) {
    throw std::invalid_argument("Ivalid data on getLeastSpanningTree function");
  }
  std::set<int> visit, unvisit;
  S21Matrix<double> rez(graph.getSize(), graph.getSize());
  for (int i = 0; i < graph.getSize(); ++i) {
    unvisit.insert(i);
  }
  visit.insert(0);
  unvisit.erase(0);
  while (!unvisit.empty()) {
    double min = __DBL_MAX__;
    std::pair<int, int> pos;
    std::set<int>::iterator visits = visit.begin();
    while (visits != visit.end()) {
      for (int i = 0; i < graph.getSize(); ++i) {
        if (visit.count(i) == 0 && graph(*visits, i) > 0) {
          if (graph(*visits, i) < min) {
            min = graph(*visits, i);
            pos = {*visits, i};
          }
        }
      }
      ++visits;
    }
    visit.insert(pos.second);
    unvisit.erase(pos.second);
    rez(pos.first, pos.second) = min;
  }
  return rez;
}

TsmResult GraphAlgorithms::solveTravelingSalesmanProblem(Graph &graph) {
  TsmResult rez = helperSolveTravel(graph);
  int count = 0;
  while (count++ < 5) {
    TsmResult tmp = helperSolveTravel(graph);
    if (tmp.distance < rez.distance) rez = tmp;
  }
  for (int i = 0; i < (int)rez.visit.size(); ++i) {
    rez.visit[i] += 1;
  }

  return rez;
}

int GraphAlgorithms::randomVal(int val) {
  std::random_device rd;
  std::default_random_engine eng(rd());
  std::uniform_int_distribution<int> distr(0, val);
  return distr(eng);
}

TsmResult GraphAlgorithms::helperSolveTravel(Graph &graph) {
  if (graph.getSize() < 2) {
    throw std::invalid_argument(
        "Error on solveTravelingSalesmanProblem. Graph < 2");
  }
  S21Matrix<double> countFeromn(graph.getSize(), graph.getSize());
  countFeromn.fillMatrix(0.2);
  S21Matrix<double> deltaTau(graph.getSize(), graph.getSize());
  std::vector<int> path;
  double curDistance = 0.0;
  const int CYCLE_STEPS = 40, ANT_COUNTS = 200;
  const double VAPE = 0.5;
  int count = 0;
  while (count < CYCLE_STEPS) {
    int stepsAnt = 0;
    while (++stepsAnt < ANT_COUNTS) {
      int currVertex = randomVal(graph.getSize() - 1);
      deltaTau += oneAntPath(graph, countFeromn, currVertex);
    }
    countFeromn *= VAPE;
    countFeromn += deltaTau;
    deltaTau.fillMatrix(0.0);
    ++count;
    if (count == CYCLE_STEPS) {
      double newDistance =
          checkDistance(graph, curDistance, countFeromn, &path);
      if (newDistance < 0.001) {
        count = 0;
      } else {
        curDistance = newDistance;
      }
    }
  }
  TsmResult rez{path, curDistance};
  return rez;
}

S21Matrix<double> GraphAlgorithms::oneAntPath(
    Graph &graph, const S21Matrix<double> &countFeromn, int curVertex) {
  const double Q = 10.0;
  S21Matrix<double> deltaTau(graph.getSize(), graph.getSize());
  Ant ant(curVertex);
  while ((int)ant.tabu.size() < graph.getSize()) {
    ant.tabu.insert(curVertex);
    double sumWish = 0.0;
    std::vector<double> probability =
        calculateProbability(graph, curVertex, &sumWish, countFeromn, ant);
    int tmpVertex = curVertex;
    double step = ant.rand();
    for (int n = 0; n < (int)probability.size(); ++n) {
      if (ant.tabu.count(n) == 0) {
        step -= probability[n];
        if (step < 0) {
          curVertex = n;
          break;
        }
      }
    }
    if (curVertex == tmpVertex && (int)ant.tabu.size() < graph.getSize()) {
      throw std::invalid_argument("Error on oneAntPath. Graph is incomplete");
    }
    if (graph(tmpVertex, curVertex) > 0) {
      deltaTau(tmpVertex, curVertex) += Q / graph(tmpVertex, curVertex);
    }
  }
  return deltaTau;
}

std::vector<double> GraphAlgorithms::calculateProbability(
    Graph &graph, int curVertex, double *sumWish,
    const S21Matrix<double> &countFeromn, Ant ant) {
  const int ALPHA = 5, BETTA = 2;
  std::vector<double> accesableWish;
  for (int j = 0; j < graph.getSize(); ++j) {
    double n = 0.0;
    double tau = 0.0;
    if (graph(curVertex, j) > 0 && ant.tabu.count(j) == 0) {
      n = 1.0 / graph(curVertex, j);
      tau = countFeromn(curVertex, j);
      *sumWish += (pow(tau, ALPHA) * pow(n, BETTA));
    }
    accesableWish.push_back((pow(tau, ALPHA) * pow(n, BETTA)));
  }
  std::vector<double> probability;
  for (int k = 0; k < graph.getSize(); ++k) {
    if (graph(curVertex, k) > 0.0 && ant.tabu.count(k) == 0) {
      probability.push_back((accesableWish[k] / *sumWish));
    } else {
      probability.push_back(0.0);
    }
  }

  return probability;
}

double GraphAlgorithms::checkDistance(Graph &graph, const double &curDistance,
                                      const S21Matrix<double> &countFeromn,
                                      std::vector<int> *path) {
  int iter = 0;
  std::set<int> check;
  double newDistance = 0.0;
  path->push_back(iter);
  check.insert(iter);
  while ((int)path->size() < graph.getSize()) {
    int tmpIter = 0;
    double max = 0.0;
    for (int z = 0; z < graph.getSize(); ++z) {
      if (max < countFeromn(iter, z) && check.count(z) == 0) {
        max = countFeromn(iter, z);
        tmpIter = z;
      }
    }
    newDistance += graph(iter, tmpIter);
    iter = tmpIter;
    path->push_back(iter);
    check.insert(iter);
  }
  if (fabs(newDistance - curDistance) < 0.9) {
    newDistance = 0.00000;
    path->clear();
  } else {
    int tmp = path->operator[](path->size() - 1);
    newDistance += graph(tmp, 0);
    path->push_back(0);
  }
  return newDistance;
}

TsmResult GraphAlgorithms::annealingAlgorithm(Graph &graph) {
  const double SPEED = 0.85;
  double currentTemp = 400.0;
  const double TEMP_MIN = 0.1;
  const int CONS_LEN = graph.getSize() * 200;  // константа всего пути
  TsmResult buf{{}, 0};
  TsmResult rez{{}, 0};
  int startRandomVertex = randomVal(graph.getSize() - 1);
  for (int i = 0; i < graph.getSize(); ++i) {
    rez.visit.push_back(startRandomVertex);
    startRandomVertex = (startRandomVertex + 1) % graph.getSize();
  }
  rez.distance = distanceLen(rez.visit, graph);
  buf = rez;
  while (currentTemp > TEMP_MIN) {
    for (int i = 0; i < CONS_LEN; i++) {
      std::swap(buf.visit[randomVal(graph.getSize() - 1)],
                buf.visit[randomVal(graph.getSize() - 1)]);
      buf.distance = distanceLen(buf.visit, graph);
      if (rez.distance > buf.distance ||
          ((int)(exp((rez.distance - buf.distance) / currentTemp) * 100) >
           (randomVal(100)))) {
        rez = buf;
      } else {
        buf = rez;
      }
    }
    currentTemp *= SPEED;
  }
  return rez;
}

double GraphAlgorithms::distanceLen(const std::vector<int> &path,
                                    Graph &graph) {
  double dis = 0;
  for (int i = 1; i < graph.getSize(); i++) {
    dis += graph(path[i - 1], path[i]);
  }
  dis += graph(path[graph.getSize() - 1], path[0]);
  return dis;
}

TsmResult GraphAlgorithms::SolveTspBranchAndBound(Graph &graph) {
  TsmResult result;
  const int kCities = graph.getSize();
  int travel_path[kCities] = {0}, pos = 0;
  int been_there[kCities] = {0};
  been_there[0] = 1;
  S21Matrix<double> reduced_matrix(*graph.getMatrix());
  double total_cost = ReduceMatrix(&reduced_matrix);
  BranchAndBound(reduced_matrix, total_cost, 0, kCities, been_there,
                 travel_path, pos);
  double path_length = 0;
  int prevIndex = 1;
  result.visit.push_back(1);
  for (int i = 1; i < kCities; ++i) {
    path_length += graph(prevIndex - 1, travel_path[i] - 1);
    prevIndex = travel_path[i];
    result.visit.push_back(travel_path[i]);
  }
  result.visit.push_back(1);
  path_length += graph(prevIndex - 1, 0);
  result.distance = path_length;
  return result;
}

const S21Matrix<double> GraphAlgorithms::matrixTransform(
    S21Matrix<double> new_matrix, S21Matrix<double> matrix,
    const int &where_am_i, const int &i, std::vector<double> *row,
    std::vector<double> *col, double *travel_cost, int step) {
  S21Matrix<double> out_matrix(new_matrix.get_rows(), new_matrix.get_cols());
  for (int k = 0; k < new_matrix.get_rows(); k++) {
    if (step == 2) {
      row->operator[](k) = new_matrix(k, 0);
    } else if (step == 4) {
      col->operator[](k) = new_matrix(0, k);
    }
    for (int l = 0; l < new_matrix.get_cols(); l++) {
      if (step == 1) {
        if (k == where_am_i || l == i || (l == where_am_i && k == i))
          out_matrix(k, l) = INFINITY;
        else
          out_matrix(k, l) = matrix(k, l);
      } else if (step == 2) {
        if (std::isinf(new_matrix(k, l)) != true) {
          if (new_matrix(k, l) < row->operator[](k))
            row->operator[](k) = new_matrix(k, l);
        }
      } else if (step == 3) {
        if (std::isinf(new_matrix(k, l)) != true) {
          out_matrix(k, l) = new_matrix(k, l) - row->operator[](k);
        } else {
          out_matrix(k, l) = new_matrix(k, l);
        }
      } else if (step == 4) {
        if (std::isinf(new_matrix(l, k)) != true) {
          if (new_matrix(l, k) < col->operator[](k))
            col->operator[](k) = new_matrix(l, k);
        }
      } else if (step == 5) {
        if (std::isinf(new_matrix(k, l)) != true) {
          out_matrix(k, l) = new_matrix(k, l) - col->operator[](k);
        } else {
          out_matrix(k, l) = new_matrix(k, l);
        }
      }
    }
    if (step == 2) {
      if (std::isinf(row->operator[](k)) == true) row->operator[](k) = 0;
      *travel_cost += row->operator[](k);
    } else if (step == 4) {
      if (std::isinf(col->operator[](k)) == true) col->operator[](k) = 0;
      *travel_cost += col->operator[](k);
    }
  }
  return out_matrix;
}

double GraphAlgorithms::BranchAndBound(S21Matrix<double> matrix,
                                       double total_cost, int where_am_i,
                                       const int kCities, int *been_there,
                                       int *travel_path, int pos) {
  double least_cost_bound = INFINITY;
  double travel_cost = 0;
  S21Matrix<double> new_matrix(matrix);
  std::vector<double> col(matrix.get_rows()), row(matrix.get_rows());
  int keep_node;
  S21Matrix<double> keep_matrix(kCities, kCities);
  been_there[where_am_i] = 1;
  int recurtion_times = 0;

  travel_path[pos] = where_am_i + 1;
  pos++;

  for (int i = 0; i < kCities; i++) {
    if (been_there[i] != 1) {
      recurtion_times++;
      new_matrix = matrixTransform(new_matrix, matrix, where_am_i, i, &row,
                                   &col, &travel_cost, 1);
      matrixTransform(new_matrix, matrix, where_am_i, i, &row, &col,
                      &travel_cost, 2);
      if (travel_cost > 0)
        new_matrix = matrixTransform(new_matrix, matrix, where_am_i, i, &row,
                                     &col, &travel_cost, 3);

      matrixTransform(new_matrix, matrix, where_am_i, i, &row, &col,
                      &travel_cost, 4);
      if (travel_cost > 0)
        new_matrix = matrixTransform(new_matrix, matrix, where_am_i, i, &row,
                                     &col, &travel_cost, 5);

      if (least_cost_bound > total_cost + matrix(where_am_i, i) + travel_cost) {
        least_cost_bound = total_cost + matrix(where_am_i, i) + travel_cost;
        keep_node = i;
        keep_matrix = new_matrix;
      }
    }
    travel_cost = 0;
  }
  if (recurtion_times < 2) {
    travel_path[pos] = keep_node + 1;
    return least_cost_bound;
  } else {
    total_cost = BranchAndBound(keep_matrix, least_cost_bound, keep_node,
                                kCities, been_there, travel_path, pos);
  }
  return total_cost;
}

double GraphAlgorithms::ReduceMatrix(S21Matrix<double> *matrix) {
  int kCities = matrix->get_rows();
  std::vector<double> minimalInRows;
  std::vector<double> minimalInCols;
  double minimal = INFINITY;
  double total_cost = 0;
  for (int i = 0; i < kCities; ++i) {
    minimal = INFINITY;
    for (int j = 0; j < kCities; ++j) {
      if (i == j) matrix->operator()(i, j) = INFINITY;
      if (matrix->operator()(i, j) < minimal) {
        minimal = matrix->operator()(i, j);
      }
    }
    total_cost += minimal;
    minimalInRows.push_back(minimal);
  }
  for (int i = 0; i < kCities; ++i) {
    for (int j = 0; j < kCities; ++j) {
      if (std::isinf(matrix->operator()(i, j)) == true) {
        continue;
      } else {
        if (std::isinf(minimalInRows[i]) == false)
          matrix->operator()(i, j) =
              matrix->operator()(i, j) - minimalInRows[i];
      }
    }
  }
  for (int i = 0; i < kCities; ++i) {
    minimal = INFINITY;
    for (int j = 0; j < kCities; ++j) {
      if (matrix->operator()(j, i) < minimal) {
        minimal = matrix->operator()(j, i);
      }
    }
    total_cost += minimal;
    minimalInCols.push_back(minimal);
  }
  for (int i = 0; i < kCities; ++i) {
    for (int j = 0; j < kCities; ++j) {
      if (std::isinf(matrix->operator()(j, i)) == true) {
        continue;
      } else {
        if (std::isinf(minimalInCols[i]) == false)
          matrix->operator()(j, i) =
              matrix->operator()(j, i) - minimalInCols[i];
      }
    }
  }
  return total_cost;
}
