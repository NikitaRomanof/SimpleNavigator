#include "console_interface.h"

#include <time.h>

void ConsoleInterface::start() {
  std::cout << style1 << "\n*******************************************" << end
            << std::endl;
  std::cout << style2 << "                SIMPLE NAVIGATOR             " << end
            << std::endl;
  std::cout << style1 << "*******************************************" << end
            << std::endl;
  read();
}

void ConsoleInterface::read() {
  while (true) {
    std::cout << style3 << " enter file name           " << end << std::endl;
    std::cout << std::endl;
    std::string str;
    std::cin >> str;
    try {
      graph_.loadGraphFromFile("../datasets/" + str);
      break;
    } catch (std::exception &e) {
      std::cerr << style4 << "\nError: " << e.what() << end << "\n\n";
    }
  }
  menu();
}

void ConsoleInterface::menu() {
  std::cout << style1 << "*********************************************" << end
            << std::endl;
  std::cout << style2 << "                   MENU                      " << end
            << std::endl;
  std::cout << style1 << "*********************************************" << end1
            << std::endl;
  std::cout << style3
            << " 1. ----> DFS                                         " << end
            << std::endl;
  std::cout << style3
            << " 2. ----> BFS                                         " << end
            << std::endl;
  std::cout << style3
            << " 3. ----> Dijkstra's algorithm                        " << end
            << std::endl;
  std::cout << style3
            << " 4. ----> Floyd-Warshall algorithm                    " << end
            << std::endl;
  std::cout << style3
            << " 5. ----> Prim's algorithm                            " << end
            << std::endl;
  std::cout << style3
            << " 6. ----> Ant algorithm                               " << end
            << std::endl;
  std::cout
      << style3
      << " 7. ----> Ant algorithm vs annealingAlgorithm vs branchAndBound "
      << end << std::endl;
  std::cout << style3
            << " 8. ----> Save graph to .dot                          " << end
            << std::endl;
  std::cout << style3
            << " 9. ----> Upload new file                             " << end
            << std::endl;
  std::cout << style3
            << " 0. ----> Exit from the program                       " << end1
            << std::endl;
  int choice;
  std::cin >> choice;
  if (choice == 1) {
    DFS();
  } else if (choice == 2) {
    BFS();
  } else if (choice == 3) {
    dijkstra();
  } else if (choice == 4) {
    floydWarshallAlgo();
  } else if (choice == 5) {
    primAlgorithmAlgo();
  } else if (choice == 6) {
    ant();
  } else if (choice == 7) {
    antVsOtherAlgo();
  } else if (choice == 8) {
    exportDot();
  } else if (choice == 9) {
    read();
  } else if (choice == 0) {
    return;
  } else {
    std::cout << style4 << "\nError: incorrect input, try again " << end
              << std::endl;
  }
}

void ConsoleInterface::DFS() {
  std::cout << style4 << "\nenter vertex number from 1 to " << graph_.getSize()
            << " to start search " << end1 << std::endl;
  int vertex;
  std::cin >> vertex;
  try {
    std::vector<int> vecWay;
    vecWay = algo_.depthFirstSearch(graph_, vertex);
    printVector(vecWay);
    nextChange();
  } catch (std::exception &e) {
    std::cerr << style4 << "\nError: " << e.what() << end << "\n";
    DFS();
  }
}

void ConsoleInterface::BFS() {
  std::cout << style4 << "\nenter vertex number from 1 to " << graph_.getSize()
            << " to start search " << end1 << std::endl;
  int vertex;
  std::cin >> vertex;
  try {
    std::vector<int> vecWay;
    vecWay = algo_.breadthFirstSearch(graph_, vertex);
    printVector(vecWay);
    nextChange();
  } catch (std::exception &e) {
    std::cerr << style4 << "\nError: " << e.what() << end << "\n";
    DFS();
  }
}

void ConsoleInterface::dijkstra() {
  std::cout
      << style4
      << "\nenter the vertex number of the beginning of the path from 1 to "
      << graph_.getSize() << " to start search " << end1 << std::endl;
  int start;
  std::cin >> start;
  std::cout << style4
            << "\nenter the vertex number of the end of the path from 1 to "
            << graph_.getSize() << " to start search " << end1 << std::endl;
  int finish;
  std::cin >> finish;
  try {
    int distance = algo_.getShortestPathBetweenVertices(graph_, start, finish);
    std::vector<double> aa = algo_.algoD(graph_, start, finish);
    printVector(aa);
    std::cout << style4 << "\ndistance: " << distance << end1 << std::endl;
    nextChange();
  } catch (std::exception &e) {
    std::cerr << style4 << "\nError: " << e.what() << end << "\n";
    dijkstra();
  }
}

void ConsoleInterface::floydWarshallAlgo() {
  try {
    S21Matrix<double> rez = algo_.getShortestPathsBetweenAllVertices(graph_);
    printMatrix(rez);
    nextChange();
  } catch (std::exception &e) {
    std::cerr << style4 << "\nError: " << e.what() << end1 << "\n";
    nextChange();
  }
}

void ConsoleInterface::primAlgorithmAlgo() {
  try {
    S21Matrix<double> rez = algo_.getLeastSpanningTree(graph_);
    printMatrix(rez);
    nextChange();
  } catch (std::exception &e) {
    std::cerr << style4 << "\nError: " << e.what() << end1 << "\n";
    nextChange();
  }
}

void ConsoleInterface::ant() {
  try {
    TsmResult rez;
    rez = algo_.solveTravelingSalesmanProblem(graph_);
    printVector(rez.visit);
    std::cout << style4 << "Distance: " << rez.distance << end1 << std::endl;
    nextChange();
  } catch (std::exception &e) {
    std::cerr << style4 << "\nError: " << e.what() << end1 << "\n";
    nextChange();
  }
}

void ConsoleInterface::exportDot() {
  std::cout << "\n"
            << style3 << " enter file name to save              " << end1
            << std::endl;
  std::string filename;
  std::cin >> filename;
  filename = "../datasets/" + filename;
  try {
    graph_.exportGraphToDot(filename);
    std::cout << "\n"
              << style3 << " file successfully saved in Graphs directory "
              << end1 << std::endl;
    std::string command = "dot -Tpng tmp.dot -o";
    command += filename;
    const char *c_com = command.c_str();
    system(c_com);
    nextChange();
  } catch (std::exception &e) {
    std::cerr << style4 << "\nError: " << e.what() << end1 << "\n";
    nextChange();
  }
}

void ConsoleInterface::nextChange() {
  std::cout << style3 << " ----> to the menu press 1               " << end
            << std::endl;
  std::cout << style3 << " ----> to exit the program press 0                 "
            << end1 << std::endl;
  int choice;
  std::cin >> choice;
  if (choice == 1) {
    menu();
  } else if (choice == 0) {
    return;
  } else {
    std::cout << style4 << "\nError: incorrect input, try again " << end
              << std::endl;
    nextChange();
  }
}

template <typename T>
void ConsoleInterface::printVector(std::vector<T> path) {
  for (auto it : path) {
    std::cout << it << " ";
  }
  std::cout << "\n" << end1;
}

template <typename T>
void ConsoleInterface::printMatrix(S21Matrix<T> res) {
  for (int i = 0; i < res.get_rows(); ++i) {
    for (int j = 0; j < res.get_cols(); ++j) {
      std::cout << res(i, j);
      if (res(i, j) < 10)
        std::cout << "   ";
      else if (res(i, j) >= 10 && res(i, j) < 100)
        std::cout << "  ";
      else
        std::cout << " ";
    }
    std::cout << std::endl;
  }
  std::cout << "\n" << end;
}

void ConsoleInterface::antVsOtherAlgo() {
  try {
    std::cout << style4
              << "\nenter how many times to keep track of the time (1 - 1000)"
              << end1 << std::endl;
    int N;
    std::cin >> N;
    if (N > 0 && N < 1001) {
      double timeAnt = 0.0;
      double timeAnneAlingAlgorithm = 0.0;
      double timeBranchAndBound = 0.0;
      TsmResult antRez;
      TsmResult anneAlingRez;
      TsmResult branchAndBoundRez;
      clock_t start, end;

      start = clock();
      for (int i = 0; i < N; ++i) {
        if (i == 0) {
          antRez = algo_.helperSolveTravel(graph_);
        } else {
          TsmResult tmp = algo_.helperSolveTravel(graph_);
          if (tmp.distance < antRez.distance) {
            antRez = tmp;
          }
        }
      }
      end = clock();
      timeAnt = ((double)end - start) / ((double)CLOCKS_PER_SEC);

      start = clock();
      for (int i = 0; i < N; ++i) {
        if (i == 0) {
          anneAlingRez = algo_.annealingAlgorithm(graph_);
        } else {
          TsmResult tmp = algo_.annealingAlgorithm(graph_);
          if (tmp.distance < anneAlingRez.distance) {
            anneAlingRez = tmp;
          }
        }
      }
      end = clock();
      timeAnneAlingAlgorithm = ((double)end - start) / ((double)CLOCKS_PER_SEC);

      start = clock();
      for (int i = 0; i < N; ++i) {
        if (i == 0) {
          branchAndBoundRez = algo_.SolveTspBranchAndBound(graph_);
        } else {
          TsmResult tmp = algo_.SolveTspBranchAndBound(graph_);
          if (tmp.distance < branchAndBoundRez.distance) {
            branchAndBoundRez = tmp;
          }
        }
      }
      end = clock();
      timeBranchAndBound = ((double)end - start) / ((double)CLOCKS_PER_SEC);

      std::cout << style4 << "1) Ant Algoritm\n";
      std::cout << style4 << "REPS: " << N << std::endl;
      std::cout << style4 << "TIME: " << timeAnt
                << "\nDISTANCE: " << antRez.distance << end1 << std::endl;

      std::cout << style4 << "2) Annealing Algoritm\n";
      std::cout << style4 << "REPS: " << N << std::endl;
      std::cout << style4 << "TIME: " << timeAnneAlingAlgorithm
                << "\nDISTANCE: " << anneAlingRez.distance << end1 << std::endl;

      std::cout << style4 << "3) Branch and Bound Algoritm\n";
      std::cout << style4 << "REPS: " << N << std::endl;
      std::cout << style4 << "TIME: " << timeBranchAndBound
                << "\nDISTANCE: " << branchAndBoundRez.distance << end1
                << std::endl;

    } else {
      std::cout << style4 << "\nINCORRECT COUNT OF REPS!" << end1 << std::endl;
    }
    nextChange();
  } catch (std::exception &exceptionText) {
    std::cerr << style4 << "\nError: " << exceptionText.what() << end1 << "\n";
    nextChange();
  }
}
