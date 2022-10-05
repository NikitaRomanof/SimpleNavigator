#pragma once
#include <iostream>

namespace s21 {
template <class T>
class Node {
 public:
  explicit Node(const T &value);
  Node(Node<T> const &other);
  Node(Node<T> &&other);
  ~Node();
  Node<T> &operator=(Node<T> const &other);
  Node<T> &operator=(Node<T> &&other);

  const T &getData();
  void setPastElement(Node *element);
  Node *getPastElement();

 private:
  Node<T> *_pastElement;
  T _data;
  void clean();
};
}  //  namespace s21

#include "s21_node.inl"
