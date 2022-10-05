#pragma once

#include <iostream>
#include <limits>

#include "s21_node.h"

namespace s21 {
template <class T>
class queue {
 public:
  queue();
  queue(const std::initializer_list<T> &items);
  queue(const queue &s);
  queue(queue &&s);
  queue<T> &operator=(const queue<T> &s);
  queue<T> &operator=(queue<T> &&s);
  ~queue();

  const T &front();
  const T &back();
  const T &peek();

  bool empty();
  size_t size();
  void push(const T &value);
  T pop();
  void swap(queue<T> &other);
  template <typename... Args>
  void emplace_back(Args &&...args);
  void reverse();
  size_t max_size();
  bool contains(const T &val);

 private:
  size_t _size;
  Node<T> *_head;
  Node<T> *_tail;

  void setRoot(Node<T> *element);
  void setTail(Node<T> *element);
  void clean();
};
}  //  namespace s21

#include "s21_queue.inl"
