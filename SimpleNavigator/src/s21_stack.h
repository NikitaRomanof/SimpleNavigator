#pragma once

#include <iostream>
#include <limits>

#include "s21_node.h"

namespace s21 {
template <class T>
class stack {
 public:
  stack();
  stack(std::initializer_list<T> const &items);
  stack(const stack &s);
  stack(stack &&s);
  stack<T> &operator=(stack<T> &&s);
  stack<T> &operator=(stack<T> const &s);
  ~stack();

  const T &peek();
  bool empty();
  size_t size();
  void push(const T &value);
  T pop();
  void swap(stack<T> &other);

  template <typename... Args>
  void emplace_front(Args &&...args);
  void reverse();
  bool contains(T const &value);

 private:
  void setRoot(Node<T> *element);
  size_t max_size();
  void clean();

  size_t _size;
  Node<T> *_head;
};
}  //  namespace s21

#include "s21_stack.inl"
