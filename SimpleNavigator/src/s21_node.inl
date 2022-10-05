namespace s21 {
//  Construct
template <typename T>
Node<T>::Node(const T &value) : _data(value), _pastElement(nullptr) {}

template <typename T>
Node<T>::Node(Node<T> const &other) {
  *this = other;
}

template <typename T>
Node<T>::Node(Node<T> &&other) {
  *this = std::move(other);
}

//  Destruct
template <typename T>
Node<T>::~Node() {
  clean();
}

//  Operator
template <typename T>
Node<T> &Node<T>::operator=(Node<T> const &other) {
  if (this != &other) {
    clean();
    if (_pastElement != nullptr) {
      delete _pastElement;
      _pastElement = nullptr;
    }
    _pastElement = other._pastElement;
    _data = other._data;
  }
  return *this;
}
template <typename T>
Node<T> &Node<T>::operator=(Node<T> &&other) {
  if (this != &other) {
    clean();
    if (_pastElement != nullptr) {
      delete _pastElement;
      _pastElement = nullptr;
    }
    std::swap(_pastElement, other._pastElement);
    std::swap(_data, other._data);
  }
  return *this;
}

//  Getes
template <class T>
Node<T> *Node<T>::getPastElement() {
  return _pastElement;
}

template <typename T>
void Node<T>::setPastElement(Node *element) {
  this->_pastElement = element;
}

template <class T>
const T &Node<T>::getData() {
  return this->_data;
}

// Help Foo
template <typename T>
void Node<T>::clean() {
  _data = T();
}
}  //  namespace s21
