namespace s21 {
// Construct
template <typename T>
stack<T>::stack() : _size(0), _head(nullptr) {}

template <typename T>
stack<T>::stack(std::initializer_list<T> const &items)
    : _size(0), _head(nullptr) {
  for (auto it = items.begin(); it != items.end(); it++) {
    this->push(*it);
  }
}

template <typename T>
stack<T>::stack(stack const &s) : _size(0), _head(nullptr) {
  *this = s;
}

template <typename T>
stack<T>::stack(stack &&s) : _size(0), _head(nullptr) {
  *this = std::move(s);
}

// Operator
template <class T>
stack<T> &stack<T>::operator=(stack<T> const &s) {
  if (this != &s && s._size) {
    Node<T> *bufferNode = s._head;
    while (bufferNode) {
      push(bufferNode->getData());
      bufferNode = bufferNode->getPastElement();
    }
    reverse();
  }

  return *this;
}

template <class T>
stack<T> &stack<T>::operator=(stack<T> &&s) {
  if (this != &s) {
    clean();
    std::swap(_size, s._size);
    std::swap(_head, s._head);
  }
  return *this;
}

//  Destructor
template <typename T>
stack<T>::~stack() {
  clean();
}

//  Methods
template <typename T>
void stack<T>::push(const T &value) {
  if (_size > max_size()) {
    throw std::out_of_range("Can`t push! Stack is full");
  }

  Node<T> *bufferNode = new Node<T>(value);

  if (_size) {
    bufferNode->setPastElement(_head);
  }

  setRoot(bufferNode);
  ++_size;
}

template <typename T>
const T &stack<T>::peek() {
  if (!_size) {
    throw std::out_of_range("Stack is empty");
  }

  return this->_head->Node<T>::getData();
}

template <typename T>
T stack<T>::pop() {
  if (!_size) {
    throw std::out_of_range("Can`t pop! Stack is empty");
  }
  T tmp = _head->getData();
  Node<T> *buffer = _head->getPastElement();
  delete _head;
  _head = nullptr;
  this->setRoot(buffer);
  _size--;
  return tmp;
}

template <typename T>
bool stack<T>::empty() {
  return _size == 0;
}

template <typename T>
size_t stack<T>::size() {
  return _size;
}

template <class T>
void stack<T>::setRoot(Node<T> *element) {
  this->_head = element;
}

template <class T>
void stack<T>::swap(stack<T> &other) {
  if (this != &other) {
    std::swap(_head, other._head);
    std::swap(other._size, _size);
  }
}

template <typename T>
size_t stack<T>::max_size() {
  return std::numeric_limits<size_t>::max() / (sizeof(Node<T>));
}

template <typename T>
bool stack<T>::contains(T const &value) {
  Node<T> *buffer = _head;
  while (buffer != nullptr) {
    if (buffer->getData() == value) return true;
    buffer = buffer->getPastElement();
  }
  return false;
}

//  Overload

template <class T>
template <typename... Args>
void stack<T>::emplace_front(Args &&...args) {
  const size_t sizeArgs = sizeof...(Args);
  T argsArray[sizeArgs] = {args...};

  for (size_t i = 0; i < sizeArgs; ++i) {
    push(argsArray[i]);
  }
}

// Help Foo
template <typename T>
void stack<T>::clean() {
  Node<T> *buffer = _head;
  while (_head != nullptr) {
    buffer = buffer->getPastElement();
    delete _head;
    _head = nullptr;
    this->setRoot(buffer);
  }
  _size = 0;
}

template <typename T>
void stack<T>::reverse() {
  if (_size > 0) {
    stack<T> tmp;
    Node<T> *buf = _head;
    while (buf != nullptr) {
      tmp.push(buf->getData());
      buf = buf->getPastElement();
    }
    *this = std::move(tmp);
  }
}

}  //  namespace s21
