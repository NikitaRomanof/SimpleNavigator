namespace s21 {

// Construct
template <typename T>
queue<T>::queue() : _size(0), _head(nullptr), _tail(nullptr) {}

template <typename T>
queue<T>::queue(const std::initializer_list<T> &items)
    : _size(0), _head(nullptr), _tail(nullptr) {
  for (auto it = items.begin(); it != items.end(); it++) {
    this->push(*it);
  }
}

template <typename T>
queue<T>::queue(const queue &s) : _size(0), _head(nullptr), _tail(nullptr) {
  *this = s;
}

template <typename T>
queue<T>::queue(queue &&s) : _size(0), _head(nullptr), _tail(nullptr) {
  *this = std::move(s);
}

//  Operator
template <class T>
queue<T> &queue<T>::operator=(const queue<T> &s) {
  if (this != &s && s._size > 0) {
    clean();
    Node<T> *bufferNode = s._head;
    while (bufferNode != nullptr) {
      push(bufferNode->getData());
      bufferNode = bufferNode->getPastElement();
    }
    _size = s._size;
  }
  return *this;
}

template <class T>
queue<T> &queue<T>::operator=(queue<T> &&s) {
  if (this != &s) {
    clean();
    swap(s);
  }
  return *this;
}

//  Destructor
template <typename T>
queue<T>::~queue() {
  clean();
}

//  Foonction
template <typename T>
void queue<T>::push(const T &value) {
  if (_size > max_size()) {
    throw std::out_of_range("Can`t push! Queue is full");
  }
  Node<T> *bufferNode = new Node<T>(value);
  if (_size) {
    _tail->setPastElement(bufferNode);
  }
  setTail(bufferNode);
  if (!_size) {
    setRoot(bufferNode);
  }
  _size++;
}

template <typename T>
const T &queue<T>::front() {
  if (!_size) {
    throw std::out_of_range("Queue is empty");
  }

  return _head->getData();
}

template <typename T>
const T &queue<T>::back() {
  if (!_size) {
    throw std::out_of_range("Queue is empty");
  }

  return this->_tail->getData();
}

template <typename T>
T queue<T>::pop() {
  if (!_size) {
    throw std::out_of_range("Queue is empty");
  }

  Node<T> *buffer = _head->getPastElement();
  T tmp = _head->getData();
  delete _head;
  _head = nullptr;
  this->setRoot(buffer);
  _size--;
  return tmp;
}

template <typename T>
const T &queue<T>::peek() {
  return front();
}

template <typename T>
bool queue<T>::empty() {
  return _size == 0;
}

template <typename T>
size_t queue<T>::size() {
  return _size;
}

template <class T>
void queue<T>::setRoot(Node<T> *element) {
  this->_head = element;
}

template <class T>
void queue<T>::setTail(Node<T> *element) {
  this->_tail = element;
}

template <class T>
void queue<T>::swap(queue<T> &other) {
  if (this != &other) {
    std::swap(_head, other._head);
    std::swap(_tail, other._tail);
    std::swap(_size, other._size);
  }
}

template <class T>
template <typename... Args>
void queue<T>::emplace_back(Args &&...args) {
  const size_t sizeArgs = sizeof...(Args);
  T argsArray[sizeArgs] = {args...};

  for (size_t i = 0; i < sizeArgs; ++i) {
    push(argsArray[i]);
  }
}

template <typename T>
size_t queue<T>::max_size() {
  return (std::numeric_limits<size_t>::max() / (sizeof(Node<T>))) / 2;
}

// Help Foo
template <typename T>
void queue<T>::reverse() {
  if (_size) {
    _tail = _head;
    Node<T> *previousNode = nullptr, *currentNode = nullptr,
            *nextNode = nullptr;
    currentNode = _head->getPastElement();
    while (currentNode != nullptr) {
      nextNode = currentNode->getPastElement();
      currentNode->setPastElement(previousNode);
      previousNode = currentNode;
      currentNode = nextNode;
    }

    _head = previousNode;
  }
}

template <typename T>
void queue<T>::clean() {
  while (_head != nullptr) {
    Node<T> *buffer = _head->getPastElement();
    delete _head;
    _head = nullptr;
    this->setRoot(buffer);
    _size--;
  }
}

template <typename T>
bool queue<T>::contains(const T &val) {
  Node<T> *buffer = _head;
  while (buffer != nullptr) {
    if (buffer->getData() == val) return true;
    buffer = buffer->getPastElement();
  }
  return false;
}

}  //  namespace s21
