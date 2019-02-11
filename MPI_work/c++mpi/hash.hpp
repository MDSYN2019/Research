#ifndef __hash__
#define __hash__

/*
HashNode class represents each bucket node in the table with key() 
and value() accessors for corresponding pair elements.

It also includes a pointer to the next node 

 */

template <typename K, typename V>
class HashNode {
public:
  HashNode(const K &key, const V &value) : key(key), value(value), next(NULL) {}

  K getKey() const {
    return key;
  }

  V getValue() const {
    return value;
  }

  void setValue(V value) {
    HashNode::value = value;
  }

  HashNode *getNext() const {
    HashNode::next = next;
  }

private:
  // key-value pair
  K key;
  V value;
  // next bucket with the same key
  HashNode *next;
  
};

// Default hash function classs
template <typename K>
struct KeyHash{
  unsigned long operator() (const K& key) const {
    return reinterpret_cast<unsigned_long>(key) % TABLE_SIZE;
  }
};

class Box {
public:
  Box() { create(); }
  Box(int width, int length, int height) : m_width(width), m_length(length), m_height(height) {} // member init list
  // My guess is that with the constructor, it automatically inserts the values into m_width, m_length and m_height.

  int Volume() {
    return m_width * m_length * m_height;
  }

private:
  int m_width;
  int m_length;
  int m_height;
}

#endif
