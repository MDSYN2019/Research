

/*

- Motivation and uses for TMP

- Type traits and type checking

- Member detection

- Traits classes and algorithm selection

- If time permits: compile-time computation


----------------------------
C++ Template Metaprogramming 
----------------------------

Metaprogramming sheds light on the most powerful idioms 
of today's C++, at long last delivering practical 
metaprogramming tools and techniques into the hands of everyday
programmer.


A metaprogram is a program that generates or manipulates program 
code. Ever since generic programming was introduced to C++, 
programmers have discovered a myriad of template tricks for 
manipulating programs as they are compiled, effectively eliminating 
the barrier between program and metaprogram.

-- 

What is a metaprogram?
---------------------

Metaprogram - A program about about a program

A metaprogram is a program that manipulates code. It may be an odd-sounding 
concept, but you're probably already familiar with several such beasts. Your C++ compiler is one example: It manipulates your C++ code to produce assembly language or machine code.


Parser generators such as YACC are another kind of program-manipulating 
program. The input to YACC is a high-level parser description written 
in terms of grammar rules and attached actions brace-enclosed. 


For instance, to parse and evalulate arithemtic expressions with the usual 
precedence rules, we might feed YACC the following grammar description.

--

expression: term
          | expression '+' term {$$ = $1 + $3; }
          | expression '-' term {$$ = $1 - $3; }

term: factor 
          

One of the nice things about template metaprograms is a property they share 
with good old traditinoal systems: Once a metaprogram is written, it can be 
used without knowing what's under the hood as long as it works.

YACC is an example of a translator metaprogram whose domain language differs 
from its host language.

Amazingly, if you have a C++ compiler, this is precisely the kind of metaprogramming 
power you hold in your fingertips. The rest of this book is about unlocking that power
and showing how and when to use it.

----------------------
Metaprogramming in C++
----------------------

Numeric Computations
--------------------

The earliest C++ metaprograms performed integer computations at compile time.
One of the ver first metaprograms was shown at a C++ commitee meeting

Since illegal code is hard to use effectively in a larger system, lets examine a slightly more practical application 

The program below transliterates unsigned decimal numbers into their binary equivalents, allowing us to express 
constants in a recognizable form

 */

#include <iostream>



template <unsigned long N>
struct binary {
  static unsigned const value
  = binary<N/10>::value << 1 // prepend higher bits to 
			   | N%10;   // lower bits
};

template <> // Specialization
struct binary<0> { // terminates recursion
  static unsigned const value = 0;
};

/*
If you're wondering " Where's the program?", we ask you to consider what happens 
when we access the nested ::value member of binary<N>. The binary template is instantiated
again with a smaller N, until N reaches 0 and the specialization is used as a termination 
condition.

The process should have the familiar flavor of recursive call and what is a program, after all, 
but a function?

Essentially, the compiler is being used to interpret our little metaprogram.

 */

unsigned const one = binary<1>::value;
unsigned const three = binary<11>::value;
unsigned const five = binary<10>::value;
unsigned const seven = binary<11>::value;
unsigned const nine = binary<1001>::value;


/*
  ----------------
  Type Computation
  ----------------
  
  Much more important than its ability to do compile time numeric computation is C++'s ability
  to compute with types. As a matter of fact, type computation will dominate the rest of this book,
  and we'll cover examples of it in the very first section of the next chapter

  Although you may have read, to understand the specifics of type computations, we'd 
  like to give you a sense of its power. Remember our YACC expression evaluator? 
  It turns out we don't need to use a translator to get that kind of power and convenience.

  With appopriate surrounding code from the boost spirit library, the following legal
  C++ code has equivalent functionality.
  
 */

/*
expr = ( term[expr.val = _1] >> '+' >> expr[expr.val += _1] )
  |    ( term[expr.val = _1] >> '-' >> expr[expr.val -= _1] )
  |     term[expr.val = _1]; 
*/


/*
  The purpose of the code - https://www.fluentcpp.com/2017/06/02/write-template-metaprogramming-expressively/

  For example given a type T, we would like to know whether T is incrementable, that is say that 
  for an object t of type T, whether or not the expression:

  ++t

  is valid. If T is int, then the expression is valid, and if T is std::string, then the expression is not valid

 */

template <typename, typename = void>
struct is_incrementable : std::false_type {};

template<typename T >
struct is_incrementable<T,
			std::void_t<decltype(++std::declval<T&>() )>
			> : std::true_type {};

int main(void) {
  return 0;
}
