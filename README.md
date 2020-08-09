# Prosas march

## Intro
"*Prosas march*" stands for **pro**bability **s**olver of **a**bsorbing **s**tates in **M**arkov **c**hain. This is my solution for an online programming game. 
As the title ~~clearly~~ indicates, the simple code in this repository provides a method to compute probabilities for the absorbing states in a given Markov chain. This solution will process matrices of fractions in order to avoid float approximations, for this reason it doesn't use external libraries such as numpy.


## Overview
The method `prosasmarch.solve(n)` gets a list of list `m` as input, representing the transition matrix of a Markov Chain, and returns an array containing the probability of each absorbing state or None if there is no abs. state. If there are K absorbing states, the returned array will contain K+1 integer numbers, where the first K items are the numerators of the abs. states' probabilities, and the last one is the denominator


Things to keep in mind

* An absorbing state with transition prob. of 1.0 towards itself is equivalent to a terminal state.
* The first state (first row of the input) is considered as the starting one.

## Examples

Some simple examples. Also, try the `demo.py` script.

### Example 1
Let's consider the following Markov chain:

![markov1](https://github.com/ianovir/Prosas_march/blob/master/pics/p1.JPG)

The state s0 is considered the starting state, while the absorbing states are: **s2**, **s3** and **s4**.

To solve this chain let's call:

```python
prosasmarch.solve([[0, 2, 1, 0, 0], [0, 0, 0, 3, 4], [0, 0, 0, 0, 0], [0, 0, 0, 0,0], [0, 0, 0, 0, 0]])
```

The output will be:

```
[7, 6, 8, 21]
```

which means that the probabilities for each absorbing state to be reached from **s0** are: 7/21 for **s2**, 6/21 for **s3** and 8/21 for **s4**.


### Example 2
Let's consider the following Markov chain:

![markov2](https://github.com/ianovir/Prosas_march/blob/master/pics/p2.JPG)

The state s0 is the starting state, while the absorbing states are: **s2**, **s3**, **s4** and **s5**. Please, note that **s2** is an **unreachable** state.

To solve this chain let's call:

```python
prosasmarch.solve([[0, 1, 0, 0, 0, 1], [4, 0, 0, 3, 2, 0], [0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0]])
```

The output will be:

```
[0, 3, 2, 9, 14]
```

which means that the probabilities for each absorbing state to be reached from **s0** are: 0 for **s2** (unreachable), 3/14 for **s3**, 2/14 for **s4** and 9/14 for **s5**.


# Copyright
Copyright(c) 2020 Sebastiano Campisi - [ianovir.com](https://ianovir.com). 
Read the LICENSE file for more details.
