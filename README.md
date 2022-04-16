# pyHEDGES
Here is the rough implementation of [HEDGES](https://github.com/whpress/hedges) in the Python platform.
 
In order to accept arbitrary constraints, this rough implementation is two different from the original design:

* mapping pattern will float 
with the current available nucleotides (called #C^*_i in the 
[paper](https://www.pnas.org/doi/abs/10.1073/pnas.2004821117)).
* 6 child hypotheses are expanded to more child hypotheses 
to suppose the unknown number of bits in the current message bit.