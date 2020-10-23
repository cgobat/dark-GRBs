import numpy as np

class AsymmetricUncertainty:
    """
    Class for handling and naive propagation of asymmetric uncertainties.
    Author: Caden Gobat, George Washington University

    Parameters
    ----------
    nominal : numeric
        the nominal value of the represented quantity
    pos_err : numeric
        the plus error on the value, as in value = nominal (+pos_err, -neg_err)
    neg_err : numeric
        the minus error on the value, as in value = nominal (+pos_err, -neg_err)

    Returns
    -------
    string
        a value in a string

    Raises
    ------
    ValueError
        if you don't pass numeric values for the arguments
    TypeError
        if you don't pass numeric values for the arguments
    """
    
    def __init__(self, nominal, pos_err, neg_err):
        self.value = float(nominal)
        self.plus = float(pos_err)
        self.minus = float(neg_err)
        
    def __str__(self):
        if np.isclose(self.plus, self.minus):
            return "{} ± {}".format(self.value,self.plus)
        else:
            return "{} (+{}, -{})".format(self.value,self.plus,self.minus)
        
    def _repr_latex_(self):
        if np.isclose(self.plus, self.minus):
            return "$%f \pm %f$" %(self.value,self.plus)
        else:
            return "$%f_{-%f}^{+%f}$" %(self.value,self.minus,self.plus)
        
    def __add__(self,other):
        if isinstance(other,AsymmetricUncertainty):
            pass
        else:
            other = AsymmetricUncertainty(other,0,0)
        maximum = (self.value+self.plus) + (other.value+other.plus)
        minimum = (self.value-self.minus) + (other.value-other.minus)
        result = self.value + other.value
        pos = maximum-result
        neg = result-minimum
        return AsymmetricUncertainty(result,pos,neg)
    
    def __sub__(self,other):
        if isinstance(other,AsymmetricUncertainty):
            pass
        else:
            other = AsymmetricUncertainty(other,0,0)
        maximum = (self.value+self.plus) - (other.value-other.minus)
        minimum = (self.value-self.minus) - (other.value+other.plus)
        result = self.value - other.value
        pos = maximum-result
        neg = result-minimum
        return AsymmetricUncertainty(result,pos,neg)
    
    def __mul__(self,other):
        if isinstance(other,AsymmetricUncertainty):
            pass
        else:
            other = AsymmetricUncertainty(other,0,0)
        maximum = (self.value+self.plus) * (other.value+other.plus)
        minimum = (self.value-self.minus) * (other.value-other.minus)
        result = self.value * other.value
        pos = maximum-result
        neg = result-minimum
        return AsymmetricUncertainty(result,pos,neg)
    
    def __truediv__(self,other):
        if isinstance(other,AsymmetricUncertainty):
            pass
        else:
            other = AsymmetricUncertainty(other,0,0)
        maximum = (self.value+self.plus) / (other.value-other.minus)
        minimum = (self.value-self.minus) / (other.value+other.plus)
        result = self.value / other.value
        pos = maximum-result
        neg = result-minimum
        return AsymmetricUncertainty(result,pos,neg)
    
    def __rpow__(self,other):
        return other**self.value